#pragma once

#include <cassert>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "util.hpp"

namespace sshash {

/*
    Winner-tree-based external-merge iterator over N sorted runs on disk.

    Each run is read with a small buffered std::ifstream (no mmap) so that
    process RSS stays bounded by `num_files_to_merge * buffer_records *
    sizeof(T)` regardless of total run size. Values are surfaced as
    `T const&` from each stream's in-RAM buffer; the merge logic compares
    those values directly instead of pointers into mmap'd memory.

    Required of T:
      - copy-constructible / move-constructible,
      - `static T T::max()` returning a strict upper bound (used as the
        sentinel for exhausted streams in the winner tree),
      - `bool operator<(T, T)`.
*/
template <typename T>
struct file_merging_iterator  //
{
    static constexpr uint64_t default_buffer_records = uint64_t(1) << 12;  // 4096 records
    const uint64_t scan_threshold = 16;

    template <typename FileNamesIterator>
    file_merging_iterator(FileNamesIterator file_names_iterator, uint64_t num_files_to_merge,
                          uint64_t buffer_records = default_buffer_records)  //
    {
        if (num_files_to_merge == 0) {
            m_num_files_to_merge = 0;
            return;
        }

        m_streams.reserve(num_files_to_merge);
        for (uint64_t i = 0; i != num_files_to_merge; ++i, ++file_names_iterator) {
            m_streams.emplace_back();
            m_streams.back().open(*file_names_iterator, buffer_records);
        }

        m_num_files_to_merge = num_files_to_merge;
        m_min_idx = 0;
        if (m_streams.size() <= scan_threshold) {
            compute_min();
        } else {
            /* build a winner tree (same shape as before, but the leaves
               index into m_streams instead of carrying raw pointers). */
            uint64_t n = num_files_to_merge;
            uint64_t m = 2 * n - 1;
            m_size = n;
            m_tree.resize(m);
            m_begin = (1ULL << uint64_t(std::ceil(std::log2(n)))) - 1;
            uint64_t i = 0;
            for (; m_begin + i != m; ++i) m_tree[m_begin + i] = i;
            for (uint64_t j = 0; i != n; ++i, ++j) m_tree[n - 1 + j] = i;
            m_min_idx = build(0);
        }
    }

    bool has_next() { return m_num_files_to_merge != 0; }
    void next() { update(); }
    T operator*() const { return m_streams[m_min_idx].current(); }

    void close() {
        for (auto& s : m_streams) s.close();
        m_streams.clear();
        m_streams.shrink_to_fit();
        m_tree.clear();
        m_tree.shrink_to_fit();
    }

private:
    /*
        A buffered, forward-only reader over a single run file. Reads in
        chunks of `m_buf.size()` records via std::ifstream and presents a
        T-by-reference current-value interface.
    */
    struct buffered_stream {
        buffered_stream() = default;
        buffered_stream(buffered_stream const&) = delete;
        buffered_stream& operator=(buffered_stream const&) = delete;
        buffered_stream(buffered_stream&&) = default;
        buffered_stream& operator=(buffered_stream&&) = default;

        void open(std::string const& filename, uint64_t buffer_records) {
            m_buf.resize(std::max<uint64_t>(1, buffer_records));
            m_in.open(filename, std::ifstream::binary);
            if (!m_in.is_open()) {
                throw std::runtime_error("cannot open run file '" + filename + "'");
            }
            m_pos = 0;
            m_size = 0;
            m_eof = false;
            refill();
        }

        void close() {
            if (m_in.is_open()) m_in.close();
            m_buf.clear();
            m_buf.shrink_to_fit();
            m_pos = 0;
            m_size = 0;
            m_eof = true;
        }

        bool empty() const { return m_pos >= m_size; }

        T const& current() const {
            assert(!empty());
            return m_buf[m_pos];
        }

        void advance() {
            assert(!empty());
            ++m_pos;
            if (m_pos >= m_size && !m_eof) refill();
        }

    private:
        std::ifstream m_in;
        std::vector<T> m_buf;
        uint64_t m_pos = 0;
        uint64_t m_size = 0;
        bool m_eof = true;

        void refill() {
            m_pos = 0;
            m_in.read(reinterpret_cast<char*>(m_buf.data()),
                      static_cast<std::streamsize>(m_buf.size() * sizeof(T)));
            const std::streamsize got = m_in.gcount();
            m_size = static_cast<uint64_t>(got) / sizeof(T);
            if (m_size == 0) m_eof = true;
        }
    };

    std::vector<buffered_stream> m_streams;
    std::vector<uint32_t> m_tree;

    uint64_t m_begin = 0, m_size = 0;
    uint64_t m_min_idx = 0, m_num_files_to_merge = 0;

    void update() {
        if (m_streams.size() <= scan_threshold) {
            auto& s = m_streams[m_min_idx];
            s.advance();
            if (s.empty()) {
                m_streams.erase(m_streams.begin() + m_min_idx);
                m_min_idx = 0;
                --m_num_files_to_merge;
                if (m_num_files_to_merge == 0) return;
            }
            compute_min();
        } else {  // winner tree
            m_min_idx = m_tree[0];
            assert(m_min_idx < m_streams.size());
            auto& s = m_streams[m_min_idx];
            s.advance();
            uint64_t p = m_begin + m_min_idx;
            p -= (p >= m_tree.size()) * m_size;  // p is the leaf index
            if (s.empty()) {
                m_tree[p] = uint32_t(-1);
                --m_num_files_to_merge;
            }
            const T inf = T::max();
            while (p) {
                uint64_t is_r_child = (p & 1) == 0;
                uint32_t l = m_tree[p - is_r_child];
                uint32_t r = m_tree[p + 1 - is_r_child];
                T const& vl = (l == uint32_t(-1)) ? inf : m_streams[l].current();
                T const& vr = (r == uint32_t(-1)) ? inf : m_streams[r].current();
                uint32_t i = (vl < vr) ? l : r;
                uint64_t parent = (p - 1) / 2;
                m_tree[parent] = i;
                p = parent;
            }
            m_min_idx = m_tree[0];
        }
    }

    uint32_t build(uint32_t p) {
        if (p >= m_tree.size()) return uint32_t(-1);
        if (p >= m_size - 1) return m_tree[p];  // leaf
        uint32_t l = build(2 * p + 1);
        uint32_t r = build(2 * p + 2);
        const T inf = T::max();
        T const& vl = (l == uint32_t(-1)) ? inf : m_streams[l].current();
        T const& vr = (r == uint32_t(-1)) ? inf : m_streams[r].current();
        uint32_t i = (vl < vr) ? l : r;
        m_tree[p] = i;
        return i;
    }

    void compute_min() {
        m_min_idx = 0;
        T min_val = m_streams.front().current();
        for (uint64_t i = 1; i != m_streams.size(); ++i) {
            assert(!m_streams[i].empty());
            T const& val = m_streams[i].current();
            if (val < min_val) {
                min_val = val;
                m_min_idx = i;
            }
        }
    }
};

}  // namespace sshash
