#pragma once

namespace sshash {

template <typename T>
struct file_merging_iterator  //
{
    const uint64_t scan_threshold = 16;

    template <typename FileNamesIterator>
    file_merging_iterator(FileNamesIterator file_names_iterator, uint64_t num_files_to_merge)
        : m_mm_files(num_files_to_merge)  //
    {
        if (num_files_to_merge == 0) return;

        /* open files and create the input iterators */
        m_iterators.reserve(num_files_to_merge);
        for (uint64_t i = 0; i != num_files_to_merge; ++i, ++file_names_iterator) {
            m_mm_files[i].open(*file_names_iterator, mm::advice::sequential);
            m_iterators.push_back(
                {m_mm_files[i].data(), m_mm_files[i].data() + m_mm_files[i].size()});
        }

        m_num_files_to_merge = num_files_to_merge;
        m_min_idx = 0;
        if (m_iterators.size() <= scan_threshold) {
            compute_min();
        } else {
            /* build a looser tree */
            uint64_t n = num_files_to_merge;
            uint64_t m = 2 * n - 1;
            m_size = n;
            m_tree.resize(m);
            m_begin = (1ULL << uint64_t(std::ceil(std::log2(n)))) - 1;
            uint64_t i = 0;
            for (; m_begin + i != m; ++i) m_tree[m_begin + i] = i;
            for (uint64_t j = 0; i != n; ++i, ++j) m_tree[n - 1 + j] = i;
            build(0);
            m_min_idx = m_tree[0];
        }
    }

    bool has_next() { return m_num_files_to_merge != 0; }
    void next() { update(); }
    T operator*() const { return *(m_iterators[m_min_idx].begin); }

    void close() {
        for (auto& mm_file : m_mm_files) mm_file.close();
        m_iterators.clear();
        m_mm_files.clear();
        m_tree.clear();
    }

private:
    struct pointer_pair {
        T const* begin;
        T const* end;
    };
    std::vector<pointer_pair> m_iterators;
    std::vector<mm::file_source<T const>> m_mm_files;
    std::vector<uint32_t> m_tree;

    uint64_t m_begin, m_size;
    uint64_t m_min_idx, m_num_files_to_merge;

    void update() {
        if (m_iterators.size() <= scan_threshold) {  // compute min with a linear scan
            auto& it = m_iterators[m_min_idx];
            it.begin += 1;
            if (it.begin == it.end) {
                m_iterators.erase(m_iterators.begin() + m_min_idx);
                m_min_idx = 0;
                --m_num_files_to_merge;
                if (m_num_files_to_merge == 0) return;
            }
            compute_min();
        } else {  // update the looser tree
            m_min_idx = m_tree[0];
            assert(m_min_idx < m_iterators.size());
            auto& it = m_iterators[m_min_idx];
            it.begin += 1;
            uint64_t p = m_begin + m_min_idx;
            p -= (p >= m_tree.size()) * m_size;  // p is the index of the leaf
            if (it.begin == it.end) {
                m_tree[p] = uint32_t(-1);
                --m_num_files_to_merge;
            }
            while (p) {
                uint64_t is_right_child = (p & 1) == 0;
                uint32_t i = 0;
                uint32_t l = m_tree[p - is_right_child];
                uint32_t r = m_tree[p + 1 - is_right_child];
                if (l == uint32_t(-1)) {
                    i = r;
                } else if (r == uint32_t(-1)) {
                    i = l;
                } else {
                    i = *(m_iterators[l].begin) < *(m_iterators[r].begin) ? l : r;
                }
                uint64_t parent = (p - 1) / 2;
                m_tree[parent] = i;
                p = parent;
            }
            m_min_idx = m_tree[0];
        }
    };

    uint32_t build(uint32_t p) {
        if (p >= m_tree.size()) return uint32_t(-1);
        if (p >= m_size - 1) return m_tree[p];  // leaf
        uint32_t l = build(2 * p + 1);
        uint32_t r = build(2 * p + 2);
        uint32_t i = 0;
        if (l == uint32_t(-1)) {
            i = r;
        } else if (r == uint32_t(-1)) {
            i = l;
        } else {
            i = *(m_iterators[l].begin) < *(m_iterators[r].begin) ? l : r;
        }
        m_tree[p] = i;
        return i;
    }

    void compute_min() {
        m_min_idx = 0;
        auto min_val = *m_iterators.front().begin;
        for (uint64_t i = 1; i != m_iterators.size(); ++i) {
            assert(m_iterators[i].begin != m_iterators[i].end);
            auto val = *m_iterators[i].begin;
            if (val < min_val) {
                min_val = val;
                m_min_idx = i;
            }
        }
    }
};

}  // namespace sshash