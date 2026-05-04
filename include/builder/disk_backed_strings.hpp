#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "external/pthash/external/bits/include/bit_vector.hpp"

namespace sshash {

/*
    Disk-backed storage for the SSHash `strings` bit-vector.

    During step 1 (encode_strings) bits are appended via `append_bits`.
    Internally only the trailing words are kept in RAM (a small "write
    window"); completed words are flushed to a tmp file. RAM usage of the
    writer is bounded by the window size, independently of the total
    bit-vector size.

    After `freeze()`, callers create one or more `reader`s. Each reader owns
    an ifstream and a small in-RAM read window, and supports
    forward-monotonic `get_word64(bit_pos)` reads. The reader matches the
    interface that `kmer_iterator<Kmer, BitVector>` expects.

    For the standard dictionary save path, `load_into(bits::bit_vector&)`
    materializes the full bit-vector in RAM (peaks briefly at strings size).
*/
struct disk_backed_strings {
    static constexpr uint64_t default_writer_buffer_words = uint64_t(1) << 16;  // 512 KiB
    static constexpr uint64_t default_reader_window_words = uint64_t(1) << 16;  // 512 KiB

    disk_backed_strings()
        : m_num_bits(0)
        , m_writer_buffer_words(default_writer_buffer_words)
        , m_words_on_disk(0)
        , m_frozen(false) {}

    disk_backed_strings(disk_backed_strings const&) = delete;
    disk_backed_strings& operator=(disk_backed_strings const&) = delete;

    /* Open `filename` for writing; truncates any existing contents. */
    void open_for_writing(std::string const& filename,
                          uint64_t writer_buffer_words = default_writer_buffer_words) {
        m_filename = filename;
        m_writer_buffer_words = std::max<uint64_t>(2, writer_buffer_words);
        m_num_bits = 0;
        m_words_on_disk = 0;
        m_frozen = false;
        m_buf.clear();
        m_buf.reserve(m_writer_buffer_words);
        m_writer.open(m_filename, std::ofstream::binary | std::ofstream::trunc);
        if (!m_writer.is_open()) {
            throw std::runtime_error("cannot open strings tmp file '" + m_filename + "'");
        }
    }

    /* No-op: kept for source-compatibility with bits::bit_vector::builder. */
    void reserve(uint64_t /*num_bits*/) {}

    /* Append `len` bits (`len` <= 64) from `bits`. Same semantics as
       bits::bit_vector::builder::append_bits. */
    void append_bits(uint64_t bits, uint64_t len) {
        assert(len <= 64);
        assert(len == 64 || (bits >> len) == 0);
        if (!len) return;
        const uint64_t pos_in_word = m_num_bits & 63;
        m_num_bits += len;
        if (pos_in_word == 0) {
            m_buf.push_back(bits);
        } else {
            m_buf.back() |= bits << pos_in_word;
            if (len > 64 - pos_in_word) m_buf.push_back(bits >> (64 - pos_in_word));
        }
        if (m_buf.size() > m_writer_buffer_words) flush_completed_words();
    }

    /* Flush any remaining buffered words and close the writer. After this,
       the file is ready for `make_reader()` and `load_into()`. */
    void freeze() {
        if (m_frozen) return;
        if (!m_buf.empty()) {
            m_writer.write(reinterpret_cast<char const*>(m_buf.data()),
                           m_buf.size() * sizeof(uint64_t));
            m_words_on_disk += m_buf.size();
            m_buf.clear();
            m_buf.shrink_to_fit();
        }
        m_writer.close();
        m_frozen = true;
    }

    uint64_t num_bits() const { return m_num_bits; }
    std::string const& filename() const { return m_filename; }
    bool frozen() const { return m_frozen; }

    /*
        Forward-monotonic reader over the strings file.

        `get_word64(bit_pos)` returns the 64-bit word starting at bit
        position `bit_pos`. Successive calls must satisfy a forward-monotonic
        access pattern in word units (calling code may seek forward via
        `at()`-style calls in `kmer_iterator`, but never backward). Reads
        past the end-of-file are returned as zero (matches the sentinel
        zero-padding the SSHash builder writes at the tail of `strings`).
    */
    struct reader {
        reader() = default;
        reader(reader&& other) noexcept { move_from(std::move(other)); }
        reader& operator=(reader&& other) noexcept {
            if (this != &other) {
                close();
                move_from(std::move(other));
            }
            return *this;
        }
        reader(reader const&) = delete;
        reader& operator=(reader const&) = delete;
        ~reader() { close(); }

        void open(std::string const& filename, uint64_t num_bits,
                  uint64_t window_capacity_words = default_reader_window_words) {
            m_num_bits = num_bits;
            m_total_words = (num_bits + 63) / 64;
            m_window_capacity = std::max<uint64_t>(2, window_capacity_words);
            m_window.assign(m_window_capacity, 0);
            m_window_size = 0;
            m_window_start_word = 0;
            m_in.open(filename, std::ifstream::binary);
            if (!m_in.is_open()) {
                throw std::runtime_error("cannot open strings tmp file '" + filename + "'");
            }
            seek_window_to(0);
        }

        bool is_open() const { return m_in.is_open(); }

        void close() {
            if (m_in.is_open()) m_in.close();
            m_window.clear();
            m_window.shrink_to_fit();
            m_window_size = 0;
            m_window_start_word = 0;
        }

        uint64_t num_bits() const { return m_num_bits; }

        uint64_t get_word64(uint64_t bit_pos) const {
            const uint64_t block = bit_pos >> 6;
            const uint64_t shift = bit_pos & 63;
            ensure_window_covers(block);
            uint64_t a = (block >= m_window_start_word &&
                          block < m_window_start_word + m_window_size)
                             ? m_window[block - m_window_start_word]
                             : uint64_t(0);
            uint64_t word = a >> shift;
            if (shift) {
                const uint64_t next = block + 1;
                uint64_t b = (next >= m_window_start_word &&
                              next < m_window_start_word + m_window_size)
                                 ? m_window[next - m_window_start_word]
                                 : uint64_t(0);
                word |= b << (64 - shift);
            }
            return word;
        }

    private:
        mutable std::ifstream m_in;
        uint64_t m_num_bits = 0;
        uint64_t m_total_words = 0;
        mutable std::vector<uint64_t> m_window;
        uint64_t m_window_capacity = 0;
        mutable uint64_t m_window_size = 0;
        mutable uint64_t m_window_start_word = 0;

        void seek_window_to(uint64_t target_word) const {
            m_window_start_word = target_word;
            if (target_word >= m_total_words) {
                m_window_size = 0;
                return;
            }
            m_in.clear();  // clear any prior eof
            m_in.seekg(static_cast<std::streamoff>(target_word * sizeof(uint64_t)),
                       std::ios::beg);
            const uint64_t to_read = std::min(m_window_capacity, m_total_words - target_word);
            m_in.read(reinterpret_cast<char*>(m_window.data()),
                      static_cast<std::streamsize>(to_read * sizeof(uint64_t)));
            const std::streamsize nread = m_in.gcount();
            m_window_size = static_cast<uint64_t>(nread) / sizeof(uint64_t);
        }

        void ensure_window_covers(uint64_t block) const {
            // We may need both `block` and `block + 1` (for cross-word shifts).
            // The window covers [m_window_start_word, m_window_start_word + m_window_size).
            const uint64_t need_end = block + 2;  // exclusive
            if (block >= m_window_start_word && need_end <= m_window_start_word + m_window_size) {
                return;
            }
            // Slide forward (backward seeks are not supported).
            seek_window_to(block);
        }

        void move_from(reader&& other) {
            m_in = std::move(other.m_in);
            m_num_bits = other.m_num_bits;
            m_total_words = other.m_total_words;
            m_window = std::move(other.m_window);
            m_window_capacity = other.m_window_capacity;
            m_window_size = other.m_window_size;
            m_window_start_word = other.m_window_start_word;
            other.m_num_bits = 0;
            other.m_total_words = 0;
            other.m_window_capacity = 0;
            other.m_window_size = 0;
            other.m_window_start_word = 0;
        }
    };

    /* Create a new reader over the frozen file. */
    reader make_reader(uint64_t window_capacity_words = default_reader_window_words) const {
        if (!m_frozen) {
            throw std::runtime_error("disk_backed_strings: must freeze() before make_reader()");
        }
        reader r;
        r.open(m_filename, m_num_bits, window_capacity_words);
        return r;
    }

    /*
        Stream the strings to `os` in the same byte format that
        `essentials::generic_saver::visit(bits::bit_vector const&)` would
        produce — i.e.,
            uint64_t m_num_bits;
            size_t   n;            // number of 64-bit words
            uint64_t m_data[n];
        — without ever materializing the full bit-vector in RAM. The bytes
        are read from the tmp file in fixed-size chunks.

        This relies on `bits::bit_vector::visit_impl` writing exactly two
        fields (`m_num_bits` and the `m_data` owning_span<uint64_t>) and on
        `generic_saver::visit_seq` writing `size_t n` followed by the raw
        `n * sizeof(uint64_t)` bytes. If `bits::bit_vector` ever changes its
        on-disk representation, this method must be updated to match.
    */
    void save_to(std::ostream& os) const {
        if (!m_frozen) {
            throw std::runtime_error("disk_backed_strings: must freeze() before save_to()");
        }
        const uint64_t num_bits = m_num_bits;
        os.write(reinterpret_cast<char const*>(&num_bits), sizeof(uint64_t));
        const uint64_t total_words = (num_bits + 63) / 64;
        const std::size_t n = static_cast<std::size_t>(total_words);
        os.write(reinterpret_cast<char const*>(&n), sizeof(std::size_t));
        if (total_words == 0) return;
        std::ifstream in(m_filename, std::ifstream::binary);
        if (!in.is_open()) {
            throw std::runtime_error("cannot open strings tmp file '" + m_filename + "'");
        }
        std::vector<char> buffer(uint64_t(64) << 10);  // 64 KiB
        uint64_t bytes_remaining = total_words * sizeof(uint64_t);
        while (bytes_remaining > 0) {
            const std::streamsize chunk = static_cast<std::streamsize>(
                std::min<uint64_t>(buffer.size(), bytes_remaining));
            in.read(buffer.data(), chunk);
            const std::streamsize got = in.gcount();
            if (got <= 0) {
                throw std::runtime_error("unexpected EOF in strings tmp file '" + m_filename + "'");
            }
            os.write(buffer.data(), got);
            bytes_remaining -= static_cast<uint64_t>(got);
        }
        in.close();
    }

    /*
        Materialize the full bit-vector in RAM. This briefly peaks at the
        bit-vector size and is used immediately before `essentials::save`.
    */
    void load_into(bits::bit_vector& bv) const {
        if (!m_frozen) {
            throw std::runtime_error("disk_backed_strings: must freeze() before load_into()");
        }
        bits::bit_vector::builder b(m_num_bits);
        auto& data_vec = b.data();
        const uint64_t total_words = (m_num_bits + 63) / 64;
        if (total_words > 0) {
            std::ifstream in(m_filename, std::ifstream::binary);
            if (!in.is_open()) {
                throw std::runtime_error("cannot open strings tmp file '" + m_filename + "'");
            }
            in.read(reinterpret_cast<char*>(data_vec.data()),
                    static_cast<std::streamsize>(total_words * sizeof(uint64_t)));
            in.close();
        }
        b.build(bv);
    }

    /* Delete the on-disk strings file. */
    void remove_file() {
        if (!m_filename.empty()) std::remove(m_filename.c_str());
    }

private:
    std::string m_filename;
    std::ofstream m_writer;
    uint64_t m_num_bits;
    std::vector<uint64_t> m_buf;
    uint64_t m_writer_buffer_words;
    uint64_t m_words_on_disk;
    bool m_frozen;

    void flush_completed_words() {
        if (m_buf.size() < 2) return;
        const uint64_t to_flush = m_buf.size() - 1;  // keep last (possibly partial) word
        m_writer.write(reinterpret_cast<char const*>(m_buf.data()),
                       static_cast<std::streamsize>(to_flush * sizeof(uint64_t)));
        m_words_on_disk += to_flush;
        m_buf[0] = m_buf.back();
        m_buf.resize(1);
    }
};

}  // namespace sshash
