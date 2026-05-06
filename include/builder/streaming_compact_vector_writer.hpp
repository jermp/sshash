#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

namespace sshash {

/*
    Streams a `bits::compact_vector` to disk one entry at a time, accepting
    `set(index, value)` calls in monotonically non-decreasing index order
    (gaps are filled with zero, matching the default-zero semantics of an
    in-RAM compact_vector::builder).

    The on-disk byte layout matches `bits::compact_vector::visit_impl`:
        uint64_t   m_size
        uint64_t   m_width
        uint64_t   m_mask
        size_t     n              (= ceil(m_size * m_width / 64))
        uint64_t   m_data[n]      (little-endian bit-packed)

    Total RAM footprint: a 2-word rolling window plus the std::ofstream's
    own buffer. Independently of the number of entries.
*/
struct streaming_compact_vector_writer {
    streaming_compact_vector_writer() = default;
    streaming_compact_vector_writer(streaming_compact_vector_writer const&) = delete;
    streaming_compact_vector_writer& operator=(streaming_compact_vector_writer const&) = delete;

    void open(std::string const& filename, uint64_t num_entries, uint64_t width) {
        if (width == 0)
            throw std::runtime_error("streaming_compact_vector_writer: width must be > 0");
        if (width > 64)
            throw std::runtime_error("streaming_compact_vector_writer: width must be <= 64");
        m_filename = filename;
        m_num_entries = num_entries;
        m_width = width;
        /* Match `bits::compact_vector::builder`'s data layout, which
           allocates `words_for(size*width) + 1` words: the trailing
           padding word allows in-RAM `set_bits` to write across word
           boundaries without bounds checking and is part of the
           serialized `m_data` owning_span. */
        const uint64_t packed_words = (num_entries == 0) ? 0 : (num_entries * width + 63) / 64;
        m_total_words = packed_words + 1;
        m_words_written = 0;
        m_buf[0] = 0;
        m_buf[1] = 0;
        m_have_last_index = false;

        m_out.open(filename, std::ofstream::binary | std::ofstream::trunc);
        if (!m_out.is_open()) {
            throw std::runtime_error("cannot open compact_vector tmp file '" + filename + "'");
        }

        /* Header (matches bits::compact_vector::visit_impl). */
        write_pod(m_num_entries);
        write_pod(m_width);
        const uint64_t mask = (m_width == 64) ? uint64_t(-1) : ((uint64_t(1) << m_width) - 1);
        write_pod(mask);
        const std::size_t n = static_cast<std::size_t>(m_total_words);
        write_pod(n);
    }

    /* Write a value at position `index`. Successive calls must satisfy
       `index >= previous_index`; gaps are filled with zero. */
    void set(uint64_t index, uint64_t value) {
        if (m_have_last_index) { assert(index >= m_last_index); }
        m_have_last_index = true;
        m_last_index = index;

        const uint64_t bit_offset = index * m_width;
        const uint64_t word_index = bit_offset / 64;
        const uint64_t bit_in_word = bit_offset % 64;

        /* Slide the 2-word window forward to cover word_index. Words below
           are now finalized; emit them. */
        while (m_words_written < word_index) {
            write_word(m_buf[0]);
            m_buf[0] = m_buf[1];
            m_buf[1] = 0;
            ++m_words_written;
        }

        /* OR `value` (m_width low bits of it) into the window starting at
           bit_in_word of m_buf[0]; overflow goes into m_buf[1]. */
        const uint64_t fits_in_word_0 = 64 - bit_in_word;
        if (m_width <= fits_in_word_0) {
            if (m_width == 64) {
                /* bit_in_word must be 0 here */
                m_buf[0] = value;
            } else {
                m_buf[0] |= value << bit_in_word;
            }
        } else {
            m_buf[0] |= value << bit_in_word;
            m_buf[1] |= value >> fits_in_word_0;
        }
    }

    /* Flush remaining buffered words and close the file. */
    void finalize() {
        while (m_words_written < m_total_words) {
            write_word(m_buf[0]);
            m_buf[0] = m_buf[1];
            m_buf[1] = 0;
            ++m_words_written;
        }
        if (m_out.is_open()) m_out.close();
    }

    std::string const& filename() const { return m_filename; }
    uint64_t num_entries() const { return m_num_entries; }
    uint64_t width() const { return m_width; }

    void remove_file() {
        if (m_out.is_open()) m_out.close();
        if (!m_filename.empty()) std::remove(m_filename.c_str());
    }

private:
    std::string m_filename;
    std::ofstream m_out;
    uint64_t m_num_entries = 0;
    uint64_t m_width = 0;
    uint64_t m_total_words = 0;
    uint64_t m_words_written = 0;
    uint64_t m_buf[2] = {0, 0};
    uint64_t m_last_index = 0;
    bool m_have_last_index = false;

    template <typename T>
    void write_pod(T const& v) {
        m_out.write(reinterpret_cast<char const*>(&v), sizeof(T));
    }
    void write_word(uint64_t w) {
        m_out.write(reinterpret_cast<char const*>(&w), sizeof(uint64_t));
    }
};

}  // namespace sshash
