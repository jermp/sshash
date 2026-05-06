#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace sshash {

/*
    A small buffered, forward-only reader of fixed-size records over a
    binary file. Records are read in fixed-capacity chunks (~`buffer_records
    * sizeof(T)` bytes of RAM) so the per-instance footprint is bounded
    independently of the file size.

    Used as the underlying primitive by all of SSHash's builder readers
    over on-disk record files (minimizer tuples, kmer requests, kmer
    records, offset values, sorted-run records, etc.). The class is
    move-only; for callers that need a copyable forward iterator (e.g.
    pthash's `build_in_external_memory`, which takes an iterator by
    value), wrap an instance in a `std::shared_ptr`.

    Usage:
        buffered_record_stream<T> s;
        s.open(filename);
        while (!s.empty()) {
            consume(s.current());
            s.advance();
        }
        s.close();
*/
template <typename T>
struct buffered_record_stream {
    static constexpr uint64_t default_buffer_records = 4096;

    buffered_record_stream() = default;
    buffered_record_stream(buffered_record_stream const&) = delete;
    buffered_record_stream& operator=(buffered_record_stream const&) = delete;
    buffered_record_stream(buffered_record_stream&&) = default;
    buffered_record_stream& operator=(buffered_record_stream&&) = default;

    /* Open `filename` for forward reading; optionally seek to byte
       `start_byte` before priming the read window. */
    void open(std::string const& filename,
              uint64_t buffer_records = default_buffer_records,
              std::streamoff start_byte = 0) {
        m_buf.resize(std::max<uint64_t>(1, buffer_records));
        m_in.open(filename, std::ifstream::binary);
        if (!m_in.is_open()) {
            throw std::runtime_error("cannot open file '" + filename + "'");
        }
        if (start_byte != 0) m_in.seekg(start_byte, std::ios::beg);
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

    bool is_open() const { return m_in.is_open(); }

    /* True iff there are no more records in the stream. */
    bool empty() const { return m_pos >= m_size; }

    /* Reference to the current record. Valid until the next `advance()`. */
    T const& current() const {
        assert(!empty());
        return m_buf[m_pos];
    }

    /* Move to the next record; refills the buffer from disk on demand. */
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

}  // namespace sshash
