#pragma once

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "include/builder/buffered_record_stream.hpp"
#include "include/offsets.hpp"

namespace sshash {

/*
    A disk-backed drop-in for `Offsets::builder` that spills offset values
    to a tmp file as they're appended, keeping only a small in-RAM write
    buffer. RAM usage of the builder is bounded by the buffer size,
    independently of the number of strings.

    Interface mirrors `decoded_offsets::builder` / `encoded_offsets::builder`
    enough for the SSHash build path:

        reserve(n)                      no-op (kept for source compat)
        push_back(val)                  append to file (buffered)
        front() / back() / size()       O(1), tracked separately
        num_bytes()                     in-RAM footprint of the builder
        set_num_bits(nb)                stash bit-width metadata
        num_bits_per_offset()           dispatched per Offsets type
        encode(offset, begin, sid)      pure, dispatched per Offsets type
        build(target)                   stream from disk into the final
                                        compact / endpoints sequence
        remove_file()                   cleanup

    Random-access `operator[]` is *not* supported. Callers that need to
    walk a contiguous range of offsets must use a `reader`, which provides
    forward-sequential reads via a small in-RAM buffer.
*/
template <typename Offsets>
struct disk_backed_offsets_builder {
    static_assert(std::is_same_v<Offsets, decoded_offsets> ||
                      std::is_same_v<Offsets, encoded_offsets>,
                  "disk_backed_offsets_builder supports decoded_offsets and encoded_offsets");

    static constexpr uint64_t default_writer_buffer_records = uint64_t(1) << 12;  // 32 KiB
    static constexpr uint64_t default_reader_buffer_records = uint64_t(1) << 12;  // 32 KiB

    disk_backed_offsets_builder() = default;
    disk_backed_offsets_builder(disk_backed_offsets_builder const&) = delete;
    disk_backed_offsets_builder& operator=(disk_backed_offsets_builder const&) = delete;

    void open_for_writing(std::string const& filename,
                          uint64_t writer_buffer_records = default_writer_buffer_records) {
        m_filename = filename;
        m_writer_buffer_capacity = std::max<uint64_t>(1, writer_buffer_records);
        m_buf.clear();
        m_buf.reserve(m_writer_buffer_capacity);
        m_size = 0;
        m_front = 0;
        m_back = 0;
        m_have_front = false;
        m_frozen = false;
        m_writer.open(m_filename, std::ofstream::binary | std::ofstream::trunc);
        if (!m_writer.is_open()) {
            throw std::runtime_error("cannot open offsets tmp file '" + m_filename + "'");
        }
    }

    /* No-op: kept for source-compatibility with the in-RAM builder. */
    void reserve(uint64_t /*n*/) {}

    void push_back(uint64_t val) {
        if (!m_have_front) {
            m_front = val;
            m_have_front = true;
        }
        m_back = val;
        m_buf.push_back(val);
        ++m_size;
        if (m_buf.size() >= m_writer_buffer_capacity) flush_buffer();
    }

    /* Finish writing: flush the in-RAM buffer and close the writer. */
    void freeze() {
        if (m_frozen) return;
        flush_buffer();
        if (m_writer.is_open()) m_writer.close();
        m_frozen = true;
    }

    uint64_t size() const { return m_size; }
    uint64_t front() const { return m_front; }
    uint64_t back() const { return m_back; }
    std::string const& filename() const { return m_filename; }

    /* In-RAM footprint of the builder (excluding the on-disk file). */
    uint64_t num_bytes() const { return sizeof(m_nb) + m_buf.capacity() * sizeof(uint64_t); }

    void set_num_bits(num_bits nb) { m_nb = nb; }

    uint64_t num_bits_per_offset() const {
        if constexpr (std::is_same_v<Offsets, decoded_offsets>) {
            return m_nb.per_absolute_offset;
        } else {
            return m_nb.per_string_id + m_nb.per_relative_offset;
        }
    }

    /* Pure: matches `decoded_offsets::builder::encode` /
       `encoded_offsets::builder::encode`. Safe to call concurrently from
       multiple threads (depends only on m_nb and m_size, both of which
       are stable while compute_minimizer_tuples runs). */
    uint64_t encode(uint64_t offset, uint64_t begin, uint64_t string_id) const {
        if constexpr (std::is_same_v<Offsets, decoded_offsets>) {
            (void)begin;
            (void)string_id;
            return offset;
        } else {
            assert(string_id < m_size);
            assert(offset >= begin);
            assert((offset - begin) < (uint64_t(1) << m_nb.per_relative_offset));
            uint64_t relative_offset = offset - begin;
            return (string_id << m_nb.per_relative_offset) + relative_offset;
        }
    }

    /*
        Forward-sequential reader over the offsets file. Each thread in
        compute_minimizer_tuples should construct one for its assigned
        index range; per-thread RAM footprint is the buffer size only.
        Built on top of `buffered_record_stream<uint64_t>`.
    */
    struct reader {
        reader() = default;

        /* Open the file and seek so that the next `next()` call returns
           `*(values + start_index)`. */
        void open(std::string const& filename, uint64_t start_index,
                  uint64_t buffer_records = default_reader_buffer_records) {
            m_stream.open(filename, buffer_records,
                          static_cast<std::streamoff>(start_index * sizeof(uint64_t)));
        }

        void close() { m_stream.close(); }

        /* Return the next offset and advance. Caller must ensure they
           don't read past the end of the file. */
        uint64_t next() {
            if (m_stream.empty()) {
                throw std::runtime_error("disk_backed_offsets_builder: read past end of file");
            }
            const uint64_t v = m_stream.current();
            m_stream.advance();
            return v;
        }

    private:
        buffered_record_stream<uint64_t> m_stream;
    };

    /* Construct a reader positioned at `start_index`. Requires freeze(). */
    reader make_reader(uint64_t start_index,
                       uint64_t buffer_records = default_reader_buffer_records) const {
        if (!m_frozen) {
            throw std::runtime_error(
                "disk_backed_offsets_builder: must freeze() before make_reader()");
        }
        reader r;
        r.open(m_filename, start_index, buffer_records);
        return r;
    }

    /*
        A copyable forward iterator over the entire offsets file, suitable
        for the `Iterator`-template `encode` / `build` calls in
        `bits::endpoints_sequence` and `bits::compact_vector`. Wraps a
        shared_ptr<buffered_record_stream<uint64_t>> so the iterator is
        copyable; copies share the underlying stream state, which is what
        those APIs expect.
    */
    struct full_iterator {
        using iterator_category = std::forward_iterator_tag;
        using value_type = uint64_t;
        using difference_type = std::ptrdiff_t;
        using reference = uint64_t const&;
        using pointer = uint64_t const*;

        full_iterator() = default;

        void open(std::string const& filename,
                  uint64_t buffer_records = default_reader_buffer_records) {
            m_stream = std::make_shared<buffered_record_stream<uint64_t>>();
            m_stream->open(filename, buffer_records);
        }

        uint64_t operator*() const {
            assert(m_stream);
            return m_stream->current();
        }
        full_iterator& operator++() {
            assert(m_stream);
            m_stream->advance();
            return *this;
        }

    private:
        std::shared_ptr<buffered_record_stream<uint64_t>> m_stream;
    };

    /*
        Stream the offsets file into the target Offsets structure (mirrors
        the in-RAM builder's `build`). After return, the file is removed
        and `size()` resets to 0 to match the in-RAM builder, which clears
        its m_v in `build`.
    */
    void build(Offsets& target) {
        if (!m_frozen) freeze();
        if (m_size == 0) {
            remove_file();
            reset_state();
            return;
        }

        if constexpr (std::is_same_v<Offsets, decoded_offsets>) {
            full_iterator it;
            it.open(m_filename);
            target.m_seq.encode(it, m_size, m_back);
        } else {
            full_iterator it;
            it.open(m_filename);
            target.m_seq.build(it, m_size, m_nb.per_absolute_offset);
            target.m_num_bits_per_relative_offset = m_nb.per_relative_offset;
        }

        remove_file();
        reset_state();
    }

    /* Remove the on-disk tmp file (if any). */
    void remove_file() {
        if (m_writer.is_open()) m_writer.close();
        if (!m_filename.empty()) std::remove(m_filename.c_str());
    }

private:
    std::string m_filename;
    std::ofstream m_writer;
    std::vector<uint64_t> m_buf;
    uint64_t m_writer_buffer_capacity = default_writer_buffer_records;
    uint64_t m_size = 0;
    uint64_t m_front = 0;
    uint64_t m_back = 0;
    bool m_have_front = false;
    bool m_frozen = false;
    num_bits m_nb;

    void flush_buffer() {
        if (m_buf.empty()) return;
        m_writer.write(reinterpret_cast<char const*>(m_buf.data()),
                       static_cast<std::streamsize>(m_buf.size() * sizeof(uint64_t)));
        m_buf.clear();
    }

    void reset_state() {
        m_size = 0;
        m_buf.clear();
        m_buf.shrink_to_fit();
        m_have_front = false;
        m_front = 0;
        m_back = 0;
        m_frozen = false;
    }
};

}  // namespace sshash
