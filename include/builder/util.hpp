#pragma once

#include <atomic>
#include <fstream>
#include <memory>
#include <vector>

#include "buffered_record_stream.hpp"
#include "file_merging_iterator.hpp"
#include "parallel_sort.hpp"

namespace sshash {

struct parse_runtime_error : public std::runtime_error {
    parse_runtime_error() : std::runtime_error("did you provide an input file with weights?") {}
};

[[maybe_unused]] static void expect(char got, char expected) {
    if (got != expected) {
        std::cout << "got '" << got << "' but expected '" << expected << "'" << std::endl;
        throw parse_runtime_error();
    }
}

typedef uint8_t num_kmers_in_super_kmer_uint_type;

#pragma pack(push, 2)
struct minimizer_tuple {
    minimizer_tuple() {}
    minimizer_tuple(minimizer_info mini_info, uint64_t num_kmers_in_super_kmer)
        : minimizer(mini_info.minimizer)
        , pos_in_seq(mini_info.pos_in_seq)
        , pos_in_kmer(mini_info.pos_in_kmer)
        , num_kmers_in_super_kmer(num_kmers_in_super_kmer) {}

    bool operator>(minimizer_tuple other) const {
        if (minimizer != other.minimizer) return minimizer > other.minimizer;
        return pos_in_seq > other.pos_in_seq;
    }
    bool operator<(minimizer_tuple other) const {
        if (minimizer != other.minimizer) return minimizer < other.minimizer;
        return pos_in_seq < other.pos_in_seq;
    }

    static minimizer_tuple max() {
        return {{uint64_t(-1), uint64_t(-1), uint64_t(-1)}, uint64_t(-1)};
    }

    uint64_t minimizer;
    uint64_t pos_in_seq;
    num_kmers_in_super_kmer_uint_type pos_in_kmer;
    num_kmers_in_super_kmer_uint_type num_kmers_in_super_kmer;
};
#pragma pack(pop)

inline std::ostream& operator<<(std::ostream& os, minimizer_tuple const& mt) {
    os << "minimizer = " << mt.minimizer << std::endl;
    os << "pos_in_seq = " << mt.pos_in_seq << std::endl;
    os << "pos_in_kmer = " << int(mt.pos_in_kmer) << std::endl;
    os << "num_kmers_in_super_kmer = " << int(mt.num_kmers_in_super_kmer) << std::endl;
    return os;
}

/*
    Streaming forward iterator over a sorted minimizers tmp file that
    yields each distinct `minimizer` value exactly once (i.e., one value
    per bucket), built on top of `buffered_record_stream<minimizer_tuple>`
    so RAM usage is constant.

    Copyable: pthash's `build_in_external_memory` takes the iterator by
    value, so the underlying buffered stream is held via shared_ptr.
    Copies share the stream state; pthash's local copy advances the
    shared stream, and the original at the call site is unused after the
    build returns.
*/
struct streaming_minimizers_iterator {
    using iterator_category = std::forward_iterator_tag;
    using value_type = uint64_t;
    using difference_type = std::ptrdiff_t;
    using reference = uint64_t const&;
    using pointer = uint64_t const*;

    streaming_minimizers_iterator() = default;

    void open(std::string const& filename) {
        m_stream = std::make_shared<buffered_record_stream<minimizer_tuple>>();
        m_stream->open(filename);
        m_current = m_stream->empty() ? uint64_t(-1) : m_stream->current().minimizer;
    }

    void close() {
        if (m_stream) m_stream->close();
        m_stream.reset();
    }

    uint64_t operator*() const { return m_current; }
    streaming_minimizers_iterator& operator++() {
        advance_to_next_minimizer();
        return *this;
    }

private:
    std::shared_ptr<buffered_record_stream<minimizer_tuple>> m_stream;
    uint64_t m_current = uint64_t(-1);

    void advance_to_next_minimizer() {
        const uint64_t prev = m_current;
        while (!m_stream->empty()) {
            m_stream->advance();
            if (m_stream->empty()) return;  // m_current holds last value
            if (m_stream->current().minimizer != prev) {
                m_current = m_stream->current().minimizer;
                return;
            }
        }
    }
};

/*
    Streaming reader over a minimizers tmp file. Reads minimizer_tuple
    records via std::ifstream (no mmap), and groups consecutive tuples by
    minimizer into "buckets" with bounded RAM (~ one bucket at a time plus
    one record of lookahead).

    The caller passes a vector to receive the bucket's tuples; for typical
    inputs this peaks at max_bucket_size * sizeof(minimizer_tuple).
*/
struct streaming_minimizer_bucket_reader {
    void open(std::string const& filename) { m_stream.open(filename); }

    void close() { m_stream.close(); }

    bool has_next_bucket() const { return !m_stream.empty(); }

    /* Read the next bucket into `bucket_out` (cleared first). All tuples in
       a bucket share the same minimizer. Returns the bucket's minimizer. */
    uint64_t next_bucket(std::vector<minimizer_tuple>& bucket_out) {
        bucket_out.clear();
        assert(has_next_bucket());
        const uint64_t mm = m_stream.current().minimizer;
        do {
            bucket_out.push_back(m_stream.current());
            m_stream.advance();
        } while (!m_stream.empty() && m_stream.current().minimizer == mm);
        return mm;
    }

private:
    buffered_record_stream<minimizer_tuple> m_stream;
};

struct minimizers_tuples {
    minimizers_tuples() {}
    minimizers_tuples(build_configuration const& build_config)
        : m_num_minimizers(0)
        , m_num_minimizer_positions(0)
        , m_num_super_kmers(0)
        , m_run_identifier(pthash::clock_type::now().time_since_epoch().count())
        , m_build_config(build_config)  //
    {
        init();
    }

    void init() { m_num_files_to_merge = 0; }

    void sort_and_flush(std::vector<minimizer_tuple>& buffer) {
        parallel_sort(buffer, m_build_config.num_threads,
                      [](minimizer_tuple const& x, minimizer_tuple const& y) {
                          return (x.minimizer < y.minimizer) or
                                 (x.minimizer == y.minimizer and x.pos_in_seq < y.pos_in_seq);
                      });
        uint64_t id = m_num_files_to_merge.fetch_add(1);
        auto tmp_output_filename = get_tmp_output_filename(id);
        if (m_build_config.verbose) {
            std::cout << "saving to file '" << tmp_output_filename << "'..." << std::endl;
        }
        std::ofstream out(tmp_output_filename.c_str(), std::ofstream::binary);
        if (!out.is_open()) throw std::runtime_error("cannot open file");
        out.write(reinterpret_cast<char const*>(buffer.data()),
                  buffer.size() * sizeof(minimizer_tuple));
        out.close();
        buffer.clear();
    }

    std::string get_minimizers_filename() const {
        assert(m_num_files_to_merge > 0);
        std::stringstream filename;
        filename << m_build_config.tmp_dirname << "/sshash.tmp.run_" << m_run_identifier
                 << ".minimizers.bin";
        return filename.str();
    }

    struct files_name_iterator {
        files_name_iterator(minimizers_tuples const* ptr) : m_id(0), m_ptr(ptr) {}

        std::string operator*() { return m_ptr->get_tmp_output_filename(m_id); }
        void operator++() { ++m_id; }

    private:
        uint64_t m_id;
        minimizers_tuples const* m_ptr;
    };

    files_name_iterator files_name_iterator_begin() { return files_name_iterator(this); }

    void merge() {
        if (m_num_files_to_merge == 0) return;

        if (m_num_files_to_merge == 1) {
            std::rename(get_tmp_output_filename(0).c_str(), get_minimizers_filename().c_str());
            if (m_num_minimizers != 0) return;

            assert(m_num_minimizers == 0);
            assert(m_num_minimizer_positions == 0);
            assert(m_num_super_kmers == 0);

            /* Single-pass count via streaming ifstream (no mmap). */
            streaming_minimizer_bucket_reader reader;
            reader.open(get_minimizers_filename());
            std::vector<minimizer_tuple> bucket_buf;
            while (reader.has_next_bucket()) {
                reader.next_bucket(bucket_buf);
                uint64_t bucket_size = 0;
                {
                    uint64_t prev = constants::invalid_uint64;
                    for (auto const& mt : bucket_buf) {
                        if (mt.pos_in_seq != prev) {
                            ++bucket_size;
                            prev = mt.pos_in_seq;
                        }
                    }
                }
                m_num_minimizers += 1;
                m_num_minimizer_positions += bucket_size;
                m_num_super_kmers += bucket_buf.size();
            }
            reader.close();
            return;
        }

        assert(m_num_files_to_merge > 1);
        file_merging_iterator<minimizer_tuple> fm_iterator(files_name_iterator_begin(),
                                                           m_num_files_to_merge);

        std::ofstream out(get_minimizers_filename().c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open file");

        if (m_build_config.verbose) {
            std::cout << "saving to file '" << get_minimizers_filename() << "'" << std::endl;
        }

        m_num_minimizers = 0;
        m_num_minimizer_positions = 0;
        m_num_super_kmers = 0;
        uint64_t prev_minimizer = constants::invalid_uint64;
        uint64_t prev_pos_in_seq = constants::invalid_uint64;
        while (fm_iterator.has_next()) {
            minimizer_tuple mt = *fm_iterator;
            if (mt.minimizer != prev_minimizer) {
                prev_minimizer = mt.minimizer;
                ++m_num_minimizers;
                ++m_num_minimizer_positions;
            } else {
                if (mt.pos_in_seq != prev_pos_in_seq) ++m_num_minimizer_positions;
            }
            out.write(reinterpret_cast<char const*>(&mt), sizeof(minimizer_tuple));
            prev_pos_in_seq = mt.pos_in_seq;
            ++m_num_super_kmers;
            if (m_build_config.verbose and m_num_super_kmers % 100'000'000 == 0) {
                std::cout << "processed " << m_num_super_kmers << " minimizer tuples" << std::endl;
            }
            fm_iterator.next();
        }

        out.close();
        fm_iterator.close();

        /* remove tmp files */
        for (uint64_t i = 0; i != m_num_files_to_merge; ++i) {
            auto tmp_output_filename = get_tmp_output_filename(i);
            std::remove(tmp_output_filename.c_str());
        }
    }

    uint64_t num_files_to_merge() const { return m_num_files_to_merge; }
    uint64_t num_minimizers() const { return m_num_minimizers; }
    uint64_t num_minimizer_positions() const { return m_num_minimizer_positions; }
    uint64_t num_super_kmers() const { return m_num_super_kmers; }

    void remove_tmp_file() { std::remove(get_minimizers_filename().c_str()); }

private:
    std::atomic<uint64_t> m_num_files_to_merge;
    uint64_t m_num_minimizers;
    uint64_t m_num_minimizer_positions;
    uint64_t m_num_super_kmers;
    uint64_t m_run_identifier;
    build_configuration m_build_config;

    std::string get_tmp_output_filename(uint64_t id) const {
        std::stringstream filename;
        filename << m_build_config.tmp_dirname << "/sshash.tmp.run_" << m_run_identifier
                 << ".minimizers." << id << ".bin";
        return filename.str();
    }
};

}  // namespace sshash
