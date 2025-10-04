#pragma once

#include <vector>
#include <atomic>

#if defined(__x86_64__)
#include <immintrin.h>
#include <x86intrin.h>
#endif

#include "file_merging_iterator.hpp"
#include "parallel_sort.hpp"

namespace sshash {

[[maybe_unused]] static void print_time(double time, uint64_t num_kmers,
                                        std::string const& message) {
    std::cout << "=== " << message << " " << time / 1000000 << " [sec] ("
              << (time * 1000) / num_kmers << " [ns/kmer])" << std::endl;
}

struct parse_runtime_error : public std::runtime_error {
    parse_runtime_error() : std::runtime_error("did you provide an input file with weights?") {}
};

[[maybe_unused]] static void expect(char got, char expected) {
    if (got != expected) {
        std::cout << "got '" << got << "' but expected '" << expected << "'" << std::endl;
        throw parse_runtime_error();
    }
}

#if defined(__x86_64__)
/*
    This function takes 32 bytes and packs the two bits
    in positions 1 and 2 (from right) of each byte into
    a single 64-bit word.

    This works with the map:
    A -> 00; C -> 01; G -> 11; T -> 10.
*/
inline uint64_t pack2bits_shift1(__m256i v) {
    // shift >> 1, then mask by 3 to isolate the relevant bits
    __m256i shifted = _mm256_srli_epi16(v, 1);
    __m256i values = _mm256_and_si256(shifted, _mm256_set1_epi8(3));

    // collect bit-0 plane
    __m256i bit0 = _mm256_slli_epi16(values, 7);
    uint32_t mask0 = _mm256_movemask_epi8(bit0);

    // collect bit-1 plane
    __m256i bit1 = _mm256_slli_epi16(values, 6);
    uint32_t mask1 = _mm256_movemask_epi8(bit1);

    // interleave into the 64-bit result
    uint64_t even = _pdep_u64(mask0, 0x5555555555555555ULL);  // 010101...
    uint64_t odd = _pdep_u64(mask1, 0xAAAAAAAAAAAAAAAAULL);   // 101010...
    return even | odd;
}
#endif

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

struct bucket_type {
    bucket_type(minimizer_tuple const* begin, minimizer_tuple const* end)
        : m_begin(begin)
        , m_end(end)
        , m_num_super_kmers(std::distance(begin, end))
        , m_num_minimizer_positions(0)  //
    {
        uint64_t prev_pos_in_seq = constants::invalid_uint64;
        while (begin != end) {
            uint64_t pos_in_seq = (*begin).pos_in_seq;
            if (pos_in_seq != prev_pos_in_seq) {
                ++m_num_minimizer_positions;
                prev_pos_in_seq = pos_in_seq;
            }
            ++begin;
        }
        assert(m_num_minimizer_positions <= m_num_super_kmers);
    }

    struct iterator {
        iterator(minimizer_tuple const* begin) : m_begin(begin) {}

        inline minimizer_tuple operator*() const { return *m_begin; }
        inline void operator++() { ++m_begin; }
        bool operator==(iterator const& other) const { return m_begin == other.m_begin; }
        bool operator!=(iterator const& other) const { return !(*this == other); }

    private:
        minimizer_tuple const* m_begin;
    };

    iterator begin() const { return iterator(m_begin); }
    iterator end() const { return iterator(m_end); }

    /*
        When a canonical index is built (option `--canonical`),
        a minimizer offset can correspond to more than one super-kmer.
        A super-kmer is uniquely identified by the couple
          (minimizer offset, position of minimizer in the first kmer of the super-kmer).
        These two components, together, give the
        starting position of a super-kmer in the sequence.

        So the method size() returns the number of minimizer
        positions which is <= the number of superkmers.
    */
    uint64_t num_super_kmers() const { return m_num_super_kmers; }
    uint64_t size() const { return m_num_minimizer_positions; }

    minimizer_tuple const* begin_ptr() const { return m_begin; }
    minimizer_tuple const* end_ptr() const { return m_end; }

private:
    minimizer_tuple const* m_begin;
    minimizer_tuple const* m_end;
    uint64_t m_num_super_kmers;
    uint64_t m_num_minimizer_positions;
};

/*
    Iterate over the "bucket" of a minimizer, i.e.,
    the sorted list of minimizer tuples
    (minimizer, pos_in_seq, pos_in_kmer, num_kmers_in_superkmer).
*/
struct minimizers_tuples_iterator {
    typedef minimizer_tuple value_type;
    using iterator_category = std::forward_iterator_tag;

    minimizers_tuples_iterator(minimizer_tuple const* begin, minimizer_tuple const* end)
        : m_bucket_begin(begin), m_bucket_end(begin), m_end(end) {
        m_bucket_end = next_begin();
    }

    inline uint64_t minimizer() const { return (*m_bucket_begin).minimizer; }
    inline uint64_t operator*() const { return minimizer(); }
    inline void next() {
        m_bucket_begin = m_bucket_end;
        m_bucket_end = next_begin();
    }
    inline void operator++() { next(); }
    bool has_next() const { return m_bucket_begin != m_end; }
    bucket_type bucket() const { return bucket_type(m_bucket_begin, m_bucket_end); }

private:
    minimizer_tuple const* m_bucket_begin;
    minimizer_tuple const* m_bucket_end;
    minimizer_tuple const* m_end;

    minimizer_tuple const* next_begin() {
        minimizer_tuple const* begin = m_bucket_begin;
        uint64_t prev_minimizer = (*begin).minimizer;
        while (begin != m_end) {
            ++begin;
            uint64_t curr_minimizer = (*begin).minimizer;
            if (curr_minimizer != prev_minimizer) break;
        }
        return begin;
    }
};

struct minimizers_tuples {
    minimizers_tuples() {}
    minimizers_tuples(build_configuration const& build_config)
        : m_num_minimizers(0)
        , m_num_minimizer_positions(0)
        , m_num_super_kmers(0)
        , m_run_identifier(pthash::clock_type::now().time_since_epoch().count())
        , m_num_threads(build_config.num_threads)
        , m_tmp_dirname(build_config.tmp_dirname)  //
    {
        init();
    }

    void init() { m_num_files_to_merge = 0; }

    void sort_and_flush(std::vector<minimizer_tuple>& buffer) {
        parallel_sort(buffer, m_num_threads,
                      [](minimizer_tuple const& x, minimizer_tuple const& y) {
                          return (x.minimizer < y.minimizer) or
                                 (x.minimizer == y.minimizer and x.pos_in_seq < y.pos_in_seq);
                      });
        uint64_t id = m_num_files_to_merge.fetch_add(1);
        auto tmp_output_filename = get_tmp_output_filename(id);
        std::cout << "saving to file '" << tmp_output_filename << "'..." << std::endl;
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
        filename << m_tmp_dirname << "/sshash.tmp.run_" << m_run_identifier << ".minimizers.bin";
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
            mm::file_source<minimizer_tuple> input(get_minimizers_filename(),
                                                   mm::advice::sequential);
            for (minimizers_tuples_iterator it(input.data(), input.data() + input.size());
                 it.has_next(); it.next())  //
            {
                auto bucket = it.bucket();
                m_num_minimizer_positions += bucket.size();
                m_num_super_kmers += bucket.num_super_kmers();
                ++m_num_minimizers;
            }
            input.close();
            return;
        }

        std::cout << " == files to merge = " << m_num_files_to_merge << std::endl;

        assert(m_num_files_to_merge > 1);

        file_merging_iterator<minimizer_tuple> fm_iterator(files_name_iterator_begin(),
                                                           m_num_files_to_merge);

        std::cout << "saving tuples to '" << get_minimizers_filename() << "'" << std::endl;
        std::ofstream out(get_minimizers_filename().c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open file");

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
            if (m_num_super_kmers % 50000000 == 0) {
                std::cout << "num_super_kmers = " << m_num_super_kmers << std::endl;
            }
            fm_iterator.next();
        }

        out.close();
        fm_iterator.close();

        std::cout << "num_minimizers = " << m_num_minimizers << std::endl;
        std::cout << "num_minimizer_positions = " << m_num_minimizer_positions << std::endl;
        std::cout << "num_super_kmers = " << m_num_super_kmers << std::endl;

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
    uint64_t m_num_threads;
    std::string m_tmp_dirname;

    std::string get_tmp_output_filename(uint64_t id) const {
        std::stringstream filename;
        filename << m_tmp_dirname << "/sshash.tmp.run_" << m_run_identifier << ".minimizers." << id
                 << ".bin";
        return filename.str();
    }
};

}  // namespace sshash