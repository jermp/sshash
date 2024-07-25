#pragma once

#include <array>
#include <cassert>
#include <fstream>
#include <thread>
#include <cmath>  // for std::ceil on linux

#include "hash_util.hpp"

namespace sshash {

struct streaming_query_report {
    streaming_query_report()
        : num_kmers(0), num_positive_kmers(0), num_searches(0), num_extensions(0) {}
    uint64_t num_kmers;
    uint64_t num_positive_kmers;
    uint64_t num_searches;
    uint64_t num_extensions;
};

struct lookup_result {
    lookup_result()
        : kmer_id(constants::invalid_uint64)
        , kmer_id_in_contig(constants::invalid_uint64)
        , kmer_orientation(constants::forward_orientation)
        , contig_id(constants::invalid_uint64)
        , contig_size(constants::invalid_uint64) {}
    uint64_t kmer_id;            // "absolute" kmer-id
    uint64_t kmer_id_in_contig;  // "relative" kmer-id: 0 <= kmer_id_in_contig < contig_size
    uint64_t kmer_orientation;
    uint64_t contig_id;
    uint64_t contig_size;
};

template <class kmer_t>
struct neighbourhood {
    std::array<lookup_result, kmer_t::alphabet_size> forward;
    std::array<lookup_result, kmer_t::alphabet_size> backward;
};

[[maybe_unused]] static bool equal_lookup_result(lookup_result expected, lookup_result got) {
    bool good = true;
    if (expected.kmer_id != got.kmer_id) {
        std::cout << "expected kmer_id " << expected.kmer_id << " but got " << got.kmer_id
                  << std::endl;
        good = false;
    }
    if (expected.kmer_id_in_contig != got.kmer_id_in_contig) {
        std::cout << "expected kmer_id_in_contig " << expected.kmer_id_in_contig << " but got "
                  << got.kmer_id_in_contig << std::endl;
        good = false;
    }
    if (got.kmer_id != constants::invalid_uint64 and
        expected.kmer_orientation != got.kmer_orientation) {
        std::cout << "expected kmer_orientation " << expected.kmer_orientation << " but got "
                  << got.kmer_orientation << std::endl;
        good = false;
    }
    if (expected.contig_id != got.contig_id) {
        std::cout << "expected contig_id " << expected.contig_id << " but got " << got.contig_id
                  << std::endl;
        good = false;
    }
    if (expected.contig_size != got.contig_size) {
        std::cout << "expected contig_size " << expected.contig_size << " but got "
                  << got.contig_size << std::endl;
        good = false;
    }
    return good;
}

struct build_configuration {
    build_configuration()
        : k(31)
        , m(17)
        , seed(constants::seed)
        , num_threads(std::thread::hardware_concurrency())

        , l(constants::min_l)
        , c(constants::c)

        , canonical_parsing(false)
        , weighted(false)
        , verbose(true)

        , tmp_dirname(constants::default_tmp_dirname) {}

    uint64_t k;  // kmer size
    uint64_t m;  // minimizer size
    uint64_t seed;
    uint64_t num_threads;

    uint64_t l;  // drive dictionary trade-off
    double c;    // drive PTHash trade-off

    bool canonical_parsing;
    bool weighted;
    bool verbose;

    std::string tmp_dirname;

    void print() const {
        std::cout << "k = " << k << ", m = " << m << ", seed = " << seed << ", l = " << l
                  << ", c = " << c
                  << ", canonical_parsing = " << (canonical_parsing ? "true" : "false")
                  << ", weighted = " << (weighted ? "true" : "false") << std::endl;
    }
};

namespace util {

static uint64_t get_seed_for_hash_function(build_configuration const& build_config) {
    static const uint64_t my_favourite_seed = 1234567890;
    return build_config.seed != my_favourite_seed ? my_favourite_seed : ~my_favourite_seed;
}

/* return the position of the most significant bit */
static inline uint32_t msb(uint32_t x) {
    assert(x > 0);
    return 31 - __builtin_clz(x);
}

static inline uint32_t ceil_log2_uint32(uint32_t x) { return (x > 1) ? msb(x - 1) + 1 : 0; }

[[maybe_unused]] static bool ends_with(std::string const& str, std::string const& pattern) {
    if (pattern.size() > str.size()) return false;
    return std::equal(pattern.begin(), pattern.end(), str.end() - pattern.size());
}

template <class kmer_t>
[[maybe_unused]] static kmer_t string_to_uint_kmer(char const* str, uint64_t k) {
    assert(k <= kmer_t::max_k);
    kmer_t x = 0;
    for (int i = k - 1; i >= 0; i--) { x.append_char(kmer_t::char_to_uint(str[i])); }
    return x;
}

template <class kmer_t>
static void uint_kmer_to_string(kmer_t x, char* str, uint64_t k) {
    assert(k <= kmer_t::max_k);
    for (uint64_t i = 0; i != k; ++i) { str[i] = kmer_t::uint64_to_char(x.pop_char()); }
}

template <class kmer_t>
[[maybe_unused]] static std::string uint_kmer_to_string(kmer_t x, uint64_t k) {
    assert(k <= kmer_t::max_k);
    std::string str;
    str.resize(k);
    uint_kmer_to_string(x, str.data(), k);
    return str;
}

template <class kmer_t>
[[maybe_unused]] static bool is_valid(char const* str, uint64_t size) {
    for (uint64_t i = 0; i != size; ++i) {
        if (!kmer_t::is_valid(str[i])) { return false; }
    }
    return true;
}

/*
    This implements the random minimizer.
*/
template <class kmer_t, typename Hasher = murmurhash2_64>
uint64_t compute_minimizer(kmer_t kmer, const uint64_t k, const uint64_t m, const uint64_t seed) {
    assert(m <= kmer_t::max_m);
    assert(m <= k);
    uint64_t min_hash = uint64_t(-1);
    kmer_t minimizer = kmer_t(-1);
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        kmer_t mmer = kmer;
        mmer.take_chars(m);
        uint64_t hash = Hasher::hash(uint64_t(mmer), seed);
        if (hash < min_hash) {
            min_hash = hash;
            minimizer = mmer;
        }
        kmer.drop_char();
    }
    return uint64_t(minimizer);
}

/* used in dump.cpp */
template <class kmer_t, typename Hasher = murmurhash2_64>
std::pair<kmer_t, uint64_t> compute_minimizer_pos(kmer_t kmer, uint64_t k, uint64_t m,
                                                  uint64_t seed) {
    assert(m <= kmer_t::max_m);
    assert(m <= k);
    uint64_t min_hash = uint64_t(-1);
    kmer_t minimizer = kmer_t(-1);
    uint64_t pos = 0;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        kmer_t mmer = kmer;
        mmer.take_chars(m);
        uint64_t hash = Hasher::hash(uint64_t(mmer), seed);
        if (hash < min_hash) {
            min_hash = hash;
            minimizer = mmer;
            pos = i;
        }
        kmer.drop_char();
    }
    return {minimizer, pos};
}

}  // namespace util

// taken from tlx
static inline std::istream& appendline(std::istream& is, std::string& str, char delim = '\n') {
    size_t size = str.size();
    size_t capacity = str.capacity();
    std::streamsize rest = capacity - size;

    if (rest == 0) {
        // if rest is zero, already expand string
        capacity = std::max(static_cast<size_t>(8), capacity * 2);
        rest = capacity - size;
    }

    // give getline access to all of capacity
    str.resize(capacity);

    // get until delim or rest is filled
    is.getline(const_cast<char*>(str.data()) + size, rest, delim);

    // gcount includes the delimiter
    size_t new_size = size + is.gcount();

    // is failbit set?
    if (!is) {
        // if string ran out of space, expand, and retry
        if (is.gcount() + 1 == rest) {
            is.clear();
            str.resize(new_size);
            str.reserve(capacity * 2);
            return appendline(is, str, delim);
        }
        // else fall through and deliver error
    } else if (!is.eof()) {
        // subtract delimiter
        --new_size;
    }

    // resize string to fit its contents
    str.resize(new_size);
    return is;
}

struct buffered_lines_iterator {
    static const uint64_t BUFFER_SIZE = 1024;

    buffered_lines_iterator(std::istream& is, uint64_t buffer_size = BUFFER_SIZE)
        : m_is(is), m_buffer_size(buffer_size), m_read_chars(0) {}

    bool fill_buffer(std::string& buffer,
                     bool force = false /* force reading of m_buffer_size characters */
    ) {
        bool empty_line_was_read = false;
        uint64_t size = buffer.size();
        uint64_t target_size = size + m_buffer_size;
        if (force) target_size += m_buffer_size;

        buffer.resize(target_size);

        char* ptr = buffer.data() + size;
        while (size != target_size) {
            // read until '\n' or rest is filled
            uint64_t rest = target_size - size;
            m_is.getline(ptr, rest, '\n');
            uint64_t read_chars = m_is.gcount();
            m_read_chars += read_chars;

            if (!m_is) {
                if (read_chars + 1 == rest) {  // '\n' not found
                    m_is.clear();
                    size += read_chars;
                    break;
                }
            } else if (!eof()) {
                assert(read_chars > 0);
                --read_chars;  // discard the delimiter
            }

            if (read_chars == 0) {  // empty line was read
                empty_line_was_read = true;
                break;
            }

            size += read_chars;
            ptr += read_chars;
        }

        buffer.resize(size);
        return empty_line_was_read;
    }

    bool eof() const { return m_is.eof(); }

    uint64_t read_chars() const { return m_read_chars; }

private:
    std::istream& m_is;
    uint64_t m_buffer_size;
    uint64_t m_read_chars;
};

}  // namespace sshash
