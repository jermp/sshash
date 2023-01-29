#pragma once

#include <vector>
#include <cassert>
#include <fstream>
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
        , kmer_id_in_contig(constants::invalid_uint32)
        , kmer_orientation(constants::forward_orientation)
        , contig_id(constants::invalid_uint32)
        , contig_size(constants::invalid_uint32) {}
    uint64_t kmer_id;            // "absolute" kmer-id
    uint32_t kmer_id_in_contig;  // "relative" kmer-id: 0 <= kmer_id_in_contig < contig_size
    uint32_t kmer_orientation;
    uint32_t contig_id;
    uint32_t contig_size;
};

struct neighbourhood {
    /* forward */
    lookup_result forward_A;
    lookup_result forward_C;
    lookup_result forward_G;
    lookup_result forward_T;
    /* backward */
    lookup_result backward_A;
    lookup_result backward_C;
    lookup_result backward_G;
    lookup_result backward_T;
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

        , l(constants::min_l)
        , c(constants::c)

        , canonical_parsing(false)
        , weighted(false)
        , verbose(true)

        , tmp_dirname(constants::default_tmp_dirname) {}

    uint64_t k;  // kmer size
    uint64_t m;  // minimizer size
    uint64_t seed;

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

/*
char decimal  binary
 A     65     01000-00-1 -> 00
 C     67     01000-01-1 -> 01
 G     71     01000-11-1 -> 11
 T     84     01010-10-0 -> 10
*/
static kmer_t char_to_uint(char c) { return (c >> 1) & 3; }

static char uint64_to_char(uint64_t x) {
    assert(x <= 3);
    static char nucleotides[4] = {'A', 'C', 'T', 'G'};
    return nucleotides[x];
}

/*
    Traditional mapping.
*/
// uint64_t char_to_uint(char c) {
//     switch (c) {
//         case 'A':
//             return 0;
//         case 'C':
//             return 1;
//         case 'G':
//             return 2;
//         case 'T':
//             return 3;
//     }
//     assert(false);
//     return -1;
// }
// char uint64_to_char(uint64_t x) {
//     switch (x) {
//         case 0:
//             return 'A';
//         case 1:
//             return 'C';
//         case 2:
//             return 'G';
//         case 3:
//             return 'T';
//     }
//     assert(false);
//     return 0;
// }

/****************************************************************************
    The following two functions preserves the lexicographic order of k-mers,
    that is: if g and t are two k-mers and g < t lexicographically,
    then also id(g) < id(t).
*/
[[maybe_unused]] static kmer_t string_to_uint_kmer(char const* str, uint64_t k) {
    assert(k <= constants::max_k);
    kmer_t x = 0;
    for (uint64_t i = 0; i != k; ++i) {
        x <<= 2;
        x += char_to_uint(str[i]);
    }
    return x;
}
[[maybe_unused]] static void uint_kmer_to_string(kmer_t x, char* str, uint64_t k) {
    assert(k <= constants::max_k);
    for (int i = k - 1; i >= 0; --i) {
        str[i] = uint64_to_char(x & 3);
        x >>= 2;
    }
}
/****************************************************************************/

[[maybe_unused]] static std::string uint_kmer_to_string(kmer_t x, uint64_t k) {
    assert(k <= constants::max_k);
    std::string str;
    str.resize(k);
    uint_kmer_to_string(x, str.data(), k);
    return str;
}

[[maybe_unused]] static kmer_t string_to_uint_kmer_no_reverse(char const* str, uint64_t k) {
    assert(k <= constants::max_k);
    kmer_t x = 0;
    for (uint64_t i = 0; i != k; ++i) x += char_to_uint(str[i]) << (2 * i);
    return x;
}

static void uint_kmer_to_string_no_reverse(kmer_t x, char* str, uint64_t k) {
    assert(k <= constants::max_k);
    for (uint64_t i = 0; i != k; ++i) {
        str[i] = uint64_to_char(x & 3);
        x >>= 2;
    }
}

[[maybe_unused]] static std::string uint_kmer_to_string_no_reverse(kmer_t x, uint64_t k) {
    assert(k <= constants::max_k);
    std::string str;
    str.resize(k);
    uint_kmer_to_string_no_reverse(x, str.data(), k);
    return str;
}

/*
    This works with the map:
    A -> 00; C -> 01; G -> 11; T -> 10.

    Example.
    reverse_complement("ACTCACG") = CGTGAGT, in binary:
    reverse_complement("00.01.10.01.00.01.11") = 01.11.10.11.00.11.10.
*/
template <bool align>
[[maybe_unused]] static uint64_t crc(uint64_t x, uint64_t k) {
    assert(k <= 32);

    /* Complement, swap byte order */
    uint64_t res = __builtin_bswap64(x ^ 0xaaaaaaaaaaaaaaaa);

    /* Swap nuc order in bytes */
    const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
    const uint64_t c2 = 0x3333333333333333;
    res = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4);  // swap 2-nuc order in bytes
    res = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2);  // swap nuc order in 2-nuc

    /* Realign to the right */
    if constexpr (align) res >>= 64 - 2 * k;

    return res;
}

[[maybe_unused]] static kmer_t compute_reverse_complement(kmer_t x, uint64_t k) {
    assert(k <= constants::max_k);
    if constexpr (constants::uint_kmer_bits == 64) {
        return crc<true>(x, k);
    } else {
        assert(constants::uint_kmer_bits == 128);
        uint64_t low = static_cast<uint64_t>(x);
        uint64_t high = static_cast<uint64_t>(x >> 64);
        uint64_t k_low = 32;
        uint64_t k_high = k - 32;
        uint64_t shift = 128 - 2 * k;
        if (k < 32) {
            k_low = k;
            k_high = 0;
            shift = 64;
        }
        kmer_t low_rc = crc<true>(low, k_low);
        kmer_t high_rc = crc<false>(high, k_high);
        kmer_t res = (low_rc << 64) + high_rc;

        /* Realign to the right */
        res >>= shift;

        return res;
    }
}

// forward character map. A -> A, C -> C, G -> G, T -> T. rest maps to zero.
static const char canonicalize_basepair_forward_map[256] = {
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 65, 0, 67, 0, 0, 0, 71, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 84, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

// reverse character map. A -> T, C -> G, G -> C, T -> A. rest maps to zero.
static const char canonicalize_basepair_reverse_map[256] = {
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 84, 0, 71, 0, 0, 0, 67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 65, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

[[maybe_unused]] static void compute_reverse_complement(char const* input, char* output,
                                                        uint64_t size) {
    for (uint64_t i = 0; i != size; ++i) {
        int c = input[i];
        output[size - i - 1] = canonicalize_basepair_reverse_map[c];
    }
}

static inline bool is_valid(int c) { return canonicalize_basepair_forward_map[c]; }

[[maybe_unused]] static bool is_valid(char const* str, uint64_t size) {
    for (uint64_t i = 0; i != size; ++i) {
        int c = str[i];
        if (canonicalize_basepair_forward_map[c] == 0) return false;
        // if (c != 'A' and c != 'C' and c != 'G' and c != 'T') return false;
    }
    return true;
}

template <typename Hasher = murmurhash2_64>
uint64_t compute_minimizer(kmer_t kmer, uint64_t k, uint64_t m, uint64_t seed) {
    assert(m <= constants::max_m);
    assert(m <= k);
    uint64_t min_hash = uint64_t(-1);
    uint64_t minimizer = uint64_t(-1);
    kmer_t mask = (kmer_t(1) << (2 * m)) - 1;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        uint64_t mmer = static_cast<uint64_t>(kmer & mask);
        uint64_t hash = Hasher::hash(mmer, seed);
        if (hash < min_hash) {
            min_hash = hash;
            minimizer = mmer;
        }
        kmer >>= 2;
    }
    return minimizer;
}

/* used in dump.cpp */
template <typename Hasher = murmurhash2_64>
std::pair<uint64_t, uint64_t> compute_minimizer_pos(kmer_t kmer, uint64_t k, uint64_t m,
                                                    uint64_t seed) {
    assert(m <= constants::max_m);
    assert(m <= k);
    uint64_t min_hash = uint64_t(-1);
    uint64_t minimizer = uint64_t(-1);
    kmer_t mask = (kmer_t(1) << (2 * m)) - 1;
    uint64_t pos = 0;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        uint64_t mmer = static_cast<uint64_t>(kmer & mask);
        uint64_t hash = Hasher::hash(mmer, seed);
        if (hash < min_hash) {
            min_hash = hash;
            minimizer = mmer;
            pos = i;
        }
        kmer >>= 2;
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