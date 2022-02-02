#pragma once

#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>  // for std::ceil on linux

#include "../external/pthash/include/pthash.hpp"

namespace sshash {

namespace constants {
constexpr uint64_t max_k = 32;
constexpr uint64_t invalid = uint64_t(-1);
constexpr uint64_t seed = 1;
constexpr uint64_t hashcode_bits = 64;
constexpr double c = 3.0;  // for PThash
constexpr uint64_t min_l = 6;
constexpr uint64_t max_l = 12;
}  // namespace constants

typedef pthash::murmurhash2_64 base_hasher_type;

typedef pthash::single_phf<base_hasher_type,               // base hasher
                           pthash::dictionary_dictionary,  // encoder type
                           true                            // minimal output
                           >
    pthash_mphf_type;

namespace util {

void check_hash_collision_probability(uint64_t size) {
    /*
        From: https://preshing.com/20110504/hash-collision-probabilities/
        Given a universe of size U (total number of possible hash values),
        which is U = 2^b for b-bit hash codes,
        the collision probability for n keys is (approximately):
            1 - e^{-n(n-1)/(2U)}.
        For example, for U=2^32 (32-bit hash codes), this probability
        gets to 50% already for n = 77,163 keys.
        We can approximate 1-e^{-X} with X when X is sufficiently small.
        Then our collision probability is
            n(n-1)/(2U) ~ n^2/(2U).
        So it can derived that ~1.97B keys and 64-bit hash codes,
        we have a probability of collision that is ~0.1 (10%), which may not be
        so small for certain applications.
        For n = 2^30, the probability of collision is ~0.031 (3.1%).
    */
    if (constants::hashcode_bits == 64 and size > (1ULL << 30)) {
        throw std::runtime_error(
            "Using 64-bit hash codes with more than 2^30 keys can be dangerous due to "
            "collisions: use 128-bit hash codes instead.");
    }
}

}  // namespace util

struct build_configuration {
    build_configuration()
        : k(31)
        , m(17)
        , seed(constants::seed)
        , max_num_kmers(-1)

        , l(constants::min_l)
        , c(constants::c)

        , canonical_parsing(false)
        , verbose(true) {}

    uint64_t k;  // kmer size
    uint64_t m;  // minimizer size
    uint64_t seed;
    uint64_t max_num_kmers;

    uint64_t l;  // drive dictionary trade-off
    double c;    // drive PTHash trade-off

    bool canonical_parsing;
    bool verbose;

    void print() const {
        std::cout << "k = " << k << ", m = " << m << ", seed = " << seed << ", l = " << l
                  << ", c = " << c
                  << ", canonical_parsing = " << (canonical_parsing ? "true" : "false")
                  << std::endl;
    }
};

struct buckets_statistics {
    static const uint64_t max_bucket_size = 4 * 1024;
    static const uint64_t max_string_size = 256;

    buckets_statistics(uint64_t num_buckets, uint64_t num_kmers, uint64_t num_strings = 0)
        : m_num_buckets(num_buckets)
        , m_num_kmers(num_kmers)
        // , m_num_strings(num_strings)
        , m_max_num_kmers_in_string(0)
        , m_max_num_strings_in_bucket(0) {
        (void)num_strings;
        m_bucket_sizes.resize(max_bucket_size + 1, 0);
        m_total_kmers.resize(max_bucket_size + 1, 0);
        m_string_sizes.resize(max_string_size + 1, 0);
    }

    void add_num_strings_in_bucket(uint64_t num_strings_in_bucket) {
        if (num_strings_in_bucket < max_bucket_size + 1) {
            m_bucket_sizes[num_strings_in_bucket] += 1;
        }
    }

    void add_num_kmers_in_string(uint64_t num_strings_in_bucket, uint64_t num_kmers_in_string) {
        if (num_strings_in_bucket < max_bucket_size + 1) {
            m_total_kmers[num_strings_in_bucket] += num_kmers_in_string;
        }
        if (num_kmers_in_string > m_max_num_kmers_in_string) {
            m_max_num_kmers_in_string = num_kmers_in_string;
        }
        if (num_strings_in_bucket > m_max_num_strings_in_bucket) {
            m_max_num_strings_in_bucket = num_strings_in_bucket;
        }
        if (num_kmers_in_string < max_string_size + 1) m_string_sizes[num_kmers_in_string] += 1;
    }

    uint64_t num_kmers() const { return m_num_kmers; }
    uint64_t num_buckets() const { return m_num_buckets; }
    uint64_t max_num_strings_in_bucket() const { return m_max_num_strings_in_bucket; }

    void print() const {
        // full statistics
        // std::cout << " === bucket statistics === \n";
        // for (uint64_t bucket_size = 1, prev_bucket_size = 0, prev_kmers_in_buckets = 0,
        //               kmers_in_buckets = 0;
        //      bucket_size != max_bucket_size + 1; ++bucket_size) {
        //     if (m_bucket_sizes[bucket_size] > 0) {
        //         std::cout << "buckets with " << bucket_size
        //                   << " strings=" << m_bucket_sizes[bucket_size] << "("
        //                   << (m_bucket_sizes[bucket_size] * 100.0) / m_num_buckets
        //                   << "%)|total_kmers=" << m_total_kmers[bucket_size] << "("
        //                   << (m_total_kmers[bucket_size] * 100.0) / m_num_kmers << "%)"
        //                   << "|avg_num_kmers_per_bucket="
        //                   << static_cast<double>(m_total_kmers[bucket_size]) /
        //                          m_bucket_sizes[bucket_size]
        //                   << "|avg_num_kmers_per_string="
        //                   << static_cast<double>(m_total_kmers[bucket_size]) /
        //                          (m_bucket_sizes[bucket_size] * bucket_size)
        //                   << std::endl;
        //         kmers_in_buckets += m_total_kmers[bucket_size];
        //     }
        //     if (bucket_size == 4 or bucket_size == 8 or bucket_size == 16 or bucket_size == 32 or
        //         bucket_size == 64 or bucket_size == 128 or bucket_size == 256 or
        //         bucket_size == 512 or bucket_size == 1024 or bucket_size == max_bucket_size) {
        //         assert(kmers_in_buckets >= prev_kmers_in_buckets);

        //         std::cout << " *** kmers in buckets of size > " << prev_bucket_size
        //                   << " and <= " << bucket_size << ": "
        //                   << kmers_in_buckets - prev_kmers_in_buckets << "("
        //                   << (100.0 * (kmers_in_buckets - prev_kmers_in_buckets)) / m_num_kmers
        //                   << "%)" << std::endl;
        //         std::cout << " *** kmers in buckets of size <= " << bucket_size << ": "
        //                   << kmers_in_buckets << "(" << (100.0 * kmers_in_buckets) / m_num_kmers
        //                   << "%)" << std::endl;

        //         prev_bucket_size = bucket_size;
        //         prev_kmers_in_buckets = kmers_in_buckets;
        //     }
        // }

        std::cout << " === bucket statistics (less) === \n";
        for (uint64_t bucket_size = 1; bucket_size != 16 + 1; ++bucket_size) {
            if (m_bucket_sizes[bucket_size] > 0) {
                std::cout << "buckets with " << bucket_size
                          << " strings = " << (m_bucket_sizes[bucket_size] * 100.0) / m_num_buckets
                          << "%" << std::endl;
            }
        }
        std::cout << "max_num_strings_in_bucket " << m_max_num_strings_in_bucket << std::endl;

        // std::cout << " === string statistics === \n";
        // uint64_t total_strings = 0;
        // uint64_t total_kmers = 0;
        // for (uint64_t string_size = 1; string_size != max_string_size + 1; ++string_size) {
        //     if (m_string_sizes[string_size] > 0) {
        //         std::cout << "strings with " << string_size
        //                   << " kmer=" << m_string_sizes[string_size] << "("
        //                   << (m_string_sizes[string_size] * 100.0) / m_num_strings
        //                   << "%)|total_kmers=" << (string_size * m_string_sizes[string_size]) <<
        //                   "("
        //                   << (string_size * m_string_sizes[string_size] * 100.0) / m_num_kmers
        //                   << "%)" << std::endl;
        //         total_strings += m_string_sizes[string_size];
        //         total_kmers += string_size * m_string_sizes[string_size];
        //     }
        // }
        // std::cout << "total_strings " << total_strings << "/" << m_num_strings << " ("
        //           << (total_strings * 100.0) / m_num_strings << "%)" << std::endl;
        // std::cout << "total_kmers " << total_kmers << "/" << m_num_kmers << " ("
        //           << (total_kmers * 100.0) / m_num_kmers << "%)" << std::endl;
        // std::cout << "max_num_kmers_in_string " << m_max_num_kmers_in_string << std::endl;
    }

private:
    uint64_t m_num_buckets;
    uint64_t m_num_kmers;
    // uint64_t m_num_strings;
    uint64_t m_max_num_kmers_in_string;
    uint64_t m_max_num_strings_in_bucket;
    std::vector<uint64_t> m_bucket_sizes;
    std::vector<uint64_t> m_total_kmers;
    std::vector<uint64_t> m_string_sizes;
};

namespace util {

/* return the position of the most significant bit */
inline uint32_t msb(uint32_t x) {
    assert(x > 0);
    return 31 - __builtin_clz(x);
}

inline uint32_t ceil_log2_uint32(uint32_t x) { return (x > 1) ? msb(x - 1) + 1 : 0; }

bool ends_with(std::string const& str, std::string const& pattern) {
    if (pattern.size() > str.size()) return false;
    return std::equal(pattern.begin(), pattern.end(), str.end() - pattern.size());
}

// for a sorted list of size n whose universe is u
uint64_t elias_fano_bitsize(uint64_t n, uint64_t u) {
    return n * ((u > n ? (std::ceil(std::log2(static_cast<double>(u) / n))) : 0) + 2);
}

/*
char decimal  binary
 A     65     01000-00-1 -> 00
 C     67     01000-01-1 -> 01
 G     71     01000-11-1 -> 11
 T     84     01010-10-0 -> 10
*/
uint64_t char_to_uint64(char c) { return (c >> 1) & 3; }

char uint64_to_char(uint64_t x) {
    assert(x <= 3);
    static char nucleotides[4] = {'A', 'C', 'T', 'G'};
    return nucleotides[x];
}

/*
    Traditional mapping.
*/
// uint64_t char_to_uint64(char c) {
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
uint64_t string_to_uint64(char const* str, uint64_t k) {
    assert(k <= 32);
    uint64_t x = 0;
    for (uint64_t i = 0; i != k; ++i) {
        x <<= 2;
        x += char_to_uint64(str[i]);
    }
    return x;
}
void uint64_to_string(uint64_t x, char* str, uint64_t k) {
    assert(k <= 32);
    for (int i = k - 1; i >= 0; --i) {
        str[i] = uint64_to_char(x & 3);
        x >>= 2;
    }
}
/****************************************************************************/

std::string uint64_to_string(uint64_t x, uint64_t k) {
    assert(k <= 32);
    std::string str;
    str.resize(k);
    uint64_to_string(x, str.data(), k);
    return str;
}

uint64_t string_to_uint64_no_reverse(char const* str, uint64_t k) {
    assert(k <= 32);
    uint64_t x = 0;
    for (uint64_t i = 0; i != k; ++i) x += char_to_uint64(str[i]) << (2 * i);
    return x;
}

void uint64_to_string_no_reverse(uint64_t x, char* str, uint64_t k) {
    assert(k <= 32);
    for (uint64_t i = 0; i != k; ++i) {
        str[i] = uint64_to_char(x & 3);
        x >>= 2;
    }
}

std::string uint64_to_string_no_reverse(uint64_t x, uint64_t k) {
    assert(k <= 32);
    std::string str;
    str.resize(k);
    uint64_to_string_no_reverse(x, str.data(), k);
    return str;
}

/*
    taken from Blight:
    it works with the map
    A -> 00; C -> 01; G -> 11; T -> 10
    Example:
    reverse_complement("ACTCACG") = CGTGAGT
    in binary:
    reverse_complement("00011001000111") = 01111011001110
*/
uint64_t compute_reverse_complement(uint64_t x, uint64_t size) {
    assert(size <= 32);
    // Complement, swap byte order
    uint64_t res = __builtin_bswap64(x ^ 0xaaaaaaaaaaaaaaaa);
    // Swap nuc order in bytes
    const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
    const uint64_t c2 = 0x3333333333333333;
    res = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4);  // swap 2-nuc order in bytes
    res = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2);  // swap nuc order in 2-nuc
    // Realign to the right
    res >>= 64 - 2 * size;
    return res;
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

void compute_reverse_complement(char const* input, char* output, uint64_t size) {
    for (uint64_t i = 0; i != size; ++i) {
        int c = input[i];
        output[size - i - 1] = canonicalize_basepair_reverse_map[c];
    }
}

inline bool is_valid(int c) { return canonicalize_basepair_forward_map[c]; }

bool is_valid(char const* str, uint64_t size) {
    for (uint64_t i = 0; i != size; ++i) {
        int c = str[i];
        if (canonicalize_basepair_forward_map[c] == 0) return false;
        // if (c != 'A' and c != 'C' and c != 'G' and c != 'T') return false;
    }
    return true;
}

/*
    This code is an adaptation from
    https://github.com/aappleby/smhasher/blob/master/src/MurmurHash2.cpp
        by Austin Appleby
*/
uint64_t MurmurHash2_64(void const* key, size_t len, uint64_t seed) {
    const uint64_t m = 0xc6a4a7935bd1e995ULL;
    const int r = 47;

    uint64_t h = seed ^ (len * m);

#if defined(__arm) || defined(__arm__)
    const size_t ksize = sizeof(uint64_t);
    const unsigned char* data = (const unsigned char*)key;
    const unsigned char* end = data + (std::size_t)(len / 8) * ksize;
#else
    const uint64_t* data = (const uint64_t*)key;
    const uint64_t* end = data + (len / 8);
#endif

    while (data != end) {
#if defined(__arm) || defined(__arm__)
        uint64_t k;
        memcpy(&k, data, ksize);
        data += ksize;
#else
        uint64_t k = *data++;
#endif

        k *= m;
        k ^= k >> r;
        k *= m;

        h ^= k;
        h *= m;
    }

    const unsigned char* data2 = (const unsigned char*)data;

    switch (len & 7) {
        // fall through
        case 7:
            h ^= uint64_t(data2[6]) << 48;
        // fall through
        case 6:
            h ^= uint64_t(data2[5]) << 40;
        // fall through
        case 5:
            h ^= uint64_t(data2[4]) << 32;
        // fall through
        case 4:
            h ^= uint64_t(data2[3]) << 24;
        // fall through
        case 3:
            h ^= uint64_t(data2[2]) << 16;
        // fall through
        case 2:
            h ^= uint64_t(data2[1]) << 8;
        // fall through
        case 1:
            h ^= uint64_t(data2[0]);
            h *= m;
    };

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}

struct byte_range {
    char const* begin;
    char const* end;
};

struct murmurhash2_64 {
    // generic range of bytes
    static inline uint64_t hash(byte_range range, uint64_t seed) {
        return MurmurHash2_64(range.begin, range.end - range.begin, seed);
    }

    // // specialization for std::string
    // static inline uint64_t hash(std::string const& val, uint64_t seed) {
    //     return MurmurHash2_64(val.data(), val.size(), seed);
    // }

    // specialization for uint64_t
    static inline uint64_t hash(uint64_t val, uint64_t seed) {
        return MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed);
    }
};

// return a pointer to the start of the minimizer within the term
template <typename Hasher = murmurhash2_64>
char const* compute_minimizer(char const* kmer, uint64_t k, uint64_t m, uint64_t seed) {
    assert(m <= k);
    uint64_t min_hash = uint64_t(-1);
    uint64_t min_index = 0;
    for (uint64_t i = 0; i + m <= k; ++i) {
        byte_range sub_kmer{kmer + i, kmer + i + m};
        uint64_t hash = Hasher::hash(sub_kmer, seed);
        if (hash < min_hash) {
            min_hash = hash;
            min_index = i;
        }
    }
    return kmer + min_index;
}

template <typename Hasher = murmurhash2_64>
uint64_t compute_minimizer(uint64_t kmer, uint64_t k, uint64_t m, uint64_t seed) {
    assert(m < 32);
    assert(m <= k);
    uint64_t min_hash = uint64_t(-1);
    uint64_t minimizer = uint64_t(-1);
    uint64_t mask = (uint64_t(1) << (2 * m)) - 1;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        uint64_t sub_kmer = kmer & mask;
        uint64_t hash = Hasher::hash(sub_kmer, seed);
        if (hash < min_hash) {
            min_hash = hash;
            minimizer = sub_kmer;
        }
        kmer >>= 2;
    }
    return minimizer;
}

}  // namespace util

// taken from tlx
std::istream& appendline(std::istream& is, std::string& str, char delim = '\n') {
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