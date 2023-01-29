#pragma once

#include "../external/pthash/include/pthash.hpp"
#include "constants.hpp"

namespace sshash {

struct kmers_pthash_hasher_64 {
    typedef pthash::hash64 hash_type;

    /* specialization for kmer_t */
    static inline pthash::hash64 hash(kmer_t x, uint64_t seed) {
        if constexpr (constants::uint_kmer_bits == 64) {
            return pthash::MurmurHash2_64(reinterpret_cast<char const*>(&x), sizeof(x), seed);
        } else {
            assert(constants::uint_kmer_bits == 128);
            uint64_t low = static_cast<uint64_t>(x);
            uint64_t high = static_cast<uint64_t>(x >> 64);
            uint64_t hash = pthash::MurmurHash2_64(reinterpret_cast<char const*>(&low),
                                                   sizeof(uint64_t), seed) ^
                            pthash::MurmurHash2_64(reinterpret_cast<char const*>(&high),
                                                   sizeof(uint64_t), seed);
            return hash;
        }
    }
};

typedef pthash::murmurhash2_64 minimizers_base_hasher_type;
// typedef pthash::murmurhash2_128 minimizers_base_hasher_type;

typedef kmers_pthash_hasher_64 kmers_base_hasher_type;
// typedef kmers_pthash_hasher_128 kmers_base_hasher_type;

typedef pthash::single_phf<minimizers_base_hasher_type,    // base hasher
                           pthash::dictionary_dictionary,  // encoder type
                           true                            // minimal output
                           >
    minimizers_pthash_type;

typedef pthash::single_phf<kmers_base_hasher_type,         // base hasher
                           pthash::dictionary_dictionary,  // encoder type
                           true                            // minimal output
                           >
    kmers_pthash_type;

struct murmurhash2_64 {
    /* specialization for kmer_t */
    static inline uint64_t hash(kmer_t x, uint64_t seed) {
        if constexpr (constants::uint_kmer_bits == 64) {
            return pthash::MurmurHash2_64(reinterpret_cast<char const*>(&x), sizeof(x), seed);
        } else {
            assert(constants::uint_kmer_bits == 128);
            uint64_t low = static_cast<uint64_t>(x);
            uint64_t high = static_cast<uint64_t>(x >> 64);
            uint64_t hash = pthash::MurmurHash2_64(reinterpret_cast<char const*>(&low),
                                                   sizeof(uint64_t), seed) ^
                            pthash::MurmurHash2_64(reinterpret_cast<char const*>(&high),
                                                   sizeof(uint64_t), seed);
            return hash;
        }
    }
};

namespace util {

// TODO: move in PTHash
static inline void check_hash_collision_probability(uint64_t size) {
    /*
        Adapted from: https://preshing.com/20110504/hash-collision-probabilities.
        Given a universe of size U (total number of possible hash values),
        which is U = 2^b for b-bit hash codes,
        the collision probability for n keys is (approximately):
            1 - e^{-n(n-1)/(2U)}.
        For example, for U=2^32 (32-bit hash codes), this probability
        gets to 50% already for n = 77,163 keys.
        We can approximate 1-e^{-X} with X when X is sufficiently small.
        Then our collision probability is
            n(n-1)/(2U) ~ n^2/(2U).
        So it can derived that for ~1.97B keys and 64-bit hash codes,
        the probability of collision is ~0.1 (10%), which may not be
        so small for some applications.
        For n = 2^30, the probability of collision is ~0.031 (3.1%).
    */
    if (sizeof(minimizers_base_hasher_type::hash_type) * 8 == 64 and size > (1ULL << 30)) {
        throw std::runtime_error(
            "Using 64-bit hash codes with more than 2^30 keys can be dangerous due to "
            "collisions: use 128-bit hash codes instead.");
    }
}

}  // namespace util

}  // namespace sshash