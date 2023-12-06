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
            uint64_t hash =
                pthash::MurmurHash2_64(reinterpret_cast<char const*>(&low), sizeof(low), seed) ^
                pthash::MurmurHash2_64(reinterpret_cast<char const*>(&high), sizeof(high), ~seed);
            return hash;
        }
    }
};

struct kmers_pthash_hasher_128 {
    typedef pthash::hash128 hash_type;

    /* specialization for kmer_t */
    static inline pthash::hash128 hash(kmer_t x, uint64_t seed) {
        if constexpr (constants::uint_kmer_bits == 64) {
            return {pthash::MurmurHash2_64(reinterpret_cast<char const*>(&x), sizeof(x), seed),
                    pthash::MurmurHash2_64(reinterpret_cast<char const*>(&x), sizeof(x), ~seed)};
        } else {
            assert(constants::uint_kmer_bits == 128);
            uint64_t low = static_cast<uint64_t>(x);
            uint64_t high = static_cast<uint64_t>(x >> 64);
            return {
                pthash::MurmurHash2_64(reinterpret_cast<char const*>(&low), sizeof(low), seed) ^
                    pthash::MurmurHash2_64(reinterpret_cast<char const*>(&high), sizeof(high),
                                           ~seed),
                pthash::MurmurHash2_64(reinterpret_cast<char const*>(&low), sizeof(low), seed + 1) ^
                    pthash::MurmurHash2_64(reinterpret_cast<char const*>(&high), sizeof(high),
                                           ~(seed + 1))};
        }
    }
};

// typedef pthash::murmurhash2_64 minimizers_base_hasher_type;
typedef pthash::murmurhash2_128 minimizers_base_hasher_type;

// typedef kmers_pthash_hasher_64 kmers_base_hasher_type;
typedef kmers_pthash_hasher_128 kmers_base_hasher_type;

// typedef pthash::single_phf<minimizers_base_hasher_type,    // base hasher
//                            pthash::dictionary_dictionary,  // encoder type
//                            true                            // minimal output
//                            >
//     minimizers_pthash_type;
typedef pthash::partitioned_phf<minimizers_base_hasher_type,    // base hasher
                                pthash::dictionary_dictionary,  // encoder type
                                true                            // minimal output
                                >
    minimizers_pthash_type;

// typedef pthash::single_phf<kmers_base_hasher_type,         // base hasher
//                            pthash::dictionary_dictionary,  // encoder type
//                            true                            // minimal output
//                            >
//     kmers_pthash_type;
typedef pthash::partitioned_phf<kmers_base_hasher_type,         // base hasher
                                pthash::dictionary_dictionary,  // encoder type
                                true                            // minimal output
                                >
    kmers_pthash_type;

/* used to hash m-mers and determine the minimizer of a k-mer */
struct murmurhash2_64 {
    /* specialization for uint64_t */
    static inline uint64_t hash(uint64_t x, uint64_t seed) {
        return pthash::MurmurHash2_64(reinterpret_cast<char const*>(&x), sizeof(x), seed);
    }
};

}  // namespace sshash