#pragma once

#include "external/pthash/include/pthash.hpp"
#include "constants.hpp"

namespace sshash {

template <class kmer_t>
struct kmers_pthash_hasher_64 {
    typedef pthash::hash64 hash_type;

    /* specialization for kmer_t */
    static inline pthash::hash64 hash(kmer_t x, uint64_t seed) {
        uint64_t hash = 0;
        for (int i = 0; i < kmer_t::uint_kmer_bits; i += 64) {
            uint64_t block = x.pop64();
            hash ^= pthash::MurmurHash2_64(reinterpret_cast<char const*>(&block), sizeof(block),
                                           seed + i);
        }
        return hash;
    }
};

template <class kmer_t>
struct kmers_pthash_hasher_128 {
    typedef pthash::hash128 hash_type;

    /* specialization for kmer_t */
    static inline pthash::hash128 hash(kmer_t x, uint64_t seed) {
        uint64_t hash_first = 0;
        uint64_t hash_second = 0;
        for (int i = 0; i < kmer_t::uint_kmer_bits; i += 64) {
            uint64_t block = x.pop64();
            hash_first ^= pthash::MurmurHash2_64(reinterpret_cast<char const*>(&block),
                                                 sizeof(block), seed + i);
            hash_second ^= pthash::MurmurHash2_64(reinterpret_cast<char const*>(&block),
                                                  sizeof(block), ~seed + i);
        }
        return {hash_first, hash_second};
    }
};

using minimizers_base_hasher_type = pthash::murmurhash2_128;

using minimizers_pthash_type =
    pthash::partitioned_phf<minimizers_base_hasher_type,                    // base hasher
                            pthash::skew_bucketer,                          // bucketer type
                            pthash::dictionary_dictionary,                  // encoder type
                            true,                                           // minimal output
                            pthash::pthash_search_type::xor_displacement>;  // search type>;

template <class kmer_t>
using kmers_base_hasher_type = kmers_pthash_hasher_128<kmer_t>;

template <class kmer_t>
using kmers_pthash_type =
    pthash::partitioned_phf<kmers_base_hasher_type<kmer_t>,                 // base hasher
                            pthash::skew_bucketer,                          // bucketer type
                            pthash::dictionary_dictionary,                  // encoder type
                            true,                                           // minimal output
                            pthash::pthash_search_type::xor_displacement>;  // search type>;

/* used to hash m-mers and determine the minimizer of a k-mer */
struct murmurhash2_64 {
    /* specialization for uint64_t */
    static inline uint64_t hash(uint64_t x, uint64_t seed) {
        return pthash::MurmurHash2_64(reinterpret_cast<char const*>(&x), sizeof(x), seed);
    }
};

}  // namespace sshash