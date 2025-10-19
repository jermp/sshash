#pragma once

#include "external/pthash/include/pthash.hpp"
#include "external/cityhash/cityhash.cpp"
#include "constants.hpp"

namespace sshash {

struct minimizers_city_hasher_128 {
    typedef pthash::hash128 hash_type;

    static inline pthash::hash128 hash(uint64_t const minimizer, uint64_t seed) {
        auto ret = CityMurmur(reinterpret_cast<char const*>(&minimizer),  //
                              sizeof(minimizer), {seed, ~seed});
        return {ret.first, ret.second};
    }
};

struct minimizers_xx_hasher_128 {
    typedef pthash::hash128 hash_type;

    static inline pthash::hash128 hash(uint64_t const minimizer, uint64_t seed) {
        /*
            Cannot use XXH128 directly because on some processors (e.g., AMD)
            it just does not work in Release mode, e.g., when compiling *without*
            sanitizers: -fsanitize=address -fno-omit-frame-pointer.
            We therefore rely on XXH64 that produces 64-bit hashes but does not
            seem to have any issue.
        */
        uint8_t const* begin = reinterpret_cast<uint8_t const*>(&minimizer);
        uint8_t const* end = begin + sizeof(minimizer);
        return {XXH64(begin, end - begin, seed), XXH64(begin, end - begin, ~seed)};
    }
};

// using minimizers_base_hasher_type = minimizers_xx_hasher_128;
using minimizers_base_hasher_type = minimizers_city_hasher_128;

using minimizers_pthash_type =        //
    pthash::partitioned_phf<          //
        minimizers_base_hasher_type,  // base hasher
        pthash::opt_bucketer,         // bucketer type
        pthash::compact,              // encoder type
        true                          // minimal output
        >;                            //

template <typename Kmer>
struct kmers_xx_hasher_128 {
    typedef pthash::hash128 hash_type;

    static inline pthash::hash128 hash(Kmer const x, uint64_t seed) {
        uint8_t const* begin = reinterpret_cast<uint8_t const*>(&(x.bits));
        uint8_t const* end = begin + sizeof(x.bits);
        return {XXH64(begin, end - begin, seed), XXH64(begin, end - begin, ~seed)};
    }
};

template <typename Kmer>
struct kmers_city_hasher_128 {
    typedef pthash::hash128 hash_type;

    static inline pthash::hash128 hash(Kmer const x, uint64_t seed) {
        auto ret = CityMurmur(reinterpret_cast<char const*>(&(x.bits)),  //
                              sizeof(x.bits), {seed, ~seed});
        return {ret.first, ret.second};
    }
};

// template <typename Kmer>
// using kmers_base_hasher_type = kmers_xx_hasher_128<Kmer>;

template <typename Kmer>
using kmers_base_hasher_type = kmers_city_hasher_128<Kmer>;

template <typename Kmer>
using kmers_pthash_type =              //
    pthash::partitioned_phf<           //
        kmers_base_hasher_type<Kmer>,  // base hasher
        pthash::opt_bucketer,          // bucketer type
        pthash::compact,               // encoder type
        true                           // minimal output
        >;                             //

struct mixer_64 {
    mixer_64() { seed(0); }
    mixer_64(const uint64_t seed) { this->seed(seed); }

    void seed(const uint64_t seed) { m_magic = pthash::xxhash_64::hash(seed, 0).first(); }

    /* specialization for uint64_t */
    inline uint64_t hash(uint64_t x) const { return (x * 0x517cc1b727220a95) ^ m_magic; }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visitor.visit(m_magic);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_magic);
    }

private:
    uint64_t m_magic;
};

/* used to hash m-mers of a k-mer to compute its minimizer */
using hasher_type = mixer_64;

}  // namespace sshash