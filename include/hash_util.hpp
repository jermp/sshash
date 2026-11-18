#pragma once

#include "external/pthash/include/pthash.hpp"
#include "external/phastr/include/phastr/phastr.hpp"
#include "external/cityhash/cityhash.hpp"
#include "constants.hpp"

namespace sshash {

/*
    The minimal perfect hash functions of SSHash are PHastR functions (see
    external/phastr): each key is hashed to a 64-bit signature by a hasher
    that is constructed once from a seed, which it preprocesses, so that
    queries pay only the per-key work.

    Minimizers are 64-bit integers, hashed with phastr::default_hasher (an
    inlined SplitMix64-based mix). K-mers are hashed lane by lane with
    independent mixes of the same family, combined with a xor: for the
    default 128-bit k-mers, two independent mixes.
*/

using minimizers_mphf_type = phastr::function<phastr::default_hasher>;

template <typename Kmer>
struct kmers_phastr_hasher {
    static constexpr size_t num_lanes = (sizeof(((Kmer*)nullptr)->bits) + 7) / 8;
    uint64_t prepared[num_lanes];

    kmers_phastr_hasher() : kmers_phastr_hasher(0) {}
    explicit kmers_phastr_hasher(uint64_t seed) {
        for (size_t i = 0; i != num_lanes; ++i) {
            prepared[i] = phastr::splitmix64(seed + i * 0x9e3779b97f4a7c15ULL);
        }
    }

    __attribute__((always_inline)) inline uint64_t operator()(Kmer const& x) const {
        uint64_t lanes[num_lanes] = {};
        memcpy(lanes, &(x.bits), sizeof(x.bits));
        uint64_t h = phastr::splitmix64(lanes[0] + prepared[0]);
        for (size_t i = 1; i != num_lanes; ++i) { h ^= phastr::splitmix64(lanes[i] + prepared[i]); }
        return h;
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }
    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        for (size_t i = 0; i != num_lanes; ++i) visitor.visit(t.prepared[i]);
    }
};

template <typename Kmer>
using kmers_mphf_type = phastr::function<kmers_phastr_hasher<Kmer>>;

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
