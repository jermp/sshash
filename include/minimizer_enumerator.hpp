#pragma once

#include "kmer.hpp"

namespace sshash {

/*
    "Re-scan" method.
*/
template <class kmer_t, typename Hasher = murmurhash2_64>
struct minimizer_enumerator {
    minimizer_enumerator() {}

    minimizer_enumerator(uint64_t k, uint64_t m, uint64_t seed)
        : m_k(k)
        , m_m(m)
        , m_seed(seed)
        , m_position(0)
        , m_min_value(constants::invalid_uint64)
        , m_min_position(0)
        , m_min_hash(constants::invalid_uint64) {}

    template <bool reverse = false>
    uint64_t next(kmer_t kmer, bool clear) {
        if (clear) {
            rescan<reverse>(kmer);
        } else {
            m_position += 1;
            if (m_min_position <= (m_position + m_m - 1) - m_k) {
                /* min leaves the window: re-scan to compute the new min */
                rescan<reverse>(kmer);
            } else {
                kmer_t mmer = kmer;
                if constexpr (reverse) {
                    mmer.take_chars(m_m);
                } else {
                    mmer.drop_chars(m_k - m_m);
                }
                uint64_t hash = Hasher::hash(uint64_t(mmer), m_seed);
                // uint64_t hash = mix(uint64_t(mmer));
                if (hash < m_min_hash) {
                    m_min_hash = hash;
                    m_min_value = uint64_t(mmer);
                    m_min_position = m_position;
                }
            }
        }
        assert(m_min_value == util::compute_minimizer<kmer_t>(kmer, m_k, m_m, m_seed));
        return m_min_value;
    }

private:
    uint64_t m_k, m_m, m_seed;
    uint64_t m_position;
    uint64_t m_min_value, m_min_position, m_min_hash;

    template <bool reverse = false>
    void rescan(kmer_t kmer) {
        m_min_hash = constants::invalid_uint64;
        m_min_value = constants::invalid_uint64;
        for (uint64_t i = 0; i != m_k - m_m + 1; ++i) {
            kmer_t mmer = kmer;
            if constexpr (reverse) {
                mmer.drop_chars(m_k - m_m - i);
            } else {
                kmer.drop_char();
            }
            mmer.take_chars(m_m);
            uint64_t hash = Hasher::hash(uint64_t(mmer), m_seed);
            // uint64_t hash = mix(uint64_t(mmer));
            if (hash < m_min_hash) {
                m_min_hash = hash;
                m_min_value = uint64_t(mmer);
                m_min_position = m_position;
            }
            m_position += 1;
        }
    }
};

}  // namespace sshash