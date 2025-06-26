#pragma once

#include "kmer.hpp"

namespace sshash {

/*
    "Re-scan" method.
*/
template <class kmer_t>
struct minimizer_enumerator {
    minimizer_enumerator() {}

    minimizer_enumerator(uint64_t k, uint64_t m, hasher_type const& hasher)
        : m_k(k)
        , m_m(m)
        , m_position(0)
        , m_min_value(constants::invalid_uint64)
        , m_min_position(0)
        , m_min_hash(constants::invalid_uint64)
        , m_hasher(hasher) {}

    void set_position(uint64_t position) { m_position = position; }

    template <bool reverse = false>
    uint64_t next(kmer_t kmer, bool clear) {
        _next<reverse>(kmer, clear);
        return m_min_value;
    }

    template <bool reverse = false>
    std::pair<uint64_t, uint64_t>  // (minimizer, pos of minimizer in sequence)
    next_with_pos(kmer_t kmer, bool clear) {
        _next<reverse>(kmer, clear);
        return {m_min_value, m_min_position};
    }

private:
    uint64_t m_k, m_m;
    uint64_t m_position;
    uint64_t m_min_value, m_min_position, m_min_hash;
    hasher_type m_hasher;

    template <bool reverse = false>
    void _next(kmer_t kmer, bool clear) {
        if (clear) {
            rescan<reverse>(kmer);
        } else {
            m_position += 1;
            // std::cout << "m_min_position + (m_k - m_m + 1) = " << (m_min_position + (m_k - m_m +
            // 1))
            //           << "; m_position = " << m_position << std::endl;
            if (m_min_position + (m_k - m_m + 1) <= m_position) {
                /* min leaves the window: re-scan to compute the new min */
                m_position = m_min_position + 1;
                rescan<reverse>(kmer);
            } else {
                kmer_t mmer = kmer;
                if constexpr (reverse) {
                    mmer.take_chars(m_m);
                } else {
                    mmer.drop_chars(m_k - m_m);
                }
                uint64_t hash = m_hasher.hash(uint64_t(mmer));
                if (hash < m_min_hash) {
                    m_min_hash = hash;
                    m_min_value = uint64_t(mmer);
                    m_min_position = m_position;
                }
            }
        }
        assert(m_min_value == util::compute_minimizer<kmer_t>(kmer, m_k, m_m, m_hasher));
    }

    template <bool reverse = false>
    void rescan(kmer_t kmer) {
        // std::cout << "RESCAN at position " << m_position << std::endl;
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
            uint64_t hash = m_hasher.hash(uint64_t(mmer));
            if (hash < m_min_hash) {
                m_min_hash = hash;
                m_min_value = uint64_t(mmer);
                m_min_position = m_position;
            }
            m_position += 1;
        }
        m_position -= 1;
    }
};

}  // namespace sshash