#pragma once

#include "kmer.hpp"

namespace sshash {

/*
    "Re-scan" method.
*/
template <class kmer_t>
struct minimizer_iterator {
    minimizer_iterator() {}

    minimizer_iterator(uint64_t k, uint64_t m, hasher_type const& hasher)
        : m_k(k)
        , m_m(m)
        , m_position(0)
        , m_min_position_in_kmer(0)
        , m_min_value(constants::invalid_uint64)
        , m_min_position(0)
        , m_min_hash(constants::invalid_uint64)
        , m_hasher(hasher)  //
    {
        assert(k > 0 and m <= k);
    }

    void set_position(uint64_t position) { m_position = position; }

    template <bool reverse = false>
    uint64_t next(kmer_t kmer, bool clear) {
        _next<reverse>(kmer, clear);
        return m_min_value;
    }

    template <bool reverse = false>
    minimizer_info next_advanced(kmer_t kmer, bool clear) {
        _next<reverse>(kmer, clear);
        return {m_min_value, m_min_position, m_min_position_in_kmer};
    }

private:
    uint64_t m_k, m_m;
    uint64_t m_position, m_min_position_in_kmer;
    uint64_t m_min_value, m_min_position, m_min_hash;
    hasher_type m_hasher;

    template <bool reverse = false>
    void _next(kmer_t kmer, bool clear) {
        if (clear) {
            rescan<reverse>(kmer);
        } else {
            m_position += 1;
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

                if constexpr (reverse) {
                    if (hash <= m_min_hash) {
                        m_min_hash = hash;
                        m_min_value = uint64_t(mmer);
                        m_min_position = m_position;
                        m_min_position_in_kmer = 0;
                    } else {
                        m_min_position_in_kmer += 1;
                        assert(m_min_position_in_kmer <= m_k - m_m);
                    }
                } else {
                    if (hash < m_min_hash) {
                        m_min_hash = hash;
                        m_min_value = uint64_t(mmer);
                        m_min_position = m_position;
                        m_min_position_in_kmer = m_k - m_m;
                    } else {
                        assert(m_min_position_in_kmer > 0);
                        m_min_position_in_kmer -= 1;
                    }
                }
            }
        }

        // auto got = minimizer_info{m_min_value, constants::invalid_uint64,
        // m_min_position_in_kmer}; auto expected = util::compute_minimizer<kmer_t>(kmer, m_k, m_m,
        // m_hasher); if (got != expected) {
        //     std::cout << "kmer " << util::uint_kmer_to_string<kmer_t>(kmer, m_k) << std::endl;
        //     std::cout << "expected minimizer = "
        //               << util::uint_minimizer_to_string<kmer_t>(expected.minimizer, m_m)
        //               << std::endl;
        //     std::cout << "got minimizer = "
        //               << util::uint_minimizer_to_string<kmer_t>(got.minimizer, m_m) << std::endl;
        //     std::cout << "expected pos in kmer = " << expected.position_in_kmer << std::endl;
        //     std::cout << "got pos in kmer = " << got.position_in_kmer << std::endl;
        // }

        assert((minimizer_info{m_min_value, constants::invalid_uint64, m_min_position_in_kmer}) ==
               util::compute_minimizer<kmer_t>(kmer, m_k, m_m, m_hasher));
    }

    template <bool reverse = false>
    void rescan(kmer_t kmer) {
        m_min_hash = constants::invalid_uint64;
        m_min_value = constants::invalid_uint64;
        m_min_position_in_kmer = 0;
        for (uint64_t i = 0; i != m_k - m_m + 1; ++i) {
            kmer_t mmer = kmer;
            if constexpr (reverse) {
                mmer.drop_chars(m_k - m_m - i);
            } else {
                kmer.drop_char();
            }
            mmer.take_chars(m_m);
            uint64_t hash = m_hasher.hash(uint64_t(mmer));

            if constexpr (reverse) {  // rightmost
                if (hash <= m_min_hash) {
                    m_min_hash = hash;
                    m_min_value = uint64_t(mmer);
                    m_min_position = m_position;
                    m_min_position_in_kmer = m_k - m_m - i;
                }
            } else {  // leftmost
                if (hash < m_min_hash) {
                    m_min_hash = hash;
                    m_min_value = uint64_t(mmer);
                    m_min_position = m_position;
                    m_min_position_in_kmer = i;
                }
            }

            m_position += 1;
        }
        m_position -= 1;
    }
};

}  // namespace sshash