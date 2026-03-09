#pragma once

#include "kmer.hpp"

namespace sshash {

/*
    "Re-scan" method.
*/
template <typename Kmer>
struct minimizer_iterator {
    minimizer_iterator() {}

    minimizer_iterator(uint64_t k, uint64_t m, hasher_type const& hasher, uint64_t position = 0)
        : m_k(k)
        , m_m(m)
        , m_min_value(constants::invalid_uint64)
        , m_min_hash(constants::invalid_uint64)
        , m_hasher(hasher)  //
    {
        assert(k > 0 and m <= k);
        set_position(position);
    }

    void set_position(uint64_t position) {
        m_position = position;
        reset();
    }

    void reset() {
        m_min_pos_in_kmer = 0;
        m_min_position = m_position - 1;
    }

    minimizer_info next(Kmer kmer) {
        if (m_min_pos_in_kmer == 0) {
            /* min leaves the window: re-scan to compute the new min */
            m_position = m_min_position + 1;
            rescan(kmer);
        } else {
            m_position += 1;
            Kmer mmer = kmer;
            mmer.drop_chars(m_k - m_m);
            uint64_t hash = m_hasher.hash(uint64_t(mmer));
            if (hash < m_min_hash) {
                m_min_hash = hash;
                m_min_value = uint64_t(mmer);
                m_min_position = m_position;
                m_min_pos_in_kmer = m_k - m_m;
            } else {
                assert(m_min_pos_in_kmer > 0);
                m_min_pos_in_kmer -= 1;
            }
        }

        assert(minimizer_info(m_min_value, m_min_pos_in_kmer) ==
               util::compute_minimizer<Kmer>(kmer, m_k, m_m, m_hasher));

        return {m_min_value, m_min_position, m_min_pos_in_kmer};
    }

private:
    uint64_t m_k, m_m;
    uint64_t m_position, m_min_pos_in_kmer;
    uint64_t m_min_value, m_min_position, m_min_hash;
    hasher_type m_hasher;

    void rescan(Kmer kmer) {
        m_min_hash = constants::invalid_uint64;
        m_min_value = constants::invalid_uint64;
        m_min_pos_in_kmer = 0;
        uint64_t begin = m_position;
        for (uint64_t i = 0; i != m_k - m_m + 1; ++i, ++m_position) {
            Kmer mmer = kmer;
            kmer.drop_char();
            mmer.take_chars(m_m);
            uint64_t hash = m_hasher.hash(uint64_t(mmer));
            if (hash < m_min_hash) {  // leftmost
                m_min_hash = hash;
                m_min_value = uint64_t(mmer);
                m_min_pos_in_kmer = i;
            }
        }
        m_position -= 1;
        m_min_position = begin + m_min_pos_in_kmer;
    }
};

/*
    "Re-scan" method.
*/
template <typename Kmer>
struct minimizer_iterator_rc {
    minimizer_iterator_rc() {}

    minimizer_iterator_rc(uint64_t k, uint64_t m, hasher_type const& hasher, uint64_t position = 0)
        : m_k(k)
        , m_m(m)
        , m_min_value(constants::invalid_uint64)
        , m_min_hash(constants::invalid_uint64)
        , m_hasher(hasher)  //
    {
        assert(k > 0 and m <= k);
        set_position(position);
    }

    void set_position(uint64_t position) {
        m_position = position;
        reset();
    }

    void reset() {
        m_min_pos_in_kmer = m_k - m_m;
        m_min_position = m_position - 1;
    }

    minimizer_info next(Kmer kmer) {
        if (m_min_pos_in_kmer == m_k - m_m) {
            /* min leaves the window: re-scan to compute the new min */
            m_position = m_min_position + 1;
            rescan(kmer);
        } else {
            m_position += 1;
            Kmer mmer = kmer;
            mmer.take_chars(m_m);
            uint64_t hash = m_hasher.hash(uint64_t(mmer));
            if (hash <= m_min_hash) {
                m_min_hash = hash;
                m_min_value = uint64_t(mmer);
                m_min_position = m_position;
                m_min_pos_in_kmer = 0;
            } else {
                m_min_pos_in_kmer += 1;
                assert(m_min_pos_in_kmer <= m_k - m_m);
            }
        }

        assert(minimizer_info(m_min_value, m_min_pos_in_kmer) ==
               util::compute_minimizer<Kmer>(kmer, m_k, m_m, m_hasher));

        return {m_min_value, m_min_position, m_min_pos_in_kmer};
    }

private:
    uint64_t m_k, m_m;
    uint64_t m_position, m_min_pos_in_kmer;
    uint64_t m_min_value, m_min_position, m_min_hash;
    hasher_type m_hasher;

    void rescan(Kmer kmer) {
        m_min_hash = constants::invalid_uint64;
        m_min_value = constants::invalid_uint64;
        m_min_pos_in_kmer = 0;
        uint64_t begin = m_position;
        for (int64_t i = m_k - m_m; i >= 0; --i, ++m_position) {
            Kmer mmer = kmer;
            mmer.drop_chars(i);
            mmer.take_chars(m_m);
            uint64_t hash = m_hasher.hash(uint64_t(mmer));
            if (hash <= m_min_hash) {  // rightmost
                m_min_hash = hash;
                m_min_value = uint64_t(mmer);
                m_min_pos_in_kmer = i;
            }
        }
        m_position -= 1;
        m_min_position = begin + (m_k - m_min_pos_in_kmer - m_m);
    }
};

}  // namespace sshash