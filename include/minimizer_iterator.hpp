#pragma once

#include "kmer.hpp"
#include "util.hpp"

namespace sshash {

/*
    "Re-scan" method: computes the minimizer of each kmer of a sequence,
    sliding one character at a time.

    See `util::compute_minimizer` for the definition of the minimizer; this
    iterator computes exactly the same thing incrementally, which the assertion
    at the end of `next` checks.

    The only extra state compared to a plain forward minimizer is `m_num_mins`,
    the number of loci of the current window attaining the minimum hash: a tie
    cannot be broken by position alone without breaking mirror-equivariance, so
    when there is one we have to look at the kmer's own orientation.

    The kmer is reverse-complemented once per call and that single value serves
    every locus, via rc(x_i) = rc(x)_{k-m-i}, so no m-mer is ever
    reverse-complemented on its own.
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
        m_num_mins = 0;
    }

    minimizer_info next(Kmer kmer) {
        Kmer kmer_rc = kmer;
        kmer_rc.reverse_complement_inplace(m_k);

        if (m_min_pos_in_kmer == 0) {
            /* min leaves the window: re-scan to compute the new min */
            m_position = m_min_position + 1;
            rescan(kmer, kmer_rc);
        } else {
            m_position += 1;
            Kmer window = kmer;
            window.drop_chars(m_k - m_m);
            uint64_t value = util::canonical_mmer_at<Kmer>(window, kmer_rc, m_k, m_m, m_k - m_m);
            uint64_t hash = m_hasher.hash(value);
            if (hash < m_min_hash) {
                m_min_hash = hash;
                m_min_value = value;
                m_min_position = m_position;
                m_min_pos_in_kmer = m_k - m_m;
                m_num_mins = 1;
            } else {
                /* Only the leftmost minimum can leave the window without a
                   re-scan, so the count stays valid across the slide. */
                assert(m_min_pos_in_kmer > 0);
                m_min_pos_in_kmer -= 1;
                if (hash == m_min_hash) m_num_mins += 1;
            }
        }

        minimizer_info mini_info(m_min_value, m_min_position, m_min_pos_in_kmer);
        if (m_num_mins > 1) break_tie(kmer, kmer_rc, mini_info);

        assert(minimizer_info(mini_info.minimizer, mini_info.pos_in_kmer) ==
               util::compute_minimizer<Kmer>(kmer, kmer_rc, m_k, m_m, m_hasher));

        return mini_info;
    }

private:
    uint64_t m_k, m_m;
    uint64_t m_position, m_min_pos_in_kmer;
    uint64_t m_min_value, m_min_position, m_min_hash;
    uint64_t m_num_mins;
    hasher_type m_hasher;

    void rescan(Kmer kmer, Kmer const& kmer_rc) {
        const uint64_t begin = m_position;

        /* first locus, peeled off the loop: see `util::compute_minimizer` */
        {
            m_min_value = util::canonical_mmer_at<Kmer>(kmer, kmer_rc, m_k, m_m, 0);
            kmer.drop_char();
            m_min_hash = m_hasher.hash(m_min_value);
            m_min_pos_in_kmer = 0;
            m_num_mins = 1;
            ++m_position;
        }

        for (uint64_t i = 1; i != m_k - m_m + 1; ++i, ++m_position) {
            uint64_t value = util::canonical_mmer_at<Kmer>(kmer, kmer_rc, m_k, m_m, i);
            kmer.drop_char();
            uint64_t hash = m_hasher.hash(value);
            if (hash < m_min_hash) {  // leftmost
                m_min_hash = hash;
                m_min_value = value;
                m_min_pos_in_kmer = i;
                m_num_mins = 1;
            } else if (hash == m_min_hash) {
                m_num_mins += 1;
            }
        }

        m_position -= 1;
        m_min_position = begin + m_min_pos_in_kmer;
    }

    /*
        Two or more loci of the window attain the minimum hash (all carrying the
        same class, so only the sampled position is at stake). Apply the
        centre-closest tie-break of `util::compute_minimizer`: take the tied
        locus closest to the window centre, and between two equally close ones
        -- mirror images i and k-m-i -- the smaller index if kmer <= rc(kmer),
        the larger otherwise. Forward and mirror-equivariant. It runs on ~1e-4
        of the windows, so the rescan below costs nothing overall.

        Tied loci cannot precede the leftmost minimum, so the scan starts there.
    */
    void break_tie(Kmer kmer, Kmer const& kmer_rc, minimizer_info& mini_info) const {
        assert(m_num_mins > 1);

        const uint64_t window_begin = m_min_position - m_min_pos_in_kmer;
        const uint64_t two_c = m_k - m_m;
        uint64_t best_dist = constants::invalid_uint64;
        uint64_t lo = m_min_pos_in_kmer;
        uint64_t hi = m_min_pos_in_kmer;

        Kmer window = kmer;
        window.drop_chars(m_min_pos_in_kmer);
        for (uint64_t i = m_min_pos_in_kmer; i != m_k - m_m + 1; ++i) {
            const uint64_t v = util::canonical_mmer_at<Kmer>(window, kmer_rc, m_k, m_m, i);
            if (m_hasher.hash(v) == m_min_hash) {
                assert(v == m_min_value);
                const uint64_t dist = 2 * i > two_c ? 2 * i - two_c : two_c - 2 * i;
                if (dist < best_dist) {
                    best_dist = dist;
                    lo = i;
                    hi = i;
                } else if (dist == best_dist) {
                    hi = i;
                }
            }
            window.drop_char();
        }

        const uint64_t chosen = (lo == hi or !(kmer_rc < kmer)) ? lo : hi;
        mini_info = minimizer_info(m_min_value, window_begin + chosen, chosen);
    }
};

}  // namespace sshash
