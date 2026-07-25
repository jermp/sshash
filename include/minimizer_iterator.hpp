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
        if (m_min_pos_in_kmer == 0) {
            /* min leaves the window: re-scan to compute the new min */
            m_position = m_min_position + 1;
            rescan(kmer);
        } else {
            m_position += 1;
            Kmer mmer = kmer;
            mmer.drop_chars(m_k - m_m);
            uint64_t value = util::canonical_mmer<Kmer>(uint64_t(mmer), m_m);
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
        if (m_num_mins > 1) break_tie(kmer, mini_info);

        assert(minimizer_info(mini_info.minimizer, mini_info.pos_in_kmer) ==
               util::compute_minimizer<Kmer>(kmer, m_k, m_m, m_hasher));

        return mini_info;
    }

private:
    uint64_t m_k, m_m;
    uint64_t m_position, m_min_pos_in_kmer;
    uint64_t m_min_value, m_min_position, m_min_hash;
    uint64_t m_num_mins;
    hasher_type m_hasher;

    void rescan(Kmer kmer) {
        const uint64_t begin = m_position;

        /* first locus, peeled off the loop: see `util::compute_minimizer` */
        {
            Kmer mmer = kmer;
            kmer.drop_char();
            mmer.take_chars(m_m);
            m_min_value = util::canonical_mmer<Kmer>(uint64_t(mmer), m_m);
            m_min_hash = m_hasher.hash(m_min_value);
            m_min_pos_in_kmer = 0;
            m_num_mins = 1;
            ++m_position;
        }

        for (uint64_t i = 1; i != m_k - m_m + 1; ++i, ++m_position) {
            Kmer mmer = kmer;
            kmer.drop_char();
            mmer.take_chars(m_m);
            uint64_t value = util::canonical_mmer<Kmer>(uint64_t(mmer), m_m);
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
        Two or more loci of the window attain the minimum hash, and the leftmost
        of them is the one currently held. Leftmost is the right answer when the
        kmer is canonical; otherwise the canonical frame is rc(kmer), whose
        leftmost tied locus is this window's rightmost one.

        This is the only place where the whole kmer, rather than just its m-mers,
        has to be reverse-complemented. It runs on ~1e-5 of the windows.
    */
    void break_tie(Kmer kmer, minimizer_info& mini_info) const {
        assert(m_num_mins > 1);
        if (util::is_canonical<Kmer>(kmer, m_k)) return;

        const uint64_t window_begin = m_min_position - m_min_pos_in_kmer;
        uint64_t pos_in_kmer = m_min_pos_in_kmer;
        uint64_t value = m_min_value;

        Kmer window = kmer;
        window.drop_chars(m_min_pos_in_kmer + 1);
        for (uint64_t i = m_min_pos_in_kmer + 1; i != m_k - m_m + 1; ++i) {
            Kmer mmer = window;
            mmer.take_chars(m_m);
            uint64_t v = util::canonical_mmer<Kmer>(uint64_t(mmer), m_m);
            if (m_hasher.hash(v) == m_min_hash) {  // rightmost
                pos_in_kmer = i;
                value = v;
            }
            window.drop_char();
        }

        mini_info = minimizer_info(value, window_begin + pos_in_kmer, pos_in_kmer);
    }
};

}  // namespace sshash
