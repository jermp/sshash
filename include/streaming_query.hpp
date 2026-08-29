#pragma once

#include "include/dictionary.hpp"
#include "include/minimizer_iterator.hpp"
#include "include/util.hpp"

namespace sshash {

template <typename Dict>
struct streaming_query  //
{
    using kmer_t = typename Dict::kmer_type;

    streaming_query(Dict const* dict)

        : m_dict(dict)

        , m_start(true)
        , m_kmer(constants::invalid_uint64)
        , m_kmer_rc(constants::invalid_uint64)
        , m_k(dict->m_k)
        , m_m(dict->m_m)

        , m_minimizer_it(dict->m_k, dict->m_m, dict->m_hasher)
        , m_curr_mini_info()

        , m_it(dict->m_spss.strings, m_k)
        , m_remaining_string_bases(0)

        , m_last_seed_mini_info()
        , m_last_seed_positive(false)

        , m_num_searches(0)
        , m_num_extensions(0)
        , m_num_invalid(0)
        , m_num_negative(0)
        , m_num_skipped_singleton_lookups(0)
        , m_num_bucket_cache_hits(0)

    {}

    void reset() {
        m_start = true;
        m_remaining_string_bases = 0;
        m_res = lookup_result();
        m_minimizer_it.reset();
        m_last_seed_mini_info = minimizer_info();
        m_last_seed_positive = false;
        m_bucket_cache.size = 0;
    }

    lookup_result lookup(char const* kmer)  //
    {
        /* 1. validation */
        bool is_valid =
            m_start ? util::is_valid<kmer_t>(kmer, m_k) : kmer_t::is_valid(kmer[m_k - 1]);
        if (!is_valid) {
            m_num_invalid += 1;
            reset();
            return m_res;
        }

        /* 2. compute uint kmer and kmer_rc from input char kmer and minimizers */
        if (!m_start) {
            m_kmer.drop_char();
            m_kmer.set(m_k - 1, kmer_t::char_to_uint(kmer[m_k - 1]));
            assert(m_kmer == util::string_to_uint_kmer<kmer_t>(kmer, m_k));
            m_kmer_rc.pad_char();
            m_kmer_rc.set(0, kmer_t::char_to_uint(
                                 kmer_t::canonicalize_basepair_reverse_map[int(kmer[m_k - 1])]));
            m_kmer_rc.take(m_k * kmer_t::bits_per_char);
        } else {
            m_kmer = util::string_to_uint_kmer<kmer_t>(kmer, m_k);
            m_kmer_rc = m_kmer;
            m_kmer_rc.reverse_complement_inplace(m_k);
        }

        m_curr_mini_info = m_minimizer_it.next(m_kmer);

        /* 3. compute result */
        if (m_remaining_string_bases == 0) {
            seed();
        } else {
            auto expected_kmer = (m_res.kmer_orientation == constants::forward_orientation)
                                     ? (m_it.next(), m_it.get())
                                     : (m_it.next_reverse(), m_it.get_reverse());
            if ((expected_kmer == m_kmer) or (expected_kmer == m_kmer_rc)) {
                ++m_num_extensions;
                m_res.kmer_id += m_res.kmer_orientation;
                m_res.kmer_id_in_string += m_res.kmer_orientation;
                m_remaining_string_bases -= 1;
            } else {
                seed();
            }
        }

        /* 4. update state */
        m_start = false;

        assert(equal_lookup_result(m_dict->lookup(kmer), m_res));
        return m_res;
    }

    uint64_t num_searches() const { return m_num_searches; }
    uint64_t num_extensions() const { return m_num_extensions; }
    uint64_t num_positive_lookups() const { return num_searches() + num_extensions(); }
    uint64_t num_negative_lookups() const { return m_num_negative; }
    uint64_t num_invalid_lookups() const { return m_num_invalid; }
    uint64_t num_skipped_singleton_lookups() const { return m_num_skipped_singleton_lookups; }
    uint64_t num_bucket_cache_hits() const { return m_num_bucket_cache_hits; }

private:
    Dict const* m_dict;

    /* result */
    lookup_result m_res;

    /* kmer state */
    bool m_start;
    kmer_t m_kmer, m_kmer_rc;
    uint64_t m_k, m_m;

    /* minimizer state */
    minimizer_iterator<kmer_t> m_minimizer_it;
    minimizer_info m_curr_mini_info;

    /* string state */
    kmer_iterator<kmer_t, bits::bit_vector> m_it;
    uint64_t m_remaining_string_bases;

    /*
        Last real seed state, for the shortcuts of `seed`: the minimizer of
        the last `m_dict->lookup` call (with, for the singleton shortcut, the stream
        position of its occurrence and whether a positive match anchors it),
        and the cached decoded locate set of that minimizer.
    */
    minimizer_info m_last_seed_mini_info;
    bool m_last_seed_positive;

    /* performance counts */
    uint64_t m_num_searches;
    uint64_t m_num_extensions;
    uint64_t m_num_invalid;
    uint64_t m_num_negative;
    uint64_t m_num_skipped_singleton_lookups;
    uint64_t m_num_bucket_cache_hits;

    /* large and cold: keep it last, away from the per-kmer members above */
    typename Dict::spss_type::bucket_cache m_bucket_cache;

    /* Whether the minimizer of the last seed is its own reverse complement
       (possible for even m only). Rarely needed, so computed on demand. */
    bool self_rc_minimizer() const {
        if constexpr (kmer_t::has_reverse_complement) {
            return kmer_t::reverse_complement_mmer(m_last_seed_mini_info.minimizer, m_m) ==
                   m_last_seed_mini_info.minimizer;
        }
        return false;
    }

    /* Outlined: `seed` inlines into `lookup` at two call sites, and letting
       its body (grown by the shortcuts) bloat the extension fast path measurably
       slows down streams that rarely seed. */
    __attribute__((noinline)) void seed()  //
    {
        m_remaining_string_bases = 0;

        if (m_curr_mini_info.minimizer == m_last_seed_mini_info.minimizer)  //
        {
            /* The minimizer was absent at the last seed: any kmer having the
               same minimizer is surely absent as well. */
            if (m_res.minimizer_found == false) {
                assert(m_res.kmer_id == constants::invalid_uint64);
                m_num_negative += 1;
                return;
            }

            /*
                The minimizer is present but the kmer may not be -- the most
                common negative in a positive stream, e.g., a sequencing error
                inside the kmer. Two shortcuts avoid the full lookup.

                1. The last seed matched via a singleton bucket, at position L
                   say, and the current kmer still carries the very same
                   minimizer occurrence (same `pos_in_seq`: the scheme is
                   forward, so a sampled position, once abandoned, is never
                   re-selected). Then the only locus the kmer could occupy is
                   L - pos_in_kmer, which the just-failed extension step
                   already compared and rejected (or it lies beyond the string
                   boundary): negative, with no text access at all. The mirror
                   locus L - (k - m - pos_in_kmer) hosts the reverse
                   complement of the minimizer occurrence, not the occurrence
                   itself, so it is excluded too -- unless the minimizer is
                   its own reverse complement (possible for even m only), in
                   which case we fall through to the verification below.
                   Only a singleton bucket admits this free skip: a larger
                   bucket has candidate loci the failed extension never
                   examined, so those must be verified against the text.
            */
            if (m_last_seed_positive and m_bucket_cache.size == 1 and
                m_curr_mini_info.pos_in_seq == m_last_seed_mini_info.pos_in_seq and
                !self_rc_minimizer())  //
            {
                m_res = lookup_result();
                m_num_negative += 1;
                m_num_skipped_singleton_lookups += 1;
                return;
            }

            /*  2. The decoded locate set of the bucket is cached: verify the
                   candidate loci directly against the text, with no codeword
                   lookup nor offset decoding. */
            if (m_bucket_cache.size != 0) {
                m_res = m_dict->m_spss.lookup_from_positions(  //
                    m_bucket_cache.positions.data(), m_bucket_cache.size, m_kmer, m_kmer_rc,
                    m_curr_mini_info);
                m_num_bucket_cache_hits += 1;
                if (m_res.kmer_id == constants::invalid_uint64) {
                    m_num_negative += 1;
                    return;
                }
                /* keep the singleton shortcut anchored to the latest match */
                m_last_seed_mini_info = m_curr_mini_info;
                m_last_seed_positive = true;
                m_num_searches += 1;
                begin_extension();
                return;
            }

            /* heavy bucket: fall through to the full lookup */
        }

        m_res = m_dict->lookup(m_kmer, m_kmer_rc, m_curr_mini_info, &m_bucket_cache);
        m_last_seed_mini_info = m_curr_mini_info;
        m_last_seed_positive = m_res.kmer_id != constants::invalid_uint64;

        if (m_res.kmer_id == constants::invalid_uint64) {
            m_num_negative += 1;
            return;
        }

        assert(m_res.minimizer_found == true);
        m_num_searches += 1;
        begin_extension();
    }

    void begin_extension() {
        uint64_t kmer_offset = 2 * (m_res.kmer_id + m_res.string_id * (m_k - 1));
        m_remaining_string_bases =
            (m_res.string_end - m_res.string_begin - m_k) - m_res.kmer_id_in_string;
        if (m_res.kmer_orientation == constants::backward_orientation) {
            kmer_offset += 2 * m_k;
            m_remaining_string_bases = m_res.kmer_id_in_string;
        }
        m_it.at(kmer_offset);
    }
};

}  // namespace sshash
