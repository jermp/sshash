#pragma once

#include "include/dictionary.hpp"
#include "include/minimizer_iterator.hpp"
#include "include/util.hpp"

namespace sshash {

/*
    Streaming membership/location queries, one kmer per call, each kmer
    sharing its first k-1 characters with the previous one.

    The query is a two-mode machine.

    - SEEDING mode: the kmer, its reverse complement, and the minimizer of
      the current window are maintained incrementally, one O(1) update per
      call, exactly as a plain streaming query would. Lookups are resolved
      either by extending the previous match along its string or by `seed()`.
      This mode makes the k negatives of an error window and alien stretches
      pay one m-mer hash per base, not one hash per window locus.

    - EXTENDING mode, entered after w = k - m + 1 consecutive successful
      extensions: no kmer, no reverse complement, no minimizer, no hashing.
      The only per-base work is comparing the next character of the query
      against the next expected character of the string, read from a
      word-wide buffer of upcoming string characters (`m_buf`, refilled every
      64 / bits_per_char bases with one word read; for backward extensions
      the word is reverse-complemented at refill time, so the per-base
      compare is orientation-independent). On a mismatch or exhaustion, all
      state is rebuilt from scratch from the current window -- the `kmer`
      pointer holds all k characters for the duration of the call -- at a
      cost of O(k) once per seed: `string_to_uint_kmer`, one reverse
      complement, and one re-scan of the window's m-mers through
      `minimizer_iterator` (the single source of truth for hashing and
      tie-breaking), positioned with `set_position` so that `pos_in_seq`
      is exactly what the incremental iterator would have produced.

    Since seeds number about one per string junction plus one per error,
    extension-dominated streams run almost entirely in EXTENDING mode.
*/
template <typename Dict>
struct streaming_query  //
{
    using kmer_t = typename Dict::kmer_type;

    static constexpr uint64_t bits_per_char = kmer_t::bits_per_char;
    static constexpr uint64_t chars_per_word = 64 / bits_per_char;
    static constexpr uint64_t char_mask = (uint64_t(1) << bits_per_char) - 1;

    streaming_query(Dict const* dict)

        : m_dict(dict)

        , m_start(true)
        , m_extending(false)
        , m_kmer(constants::invalid_uint64)
        , m_kmer_rc(constants::invalid_uint64)
        , m_k(dict->m_k)
        , m_m(dict->m_m)

        , m_minimizer_it(dict->m_k, dict->m_m, dict->m_hasher)
        , m_curr_mini_info()
        , m_window(0)

        , m_remaining_string_bases(0)
        , m_buf(0)
        , m_buf_avail(0)
        , m_buf_cursor(0)
        , m_run(0)

        , m_last_seed_mini_info()
        , m_last_seed_positive(false)

        , m_num_searches(0)
        , m_num_extensions(0)
        , m_num_invalid(0)
        , m_num_negative(0)
        , m_num_skipped_singleton_lookups(0)
        , m_num_bucket_cache_hits(0)
        , m_num_mode_switches(0)

#ifndef NDEBUG
        , m_shadow_it(dict->m_k, dict->m_m, dict->m_hasher)
#endif

    {}

    void reset() {
        /* The window counter mirrors the position accounting of the legacy
           always-incremental minimizer iterator, whose `reset()` kept
           `m_position`: after processing window W the iterator stood at
           W + (k - m), where the next re-scan would begin. */
        if (!m_start) m_window += m_k - m_m;
        m_start = true;
        m_extending = false;
        m_remaining_string_bases = 0;
        m_buf_avail = 0;
        m_run = 0;
        m_res = lookup_result();
        m_last_seed_mini_info = minimizer_info();
        m_last_seed_positive = false;
        m_bucket_cache.size = 0;
#ifndef NDEBUG
        m_shadow_it.reset();
#endif
    }

    lookup_result lookup(char const* kmer)  //
    {
        /* 1. validation */
        const bool is_valid =
            m_start ? util::is_valid<kmer_t>(kmer, m_k) : kmer_t::is_valid(kmer[m_k - 1]);
        if (!is_valid) {
            m_num_invalid += 1;
            reset();
            return m_res;
        }

        if (!m_start) m_window += 1;

#ifndef NDEBUG
        /* Shadow of the legacy always-incremental iterator: the lazy
           rebuilds below must reproduce its output exactly, `pos_in_seq`
           included. */
        m_shadow_mini_info = m_shadow_it.next(util::string_to_uint_kmer<kmer_t>(kmer, m_k));
#endif

        if (m_extending)  //
        {
            assert(!m_start);
            /* 2. fast path: one character compare, nothing else */
            if (m_remaining_string_bases != 0) {
                if (m_buf_avail == 0) fill_buffer();
                if ((m_buf & char_mask) == kmer_t::char_to_uint(kmer[m_k - 1])) {
                    extend();
                    assert(equal_lookup_result(m_dict->lookup(kmer), m_res));
                    return m_res;
                }
            }
            rebuild_and_seed(kmer);
        } else {
            lookup_seeding(kmer);
        }

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
    uint64_t num_mode_switches() const { return m_num_mode_switches; }

private:
    Dict const* m_dict;

    /* result */
    lookup_result m_res;

    /* kmer state: valid in SEEDING mode only, rebuilt from the input
       characters on the EXTENDING-mode cold path */
    bool m_start;
    bool m_extending;
    kmer_t m_kmer, m_kmer_rc;
    uint64_t m_k, m_m;

    /* minimizer state: valid in SEEDING mode only, rebuilt by a re-scan of
       the current window on the EXTENDING-mode cold path. `m_window` is the
       index of the current window in the stream of windows, maintained on
       every call so that rebuilds reproduce `pos_in_seq` exactly. */
    minimizer_iterator<kmer_t> m_minimizer_it;
    minimizer_info m_curr_mini_info;
    uint64_t m_window;

    /* string state: the expected next characters of the current string, as
       a word of codes consumed lowest-first. For forward extensions these
       are the string's characters at increasing offsets; for backward
       extensions, the complements of its characters at decreasing offsets.
       `m_buf_cursor` is the character offset of the next refill: the first
       character after (forward), or the one past the last character of
       (backward), the words still to verify. */
    uint64_t m_remaining_string_bases;
    uint64_t m_buf;
    uint64_t m_buf_avail;
    uint64_t m_buf_cursor;
    uint64_t m_run;  // consecutive successful extensions since the last seed

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
    uint64_t m_num_mode_switches;

    /* large and cold: keep it last, away from the per-kmer members above */
    typename Dict::spss_type::bucket_cache m_bucket_cache;

#ifndef NDEBUG
    minimizer_iterator<kmer_t> m_shadow_it;
    minimizer_info m_shadow_mini_info;
#endif

    kmer_t reverse_complement(kmer_t kmer) const {
        kmer.reverse_complement_inplace(m_k);
        return kmer;
    }

    /* The character offset of the current kmer in `strings`, from the result
       fields that extensions maintain (`m_res.kmer_offset` itself goes stale
       while extending). */
    uint64_t curr_kmer_offset() const {
        return m_res.kmer_id + m_res.string_id * (m_k - 1);
    }

    void extend() {
        m_buf >>= bits_per_char;
        m_buf_avail -= 1;
        m_num_extensions += 1;
        m_res.kmer_id += m_res.kmer_orientation;
        m_res.kmer_id_in_string += m_res.kmer_orientation;
        m_remaining_string_bases -= 1;
    }

    /* The whole SEEDING mode, outlined: `lookup` then holds little beyond the
       EXTENDING fast path and stays small enough to inline into the callers'
       loops, while `next` inlines here, its only hot call site. */
    __attribute__((noinline)) void lookup_seeding(char const* kmer)  //
    {
        /* 2. compute uint kmer and kmer_rc from input char kmer and minimizer */
        if (!m_start) {
            m_kmer.drop_char();
            m_kmer.set(m_k - 1, kmer_t::char_to_uint(kmer[m_k - 1]));
            m_kmer_rc.pad_char();
            const unsigned char c = kmer[m_k - 1];
            m_kmer_rc.set(0, kmer_t::char_to_uint(kmer_t::canonicalize_basepair_reverse_map[c]));
            m_kmer_rc.take(m_k * bits_per_char);
        } else {
            m_kmer = util::string_to_uint_kmer<kmer_t>(kmer, m_k);
            m_kmer_rc = m_kmer;
            m_kmer_rc.reverse_complement_inplace(m_k);
            m_minimizer_it.set_position(m_window);
            m_start = false;
        }
        m_curr_mini_info = m_minimizer_it.next(m_kmer);
        assert(m_kmer == util::string_to_uint_kmer<kmer_t>(kmer, m_k));
        assert(m_kmer_rc == reverse_complement(m_kmer));
        assert(m_curr_mini_info == m_shadow_mini_info);

        /* 3. compute result */
        if (m_remaining_string_bases != 0) {
            if (m_buf_avail == 0) fill_buffer();
            if ((m_buf & char_mask) == kmer_t::char_to_uint(kmer[m_k - 1])) {
                extend();
                if (++m_run == m_k - m_m + 1) {
                    /* the incremental state would be stale from here on */
                    m_extending = true;
                    m_num_mode_switches += 1;
                }
                return;
            }
            seed(/* from_failed_extension */ true);
        } else {
            seed(/* from_failed_extension */ false);
        }
    }

    /* The cold path of EXTENDING mode: rebuild every piece of state from the
       current window, then seed. All k characters are readable from `kmer`.
       Outlined like `seed`, to keep the extension fast path small. */
    __attribute__((noinline)) void rebuild_and_seed(char const* kmer)  //
    {
        const bool from_failed_extension = m_remaining_string_bases != 0;
        m_kmer = util::string_to_uint_kmer<kmer_t>(kmer, m_k);
        m_kmer_rc = m_kmer;
        m_kmer_rc.reverse_complement_inplace(m_k);
        m_minimizer_it.set_position(m_window);
        m_curr_mini_info = m_minimizer_it.next(m_kmer);
        assert(m_curr_mini_info == m_shadow_mini_info);
        m_extending = false;
        m_run = 0;
        m_num_mode_switches += 1;
        seed(from_failed_extension);
    }

    __attribute__((noinline)) void fill_buffer() {
        assert(m_remaining_string_bases != 0);
        auto const& strings = m_dict->m_spss.strings;
        m_buf_avail = std::min<uint64_t>(chars_per_word, m_remaining_string_bases);
        if (m_res.kmer_orientation == constants::forward_orientation) {
            assert(bits_per_char * m_buf_cursor < strings.num_bits());
            m_buf = strings.get_word64(bits_per_char * m_buf_cursor);
            m_buf_cursor += m_buf_avail;
        } else {
            if constexpr (kmer_t::has_reverse_complement) {
                static_assert(bits_per_char == 2);
                /* one word ending at the cursor, reverse-complemented:
                   code i of `m_buf` is the complement of the character at
                   offset cursor - 1 - i, which is the character the i-th
                   backward extension from here appends to the query */
                uint64_t word;
                if (m_buf_cursor >= chars_per_word) {
                    word = strings.get_word64(bits_per_char * (m_buf_cursor - chars_per_word));
                } else {
                    word = strings.get_word64(0)
                           << (bits_per_char * (chars_per_word - m_buf_cursor));
                }
                m_buf = kmer_t::crc64(word);
                m_buf_cursor -= m_buf_avail;
            } else {
                assert(false);  // unreachable: no orientations without complements
            }
        }
    }

    /* Whether the minimizer of the last seed is its own reverse complement
       (possible for even m only). Rarely needed, so computed on demand. */
    bool self_rc_minimizer() const {
        if constexpr (kmer_t::has_reverse_complement) {
            return kmer_t::reverse_complement_mmer(m_last_seed_mini_info.minimizer, m_m) ==
                   m_last_seed_mini_info.minimizer;
        }
        return false;
    }

    /* Only called right after a failed extension: whether the string, at the
       locus the extension just rejected, carries the current kmer in the
       orientation the extension did NOT test (a hairpin-like strand flip).
       The buffer compare rejects the extension's own orientation only, while
       the singleton skip of `seed` needs both rejected; when this returns
       true the skip is forfeited and the kmer is re-found by the regular
       verification below it, with the flipped orientation. */
    bool flipped_match_at_failed_locus() const {
        if constexpr (kmer_t::has_reverse_complement) {
            const uint64_t failed_offset = curr_kmer_offset() + m_res.kmer_orientation;
            auto read_kmer = util::read_kmer_at<kmer_t>(m_dict->m_spss.strings, m_k,
                                                        bits_per_char * failed_offset);
            return read_kmer == (m_res.kmer_orientation == constants::forward_orientation
                                     ? m_kmer_rc
                                     : m_kmer);
        }
        return false;
    }

    /* Outlined: `seed` inlines into `lookup` at two call sites, and letting
       its body (grown by the shortcuts) bloat the extension fast path measurably
       slows down streams that rarely seed. */
    __attribute__((noinline)) void seed(bool from_failed_extension)  //
    {
        if (m_curr_mini_info.minimizer == m_last_seed_mini_info.minimizer)  //
        {
            /* The minimizer was absent at the last seed: any kmer having the
               same minimizer is surely absent as well. */
            if (m_res.minimizer_found == false) {
                assert(m_res.kmer_id == constants::invalid_uint64);
                m_remaining_string_bases = 0;
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
                   One caveat of the word-buffered extension: its per-base
                   compare rejects a locus for the extension's own orientation
                   only, so a failed extension step must also exclude a match
                   there with the opposite orientation before the skip applies
                   (`flipped_match_at_failed_locus`, one text read, only
                   evaluated when everything else admits the skip).
            */
            if (m_last_seed_positive and m_bucket_cache.size == 1 and
                m_curr_mini_info.pos_in_seq == m_last_seed_mini_info.pos_in_seq and
                !self_rc_minimizer() and
                !(from_failed_extension and flipped_match_at_failed_locus()))  //
            {
                m_remaining_string_bases = 0;
                m_res = lookup_result();
                m_num_negative += 1;
                m_num_skipped_singleton_lookups += 1;
                return;
            }

            m_remaining_string_bases = 0;

            /*  2. The decoded locate set of the bucket is cached: verify the
                   candidate loci directly against the text, with no codeword
                   lookup nor offset decoding. */
            if (m_bucket_cache.size != 0) {
                if (m_bucket_cache.size > 1) {
                    m_dict->m_spss.prefetch_positions(m_bucket_cache.positions.data(),
                                                      m_bucket_cache.size, m_curr_mini_info);
                }
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

        m_remaining_string_bases = 0;

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
        const uint64_t kmer_offset = curr_kmer_offset();
        assert(kmer_offset == m_res.kmer_offset);
        if (m_res.kmer_orientation == constants::backward_orientation) {
            m_remaining_string_bases = m_res.kmer_id_in_string;
            m_buf_cursor = kmer_offset;
        } else {
            m_remaining_string_bases =
                (m_res.string_end - m_res.string_begin - m_k) - m_res.kmer_id_in_string;
            m_buf_cursor = kmer_offset + m_k;
        }
        m_buf_avail = 0;
        m_run = 0;
    }
};

}  // namespace sshash
