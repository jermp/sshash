#pragma once

#include "include/dictionary.hpp"
#include "include/minimizer_iterator.hpp"
#include "include/util.hpp"

namespace sshash {

/*
    Mini-batched streaming queries: one routing walk per distinct minimizer
    occurrence per batch, with the walks of a batch overlapped.

    `lookup_sequence(chars, length)` resolves all length - k + 1 kmers of a
    sequence. The sequence is cut into maximal runs of valid characters;
    every kmer spanning an invalid character is counted invalid, exactly as
    in `streaming_query`. Each run is processed in batches of
    B = min(kmers left in the run, 128) consecutive kmers:

    1. the B kmers, their reverse complements and their minimizers are
       computed incrementally, one O(1) update per base;
    2. consecutive kmers sharing the same minimizer occurrence form a group,
       and a group is resolved by a single routing walk, whatever its size;
    3. the walks of the batch are staged, so that their memory accesses
       overlap instead of serializing: first the minimizers MPHF over all
       group minimizers (phastr's `lookup_batch` hashes every key and
       prefetches the seeds of their buckets before reading any of them),
       then the control codewords of all groups (independent accesses the
       hardware overlaps), and then the walk tails, with the text of every
       singleton bucket -- whose locus is known from the codeword alone --
       prefetched before any verification starts;
    4. each group's first kmer runs the regular verification
       (`spss::lookup`), which also caches the decoded locate set of the
       bucket; the other kmers of the group verify against the cached set
       directly, with no walk of their own, and a minimizer found absent for
       the first kmer answers the whole group at no further cost.

    Compared to `streaming_query`, the per-kmer serial dependence is broken:
    the cache-missing part of a batch is `num_groups` independent walks
    instead of an alternation of extensions and seeds in which every seed
    waits for the previous result. The class keeps no cross-call state
    besides counters, and its results are identical to `dictionary::lookup`
    on every kmer.

    Counters, for compatibility with `streaming_query` reports: a positive
    kmer counts as a "search" if it is the first positive of its group's
    walk and as an "extension" otherwise; `num_bucket_cache_hits` counts the
    kmers verified from a cached locate set.
*/
template <typename Dict>
struct streaming_query_batched  //
{
    using kmer_t = typename Dict::kmer_type;

    static constexpr bool is_batched = true;
    static constexpr uint64_t max_batch_size = 128;

    streaming_query_batched(Dict const* dict)

        : m_dict(dict)
        , m_k(dict->m_k)
        , m_m(dict->m_m)

        , m_minimizer_it(dict->m_k, dict->m_m, dict->m_hasher)
        , m_window(0)

        , m_num_searches(0)
        , m_num_extensions(0)
        , m_num_invalid(0)
        , m_num_negative(0)
        , m_num_skipped_singleton_lookups(0)
        , m_num_bucket_cache_hits(0)

    {}

    /* No cross-sequence state besides counters. */
    void reset() {}

    void lookup_sequence(char const* chars, const uint64_t length) {
        lookup_sequence(chars, length, [](uint64_t, lookup_result const&) {});
    }

    /* Resolves the length - k + 1 kmers of chars[0..length) and calls
       cb(i, result) for each, i being the position of the kmer in the
       sequence, in increasing order. */
    template <typename Callback>
    void lookup_sequence(char const* chars, const uint64_t length, Callback&& cb) {
        assert(length >= m_k);
        const uint64_t num_kmers = length - m_k + 1;
        uint64_t num_valid = 0;
        uint64_t pos = 0;
        while (pos != length) {
            while (pos != length and !kmer_t::is_valid(chars[pos])) ++pos;
            const uint64_t begin = pos;
            while (pos != length and kmer_t::is_valid(chars[pos])) ++pos;
            if (pos - begin >= m_k) {
                lookup_segment(chars + begin, pos - begin, begin, cb);
                num_valid += pos - begin - m_k + 1;
            }
        }
        assert(num_valid <= num_kmers);
        m_num_invalid += num_kmers - num_valid;
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
    uint64_t m_k, m_m;

    minimizer_iterator<kmer_t> m_minimizer_it;
    uint64_t m_window;  // monotone window index, so that `set_position` at
                        // each segment keeps `pos_in_seq` values distinct

    /* the current batch */
    std::array<kmer_t, max_batch_size> m_kmers;
    std::array<kmer_t, max_batch_size> m_kmers_rc;
    std::array<minimizer_info, max_batch_size> m_infos;
    std::array<uint64_t, max_batch_size> m_minimizers;
    std::array<uint64_t, max_batch_size> m_ids;
    std::array<uint64_t, max_batch_size> m_codes;
    std::array<uint32_t, max_batch_size + 1> m_group_begin;

    typename Dict::spss_type::bucket_cache m_cache;

    /* performance counts */
    uint64_t m_num_searches;
    uint64_t m_num_extensions;
    uint64_t m_num_invalid;
    uint64_t m_num_negative;
    uint64_t m_num_skipped_singleton_lookups;
    uint64_t m_num_bucket_cache_hits;

    template <typename Callback>
    void lookup_segment(char const* chars, const uint64_t length,
                        const uint64_t offset_in_sequence, Callback& cb)  //
    {
        assert(length >= m_k);
        const uint64_t num_kmers = length - m_k + 1;
        m_minimizer_it.set_position(m_window);

        kmer_t kmer = util::string_to_uint_kmer<kmer_t>(chars, m_k);
        kmer_t kmer_rc = kmer;
        kmer_rc.reverse_complement_inplace(m_k);

        for (uint64_t batch_begin = 0; batch_begin != num_kmers;) {
            const uint64_t batch_size = std::min<uint64_t>(max_batch_size,  //
                                                           num_kmers - batch_begin);
            for (uint64_t j = 0; j != batch_size; ++j) {
                const uint64_t i = batch_begin + j;
                if (i != 0) {
                    const unsigned char c = chars[i + m_k - 1];
                    kmer.drop_char();
                    kmer.set(m_k - 1, kmer_t::char_to_uint(c));
                    kmer_rc.pad_char();
                    kmer_rc.set(0,
                                kmer_t::char_to_uint(kmer_t::canonicalize_basepair_reverse_map[c]));
                    kmer_rc.take(m_k * kmer_t::bits_per_char);
                }
                assert(kmer == util::string_to_uint_kmer<kmer_t>(chars + i, m_k));
                m_kmers[j] = kmer;
                m_kmers_rc[j] = kmer_rc;
                m_infos[j] = m_minimizer_it.next(kmer);
                m_window += 1;
            }
            lookup_batch(chars + batch_begin, batch_size, offset_in_sequence + batch_begin, cb);
            batch_begin += batch_size;
        }
    }

    template <typename Callback>
    void lookup_batch(char const* chars, const uint64_t batch_size,
                      const uint64_t kmer_index_offset, Callback& cb)  //
    {
        /* 1. group the kmers by minimizer occurrence */
        uint32_t num_groups = 0;
        for (uint64_t j = 0; j != batch_size; ++j) {
            if (j == 0 or m_infos[j].minimizer != m_infos[j - 1].minimizer or
                m_infos[j].pos_in_seq != m_infos[j - 1].pos_in_seq) {
                m_group_begin[num_groups++] = j;
            }
        }
        m_group_begin[num_groups] = batch_size;

        /* 2. the minimizers MPHF over all groups, batched: phastr hashes the
              keys and prefetches the seeds of their buckets before reading
              any of them, so the lookups of the batch overlap */
        auto const& mphf = m_dict->m_ssi.codewords.mphf;
        for (uint32_t g = 0; g != num_groups; ++g) {
            m_minimizers[g] = m_infos[m_group_begin[g]].minimizer;
        }
        constexpr uint64_t mphf_batch_size = 16;
        uint32_t g = 0;
        for (; g + mphf_batch_size <= num_groups; g += mphf_batch_size) {
            mphf.template lookup_batch<mphf_batch_size>(&m_minimizers[g], &m_ids[g]);
        }
        for (; g != num_groups; ++g) m_ids[g] = mphf(m_minimizers[g]);

        /* 3. the control codewords of all groups: prefetch every line, then
              read, so the reads of the batch overlap */
        auto const& codewords = m_dict->m_ssi.codewords.control_codewords;
        {
            auto const* bytes = reinterpret_cast<char const*>(codewords.data().data());
            const uint64_t width = codewords.width();
            for (uint32_t g = 0; g != num_groups; ++g) {
                __builtin_prefetch(bytes + ((m_ids[g] * width) >> 3));
            }
        }
        for (uint32_t g = 0; g != num_groups; ++g) m_codes[g] = codewords.access(m_ids[g]);

        /* 4. the codeword of a singleton bucket gives its text locus, and
              that of a mid-load bucket the address of its locate set:
              prefetch both before any decoding or verification runs */
        auto const& strings = m_dict->m_spss.strings;
        auto const& mid_load = m_dict->m_ssi.mid_load_buckets;
        auto const* mid_load_bytes = reinterpret_cast<char const*>(mid_load.data().data());
        const uint64_t mid_load_width = mid_load.width();
        for (uint32_t g = 0; g != num_groups; ++g) {
            const uint64_t code = m_codes[g];
            if ((code & 1) == bucket_t::SINGLETON) {
                const uint64_t minimizer_offset = code >> 1;
                const uint64_t pos_in_kmer = m_infos[m_group_begin[g]].pos_in_kmer;
                const uint64_t char_offset =
                    minimizer_offset >= pos_in_kmer ? minimizer_offset - pos_in_kmer : 0;
                __builtin_prefetch(&strings.data()[(kmer_t::bits_per_char * char_offset) >> 6]);
            } else if ((code & 3) == bucket_t::MIDLOAD) {
                constexpr uint64_t mask = (uint64_t(1) << constants::min_l) - 1;
                const uint64_t c = code >> 2;
                const uint64_t bucket_size = (c & mask) + 2;
                const uint64_t bucket_id = c >> constants::min_l;
                const uint64_t bucket_begin =
                    m_dict->m_ssi.begin_buckets_of_size[bucket_size] + bucket_id * bucket_size;
                __builtin_prefetch(mid_load_bytes + ((bucket_begin * mid_load_width) >> 3));
            }
        }

        /* 5. resolve the groups */
        for (uint32_t h = 0; h != num_groups; ++h) {
            resolve_group(m_codes[h], m_group_begin[h], m_group_begin[h + 1], chars,
                          kmer_index_offset, cb);
        }
    }

    template <typename Callback>
    void resolve_group(const uint64_t code, const uint64_t begin, const uint64_t end,
                       char const* chars, const uint64_t kmer_index_offset, Callback& cb)  //
    {
        auto const& ssi = m_dict->m_ssi;
        bool group_has_positive = false;

        /* decode the codeword, as in sparse_and_skew_index::lookup */
        uint64_t bucket_begin = 0;
        uint64_t bucket_size = 0;
        bucket_t bucket_type = bucket_t::HEAVYLOAD;
        if ((code & 1) == bucket_t::SINGLETON) {
            bucket_begin = code >> 1;
            bucket_size = 1;
            bucket_type = bucket_t::SINGLETON;
        } else if ((code & 3) == bucket_t::MIDLOAD) {
            constexpr uint64_t mask = (uint64_t(1) << constants::min_l) - 1;
            const uint64_t c = code >> 2;
            bucket_size = (c & mask) + 2;
            const uint64_t bucket_id = c >> constants::min_l;
            assert(bucket_size < ssi.begin_buckets_of_size.size());
            bucket_begin = ssi.begin_buckets_of_size[bucket_size] + bucket_id * bucket_size;
            bucket_type = bucket_t::MIDLOAD;
        } else {
            assert((code & 3) == bucket_t::HEAVYLOAD);
            /* Skew-index bucket: the walk itself is kmer-specific (the kmers
               MPHF locates each kmer), so the group shares nothing. Rare. */
            for (uint64_t j = begin; j != end; ++j) {
                auto res = m_dict->lookup(m_kmers[j], m_kmers_rc[j], m_infos[j]);
                emit(j, res, group_has_positive, chars, kmer_index_offset, cb);
            }
            return;
        }

        using bucket_iterator = typename sparse_and_skew_index<kmer_t>::bucket_iterator;

        /* the group's first kmer walks; the bucket's decoded locate set is
           cached on the way */
        bucket_iterator it(&ssi, bucket_begin, bucket_size, bucket_type);
        auto res = m_dict->m_spss.lookup(it, m_kmers[begin], m_kmers_rc[begin], m_infos[begin],
                                         &m_cache);
        emit(begin, res, group_has_positive, chars, kmer_index_offset, cb);
        if (end - begin == 1) return;

        if (res.minimizer_found == false) {
            /* The minimizer is absent (definitive for non-HEAVYLOAD buckets):
               every kmer of the group is negative at no further cost. */
            for (uint64_t j = begin + 1; j != end; ++j) {
                emit(j, lookup_result(false), group_has_positive, chars, kmer_index_offset, cb);
            }
            return;
        }

        /*
            Run extension inside the group: the kmers are consecutive, so once
            one of them is located, each following kmer is expected at the
            next text locus of its string, and one character compare against a
            word-wide buffer of upcoming text characters verifies it (the
            overlap of k - 1 characters is implied by the previous match). A
            mismatch, an exhausted string, or a lost anchor falls back to the
            regular verification -- against the cached locate set when the
            bucket fits the cache, else by redoing the (by now cache-hot)
            walk -- which also re-anchors the extension when it is positive.
        */
        extension ext;
        if (res.kmer_id != constants::invalid_uint64) anchor(ext, res);

        for (uint64_t j = begin + 1; j != end; ++j) {
            if (ext.anchored and ext.remaining != 0) {
                if (ext.avail == 0) fill_extension_buffer(ext);
                constexpr uint64_t mask = (uint64_t(1) << kmer_t::bits_per_char) - 1;
                const uint64_t code_j = kmer_t::char_to_uint(chars[j + m_k - 1]);
                if ((ext.buf & mask) == code_j) {
                    ext.buf >>= kmer_t::bits_per_char;
                    ext.avail -= 1;
                    ext.remaining -= 1;
                    res.kmer_id += res.kmer_orientation;
                    res.kmer_id_in_string += res.kmer_orientation;
                    emit(j, res, group_has_positive, chars, kmer_index_offset, cb);
                    continue;
                }
            }
            /* fall back to the regular verification */
            if (m_cache.size != 0) {
                m_num_bucket_cache_hits += 1;
                res = m_dict->m_spss.lookup_from_positions(
                    m_cache.positions.data(), m_cache.size, m_kmers[j], m_kmers_rc[j], m_infos[j]);
            } else {
                bucket_iterator jt(&ssi, bucket_begin, bucket_size, bucket_type);
                res = m_dict->m_spss.lookup(jt, m_kmers[j], m_kmers_rc[j], m_infos[j], nullptr);
            }
            ext.anchored = false;
            if (res.kmer_id != constants::invalid_uint64) anchor(ext, res);
            emit(j, res, group_has_positive, chars, kmer_index_offset, cb);
        }
    }

    /* The verified-ahead state of a run extension: a word of upcoming text
       characters (complemented, in descending text order, for backward
       matches, so the compare is orientation-independent), and the character
       offset `cursor` of its next refill: the first character after
       (forward), or one past the last character of (backward), the
       characters still to verify. */
    struct extension {
        bool anchored = false;
        bool backward = false;
        uint64_t buf = 0;
        uint64_t avail = 0;
        uint64_t cursor = 0;
        uint64_t remaining = 0;  // extension steps left in the string
    };

    void anchor(extension& ext, lookup_result const& res) const {
        assert(res.kmer_id != constants::invalid_uint64);
        const uint64_t kmer_offset = res.kmer_id + res.string_id * (m_k - 1);
        ext.anchored = true;
        ext.buf = 0;
        ext.avail = 0;
        if (res.kmer_orientation == constants::backward_orientation) {
            ext.backward = true;
            ext.remaining = res.kmer_id_in_string;
            ext.cursor = kmer_offset;
        } else {
            ext.backward = false;
            ext.remaining = (res.string_end - res.string_begin - m_k) - res.kmer_id_in_string;
            ext.cursor = kmer_offset + m_k;
        }
    }

    void fill_extension_buffer(extension& ext) const {
        assert(ext.anchored and ext.remaining != 0);
        constexpr uint64_t bits_per_char = kmer_t::bits_per_char;
        constexpr uint64_t chars_per_word = 64 / bits_per_char;
        auto const& strings = m_dict->m_spss.strings;
        ext.avail = std::min<uint64_t>(chars_per_word, ext.remaining);
        if (!ext.backward) {
            assert(bits_per_char * ext.cursor < strings.num_bits());
            ext.buf = strings.get_word64(bits_per_char * ext.cursor);
            ext.cursor += ext.avail;
        } else {
            if constexpr (kmer_t::has_reverse_complement) {
                static_assert(bits_per_char == 2);
                /* one word ending at the cursor, reverse-complemented: code i
                   of `buf` is the complement of the character at offset
                   cursor - 1 - i, the character the i-th backward extension
                   from here appends to the query */
                uint64_t word;
                if (ext.cursor >= chars_per_word) {
                    word = strings.get_word64(bits_per_char * (ext.cursor - chars_per_word));
                } else {
                    word = strings.get_word64(0)
                           << (bits_per_char * (chars_per_word - ext.cursor));
                }
                ext.buf = kmer_t::crc64(word);
                ext.cursor -= ext.avail;
            } else {
                assert(false);  // unreachable: no orientations without complements
            }
        }
    }

    template <typename Callback>
    void emit(const uint64_t j, lookup_result res, bool& group_has_positive, char const* chars,
              const uint64_t kmer_index_offset, Callback& cb)  //
    {
        if (res.kmer_id != constants::invalid_uint64) {
            if (group_has_positive) {
                m_num_extensions += 1;
            } else {
                m_num_searches += 1;
                group_has_positive = true;
            }
        } else {
            m_num_negative += 1;
        }
        assert(equal_lookup_result(m_dict->lookup(chars + j), res));
        (void)chars;
        cb(kmer_index_offset + j, res);
    }
};

}  // namespace sshash
