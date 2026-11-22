#pragma once

#include "include/dictionary.hpp"
#include "include/minimizer_iterator.hpp"
#include "include/util.hpp"

namespace sshash {

/*
    Mini-batched streaming queries: the routing walks of a block of
    sequences, overlapped.

    One-at-a-time streaming resolves a sequence as an alternation of cheap
    extensions and expensive seeds, and every seed -- a chain of dependent
    cache misses through the minimizers MPHF, the control codewords and the
    text -- must finish before the next kmer is even looked at, so on large
    indexes the query is bound by serial memory latency. This class keeps
    the walk-sparing structure (extension across minimizer changes, the
    same-occurrence memo, cached locate sets) but collects the walks of many
    sequences and resolves them together, stage by stage, so their misses
    overlap.

    Sequences are enqueued (`enqueue_sequence`) into an internal arena and
    processed in blocks of about `flush_kmers` kmers (`flush`; automatic
    when the arena fills). Per block:

    - Front end: each sequence is cut at invalid characters (kmers spanning
      one are counted invalid, as in `streaming_query`); the kmers of the
      valid segments, their reverse complements and their minimizers are
      computed incrementally into flat arrays, and consecutive kmers sharing
      a minimizer occurrence are grouped.

    - The block is covered by chains: one per segment, plus one every
      `chain_length` kmers inside long segments. A chain owns a range of
      consecutive kmers and advances through it like a streaming query:
      successful extensions (one character compare against a word-wide
      buffer of upcoming text, orientation-independent) cross group
      boundaries freely and cost no walk; a kmer that cannot be resolved
      from the chain's current state stops the chain with a walk request.

    - Passes: the pending requests of all stopped chains are resolved
      together in stages -- the minimizers MPHF over all of them through
      phastr's `lookup_batch` (which prefetches the seeds of every bucket
      before reading any), a prefetch pass then a read pass over the
      control codewords, a prefetch pass over singleton text loci and
      mid-load locate sets, then the verifications -- and the chains run
      again. Chains that stop in a non-matching region ask for the walks of
      several upcoming minimizer occurrences at once (1, 4, 16: negative
      stretches have no extension to spare walks, so speculation only
      overlaps walks that would happen anyway, while a re-anchored chain
      falls back to single requests). The number of pending chains shrinks
      geometrically with the passes.

    - A resolved walk leaves the chain a memo -- whether the minimizer is
      present, and the decoded locate set of its bucket when it fits the
      cache -- against which the kmers of the same minimizer occurrence
      that extension does not cover are verified directly, with no further
      walk: sequencing-error windows cost one walk per minimizer
      occurrence, not per kmer.

    Results are identical to `dictionary::lookup` on every kmer (debug
    assert). Counters, for report compatibility: a positive kmer resolved
    by an extension step counts as an "extension", every other positive as
    a "search"; `num_bucket_cache_hits` counts kmers verified against a
    cached locate set.
*/
template <typename Dict>
struct streaming_query_batched  //
{
    using kmer_t = typename Dict::kmer_type;

    static constexpr bool is_batched = true;
    static constexpr uint64_t flush_kmers = 1 << 15;
    static constexpr uint64_t chain_length = 1024;
    static constexpr uint64_t mphf_batch_size = 16;
    static constexpr uint32_t max_speculation = 16;
    static constexpr uint64_t char_mask = (uint64_t(1) << kmer_t::bits_per_char) - 1;

    streaming_query_batched(Dict const* dict)

        : m_dict(dict)
        , m_k(dict->m_k)
        , m_m(dict->m_m)

        , m_minimizer_it(dict->m_k, dict->m_m, dict->m_hasher)
        , m_window(0)
        , m_pending_kmers(0)

        , m_num_searches(0)
        , m_num_extensions(0)
        , m_num_invalid(0)
        , m_num_negative(0)
        , m_num_skipped_singleton_lookups(0)
        , m_num_bucket_cache_hits(0)

    {}

    /* Sequences are independent of each other; nothing to reset. */
    void reset() {}

    /* Copies the sequence into the arena; flushes when the arena is full.
       Counters are complete only after a final `flush()`. A sequence longer
       than one block is cut into overlapping pieces (by k - 1 characters,
       so every kmer appears in exactly one piece). */
    void enqueue_sequence(char const* chars, const uint64_t length) {
        assert(length >= m_k);
        uint64_t begin = 0;
        while (length - begin > flush_kmers + m_k - 1) {
            enqueue_piece(chars + begin, flush_kmers + m_k - 1);
            begin += flush_kmers;
        }
        enqueue_piece(chars + begin, length - begin);
    }

    /* One-shot: resolve the kmers of one sequence now (with whatever else
       is enqueued). */
    void lookup_sequence(char const* chars, const uint64_t length) {
        enqueue_sequence(chars, length);
        flush();
    }

    void enqueue_piece(char const* chars, const uint64_t length) {
        m_views.emplace_back(m_arena.size(), length);
        m_arena.append(chars, length);
        m_pending_kmers += length - m_k + 1;
        if (m_pending_kmers >= flush_kmers) flush();
    }

    void flush() {
        if (m_views.empty()) return;
        front_end();
        run_passes();
        m_arena.clear();
        m_views.clear();
        m_pending_kmers = 0;
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
    uint64_t m_window;  // monotone window index: `set_position` at each
                        // segment keeps `pos_in_seq` values distinct

    /* enqueued input */
    std::string m_arena;
    std::vector<std::pair<uint64_t, uint64_t>> m_views;  // (offset, length)
    uint64_t m_pending_kmers;

    /* The current block, slot = one kmer. Only the kmers themselves are
       kept per slot (the block must stay cache-resident): the reverse
       complement is recomputed at the rare verification sites, and the
       minimizer of a slot is reconstructed from its group, `pos_in_kmer`
       decreasing by one per slide within an occurrence. */
    std::vector<kmer_t> m_kmers;
    std::vector<uint32_t> m_group_of;          // slot -> group
    std::vector<uint32_t> m_group_head;        // group -> first slot
    std::vector<uint64_t> m_group_minimizer;   // group -> minimizer
    std::vector<uint32_t> m_group_head_pik;    // group -> pos_in_kmer at head

    minimizer_info slot_mini_info(const uint64_t slot) const {
        const uint32_t g = m_group_of[slot];
        assert(slot >= m_group_head[g]);
        assert(m_group_head_pik[g] >= slot - m_group_head[g]);
        return {m_group_minimizer[g], m_group_head_pik[g] - (slot - m_group_head[g])};
    }

    kmer_t slot_kmer_rc(const uint64_t slot) const {
        kmer_t rc = m_kmers[slot];
        rc.reverse_complement_inplace(m_k);
        return rc;
    }

    struct segment {
        uint64_t slot_begin, slot_end;
        char const* chars;
    };
    std::vector<segment> m_segments;

    using bucket_cache = typename Dict::spss_type::bucket_cache;

    /* The verified-ahead state of a run extension: a word of upcoming text
       characters (complemented, in descending text order, for backward
       matches, so the compare is orientation-independent), and the
       character offset `cursor` of its next refill. */
    struct extension {
        bool anchored = false;
        bool backward = false;
        uint64_t buf = 0;
        uint64_t avail = 0;
        uint64_t cursor = 0;
        uint64_t remaining = 0;  // extension steps left in the string
    };

    struct chain {
        uint64_t cursor, end;
        uint32_t segment;
        uint32_t spec;
        bool anchored_since_submit;
        extension ext;
        lookup_result cur;  // last positive result, advanced by extensions
        /* memo of the last consumed walk */
        bool memo_valid;
        uint32_t memo_group;
        bool memo_found;
        uint64_t memo_head_slot;
        lookup_result memo_head_res;
        bucket_cache memo_cache;
        /* prepaid requests of the current pass */
        uint32_t req_begin, req_count;
    };
    std::vector<chain> m_chains;
    std::vector<uint32_t> m_pending[2];  // chain ids, current and next pass

    struct request_buffer {
        std::vector<uint32_t> slot;  // the kmer the walk resolves
        std::vector<uint32_t> group;
        std::vector<uint8_t> found;
        std::vector<lookup_result> res;
        std::vector<bucket_cache> cache;
        void clear() {
            slot.clear();
            group.clear();
            found.clear();
            res.clear();
            cache.clear();
        }
        uint64_t size() const { return slot.size(); }
    };
    request_buffer m_req[2];  // current and next pass
    std::vector<uint64_t> m_minimizers, m_ids, m_codes;  // resolve scratch

    /* performance counts */
    uint64_t m_num_searches;
    uint64_t m_num_extensions;
    uint64_t m_num_invalid;
    uint64_t m_num_negative;
    uint64_t m_num_skipped_singleton_lookups;
    uint64_t m_num_bucket_cache_hits;

    /* ---------------- front end ---------------- */

    void front_end() {
        m_kmers.clear();
        m_group_of.clear();
        m_group_head.clear();
        m_group_minimizer.clear();
        m_group_head_pik.clear();
        m_segments.clear();
        m_chains.clear();

        for (auto [offset, length] : m_views) {
            char const* chars = m_arena.data() + offset;
            const uint64_t num_kmers = length - m_k + 1;
            uint64_t num_valid = 0;
            uint64_t pos = 0;
            while (pos != length) {
                while (pos != length and !kmer_t::is_valid(chars[pos])) ++pos;
                const uint64_t begin = pos;
                while (pos != length and kmer_t::is_valid(chars[pos])) ++pos;
                if (pos - begin >= m_k) {
                    fill_segment(chars + begin, pos - begin);
                    num_valid += pos - begin - m_k + 1;
                }
            }
            assert(num_valid <= num_kmers);
            m_num_invalid += num_kmers - num_valid;
        }
    }

    void fill_segment(char const* chars, const uint64_t length) {
        const uint64_t num_kmers = length - m_k + 1;
        const uint64_t slot_begin = m_kmers.size();
        m_segments.push_back({slot_begin, slot_begin + num_kmers, chars});
        const uint32_t segment_id = m_segments.size() - 1;

        m_minimizer_it.set_position(m_window);
        kmer_t kmer = util::string_to_uint_kmer<kmer_t>(chars, m_k);

        minimizer_info prev_info;
        for (uint64_t i = 0; i != num_kmers; ++i) {
            if (i != 0) {
                const unsigned char c = chars[i + m_k - 1];
                kmer.drop_char();
                kmer.set(m_k - 1, kmer_t::char_to_uint(c));
            }
            assert(kmer == util::string_to_uint_kmer<kmer_t>(chars + i, m_k));
            auto info = m_minimizer_it.next(kmer);
            m_window += 1;
            const uint64_t slot = slot_begin + i;
            if (i == 0 or info.minimizer != prev_info.minimizer or
                info.pos_in_seq != prev_info.pos_in_seq) {
                m_group_head.push_back(slot);
                m_group_minimizer.push_back(info.minimizer);
                m_group_head_pik.push_back(info.pos_in_kmer);
            }
            prev_info = info;
            m_kmers.push_back(kmer);
            m_group_of.push_back(m_group_head.size() - 1);
            assert(slot_mini_info(slot).minimizer == info.minimizer);
            assert(slot_mini_info(slot).pos_in_kmer == info.pos_in_kmer);
        }

        for (uint64_t begin = slot_begin; begin < slot_begin + num_kmers; begin += chain_length) {
            chain c;
            c.cursor = begin;
            c.end = std::min(begin + chain_length, slot_begin + num_kmers);
            c.segment = segment_id;
            c.spec = 1;
            c.anchored_since_submit = false;
            c.memo_valid = false;
            c.req_begin = c.req_count = 0;
            m_chains.push_back(c);
        }
    }

    /* ---------------- pass engine ---------------- */

    void run_passes() {
        uint32_t parity = 0;
        m_req[0].clear();
        m_req[1].clear();
        m_pending[0].clear();
        m_pending[1].clear();

        /* every chain submits its first walk */
        for (uint32_t c = 0; c != m_chains.size(); ++c) {
            submit(m_chains[c], m_group_of[m_chains[c].cursor], c, parity);
        }

        while (!m_pending[parity].empty()) {
            resolve_requests(m_req[parity]);
            const uint32_t next = parity ^ 1;
            m_req[next].clear();
            m_pending[next].clear();
            for (uint32_t c : m_pending[parity]) {
                run_chain(m_chains[c], c, parity);
            }
            parity = next;
        }
    }

    void submit(chain& c, const uint32_t group, const uint32_t chain_id, const uint32_t parity) {
        c.spec = c.anchored_since_submit ? 1 : std::min<uint32_t>(c.spec * 4, max_speculation);
        c.anchored_since_submit = false;

        auto& buf = m_req[parity];
        const uint32_t last_group = m_group_of[c.end - 1];
        const uint32_t count = std::min<uint32_t>(c.spec, last_group - group + 1);
        c.req_begin = buf.size();
        c.req_count = count;
        for (uint32_t t = 0; t != count; ++t) {
            const uint32_t g = group + t;
            const uint32_t head = t == 0 ? c.cursor : m_group_head[g];
            assert(m_group_of[head] == g);
            buf.slot.push_back(head);
            buf.group.push_back(g);
        }
        buf.found.resize(buf.size());
        buf.res.resize(buf.size());
        buf.cache.resize(buf.size());
        m_pending[parity].push_back(chain_id);
    }

    /* Resolve all pending walks, stage by stage, so that the memory
       accesses of the pass overlap across requests. */
    void resolve_requests(request_buffer& buf) {
        const uint64_t n = buf.size();
        if (n == 0) return;

        /* stage 1: the minimizers MPHF, batched */
        m_minimizers.resize(n);
        m_ids.resize(n);
        m_codes.resize(n);
        for (uint64_t r = 0; r != n; ++r) m_minimizers[r] = m_group_minimizer[buf.group[r]];
        auto const& mphf = m_dict->m_ssi.codewords.mphf;
        uint64_t r = 0;
        for (; r + mphf_batch_size <= n; r += mphf_batch_size) {
            mphf.template lookup_batch<mphf_batch_size>(&m_minimizers[r], &m_ids[r]);
        }
        for (; r != n; ++r) m_ids[r] = mphf(m_minimizers[r]);

        /* stage 2: the control codewords: prefetch all lines, then read */
        auto const& codewords = m_dict->m_ssi.codewords.control_codewords;
        {
            auto const* bytes = reinterpret_cast<char const*>(codewords.data().data());
            const uint64_t width = codewords.width();
            for (uint64_t i = 0; i != n; ++i) {
                __builtin_prefetch(bytes + ((m_ids[i] * width) >> 3));
            }
        }
        for (uint64_t i = 0; i != n; ++i) m_codes[i] = codewords.access(m_ids[i]);

        /* stage 3: the codeword of a singleton bucket gives its text locus,
           that of a mid-load bucket the address of its locate set: prefetch
           both before any decoding or verification */
        auto const& strings = m_dict->m_spss.strings;
        auto const& mid_load = m_dict->m_ssi.mid_load_buckets;
        auto const* mid_load_bytes = reinterpret_cast<char const*>(mid_load.data().data());
        const uint64_t mid_load_width = mid_load.width();
        for (uint64_t i = 0; i != n; ++i) {
            const uint64_t code = m_codes[i];
            if ((code & 1) == bucket_t::SINGLETON) {
                const uint64_t minimizer_offset = code >> 1;
                const uint64_t pos_in_kmer = slot_mini_info(buf.slot[i]).pos_in_kmer;
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

        /* stage 4: verification */
        for (uint64_t i = 0; i != n; ++i) {
            const uint64_t slot = buf.slot[i];
            const uint64_t code = m_codes[i];
            lookup_result res;
            buf.cache[i].size = 0;
            if ((code & 3) == bucket_t::HEAVYLOAD and (code & 1) != bucket_t::SINGLETON) {
                /* skew-index bucket: the walk is kmer-specific; rare */
                res = m_dict->lookup(m_kmers[slot], slot_kmer_rc(slot), slot_mini_info(slot));
            } else {
                auto it = make_bucket_iterator(code);
                res = m_dict->m_spss.lookup(it, m_kmers[slot], slot_kmer_rc(slot),
                                            slot_mini_info(slot), &buf.cache[i]);
            }
            buf.found[i] = res.minimizer_found;
            buf.res[i] = res;
        }
    }

    typename sparse_and_skew_index<kmer_t>::bucket_iterator  //
    make_bucket_iterator(const uint64_t code) const {
        auto const& ssi = m_dict->m_ssi;
        if ((code & 1) == bucket_t::SINGLETON) {
            return {&ssi, code >> 1, 1, bucket_t::SINGLETON};
        }
        assert((code & 3) == bucket_t::MIDLOAD);
        constexpr uint64_t mask = (uint64_t(1) << constants::min_l) - 1;
        const uint64_t c = code >> 2;
        const uint64_t bucket_size = (c & mask) + 2;
        const uint64_t bucket_id = c >> constants::min_l;
        const uint64_t bucket_begin =
            ssi.begin_buckets_of_size[bucket_size] + bucket_id * bucket_size;
        return {&ssi, bucket_begin, bucket_size, bucket_t::MIDLOAD};
    }

    /* Advance the chain until its range is exhausted or it stops with a new
       walk submission. */
    void run_chain(chain& c, const uint32_t chain_id, const uint32_t parity) {
        auto const& seg = m_segments[c.segment];
        auto const& buf = m_req[parity];

        while (c.cursor != c.end) {
            const uint64_t j = c.cursor;

            /* 1. extension: crosses minimizer changes, costs no walk */
            if (c.ext.anchored and c.ext.remaining != 0) {
                if (c.ext.avail == 0) fill_extension_buffer(c.ext);
                const uint64_t code_j =
                    kmer_t::char_to_uint(seg.chars[(j - seg.slot_begin) + m_k - 1]);
                if ((c.ext.buf & char_mask) == code_j) {
                    c.ext.buf >>= kmer_t::bits_per_char;
                    c.ext.avail -= 1;
                    c.ext.remaining -= 1;
                    c.cur.kmer_id += c.cur.kmer_orientation;
                    c.cur.kmer_id_in_string += c.cur.kmer_orientation;
                    m_num_extensions += 1;
                    assert(equal_lookup_result(
                        m_dict->lookup(seg.chars + (j - seg.slot_begin)), c.cur));
                    c.cursor += 1;
                    continue;
                }
            }

            const uint32_t g = m_group_of[j];

            /* 2. make the walk of group g the chain's memo */
            if (!(c.memo_valid and c.memo_group == g)) {
                if (c.req_count != 0 and g >= buf.group[c.req_begin] and
                    g < buf.group[c.req_begin] + c.req_count) {
                    const uint32_t i = c.req_begin + (g - buf.group[c.req_begin]);
                    assert(buf.group[i] == g);
                    c.memo_valid = true;
                    c.memo_group = g;
                    c.memo_found = buf.found[i];
                    c.memo_head_slot = buf.slot[i];
                    c.memo_head_res = buf.res[i];
                    c.memo_cache = buf.cache[i];
                } else {
                    submit(c, g, chain_id, parity ^ 1);
                    return;
                }
            }

            /* 3. resolve the kmer from the memo */
            lookup_result res;
            if (j == c.memo_head_slot) {
                res = c.memo_head_res;
            } else if (!c.memo_found) {
                /* absent minimizer: definitive (never set for HEAVYLOAD) */
                res = lookup_result(false);
            } else if (c.memo_cache.size != 0) {
                m_num_bucket_cache_hits += 1;
                res = m_dict->m_spss.lookup_from_positions(
                    c.memo_cache.positions.data(), c.memo_cache.size, m_kmers[j], slot_kmer_rc(j),
                    slot_mini_info(j));
            } else {
                /* heavy or over-large bucket: per-kmer walk; rare */
                res = m_dict->lookup(m_kmers[j], slot_kmer_rc(j), slot_mini_info(j));
            }

            if (res.kmer_id != constants::invalid_uint64) {
                m_num_searches += 1;
                c.cur = res;
                anchor(c.ext, res);
                c.anchored_since_submit = true;
            } else {
                m_num_negative += 1;
                c.ext.anchored = false;
            }
            assert(equal_lookup_result(m_dict->lookup(seg.chars + (j - seg.slot_begin)), res));
            c.cursor += 1;
        }
    }

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
};

}  // namespace sshash
