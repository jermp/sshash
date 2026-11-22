#pragma once

#include "include/dictionary.hpp"
#include "include/minimizer_iterator.hpp"
#include "include/util.hpp"

namespace sshash {

/*
    Mini-batched streaming queries: word-parallel extension, lazy minimizers,
    and the routing walks of a block of sequences overlapped.

    One-at-a-time streaming pays three per-kmer costs: the front end (one
    m-mer hash per base to maintain the minimizer, plus the kmer and its
    reverse complement), the extension compare, and -- at every seed -- a
    chain of dependent cache misses through the minimizers MPHF, the control
    codewords and the text, which must finish before the next kmer is even
    looked at. This class attacks all three:

    - sequences are enqueued into an arena and processed in blocks of about
      `flush_kmers` kmers; a block is cut at invalid characters into
      segments (kmers spanning an invalid character are counted invalid),
      each segment is 2-bit packed once, and segments are covered by chains,
      one per segment plus one every `chain_length` kmers;

    - a chain that holds a match extends it with a word-parallel LCP of the
      packed query against the packed text (XOR and count-trailing-zeros,
      64 / bits_per_char characters per compare; for backward matches the
      text word is reverse-complemented first, one crc64), crossing
      minimizer changes freely: a run of consecutive positive kmers costs
      O(run / word size), with no hashing and no per-kmer state at all;

    - where a match stops, the chain scans forward with an incremental
      minimizer iterator (one hash per base, only here) and requests the
      walk of the minimizer occurrence it cannot resolve -- several upcoming
      occurrences at once in non-matching stretches (1, 4, 16), since those
      have no extension to spare walks;

    - the pending requests of all chains are resolved together, stage by
      stage (the minimizers MPHF through phastr's `lookup_batch`, which
      prefetches the seeds of every bucket before reading any; a prefetch
      pass then a read pass over the control codewords; prefetches of
      singleton text loci and mid-load locate sets; decoding of the locate
      sets with the text around every candidate locus prefetched; then the
      verifications), so their misses overlap, and the chains run again:
      the number of pending chains shrinks geometrically with the passes;

    - a resolved walk covers its whole minimizer occurrence: kmers of the
      occurrence that extension does not cover are verified directly
      against the memo (absence, or the cached locate set of the bucket),
      with no walk of their own, so an error window costs one walk per
      minimizer occurrence, not per kmer.

    Results are identical to `dictionary::lookup` on every kmer (debug
    assert). Counters, for report compatibility: kmers resolved by
    extension count as "extensions", every other positive as a "search";
    `num_bucket_cache_hits` counts kmers verified against a cached locate
    set.
*/
template <typename Dict>
struct streaming_query_batched  //
{
    using kmer_t = typename Dict::kmer_type;
    using bucket_cache = typename Dict::spss_type::bucket_cache;

    static constexpr bool is_batched = true;
    static constexpr uint64_t flush_kmers = 1 << 15;
    static constexpr uint64_t chain_length = 1024;
    static constexpr uint64_t mphf_batch_size = 16;
    static constexpr uint32_t max_speculation = 16;
    static constexpr uint64_t bits_per_char = kmer_t::bits_per_char;
    static constexpr uint64_t chars_per_word = 64 / bits_per_char;

    streaming_query_batched(Dict const* dict)

        : m_dict(dict)
        , m_k(dict->m_k)
        , m_m(dict->m_m)
        , m_iterator_proto(dict->m_k, dict->m_m, dict->m_hasher)
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
    minimizer_iterator<kmer_t> m_iterator_proto;

    /* enqueued input */
    std::string m_arena;
    std::vector<std::pair<uint64_t, uint64_t>> m_views;  // (offset, length)
    uint64_t m_pending_kmers;

    void enqueue_piece(char const* chars, const uint64_t length) {
        m_views.emplace_back(m_arena.size(), length);
        m_arena.append(chars, length);
        m_pending_kmers += length - m_k + 1;
        if (m_pending_kmers >= flush_kmers) flush();
    }

    /* ---------------- block data ---------------- */

    /* The packed text of the block's valid segments: each segment starts at
       a word boundary and is followed by padding words, so any 64-bit
       window of it can be read. */
    std::vector<uint64_t> m_packed;

    struct segment {
        uint64_t slot_begin, slot_end;  // the kmers of the block it holds
        uint64_t packed_word_begin;     // its packed characters
        char const* chars;              // its raw characters (debug asserts)
    };
    std::vector<segment> m_segments;

    uint64_t packed_word(const uint64_t bit) const {
        const uint64_t block = bit >> 6;
        const uint64_t shift = bit & 63;
        uint64_t word = m_packed[block] >> shift;
        if (shift != 0) word |= m_packed[block + 1] << (64 - shift);
        return word;
    }

    /* The 64-bit window of packed query characters starting at character
       `char_offset` of the segment. */
    uint64_t query_word(segment const& seg, const uint64_t char_offset) const {
        return packed_word(64 * seg.packed_word_begin + bits_per_char * char_offset);
    }

    /* The kmer at slot `slot` of the segment, from the packed characters. */
    kmer_t slot_kmer(segment const& seg, const uint64_t slot) const {
        static_assert(kmer_t::uint_kmer_bits % 64 == 0);
        const uint64_t bit = 64 * seg.packed_word_begin +  //
                             bits_per_char * (slot - seg.slot_begin);
        kmer_t kmer = 0;
        for (int i = kmer_t::uint_kmer_bits - 64; i >= 0; i -= 64) {
            kmer.append64(packed_word(bit + i));
        }
        kmer.take(bits_per_char * m_k);
        assert(kmer == util::string_to_uint_kmer<kmer_t>(seg.chars + (slot - seg.slot_begin), m_k));
        return kmer;
    }

    kmer_t slot_kmer_rc(segment const& seg, const uint64_t slot) const {
        kmer_t rc = slot_kmer(seg, slot);
        rc.reverse_complement_inplace(m_k);
        return rc;
    }

    /* ---------------- chains ---------------- */

    /* The memo of a resolved minimizer occurrence: whether the minimizer is
       present, the result of the kmer the walk resolved, and the decoded
       locate set of the bucket when it fits the cache. It answers every
       kmer of [head_slot, extent) that extension does not cover. */
    struct occurrence {
        uint32_t head_slot, extent;
        uint32_t head_pik;  // pos_in_kmer at head_slot; it decreases by one
                            // per slide within the occurrence
        uint64_t minimizer;
        bool found;
        lookup_result head_res;
        bucket_cache cache;
    };

    struct chain {
        uint64_t cursor, end;
        uint32_t segment;
        /* match state */
        bool anchored;
        bool backward;
        uint64_t t_cursor;     // char offset in the text of the next refill
        uint64_t t_remaining;  // extension steps left in the string
        lookup_result cur;
        /* scan state: the incremental minimizer iterator has consumed the
           slots [_, scan_pos), and `frontier` is the occurrence identity
           (minimizer, pos_in_seq) of slot scan_pos - 1 */
        bool scan_valid;
        uint64_t scan_pos;
        minimizer_info frontier;
        /* the head the last scan found but did not submit */
        bool stash_valid;
        uint64_t stash_slot;
        uint32_t stash_pik;
        uint64_t stash_minimizer;
        /* speculation */
        uint32_t spec;
        bool anchored_since_submit;
        /* the memo and the prepaid requests of the current pass */
        bool memo_valid;
        occurrence memo;
        uint32_t req_next, req_end;  // indices into the pass's request buffer
    };
    std::vector<chain> m_chains;
    std::vector<minimizer_iterator<kmer_t>> m_chain_its;  // parallel to m_chains
    std::vector<uint32_t> m_pending[2];                   // chain ids, this and next pass

    struct request_buffer {
        std::vector<uint32_t> chain;
        std::vector<uint32_t> slot;    // the occurrence's head: the kmer the walk resolves
        std::vector<uint32_t> extent;  // one past the last slot of the occurrence
        std::vector<uint32_t> pik;     // pos_in_kmer at the head
        std::vector<uint64_t> minimizer;
        std::vector<uint8_t> found;
        std::vector<lookup_result> res;
        std::vector<bucket_cache> cache;
        void clear() {
            chain.clear();
            slot.clear();
            extent.clear();
            pik.clear();
            minimizer.clear();
            found.clear();
            res.clear();
            cache.clear();
        }
        uint64_t size() const { return slot.size(); }
    };
    request_buffer m_req[2];                   // this and next pass
    std::vector<uint64_t> m_ids, m_codes;      // resolve scratch
    std::vector<uint32_t> m_decoded;           // request -> decoded bucket size

    /* performance counts */
    uint64_t m_num_searches;
    uint64_t m_num_extensions;
    uint64_t m_num_invalid;
    uint64_t m_num_negative;
    uint64_t m_num_skipped_singleton_lookups;
    uint64_t m_num_bucket_cache_hits;

    /* ---------------- front end: pack and cover ---------------- */

    void front_end() {
        m_packed.clear();
        m_segments.clear();
        m_chains.clear();
        m_chain_its.clear();

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
        const uint64_t slot_begin =
            m_segments.empty() ? 0 : m_segments.back().slot_end;
        m_segments.push_back({slot_begin, slot_begin + num_kmers, m_packed.size(), chars});
        const uint32_t segment_id = m_segments.size() - 1;

        /* pack the characters, 64 / bits_per_char per word */
        uint64_t word = 0;
        uint64_t fill = 0;
        for (uint64_t i = 0; i != length; ++i) {
            word |= kmer_t::char_to_uint(chars[i]) << (bits_per_char * fill);
            if (++fill == chars_per_word) {
                m_packed.push_back(word);
                word = 0;
                fill = 0;
            }
        }
        if (fill != 0) m_packed.push_back(word);
        /* padding, so any 64-bit window up to the last character is safe,
           whatever the kmer width */
        for (int i = 0; i != 1 + int(kmer_t::uint_kmer_bits / 64); ++i) m_packed.push_back(0);

        for (uint64_t begin = slot_begin; begin < slot_begin + num_kmers; begin += chain_length) {
            chain c;
            c.cursor = begin;
            c.end = std::min(begin + chain_length, slot_begin + num_kmers);
            c.segment = segment_id;
            c.anchored = false;
            c.scan_valid = false;
            c.stash_valid = false;
            c.spec = 1;
            c.anchored_since_submit = false;
            c.memo_valid = false;
            c.req_next = c.req_end = 0;
            m_chains.push_back(c);
            m_chain_its.push_back(m_iterator_proto);
        }
    }

    /* ---------------- pass engine ---------------- */

    void run_passes() {
        uint32_t parity = 0;
        m_req[0].clear();
        m_req[1].clear();
        m_pending[0].clear();
        m_pending[1].clear();

        for (uint32_t c = 0; c != m_chains.size(); ++c) scan_and_submit(c, parity);

        while (!m_pending[parity].empty()) {
            resolve_requests(m_req[parity]);
            const uint32_t next = parity ^ 1;
            m_req[next].clear();
            m_pending[next].clear();
            for (uint32_t c : m_pending[parity]) run_chain(c, parity);
            parity = next;
        }
    }

    /* One step of the scan: hash the minimizer of slot `scan_pos` with the
       chain's incremental iterator. Returns its minimizer info. */
    minimizer_info scan_step(chain& c, minimizer_iterator<kmer_t>& it, segment const& seg) {
        assert(c.scan_valid and c.scan_pos < c.end);
        auto info = it.next(slot_kmer(seg, c.scan_pos));
        c.scan_pos += 1;
        return info;
    }

    /* The chain cannot resolve the kmer at its cursor: scan forward for up
       to `spec` minimizer occurrences and request their walks. */
    void scan_and_submit(const uint32_t chain_id, const uint32_t parity) {
        chain& c = m_chains[chain_id];
        auto& it = m_chain_its[chain_id];
        auto const& seg = m_segments[c.segment];
        assert(c.cursor < c.end);

        c.spec = c.anchored_since_submit ? 1 : std::min<uint32_t>(c.spec * 4, max_speculation);
        c.anchored_since_submit = false;

        auto& buf = m_req[parity];
        c.req_next = buf.size();

        uint64_t head_slot;
        uint32_t head_pik;
        uint64_t head_minimizer;

        if (c.stash_valid and c.scan_valid and c.stash_slot == c.cursor) {
            /* the last scan already found the occurrence at the cursor */
            head_slot = c.stash_slot;
            head_pik = c.stash_pik;
            head_minimizer = c.stash_minimizer;
        } else {
            /* (re)start scanning at the cursor */
            it.set_position(c.cursor);
            c.scan_pos = c.cursor;
            c.scan_valid = true;
            auto info = scan_step(c, it, seg);
            c.frontier = info;
            head_slot = c.cursor;
            head_pik = info.pos_in_kmer;
            head_minimizer = info.minimizer;
        }
        c.stash_valid = false;

        uint32_t count = 0;
        while (true) {
            /* scan to the end of the current occurrence */
            uint64_t extent = c.end;
            bool found_next = false;
            minimizer_info next_info;
            while (c.scan_pos != c.end) {
                auto info = scan_step(c, it, seg);
                if (info.minimizer != c.frontier.minimizer or
                    info.pos_in_seq != c.frontier.pos_in_seq) {
                    extent = c.scan_pos - 1;
                    found_next = true;
                    next_info = info;
                    break;
                }
                c.frontier = info;
            }
            buf.chain.push_back(chain_id);
            buf.slot.push_back(head_slot);
            buf.extent.push_back(extent);
            buf.pik.push_back(head_pik);
            buf.minimizer.push_back(head_minimizer);
            count += 1;
            if (!found_next) break;  // the chain's range is fully scanned
            c.frontier = next_info;
            if (count == c.spec) {
                /* stash the head we found but do not submit */
                c.stash_valid = true;
                c.stash_slot = c.scan_pos - 1;
                c.stash_pik = next_info.pos_in_kmer;
                c.stash_minimizer = next_info.minimizer;
                break;
            }
            head_slot = c.scan_pos - 1;
            head_pik = next_info.pos_in_kmer;
            head_minimizer = next_info.minimizer;
        }

        buf.found.resize(buf.size());
        buf.res.resize(buf.size());
        buf.cache.resize(buf.size());
        c.req_end = buf.size();
        m_pending[parity].push_back(chain_id);
    }

    /* Resolve all pending walks, stage by stage, so that the memory
       accesses of the pass overlap across requests. */
    void resolve_requests(request_buffer& buf) {
        const uint64_t n = buf.size();
        if (n == 0) return;

        /* stage 1: the minimizers MPHF, batched: it hashes all the keys and
           prefetches the seeds of their buckets before reading any */
        m_ids.resize(n);
        m_codes.resize(n);
        auto const& mphf = m_dict->m_ssi.codewords.mphf;
        uint64_t r = 0;
        for (; r + mphf_batch_size <= n; r += mphf_batch_size) {
            mphf.template lookup_batch<mphf_batch_size>(&buf.minimizer[r], &m_ids[r]);
        }
        for (; r != n; ++r) m_ids[r] = mphf(buf.minimizer[r]);

        /* stage 2: the control codewords: prefetch every line, then read */
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
                const uint64_t pos_in_kmer = buf.pik[i];
                const uint64_t char_offset =
                    minimizer_offset >= pos_in_kmer ? minimizer_offset - pos_in_kmer : 0;
                __builtin_prefetch(&strings.data()[(bits_per_char * char_offset) >> 6]);
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

        /* stage 4: decode the locate sets of the cacheable buckets and
           prefetch the text around every candidate locus, so the text
           accesses of the pass overlap before any verification reads them */
        m_decoded.resize(n);
        for (uint64_t i = 0; i != n; ++i) {
            m_decoded[i] = 0;
            const uint64_t code = m_codes[i];
            buf.cache[i].size = 0;
            if ((code & 3) == bucket_t::HEAVYLOAD and (code & 1) != bucket_t::SINGLETON) continue;
            auto it = make_bucket_iterator(code);
            const uint64_t size = it.size();
            if (size > bucket_cache::max_cached_size) continue;
            auto& positions = buf.cache[i].positions;
            for (uint64_t t = 0; t != size; ++t, ++it) {
                positions[t] = m_dict->m_spss.strings_offsets.decode(*it);
            }
            m_decoded[i] = size;
            const uint64_t pos_in_kmer = buf.pik[i];
            __builtin_prefetch(
                &strings.data()[(bits_per_char * positions[0].absolute_offset) >> 6]);
            for (uint64_t t = 0; t != size; ++t) {
                const uint64_t char_offset = positions[t].absolute_offset >= pos_in_kmer
                                                 ? positions[t].absolute_offset - pos_in_kmer
                                                 : 0;
                __builtin_prefetch(&strings.data()[(bits_per_char * char_offset) >> 6]);
            }
        }

        /* stage 5: verification */
        for (uint64_t i = 0; i != n; ++i) {
            auto const& seg = m_segments[m_chains[buf.chain[i]].segment];
            const uint64_t slot = buf.slot[i];
            const uint64_t code = m_codes[i];
            const minimizer_info info(buf.minimizer[i], buf.pik[i]);
            lookup_result res;
            if ((code & 3) == bucket_t::HEAVYLOAD and (code & 1) != bucket_t::SINGLETON) {
                /* skew-index bucket: the walk is kmer-specific; rare */
                res = m_dict->lookup(slot_kmer(seg, slot), slot_kmer_rc(seg, slot), info);
            } else if (const uint64_t size = m_decoded[i]; size != 0) {
                /* the decoded locate set is at hand: check the minimizer at
                   its anchor, then verify, as spss::lookup does */
                auto const& positions = buf.cache[i].positions;
                if (uint64_t read_mmer = uint64_t(util::read_kmer_at<kmer_t>(
                        strings, m_m, bits_per_char * positions[0].absolute_offset));
                    read_mmer != info.minimizer and
                    util::canonical_mmer<kmer_t>(read_mmer, m_m) != info.minimizer) {
                    res = lookup_result(false);
                } else {
                    buf.cache[i].size = size;
                    res = m_dict->m_spss.lookup_from_positions(
                        positions.data(), size, slot_kmer(seg, slot), slot_kmer_rc(seg, slot),
                        info);
                }
            } else {
                /* bucket larger than the cache */
                auto it = make_bucket_iterator(code);
                res = m_dict->m_spss.lookup(it, slot_kmer(seg, slot), slot_kmer_rc(seg, slot),
                                            info, &buf.cache[i]);
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
    void run_chain(const uint32_t chain_id, const uint32_t parity) {
        chain& c = m_chains[chain_id];
        auto const& seg = m_segments[c.segment];
        auto const& buf = m_req[parity];

        while (c.cursor != c.end) {
            /* 1. word-parallel extension: crosses minimizer changes, costs
                  no walk and no hashing */
            if (c.anchored) {
                lcp_sprint(c, seg);
                if (c.cursor == c.end) return;
                c.anchored = false;
            }

            const uint64_t j = c.cursor;

            /* 2. find the resolved occurrence covering the cursor */
            if (!(c.memo_valid and j >= c.memo.head_slot and j < c.memo.extent)) {
                while (c.req_next != c.req_end and buf.extent[c.req_next] <= j) c.req_next += 1;
                if (c.req_next != c.req_end and buf.slot[c.req_next] <= j) {
                    const uint32_t i = c.req_next;
                    c.memo_valid = true;
                    c.memo.head_slot = buf.slot[i];
                    c.memo.extent = buf.extent[i];
                    c.memo.head_pik = buf.pik[i];
                    c.memo.minimizer = buf.minimizer[i];
                    c.memo.found = buf.found[i];
                    c.memo.head_res = buf.res[i];
                    c.memo.cache = buf.cache[i];
                    c.req_next += 1;
                } else {
                    scan_and_submit(chain_id, parity ^ 1);
                    return;
                }
            }

            /* 3. resolve the kmer from the memo */
            lookup_result res;
            if (j == c.memo.head_slot) {
                res = c.memo.head_res;
            } else if (!c.memo.found) {
                /* absent minimizer: definitive (never set for HEAVYLOAD) */
                res = lookup_result(false);
            } else {
                const minimizer_info info(c.memo.minimizer,
                                          c.memo.head_pik - (j - c.memo.head_slot));
                if (c.memo.cache.size != 0) {
                    m_num_bucket_cache_hits += 1;
                    res = m_dict->m_spss.lookup_from_positions(c.memo.cache.positions.data(),
                                                               c.memo.cache.size,
                                                               slot_kmer(seg, j),
                                                               slot_kmer_rc(seg, j), info);
                } else {
                    /* heavy or over-large bucket: per-kmer walk; rare */
                    res = m_dict->lookup(slot_kmer(seg, j), slot_kmer_rc(seg, j), info);
                }
            }

            if (res.kmer_id != constants::invalid_uint64) {
                m_num_searches += 1;
                c.cur = res;
                anchor(c, res);
                c.anchored_since_submit = true;
            } else {
                m_num_negative += 1;
            }
            assert(equal_lookup_result(m_dict->lookup(seg.chars + (j - seg.slot_begin)), res));
            c.cursor += 1;
        }
    }

    void anchor(chain& c, lookup_result const& res) const {
        assert(res.kmer_id != constants::invalid_uint64);
        const uint64_t kmer_offset = res.kmer_id + res.string_id * (m_k - 1);
        c.anchored = true;
        if (res.kmer_orientation == constants::backward_orientation) {
            c.backward = true;
            c.t_remaining = res.kmer_id_in_string;
            c.t_cursor = kmer_offset;
        } else {
            c.backward = false;
            c.t_remaining = (res.string_end - res.string_begin - m_k) - res.kmer_id_in_string;
            c.t_cursor = kmer_offset + m_k;
        }
    }

    /* Extend the chain's match as far as the text agrees with the query,
       one word -- 64 / bits_per_char characters -- per compare. */
    void lcp_sprint(chain& c, segment const& seg) {
        assert(c.anchored);
        auto const& strings = m_dict->m_spss.strings;
        const uint64_t limit = std::min(c.end - c.cursor, c.t_remaining);
        const uint64_t q_char = (c.cursor - seg.slot_begin) + m_k - 1;
        uint64_t total = 0;
        while (total != limit) {
            const uint64_t qw = query_word(seg, q_char + total);
            uint64_t tw;
            if (!c.backward) {
                assert(bits_per_char * (c.t_cursor + total) < strings.num_bits());
                tw = strings.get_word64(bits_per_char * (c.t_cursor + total));
            } else {
                if constexpr (kmer_t::has_reverse_complement) {
                    static_assert(bits_per_char == 2);
                    /* the word of text ending at the cursor, reverse
                       complemented: its code i is the complement of the
                       character at offset cursor - 1 - i, the character the
                       i-th backward extension from here appends */
                    const uint64_t end = c.t_cursor - total;
                    uint64_t word;
                    if (end >= chars_per_word) {
                        word = strings.get_word64(bits_per_char * (end - chars_per_word));
                    } else {
                        word = strings.get_word64(0) << (bits_per_char * (chars_per_word - end));
                    }
                    tw = kmer_t::crc64(word);
                } else {
                    assert(false);  // unreachable: no orientations without complements
                    tw = 0;
                }
            }
            const uint64_t n = std::min(chars_per_word, limit - total);
            const uint64_t mask =
                n == chars_per_word ? ~uint64_t(0) : (uint64_t(1) << (bits_per_char * n)) - 1;
            const uint64_t diff = (qw ^ tw) & mask;
            if (diff != 0) {
                total += __builtin_ctzll(diff) / bits_per_char;
                break;
            }
            total += n;
        }

#ifndef NDEBUG
        /* re-verify every extended kmer one by one */
        {
            auto cur = c.cur;
            for (uint64_t t = 0; t != total; ++t) {
                cur.kmer_id += cur.kmer_orientation;
                cur.kmer_id_in_string += cur.kmer_orientation;
                assert(equal_lookup_result(
                    m_dict->lookup(seg.chars + (c.cursor - seg.slot_begin) + t), cur));
            }
        }
#endif

        m_num_extensions += total;
        c.cur.kmer_id += c.cur.kmer_orientation * int64_t(total);
        c.cur.kmer_id_in_string += c.cur.kmer_orientation * int64_t(total);
        c.cursor += total;
        c.t_remaining -= total;
        if (!c.backward) {
            c.t_cursor += total;
        } else {
            c.t_cursor -= total;
        }
    }
};

}  // namespace sshash
