#pragma once

#include <cstring>
#include <cstdlib>
#if defined(__x86_64__)
#include <x86intrin.h>
#endif

#include "include/dictionary.hpp"
#include "include/minimizer_iterator.hpp"
#include "include/util.hpp"

namespace sshash {

/*
    Batched streaming queries: a small batch of sequences advances in
    lockstep, and every dependent memory access of every sequence is
    prefetched one step before it is read.

    A single streaming lookup is a chain of dependent, essentially random
    memory accesses (the minimizers MPHF, a control codeword, the text),
    each costing a cache miss that leaves the rest of the memory system
    idle: one query alone cannot use the dozen or so misses a core can keep
    in flight. The engine therefore holds `num_slots` pieces of sequence --
    one piece is a whole sequence, or a chunk of `piece_kmers` kmers of a
    longer one -- and sweeps over them, advancing each piece by one bounded
    step of a state machine per sweep. Each step consumes memory that was
    prefetched when the previous step of the same piece issued it, roughly
    one full sweep earlier, and issues the prefetch for the next step. A
    finished piece is replaced from the input, so the batch stays full.

    The per-piece state machine keeps the structure that makes streaming
    SSHash fast in the first place:

    - EXTEND: a held match is extended with a word-parallel comparison of
      the 2-bit packed query against the packed text (XOR and
      count-trailing-zeros, 64 / bits_per_char characters per step; for
      backward matches the text word is reverse-complemented first, one
      crc64), crossing minimizer changes freely, with the next text word
      prefetched each step. A run of consecutive positive kmers costs
      O(run / word size) and no hashing at all.

    - SCAN: where a match stops, the minimizer of the current window is
      maintained incrementally (one hash per base, only here), and kmers of
      an already-resolved minimizer occurrence are answered from its memo:
      an absent minimizer answers them at no cost, a cached locate set by
      direct verification. A kmer of an unresolved occurrence files a walk
      request and the piece waits.

    - Walk requests are collected across the batch and started in groups
      through the minimizers MPHF's own batched lookup (which hashes every
      key and prefetches the seeds of their buckets before reading any);
      the walk then proceeds one stage per sweep -- codeword (line
      prefetched at group start), bucket decode (locate set line prefetched
      by the previous stage), verification (candidate text lines prefetched
      by the previous stage), and resolution of the matched locus to its
      string -- so the misses of concurrent walks overlap.

    Sequences are enqueued (`enqueue_sequence`) into an arena and processed
    in blocks of about `flush_kmers` kmers (`flush`; automatic when the
    arena fills); each block is cut at invalid characters into 2-bit-packed
    segments (kmers spanning an invalid character are counted invalid), and
    the pieces of all segments feed the slots. Results are identical to
    `dictionary::lookup` on every kmer (debug assert). Counters, for report
    compatibility: kmers resolved by extension count as "extensions", every
    other positive as a "search"; `num_bucket_cache_hits` counts kmers
    verified against a cached locate set.
*/
template <typename Dict>
struct streaming_query_batched  //
{
    using kmer_t = typename Dict::kmer_type;
    using bucket_cache = typename Dict::spss_type::bucket_cache;

    static constexpr bool is_batched = true;
    static constexpr uint64_t flush_kmers = 1 << 15;
    static constexpr uint64_t piece_kmers = 128;
    static constexpr uint64_t num_slots = 64;
    static constexpr uint64_t walk_group_size = 16;
    static constexpr uint64_t scan_budget = 32;
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

    {
        m_slots.resize(num_slots);
        for (auto& s : m_slots) s.it = m_iterator_proto;
        m_probe.capture(dict->m_spss.strings_offsets);
        m_timings = std::getenv("SSHASH_TIMINGS") != nullptr;
    }

    ~streaming_query_batched() { print_timings(); }

    /* Sequences are independent of each other; nothing to reset. */
    void reset() {}

    /* Copies the sequence into the arena; flushes when the arena is full.
       Counters are complete only after a final `flush()`. A sequence longer
       than one block is cut into overlapping chunks (by k - 1 characters,
       so every kmer appears in exactly one chunk). */
    void enqueue_sequence(char const* chars, const uint64_t length) {
        assert(length >= m_k);
        uint64_t begin = 0;
        while (length - begin > flush_kmers + m_k - 1) {
            enqueue_chunk(chars + begin, flush_kmers + m_k - 1);
            begin += flush_kmers;
        }
        enqueue_chunk(chars + begin, length - begin);
    }

    /* One-shot: resolve the kmers of one sequence now (with whatever else
       is enqueued). */
    void lookup_sequence(char const* chars, const uint64_t length) {
        enqueue_sequence(chars, length);
        flush();
    }

    void flush() {
        if (m_views.empty()) return;
        if (m_timings) {
            const uint64_t t0 = tsc();
            front_end();
            account(1, t0);
        } else {
            front_end();
        }
        run();
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

    void enqueue_chunk(char const* chars, const uint64_t length) {
        const uint64_t t0 = m_timings ? tsc() : 0;
        m_views.emplace_back(m_arena.size(), length);
        m_arena.append(chars, length);
        m_pending_kmers += length - m_k + 1;
        if (m_timings) account(0, t0);
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

    struct piece {
        uint32_t seg;
        uint64_t begin, end;  // kmers (block-global)
    };
    std::vector<piece> m_pieces;
    uint64_t m_next_piece;

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

    /*
        Addresses inside the string-endpoints sequence, captured through its
        own public `visit` traversal, so that the dependent reads of
        `offset_to_id` -- the hint, then the entry lines of the search --
        can be prefetched one sweep before they are read. These are pure
        prefetch hints: if the traversal shape ever changes and the capture
        comes back incomplete, the engine runs identically without them.
    */
    struct endpoints_probe {
        bool valid = false;
        char const* hints_bytes = nullptr;
        uint64_t hints_width = 0;
        uint64_t hints_size = 0;
        uint64_t const* high_words = nullptr;
        uint8_t const* low_bytes = nullptr;

        template <typename T, typename V, typename = void>
        struct has_visit : std::false_type {};
        template <typename T, typename V>
        struct has_visit<T, V,
                         std::void_t<decltype(std::declval<T const&>().visit(std::declval<V&>()))>>
            : std::true_type {};

        template <typename T>
        void visit(T const& t) {
            if constexpr (std::is_same_v<T, bits::compact_vector>) {
                if (hints_bytes == nullptr) {
                    hints_bytes = reinterpret_cast<char const*>(t.data().data());
                    hints_width = t.width();
                    hints_size = t.size();
                }
            } else if constexpr (std::is_same_v<T, bits::bit_vector>) {
                if (high_words == nullptr) high_words = t.data().data();
            } else if constexpr (std::is_same_v<T, essentials::owning_span<uint8_t>>) {
                if (low_bytes == nullptr) low_bytes = t.data();
            } else if constexpr (std::is_arithmetic_v<T> or std::is_enum_v<T>) {
                /* nothing */
            } else if constexpr (has_visit<T, endpoints_probe>::value) {
                t.visit(*this);
            }
        }

        template <typename Offsets>
        void capture(Offsets const& offsets) {
            offsets.visit(*this);
            valid = hints_bytes != nullptr and high_words != nullptr and
                    low_bytes != nullptr and hints_width > 0 and hints_width < 57;
        }
    };
    endpoints_probe m_probe;

    /* ---------------- timings (diagnostic) ----------------

       With the environment variable SSHASH_TIMINGS set, every step of the
       sweep loop is timed with the TSC and attributed to its state, and a
       share table is printed when the query object is destroyed. Off by
       default: a run with timings pays two TSC reads per step, which
       distorts the shares by a few percent but leaves their order intact. */

    static constexpr int num_timing_bins = 15;
    static constexpr char const* timing_names[num_timing_bins] = {
        "enqueue (arena copy)", "front end (validity + pack)", "scan: other",
        "extend (word LCP)",    "walk: start (MPHF group)",    "walk: codeword other",
        "walk: bucket decode",  "walk: verify",                "walk: endpoints hint",
        "walk: locate",         "sweep overhead",              "scan: iterator sync + hash",
        "scan: memo resolve",   "walk: codeword access",       "scan: heavy member walk"};
    bool m_timings = false;
    uint64_t m_cycles[num_timing_bins] = {};
    uint64_t m_steps[num_timing_bins] = {};

    static uint64_t tsc() {
#if defined(__x86_64__)
        return __rdtsc();
#else
        return 0;
#endif
    }

    void account(const int bin, const uint64_t t0) {
        m_cycles[bin] += tsc() - t0;
        m_steps[bin] += 1;
    }

    void print_timings() const {
        uint64_t cycles[num_timing_bins];
        for (int b = 0; b != num_timing_bins; ++b) cycles[b] = m_cycles[b];
        /* the sub-bins are nested inside their parents: subtract them */
        cycles[2] = cycles[2] > cycles[11] + cycles[12] ? cycles[2] - cycles[11] - cycles[12] : 0;
        cycles[12] = cycles[12] > cycles[14] ? cycles[12] - cycles[14] : 0;
        cycles[5] = cycles[5] > cycles[13] ? cycles[5] - cycles[13] : 0;
        /* the scheduler overhead is the run() time not attributed to steps */
        uint64_t inside_run = 0;
        for (int b = 2; b <= 9; ++b) inside_run += cycles[b];
        inside_run += cycles[11] + cycles[12] + cycles[13] + cycles[14];
        cycles[10] = cycles[10] > inside_run ? cycles[10] - inside_run : 0;
        uint64_t total = 0;
        for (int b = 0; b != num_timing_bins; ++b) total += cycles[b];
        if (total == 0) return;
        std::cout << "==== streaming_query_batched timings (cycles, shares):\n";
        for (int b = 0; b != num_timing_bins; ++b) {
            std::cout << "  " << timing_names[b] << " = " << cycles[b] << " ("
                      << (100.0 * cycles[b]) / total << "%), steps = " << m_steps[b] << '\n';
        }
        std::cout << "  total accounted cycles = " << total << std::endl;
    }

    /* ---------------- the state machine ---------------- */

    enum class state : uint8_t {
        idle,         // no piece
        scanning,     // maintaining the minimizer, resolving from the memo
        extending,    // holding a match, word-parallel LCP
        walk_wait,    // a walk request is filed, waiting for a group
        walk_pf,      // prefetch the control codeword line (one per step,
                      // so prefetches never burst against full fill buffers)
        walk_code,    // read the control codeword (line prefetched)
        walk_bucket,  // decode the locate set (line prefetched)
        walk_verify,  // compare the candidate loci (text prefetched)
        walk_hint,    // read the endpoints hint (line prefetched) and
                      // prefetch the entry lines of the endpoints search
        walk_locate,  // resolve the matched locus to its string
    };

    struct slot_t {
        state st = state::idle;
        uint32_t seg = 0;
        uint64_t cursor = 0, end = 0;

        /* extension */
        bool backward = false;
        uint64_t t_cursor = 0;     // char offset in the text of the next compare
        uint64_t t_remaining = 0;  // extension steps left in the string
        lookup_result cur;

        /* scan: the incremental iterator has consumed kmers [_, scan_pos),
           and `info` is the minimizer info of kmer scan_pos - 1 */
        minimizer_iterator<kmer_t> it;
        bool scan_valid = false;
        uint64_t scan_pos = 0;
        minimizer_info info;

        /* the memo of the last resolved minimizer occurrence */
        bool memo_valid = false;
        uint64_t memo_minimizer = 0;
        uint64_t memo_pos_in_seq = 0;
        uint64_t memo_head = 0;  // the kmer its walk resolved
        bool memo_found = false;
        lookup_result memo_head_res;
        bucket_cache cache;
        /* mid-load buckets larger than the cache: decoded in full here, so
           their kmers verify against positions instead of walking */
        static constexpr uint64_t max_bucket_size = (uint64_t(1) << constants::min_l) + 1;
        std::array<
            typename std::decay_t<decltype(std::declval<bucket_cache>().positions)>::value_type,
            max_bucket_size>
            big_positions;
        uint64_t big_size = 0;

        /* the walk in progress */
        uint64_t walk_id = 0;    // minimizers MPHF value
        uint64_t walk_code = 0;  // control codeword
        uint64_t walk_size = 0;  // its bucket size
        uint64_t match_offset = 0;
        bool match_backward = false;
    };
    std::vector<slot_t> m_slots;
    std::vector<uint32_t> m_walk_queue;  // slots in walk_wait
    std::vector<uint64_t> m_group_keys, m_group_ids;  // walk group scratch

    /* performance counts */
    uint64_t m_num_searches;
    uint64_t m_num_extensions;
    uint64_t m_num_invalid;
    uint64_t m_num_negative;
    uint64_t m_num_skipped_singleton_lookups;
    uint64_t m_num_bucket_cache_hits;

    /* ---------------- front end: pack and cut ---------------- */

    void front_end() {
        m_packed.clear();
        m_segments.clear();
        m_pieces.clear();
        m_next_piece = 0;

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
        const uint64_t slot_begin = m_segments.empty() ? 0 : m_segments.back().slot_end;
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

        for (uint64_t begin = slot_begin; begin < slot_begin + num_kmers; begin += piece_kmers) {
            m_pieces.push_back(
                {segment_id, begin, std::min(begin + piece_kmers, slot_begin + num_kmers)});
        }
    }

    /* ---------------- the sweep scheduler ---------------- */

    void run() {
        m_walk_queue.clear();
        uint64_t num_active = 0;
        for (auto& s : m_slots) {
            s.st = state::idle;
            if (take_piece(s)) num_active += 1;
        }

        const uint64_t t_run = m_timings ? tsc() : 0;
        while (num_active != 0) {
            bool advanced = false;
            for (auto& s : m_slots) {
                const uint64_t t0 = m_timings ? tsc() : 0;
                int bin;
                switch (s.st) {
                    case state::idle:
                    case state::walk_wait: continue;
                    case state::scanning: step_scan(s); bin = 2; break;
                    case state::extending: step_extend(s); bin = 3; break;
                    case state::walk_pf: step_walk_pf(s); bin = 5; break;
                    case state::walk_code: step_walk_code(s); bin = 5; break;
                    case state::walk_bucket: step_walk_bucket(s); bin = 6; break;
                    case state::walk_verify: step_walk_verify(s); bin = 7; break;
                    case state::walk_hint: step_walk_hint(s); bin = 8; break;
                    case state::walk_locate: step_walk_locate(s); bin = 9; break;
                    default: bin = 10; break;
                }
                if (m_timings) account(bin, t0);
                advanced = true;
                if (s.st == state::idle) {
                    if (!take_piece(s)) num_active -= 1;
                }
            }
            /* start a group of walks: when enough are waiting, or when no
               slot can advance without them */
            if (m_walk_queue.size() >= walk_group_size or (!advanced and !m_walk_queue.empty())) {
                if (m_timings) {
                    const uint64_t t0 = tsc();
                    start_walks();
                    account(4, t0);
                } else {
                    start_walks();
                }
            }
        }
        assert(m_walk_queue.empty());
        /* bin 10 accumulates the raw run() time; print_timings subtracts
           the attributed bins from it to report the scheduler overhead */
        if (m_timings) m_cycles[10] += tsc() - t_run;
    }

    bool take_piece(slot_t& s) {
        assert(s.st == state::idle);
        if (m_next_piece == m_pieces.size()) return false;
        auto const& p = m_pieces[m_next_piece++];
        s.seg = p.seg;
        s.cursor = p.begin;
        s.end = p.end;
        s.st = state::scanning;
        s.scan_valid = false;
        s.memo_valid = false;
        return true;
    }

    /* Start the walks of up to `walk_group_size` waiting slots: the
       minimizers MPHF hashes all their keys and prefetches the seeds of
       their buckets before reading any; each walk then gets its control
       codeword line prefetched and proceeds stage by stage. */
    void start_walks() {
        auto const& mphf = m_dict->m_ssi.codewords.mphf;

        while (!m_walk_queue.empty()) {
            const uint64_t n = std::min<uint64_t>(walk_group_size, m_walk_queue.size());
            m_group_keys.resize(walk_group_size);
            m_group_ids.resize(walk_group_size);
            const uint64_t first = m_walk_queue.size() - n;
            for (uint64_t i = 0; i != walk_group_size; ++i) {
                /* pad short groups by repeating the first key */
                const uint64_t q = i < n ? first + i : first;
                m_group_keys[i] = m_slots[m_walk_queue[q]].memo_minimizer;
            }
            mphf.template lookup_batch<walk_group_size>(m_group_keys.data(), m_group_ids.data());
            for (uint64_t i = 0; i != n; ++i) {
                auto& s = m_slots[m_walk_queue[first + i]];
                s.walk_id = m_group_ids[i];
                s.st = state::walk_pf;
            }
            m_walk_queue.resize(first);
            if (m_walk_queue.size() < walk_group_size) break;
        }
    }

    /* ---------------- per-sweep steps ---------------- */

    segment const& seg_of(slot_t const& s) const { return m_segments[s.seg]; }

    void finish_piece(slot_t& s) { s.st = state::idle; }

    void emit_positive(slot_t& s, lookup_result const& res) {
        m_num_searches += 1;
        s.cur = res;
        anchor(s, res);
        s.st = state::extending;
        prefetch_text(s);
        assert(equal_lookup_result(
            m_dict->lookup(seg_of(s).chars + (s.cursor - seg_of(s).slot_begin)), res));
        s.cursor += 1;
    }

    void emit_negative(slot_t& s, lookup_result const& res) {
        m_num_negative += 1;
        s.st = state::scanning;
        assert(equal_lookup_result(
            m_dict->lookup(seg_of(s).chars + (s.cursor - seg_of(s).slot_begin)), res));
        (void)res;
        s.cursor += 1;
    }

    /* SCAN: maintain the minimizer incrementally and resolve kmers of the
       memo'd occurrence; file a walk for the first kmer of any other. */
    void step_scan(slot_t& s) {
        auto const& seg = seg_of(s);
        for (uint64_t budget = scan_budget; budget != 0; --budget) {
            if (s.cursor == s.end) {
                finish_piece(s);
                return;
            }
            /* sync the iterator with the cursor */
            {
                const uint64_t t0 = m_timings ? tsc() : 0;
                if (s.scan_valid and s.scan_pos == s.cursor + 1) {
                    /* s.info is current */
                } else if (s.scan_valid and s.scan_pos == s.cursor) {
                    s.info = s.it.next(slot_kmer(seg, s.cursor));
                    s.scan_pos += 1;
                } else {
                    s.it.set_position(s.cursor);
                    s.info = s.it.next(slot_kmer(seg, s.cursor));
                    s.scan_pos = s.cursor + 1;
                    s.scan_valid = true;
                }
                if (m_timings) account(11, t0);
            }

            if (!(s.memo_valid and s.info.minimizer == s.memo_minimizer and
                  s.info.pos_in_seq == s.memo_pos_in_seq)) {
                /* an unresolved occurrence: file a walk request. The head
                   of the walk is the cursor kmer itself. */
                s.memo_valid = true;
                s.memo_minimizer = s.info.minimizer;
                s.memo_pos_in_seq = s.info.pos_in_seq;
                s.memo_head = s.cursor;
                s.st = state::walk_wait;
                m_walk_queue.push_back(&s - m_slots.data());
                return;
            }

            /* resolve from the memo */
            const uint64_t t1 = m_timings ? tsc() : 0;
            lookup_result res;
            if (!s.memo_found) {
                res = lookup_result(false);
            } else if (s.cache.size != 0) {
                m_num_bucket_cache_hits += 1;
                res = m_dict->m_spss.lookup_from_positions(s.cache.positions.data(), s.cache.size,
                                                           slot_kmer(seg, s.cursor),
                                                           slot_kmer_rc(seg, s.cursor), s.info);
            } else if (s.big_size != 0) {
                m_num_bucket_cache_hits += 1;
                res = m_dict->m_spss.lookup_from_positions(s.big_positions.data(), s.big_size,
                                                           slot_kmer(seg, s.cursor),
                                                           slot_kmer_rc(seg, s.cursor), s.info);
            } else {
                /* heavy bucket: the walk is kmer-specific */
                const uint64_t t2 = m_timings ? tsc() : 0;
                res = m_dict->lookup(slot_kmer(seg, s.cursor), slot_kmer_rc(seg, s.cursor),
                                     s.info);
                if (m_timings) account(14, t2);
            }
            if (m_timings) account(12, t1);
            if (res.kmer_id != constants::invalid_uint64) {
                emit_positive(s, res);
                return;
            }
            emit_negative(s, res);
        }
    }

    void anchor(slot_t& s, lookup_result const& res) const {
        assert(res.kmer_id != constants::invalid_uint64);
        const uint64_t kmer_offset = res.kmer_id + res.string_id * (m_k - 1);
        if (res.kmer_orientation == constants::backward_orientation) {
            s.backward = true;
            s.t_remaining = res.kmer_id_in_string;
            s.t_cursor = kmer_offset;
        } else {
            s.backward = false;
            s.t_remaining = (res.string_end - res.string_begin - m_k) - res.kmer_id_in_string;
            s.t_cursor = kmer_offset + m_k;
        }
    }

    void prefetch_text(slot_t const& s) const {
        if (s.t_remaining == 0) return;
        auto const& strings = m_dict->m_spss.strings;
        const uint64_t char_offset = s.backward
                                         ? (s.t_cursor >= chars_per_word  //
                                                ? s.t_cursor - chars_per_word
                                                : 0)
                                         : s.t_cursor;
        __builtin_prefetch(&strings.data()[(bits_per_char * char_offset) >> 6]);
    }

    /* EXTEND: one word-parallel comparison per sweep, with the next text
       word prefetched for the following sweep. */
    void step_extend(slot_t& s) {
        if (s.cursor == s.end) {
            finish_piece(s);
            return;
        }
        if (s.t_remaining == 0) {
            s.st = state::scanning;
            return;
        }
        auto const& seg = seg_of(s);
        auto const& strings = m_dict->m_spss.strings;

        const uint64_t limit =
            std::min({chars_per_word, s.end - s.cursor, s.t_remaining});
        const uint64_t qw = query_word(seg, (s.cursor - seg.slot_begin) + m_k - 1);
        uint64_t tw;
        if (!s.backward) {
            assert(bits_per_char * s.t_cursor < strings.num_bits());
            tw = strings.get_word64(bits_per_char * s.t_cursor);
        } else {
            if constexpr (kmer_t::has_reverse_complement) {
                static_assert(bits_per_char == 2);
                /* the word of text ending at the cursor, reverse
                   complemented: its code i is the complement of the
                   character at offset t_cursor - 1 - i, the character the
                   i-th backward extension from here appends */
                uint64_t word;
                if (s.t_cursor >= chars_per_word) {
                    word = strings.get_word64(bits_per_char * (s.t_cursor - chars_per_word));
                } else {
                    word = strings.get_word64(0)
                           << (bits_per_char * (chars_per_word - s.t_cursor));
                }
                tw = kmer_t::crc64(word);
            } else {
                assert(false);  // unreachable: no orientations without complements
                tw = 0;
            }
        }
        const uint64_t mask =
            limit == chars_per_word ? ~uint64_t(0) : (uint64_t(1) << (bits_per_char * limit)) - 1;
        const uint64_t diff = (qw ^ tw) & mask;
        const uint64_t matched = diff != 0 ? __builtin_ctzll(diff) / bits_per_char : limit;

#ifndef NDEBUG
        /* re-verify every extended kmer one by one */
        {
            auto cur = s.cur;
            for (uint64_t t = 0; t != matched; ++t) {
                cur.kmer_id += cur.kmer_orientation;
                cur.kmer_id_in_string += cur.kmer_orientation;
                assert(equal_lookup_result(
                    m_dict->lookup(seg.chars + (s.cursor - seg.slot_begin) + t), cur));
            }
        }
#endif

        m_num_extensions += matched;
        s.cur.kmer_id += s.cur.kmer_orientation * int64_t(matched);
        s.cur.kmer_id_in_string += s.cur.kmer_orientation * int64_t(matched);
        s.cursor += matched;
        s.t_remaining -= matched;
        if (!s.backward) {
            s.t_cursor += matched;
        } else {
            s.t_cursor -= matched;
        }

        if (diff != 0 and matched != limit) {
            /* a real mismatch within the word */
            s.st = state::scanning;
            return;
        }
        prefetch_text(s);
    }

    /* WALK, stage 0: prefetch the control codeword line. One prefetch per
       slot step: issued in a burst at group start they contend with the
       MPHF's own reads for line-fill buffers and get dropped. */
    void step_walk_pf(slot_t& s) {
        auto const& codewords = m_dict->m_ssi.codewords.control_codewords;
        auto const* bytes = reinterpret_cast<char const*>(codewords.data().data());
        __builtin_prefetch(bytes + ((s.walk_id * codewords.width()) >> 3));
        s.st = state::walk_code;
    }

    /* WALK, stage 1: read the control codeword (its line was prefetched
       by the previous step) and prefetch what the bucket kind needs. */
    void step_walk_code(slot_t& s) {
        auto const& ssi = m_dict->m_ssi;
        auto const& strings = m_dict->m_spss.strings;
        const uint64_t t0 = m_timings ? tsc() : 0;
        const uint64_t code = ssi.codewords.control_codewords.access(s.walk_id);
        if (m_timings) account(13, t0);
        s.walk_code = code;
        s.cache.size = 0;
        s.big_size = 0;

        if ((code & 1) == bucket_t::SINGLETON) {
            /* the locate set is the codeword itself: decode it now and
               prefetch the candidate text */
            s.cache.positions[0] = m_dict->m_spss.strings_offsets.decode(code >> 1);
            const uint64_t pos_in_kmer = s.info.pos_in_kmer;
            const uint64_t offset = s.cache.positions[0].absolute_offset;
            const uint64_t char_offset = offset >= pos_in_kmer ? offset - pos_in_kmer : 0;
            __builtin_prefetch(&strings.data()[(bits_per_char * char_offset) >> 6]);
            __builtin_prefetch(&strings.data()[(bits_per_char * offset) >> 6]);
            s.walk_size = 1;
            s.st = state::walk_verify;
            return;
        }
        if ((code & 3) == bucket_t::MIDLOAD) {
            /* every mid-load bucket fits the slot's big_positions array:
               decode it through the pipeline, whatever its size */
            constexpr uint64_t mask = (uint64_t(1) << constants::min_l) - 1;
            const uint64_t c = code >> 2;
            const uint64_t bucket_size = (c & mask) + 2;
            const uint64_t bucket_id = c >> constants::min_l;
            const uint64_t bucket_begin =
                ssi.begin_buckets_of_size[bucket_size] + bucket_id * bucket_size;
            auto const& mid_load = ssi.mid_load_buckets;
            auto const* bytes = reinterpret_cast<char const*>(mid_load.data().data());
            const uint64_t first_bit = bucket_begin * mid_load.width();
            const uint64_t last_bit = (bucket_begin + bucket_size) * mid_load.width();
            for (uint64_t byte = first_bit >> 3; byte <= (last_bit >> 3); byte += 64) {
                __builtin_prefetch(bytes + byte);
            }
            s.st = state::walk_bucket;
            return;
        }
        /* HEAVYLOAD: the walk is kmer-specific; rare */
        auto const& seg = seg_of(s);
        lookup_result res =
            m_dict->lookup(slot_kmer(seg, s.cursor), slot_kmer_rc(seg, s.cursor), s.info);
        s.memo_found = res.minimizer_found;
        s.memo_head_res = res;
        if (res.kmer_id != constants::invalid_uint64) {
            emit_positive(s, res);
        } else {
            emit_negative(s, res);
        }
    }

    /* WALK, stage 2 (mid-load buckets): decode the locate set (its line was
       prefetched) and prefetch the candidate text. */
    void step_walk_bucket(slot_t& s) {
        auto const& strings = m_dict->m_spss.strings;
        auto it = make_bucket_iterator(s.walk_code);
        const uint64_t size = it.size();
        assert(size >= 2 and size <= slot_t::max_bucket_size);
        const bool small = size <= bucket_cache::max_cached_size;
        auto* positions = small ? s.cache.positions.data() : s.big_positions.data();
        const uint64_t pos_in_kmer = s.info.pos_in_kmer;
        for (uint64_t t = 0; t != size; ++t, ++it) {
            positions[t] = m_dict->m_spss.strings_offsets.decode(*it);
            const uint64_t offset = positions[t].absolute_offset;
            const uint64_t char_offset = offset >= pos_in_kmer ? offset - pos_in_kmer : 0;
            __builtin_prefetch(&strings.data()[(bits_per_char * char_offset) >> 6]);
        }
        __builtin_prefetch(&strings.data()[(bits_per_char * positions[0].absolute_offset) >> 6]);
        s.walk_size = size;
        s.st = state::walk_verify;
    }

    /* the decoded locate set of the walk in progress */
    auto const* walk_positions(slot_t const& s) const {
        return s.walk_size <= bucket_cache::max_cached_size ? s.cache.positions.data()
                                                            : s.big_positions.data();
    }

    /* WALK, stage 3: the minimizer presence check and the candidate
       comparisons, against prefetched text. */
    void step_walk_verify(slot_t& s) {
        auto const& seg = seg_of(s);
        auto const& strings = m_dict->m_spss.strings;
        const uint64_t size = s.walk_size;
        auto const* positions = walk_positions(s);

        /* presence check, as spss::lookup does */
        if (uint64_t read_mmer = uint64_t(util::read_kmer_at<kmer_t>(
                strings, m_m, bits_per_char * positions[0].absolute_offset));
            read_mmer != s.info.minimizer and
            util::canonical_mmer<kmer_t>(read_mmer, m_m) != s.info.minimizer) {
            s.memo_found = false;
            emit_negative(s, lookup_result(false));
            return;
        }
        s.memo_found = true;
        if (size <= bucket_cache::max_cached_size) {
            s.cache.size = size;
        } else {
            s.big_size = size;
        }

        const kmer_t kmer = slot_kmer(seg, s.cursor);
        const kmer_t kmer_rc = slot_kmer_rc(seg, s.cursor);
        uint64_t loci[2] = {s.info.pos_in_kmer, m_k - m_m - s.info.pos_in_kmer};
        const int num_loci = kmer_t::has_reverse_complement ? 2 : 1;
        for (uint64_t t = 0; t != size; ++t) {
            for (int l = 0; l != num_loci; ++l) {
                if (positions[t].absolute_offset < loci[l]) continue;
                const uint64_t offset = positions[t].absolute_offset - loci[l];
                auto read_kmer = util::read_kmer_at<kmer_t>(strings, m_k, bits_per_char * offset);
                if (read_kmer != kmer and read_kmer != kmer_rc) continue;
                s.match_offset = offset;
                s.match_backward = read_kmer == kmer_rc;
                if (m_probe.valid and (offset >> 8) != 0) {
                    const uint64_t h_x = offset >> 8;
                    __builtin_prefetch(m_probe.hints_bytes +
                                       (((h_x - 1) * m_probe.hints_width) >> 3));
                    s.st = state::walk_hint;
                } else {
                    s.st = state::walk_locate;
                }
                return;
            }
        }
        /* the minimizer is present but the kmer is not */
        s.memo_head_res = lookup_result();
        emit_negative(s, lookup_result());
    }

    /* WALK, stage 3b: read the endpoints hint (its line was prefetched)
       and prefetch the lines where the endpoints search will start. */
    void step_walk_hint(slot_t& s) {
        const uint64_t h_x = s.match_offset >> 8;
        assert(m_probe.valid and h_x != 0);
        if (h_x - 1 < m_probe.hints_size) {
            const uint64_t bit = (h_x - 1) * m_probe.hints_width;
            uint64_t word;
            std::memcpy(&word, m_probe.hints_bytes + (bit >> 3), sizeof(word));
            const uint64_t p = (word >> (bit & 7)) &
                               ((uint64_t(1) << m_probe.hints_width) - 1);
            const uint64_t begin = p - h_x + 1;
            __builtin_prefetch(m_probe.low_bytes + begin);
            __builtin_prefetch(m_probe.high_words + (p >> 6));
        }
        s.st = state::walk_locate;
    }

    /* WALK, stage 4: resolve the matched locus to its string. */
    void step_walk_locate(slot_t& s) {
        auto const& seg = seg_of(s);
        lookup_result res;
        res.kmer_offset = s.match_offset;
        res.kmer_orientation =
            s.match_backward ? constants::backward_orientation : constants::forward_orientation;
        m_dict->m_spss.strings_offsets.offset_to_id(res, s.cache.positions[0], m_k);
        if (res.kmer_offset >= res.string_end - m_k + 1) {
            /* the matched locus spans a string boundary: rare; redo the
               verification in full, which tries the later candidates */
            res = m_dict->m_spss.lookup_from_positions(walk_positions(s), s.walk_size,
                                                        slot_kmer(seg, s.cursor),
                                                        slot_kmer_rc(seg, s.cursor), s.info);
        }
        s.memo_head_res = res;
        if (res.kmer_id != constants::invalid_uint64) {
            emit_positive(s, res);
        } else {
            emit_negative(s, res);
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
};

}  // namespace sshash
