#pragma once

#include <cstdlib>
#include <cstring>
#if defined(__x86_64__)
#include <x86intrin.h>
#endif

#include "include/dictionary.hpp"
#include "include/endpoints_probe.hpp"
#include "include/util.hpp"

namespace sshash {

/*
    Batched lookups of independent kmers: a small batch advances in
    lockstep, and every dependent memory access of every lookup is
    prefetched one step before it is read.

    A single lookup is a chain of dependent reads at essentially random
    addresses -- the minimizers MPHF, a control codeword, the locate set,
    the text, and, for positives, the string-endpoints search -- so one
    query at a time leaves almost all of the memory system's concurrent
    misses unused. The engine holds `num_slots` lookups in flight and
    sweeps over them, advancing each by one bounded stage per sweep; each
    stage consumes memory prefetched by the previous stage of the same
    lookup, roughly a sweep earlier, and issues the prefetch for the next.
    Finished slots are refilled from the input, so the batch stays full.

    The stages are those of the streaming engine's walk pipeline
    (streaming_query_batched), without any streaming state: minimizers
    MPHF in groups of `group_size` through its own batched lookup (which
    hashes every key and prefetches the seeds of their buckets before
    reading any); then, per lookup and one prefetched stage per sweep, the
    control codeword, the locate-set decode (mid-load buckets decoded in
    full, up to 2^l + 1 positions, into a slot-local array), the candidate
    text verification, the endpoints hint, and the resolution of a matched
    locus to its string. Kmers of skew buckets resolve through their own
    stages, grouped by partition through the kmers MPHF's batched lookup.

    Results are written in input order and are identical to
    `dictionary::lookup` on every kmer (debug assert).
*/
template <typename Dict>
struct lookup_query_batched  //
{
    using kmer_t = typename Dict::kmer_type;
    using bucket_cache = typename Dict::spss_type::bucket_cache;

    static constexpr uint64_t num_slots = 64;
    static constexpr uint64_t group_size = 16;
    static constexpr uint64_t bits_per_char = kmer_t::bits_per_char;

    lookup_query_batched(Dict const* dict)
        : m_dict(dict), m_k(dict->m_k), m_m(dict->m_m)  //
    {
        m_slots.resize(num_slots);
        m_probe.capture(dict->m_spss.strings_offsets);
        m_timings = std::getenv("SSHASH_TIMINGS") != nullptr;
    }

    /* Resolve kmers[0..n) into out[0..n). */
    void lookup_batch(kmer_t const* kmers, const uint64_t n, lookup_result* out) {
        struct source {
            kmer_t const* kmers;
            kmer_t get(uint64_t i) const { return kmers[i]; }
        };
        run(source{kmers}, n, out);
    }

    /* Resolve the kmers spelled by queries[0..n) (strings of k valid
       characters) into out[0..n). */
    void lookup_batch(char const* const* queries, const uint64_t n, lookup_result* out) {
        struct source {
            char const* const* queries;
            uint64_t k;
            kmer_t get(uint64_t i) const {
                char const* str = queries[i];
#if defined(__BMI2__)
                /*
                    Pack eight characters per step: for both nucleotide
                    encodings the code of a character is a function of bits
                    1-3 of that same byte, so the whole word can be shifted
                    and the two code bits of every byte extracted in one
                    pext. The tail of k % 8 characters stays scalar, so no
                    byte beyond str[k-1] is ever read.
                */
                if constexpr (kmer_t::bits_per_char == 2 and kmer_t::has_reverse_complement) {
                    static_assert(kmer_t::uint_kmer_bits % 64 == 0);
                    constexpr uint64_t num_words = kmer_t::uint_kmer_bits / 64;
                    uint64_t words[num_words] = {};
                    const uint64_t num_chunks = k / 8;
                    for (uint64_t j = 0; j != num_chunks; ++j) {
                        uint64_t v;
                        std::memcpy(&v, str + 8 * j, 8);
#ifdef SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING
                        const uint64_t y = (v >> 1) ^ (v >> 2);
#else
                        const uint64_t y = v >> 1;
#endif
                        const uint64_t codes = _pext_u64(y, 0x0303030303030303ULL);
                        words[(16 * j) >> 6] |= codes << ((16 * j) & 63);
                    }
                    for (uint64_t i = 8 * num_chunks; i != k; ++i) {
                        words[(2 * i) >> 6] |= kmer_t::char_to_uint(str[i]) << ((2 * i) & 63);
                    }
                    kmer_t x = 0;
                    std::memcpy(&x.bits, words, sizeof(x.bits));
                    assert(x == util::string_to_uint_kmer<kmer_t>(str, k));
                    return x;
                }
#endif
                return util::string_to_uint_kmer<kmer_t>(str, k);
            }
        };
        run(source{queries, m_k}, n, out);
    }

private:
    Dict const* m_dict;
    uint64_t m_k, m_m;
    endpoints_probe m_probe;

    /* ---------------- timings (diagnostic) ----------------
       With SSHASH_TIMINGS set, every step is timed with the TSC and
       attributed to its stage; a share table prints at destruction. The
       sweep-overhead bin also absorbs the cost of the TSC reads. */
    static constexpr int num_timing_bins = 10;
    static constexpr char const* timing_names[num_timing_bins] = {
        "load (kmer + minimizer)", "walk: start (MPHF group)", "walk: codeword prefetch",
        "walk: codeword",          "walk: bucket decode",      "walk: verify",
        "walk: endpoints hint",    "walk: locate",             "heavy walk stages",
        "sweep overhead"};
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
        uint64_t inside = 0;
        for (int b = 0; b != 9; ++b) inside += cycles[b];
        cycles[9] = cycles[9] > inside ? cycles[9] - inside : 0;
        uint64_t total = 0;
        for (int b = 0; b != num_timing_bins; ++b) total += cycles[b];
        if (total == 0) return;
        std::cout << "==== lookup_query_batched timings (cycles, shares):\n";
        for (int b = 0; b != num_timing_bins; ++b) {
            std::cout << "  " << timing_names[b] << " = " << cycles[b] << " ("
                      << (100.0 * cycles[b]) / total << "%), steps = " << m_steps[b] << '\n';
        }
        std::cout << "  total accounted cycles = " << total << std::endl;
    }

public:
    ~lookup_query_batched() { print_timings(); }

private:

    enum class state : uint8_t {
        idle,
        mphf_wait,    // filed, waiting for a minimizers-MPHF group
        walk_pf,      // prefetch the control codeword line
        walk_code,    // read the codeword (line prefetched)
        walk_bucket,  // decode the locate set (lines prefetched)
        walk_verify,  // compare the candidate loci (text prefetched)
        walk_hint,    // read the endpoints hint, prefetch the search entry
        walk_locate,  // resolve the matched locus to its string
        heavy_wait,   // filed, waiting for a kmers-MPHF group (skew bucket)
        heavy_pf,     // prefetch the kmer's line of the partition positions
        heavy_pos,    // read the position in the bucket (line prefetched)
        heavy_off,    // read the kmer's offset, prefetch the candidate text
    };

    struct slot_t {
        state st = state::idle;
        uint64_t out_index = 0;
        kmer_t kmer = 0, kmer_rc = 0;
        minimizer_info info;
        /* the walk in progress */
        uint64_t walk_id = 0;
        uint64_t walk_code = 0;
        uint64_t walk_size = 0;
        bool walk_heavy = false;
        uint64_t match_offset = 0;
        bool match_backward = false;
        uint8_t heavy_partition = 0;
        uint64_t heavy_begin = 0;
        kmer_t heavy_key = 0;
        /* decoded locate sets */
        static constexpr uint64_t max_bucket_size = (uint64_t(1) << constants::min_l) + 1;
        bucket_cache cache;
        std::array<
            typename std::decay_t<decltype(std::declval<bucket_cache>().positions)>::value_type,
            max_bucket_size>
            big_positions;
    };
    std::vector<slot_t> m_slots;
    std::vector<uint32_t> m_mphf_queue, m_heavy_queue;
    std::array<std::vector<uint32_t>, 8> m_heavy_parts;
    lookup_result* m_out = nullptr;

    template <typename Source>
    void run(Source&& source, const uint64_t n, lookup_result* out) {
        m_out = out;
        m_mphf_queue.clear();
        m_heavy_queue.clear();
        uint64_t next = 0;
        uint64_t active = 0;
        for (auto& s : m_slots) {
            s.st = state::idle;
            if (next != n) {
                load(s, source, next++);
                active += 1;
            }
        }

        const uint64_t t_run = m_timings ? tsc() : 0;
        while (active != 0) {
            bool advanced = false;
            for (auto& s : m_slots) {
                const uint64_t t0 = m_timings ? tsc() : 0;
                int bin;
                switch (s.st) {
                    case state::idle:
                    case state::mphf_wait:
                    case state::heavy_wait: continue;
                    case state::walk_pf: step_walk_pf(s); bin = 2; break;
                    case state::walk_code: step_walk_code(s); bin = 3; break;
                    case state::walk_bucket: step_walk_bucket(s); bin = 4; break;
                    case state::walk_verify: step_walk_verify(s); bin = 5; break;
                    case state::walk_hint: step_walk_hint(s); bin = 6; break;
                    case state::walk_locate: step_walk_locate(s); bin = 7; break;
                    case state::heavy_pf: step_heavy_pf(s); bin = 8; break;
                    case state::heavy_pos: step_heavy_pos(s); bin = 8; break;
                    case state::heavy_off: step_heavy_off(s); bin = 8; break;
                    default: bin = 9; break;
                }
                if (m_timings) account(bin, t0);
                advanced = true;
                if (s.st == state::idle) {
                    if (next != n) {
                        const uint64_t t1 = m_timings ? tsc() : 0;
                        load(s, source, next++);
                        if (m_timings) account(0, t1);
                    } else {
                        active -= 1;
                    }
                }
            }
            if (m_mphf_queue.size() >= group_size or (!advanced and !m_mphf_queue.empty())) {
                const uint64_t t0 = m_timings ? tsc() : 0;
                start_walks();
                if (m_timings) account(1, t0);
            }
            if (m_heavy_queue.size() >= group_size or (!advanced and !m_heavy_queue.empty())) {
                const uint64_t t0 = m_timings ? tsc() : 0;
                start_heavies();
                if (m_timings) account(8, t0);
            }
        }
        if (m_timings) m_cycles[9] += tsc() - t_run;
        assert(m_mphf_queue.empty());
        assert(m_heavy_queue.empty());
    }

    template <typename Source>
    void load(slot_t& s, Source const& source, const uint64_t i) {
        s.out_index = i;
        s.kmer = source.get(i);
        s.kmer_rc = s.kmer;
        s.kmer_rc.reverse_complement_inplace(m_k);
        s.info = minimizer(s.kmer, s.kmer_rc);
        s.st = state::mphf_wait;
        m_mphf_queue.push_back(&s - m_slots.data());
    }

    minimizer_info minimizer(kmer_t const& kmer, kmer_t const& kmer_rc) const {
#if defined(__AVX2__)
        if constexpr (kmer_t::has_reverse_complement and bits_per_char == 2 and
                      std::is_same_v<hasher_type, mixer_64>) {
            return minimizer_avx2(kmer, kmer_rc);
        }
#endif
        return util::compute_minimizer<kmer_t>(kmer, kmer_rc, m_k, m_m, m_dict->m_hasher);
    }

#if defined(__AVX2__)
    /*
        The minimizer of one isolated kmer, four loci per iteration: the
        m-mers of loci i..i+3 are extracted with variable shifts from the
        kmer's packed words, reverse-complemented in-vector (the crc64 of
        kmer.hpp: complement, byte swap, two nibble swaps, final shift),
        canonicalized by unsigned minimum, hashed with the mixer, and folded
        into a per-lane running minimum; a scalar pass reduces the four
        lanes and the tail loci. Any hash equality along the way raises a
        (conservative) tie flag, on which the computation falls back to the
        generic `compute_minimizer`, whose tie-breaking is the single
        source of truth; true ties occur on ~1e-4 of the windows.
    */
    minimizer_info minimizer_avx2(kmer_t const& kmer, kmer_t const& kmer_rc) const {
        static_assert(kmer_t::uint_kmer_bits % 64 == 0);
        const uint64_t w = m_k - m_m + 1;
        if (w < 8) {
            return util::compute_minimizer<kmer_t>(kmer, kmer_rc, m_k, m_m, m_dict->m_hasher);
        }
        constexpr uint64_t num_words = kmer_t::uint_kmer_bits / 64;
        alignas(32) uint64_t words[num_words + 2];
        std::memcpy(words, &kmer.bits, sizeof(kmer.bits));
        words[num_words] = 0;
        words[num_words + 1] = 0;

        const uint64_t mb = m_m * bits_per_char;
        const uint64_t mask64 = (uint64_t(1) << mb) - 1;
        const __m256i mask = _mm256_set1_epi64x(mask64);
        const __m256i mulc = _mm256_set1_epi64x(0x517cc1b727220a95ULL);
        const uint64_t magic64 = m_dict->m_hasher.hash(0);
        const __m256i magic = _mm256_set1_epi64x(magic64);
        const __m256i sign = _mm256_set1_epi64x(int64_t(uint64_t(1) << 63));
        const __m256i c64 = _mm256_set1_epi64x(64);

        /* the m-mer of locus i read on the reverse strand is the m-mer of
           kmer_rc at the mirror locus k - m - i, i.e. at bit 2*(k - m) - 2*i */
        alignas(32) uint64_t words_rc[num_words + 2];
        std::memcpy(words_rc, &kmer_rc.bits, sizeof(kmer_rc.bits));
        words_rc[num_words] = 0;
        words_rc[num_words + 1] = 0;
        const __m256i bitr_base = _mm256_set1_epi64x(int64_t(2 * (m_k - m_m)));

        __m256i min_hash = _mm256_set1_epi64x(-1);  // unsigned max
        __m256i min_value = _mm256_setzero_si256();
        __m256i min_pos = _mm256_setzero_si256();
        __m256i tie = _mm256_setzero_si256();

        auto mullo64 = [](__m256i a, __m256i b) {
            const __m256i lo = _mm256_mul_epu32(a, b);
            const __m256i ah = _mm256_srli_epi64(a, 32);
            const __m256i bh = _mm256_srli_epi64(b, 32);
            const __m256i cross =
                _mm256_add_epi64(_mm256_mul_epu32(ah, b), _mm256_mul_epu32(a, bh));
            return _mm256_add_epi64(lo, _mm256_slli_epi64(cross, 32));
        };

        uint64_t i = 0;
        for (; i + 4 <= w; i += 4) {
            const __m256i bit = _mm256_setr_epi64x(int64_t(2 * i), int64_t(2 * (i + 1)),
                                                   int64_t(2 * (i + 2)), int64_t(2 * (i + 3)));

            /* extract four m-mers at bit positions `at` from `src` */
            auto extract = [&](uint64_t const* src, const __m256i at) {
                if constexpr (num_words == 1) {
                    return _mm256_and_si256(
                        _mm256_srlv_epi64(_mm256_set1_epi64x(int64_t(src[0])), at), mask);
                } else {
                    /* per-lane word pair, stitched with variable shifts; shift
                       counts >= 64 yield zero, which handles the boundaries */
                    const __m256i in_hi = _mm256_cmpgt_epi64(at, _mm256_set1_epi64x(63));
                    const __m256i w0 =
                        _mm256_blendv_epi8(_mm256_set1_epi64x(int64_t(src[0])),
                                           _mm256_set1_epi64x(int64_t(src[1])), in_hi);
                    const __m256i w1 =
                        _mm256_blendv_epi8(_mm256_set1_epi64x(int64_t(src[1])),
                                           _mm256_set1_epi64x(int64_t(src[2])), in_hi);
                    const __m256i sh = _mm256_and_si256(at, _mm256_set1_epi64x(63));
                    return _mm256_and_si256(
                        _mm256_or_si256(_mm256_srlv_epi64(w0, sh),
                                        _mm256_sllv_epi64(w1, _mm256_sub_epi64(c64, sh))),
                        mask);
                }
            };

            const __m256i mmer = extract(words, bit);
            const __m256i r = extract(words_rc, _mm256_sub_epi64(bitr_base, bit));

            /* canonical = unsigned min(mmer, rc) */
            const __m256i rc_lt =
                _mm256_cmpgt_epi64(_mm256_xor_si256(mmer, sign), _mm256_xor_si256(r, sign));
            const __m256i value = _mm256_blendv_epi8(mmer, r, rc_lt);
            const __m256i hash = _mm256_xor_si256(mullo64(value, mulc), magic);

            /* per-lane running minimum, with equality detection */
            tie = _mm256_or_si256(tie, _mm256_cmpeq_epi64(hash, min_hash));
            const __m256i lt = _mm256_cmpgt_epi64(_mm256_xor_si256(min_hash, sign),
                                                  _mm256_xor_si256(hash, sign));
            min_hash = _mm256_blendv_epi8(min_hash, hash, lt);
            min_value = _mm256_blendv_epi8(min_value, value, lt);
            const __m256i pos = _mm256_setr_epi64x(int64_t(i), int64_t(i + 1), int64_t(i + 2),
                                                   int64_t(i + 3));
            min_pos = _mm256_blendv_epi8(min_pos, pos, lt);
        }

        if (_mm256_movemask_epi8(tie) != 0) {
            return util::compute_minimizer<kmer_t>(kmer, kmer_rc, m_k, m_m, m_dict->m_hasher);
        }

        alignas(32) uint64_t hs[4], vs[4], ps[4];
        _mm256_store_si256(reinterpret_cast<__m256i*>(hs), min_hash);
        _mm256_store_si256(reinterpret_cast<__m256i*>(vs), min_value);
        _mm256_store_si256(reinterpret_cast<__m256i*>(ps), min_pos);

        /* the tail loci, scalar, and the final reduction */
        uint64_t best_hash = hs[0], best_value = vs[0], best_pos = ps[0];
        uint64_t num_mins = 1;
        for (int l = 1; l != 4; ++l) {
            if (hs[l] < best_hash) {
                best_hash = hs[l];
                best_value = vs[l];
                best_pos = ps[l];
                num_mins = 1;
            } else if (hs[l] == best_hash) {
                num_mins += 1;
            }
        }
        for (; i != w; ++i) {
            const uint64_t bit = 2 * i;
            uint64_t mmer = words[bit >> 6] >> (bit & 63);
            if ((bit & 63) != 0) mmer |= words[(bit >> 6) + 1] << (64 - (bit & 63));
            mmer &= mask64;
            const uint64_t value =
                std::min(mmer, kmer_t::reverse_complement_mmer(mmer, m_m));
            const uint64_t hash = (value * 0x517cc1b727220a95ULL) ^ magic64;
            if (hash < best_hash) {
                best_hash = hash;
                best_value = value;
                best_pos = i;
                num_mins = 1;
            } else if (hash == best_hash) {
                num_mins += 1;
            }
        }
        if (num_mins > 1) {
            return util::compute_minimizer<kmer_t>(kmer, kmer_rc, m_k, m_m, m_dict->m_hasher);
        }
        minimizer_info info(best_value, best_pos);
        assert([&] {
            auto expected =
                util::compute_minimizer<kmer_t>(kmer, kmer_rc, m_k, m_m, m_dict->m_hasher);
            return minimizer_info(info.minimizer, info.pos_in_kmer) ==
                   minimizer_info(expected.minimizer, expected.pos_in_kmer);
        }());
        return info;
    }
#endif

    void emit(slot_t& s, lookup_result const& res) {
        assert(equal_lookup_result(m_dict->lookup(s.kmer), res));
        m_out[s.out_index] = res;
        s.st = state::idle;
    }

    /* Start a group of walks through the minimizers MPHF's batched lookup. */
    void start_walks() {
        auto const& mphf = m_dict->m_ssi.codewords.mphf;
        std::array<uint64_t, group_size> keys;
        std::array<uint64_t, group_size> ids;
        while (!m_mphf_queue.empty()) {
            const uint64_t nn = std::min<uint64_t>(group_size, m_mphf_queue.size());
            const uint64_t first = m_mphf_queue.size() - nn;
            for (uint64_t i = 0; i != group_size; ++i) {
                const uint64_t q = i < nn ? first + i : first;
                keys[i] = m_slots[m_mphf_queue[q]].info.minimizer;
            }
            mphf.template lookup_batch<group_size>(keys.data(), ids.data());
            for (uint64_t i = 0; i != nn; ++i) {
                auto& s = m_slots[m_mphf_queue[first + i]];
                s.walk_id = ids[i];
                s.st = state::walk_pf;
            }
            m_mphf_queue.resize(first);
            if (m_mphf_queue.size() < group_size) break;
        }
    }

    /* WALK, stage 0: one codeword-line prefetch per step (prefetches issued
       in bursts contend with the MPHF's reads for fill buffers and get
       dropped). */
    void step_walk_pf(slot_t& s) {
        auto const& codewords = m_dict->m_ssi.codewords.control_codewords;
        auto const* bytes = reinterpret_cast<char const*>(codewords.data().data());
        __builtin_prefetch(bytes + ((s.walk_id * codewords.width()) >> 3));
        s.st = state::walk_code;
    }

    /* WALK, stage 1: read the codeword and prefetch what the bucket kind
       needs. */
    void step_walk_code(slot_t& s) {
        auto const& ssi = m_dict->m_ssi;
        auto const& strings = m_dict->m_spss.strings;
        const uint64_t code = ssi.codewords.control_codewords.access(s.walk_id);
        s.walk_code = code;
        s.cache.size = 0;
        s.walk_heavy = false;

        if ((code & 1) == bucket_t::SINGLETON) {
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
        /* HEAVYLOAD: kmer-specific walk through the skew index */
        const uint64_t c = code >> 2;
        s.heavy_partition = c & 7;
        s.heavy_begin = c >> 3;
        if constexpr (kmer_t::has_reverse_complement) {
            s.heavy_key = std::min(s.kmer, s.kmer_rc);
        } else {
            s.heavy_key = s.kmer;
        }
        s.st = state::heavy_wait;
        m_heavy_queue.push_back(&s - m_slots.data());
    }

    /* WALK, stage 2 (mid-load buckets): decode the locate set in full and
       prefetch the candidate text. */
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

    auto const* walk_positions(slot_t const& s) const {
        return s.walk_size <= bucket_cache::max_cached_size ? s.cache.positions.data()
                                                            : s.big_positions.data();
    }

    /* WALK, stage 3: minimizer presence check, then the candidate loci,
       against prefetched text. */
    void step_walk_verify(slot_t& s) {
        auto const& strings = m_dict->m_spss.strings;
        const uint64_t size = s.walk_size;
        auto const* positions = walk_positions(s);

        /* presence check, as spss::lookup does; for a skew bucket a
           mismatch at the located position proves nothing */
        if (uint64_t read_mmer = uint64_t(util::read_kmer_at<kmer_t>(
                strings, m_m, bits_per_char * positions[0].absolute_offset));
            read_mmer != s.info.minimizer and
            util::canonical_mmer<kmer_t>(read_mmer, m_m) != s.info.minimizer) {
            emit(s, s.walk_heavy ? lookup_result() : lookup_result(false));
            return;
        }

        uint64_t loci[2] = {s.info.pos_in_kmer, m_k - m_m - s.info.pos_in_kmer};
        const int num_loci = kmer_t::has_reverse_complement ? 2 : 1;
        for (uint64_t t = 0; t != size; ++t) {
            for (int l = 0; l != num_loci; ++l) {
                if (positions[t].absolute_offset < loci[l]) continue;
                const uint64_t offset = positions[t].absolute_offset - loci[l];
                auto read_kmer = util::read_kmer_at<kmer_t>(strings, m_k, bits_per_char * offset);
                if (read_kmer != s.kmer and read_kmer != s.kmer_rc) continue;
                s.match_offset = offset;
                s.match_backward = read_kmer == s.kmer_rc;
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
        emit(s, lookup_result());
    }

    /* WALK, stage 3b: read the endpoints hint (line prefetched) and
       prefetch the lines where the endpoints search will start. */
    void step_walk_hint(slot_t& s) {
        const uint64_t h_x = s.match_offset >> 8;
        assert(m_probe.valid and h_x != 0);
        if (h_x - 1 < m_probe.hints_size) {
            const uint64_t bit = (h_x - 1) * m_probe.hints_width;
            uint64_t word;
            std::memcpy(&word, m_probe.hints_bytes + (bit >> 3), sizeof(word));
            const uint64_t p =
                (word >> (bit & 7)) & ((uint64_t(1) << m_probe.hints_width) - 1);
            const uint64_t begin = p - h_x + 1;
            __builtin_prefetch(m_probe.low_bytes + begin);
            __builtin_prefetch(m_probe.high_words + (p >> 6));
        }
        s.st = state::walk_locate;
    }

    /* WALK, stage 4: resolve the matched locus to its string. */
    void step_walk_locate(slot_t& s) {
        lookup_result res;
        res.kmer_offset = s.match_offset;
        res.kmer_orientation =
            s.match_backward ? constants::backward_orientation : constants::forward_orientation;
        m_dict->m_spss.strings_offsets.offset_to_id(res, walk_positions(s)[0], m_k);
        if (res.kmer_offset >= res.string_end - m_k + 1) {
            /* the matched locus spans a string boundary: rare; redo the
               verification in full, which tries the later candidates */
            res = m_dict->m_spss.lookup_from_positions(walk_positions(s), s.walk_size, s.kmer,
                                                       s.kmer_rc, s.info);
        }
        emit(s, res);
    }

    /* Start the filed skew-bucket walks, grouped by partition, through the
       kmers MPHF's batched lookup. */
    void start_heavies() {
        auto const& ski = m_dict->m_ssi.ski;
        for (auto& part : m_heavy_parts) part.clear();
        for (uint32_t id : m_heavy_queue) {
            m_heavy_parts[m_slots[id].heavy_partition].push_back(id);
        }
        m_heavy_queue.clear();

        std::array<kmer_t, group_size> keys;
        std::array<uint64_t, group_size> ids;
        for (uint64_t p = 0; p != m_heavy_parts.size(); ++p) {
            auto const& part = m_heavy_parts[p];
            if (part.empty()) continue;
            auto const& mphf = ski.mphfs[p];
            for (uint64_t begin = 0; begin < part.size(); begin += group_size) {
                const uint64_t nn = std::min<uint64_t>(group_size, part.size() - begin);
                for (uint64_t i = 0; i != group_size; ++i) {
                    keys[i] = m_slots[part[begin + (i < nn ? i : 0)]].heavy_key;
                }
                mphf.template lookup_batch<group_size>(keys.data(), ids.data());
                for (uint64_t i = 0; i != nn; ++i) {
                    auto& s = m_slots[part[begin + i]];
                    s.walk_id = ids[i];
                    s.st = state::heavy_pf;
                }
            }
        }
    }

    void step_heavy_pf(slot_t& s) {
        auto const& positions = m_dict->m_ssi.ski.positions[s.heavy_partition];
        auto const* bytes = reinterpret_cast<char const*>(positions.data().data());
        __builtin_prefetch(bytes + ((s.walk_id * positions.width()) >> 3));
        s.st = state::heavy_pos;
    }

    void step_heavy_pos(slot_t& s) {
        auto const& ski = m_dict->m_ssi.ski;
        const uint64_t pos_in_bucket = ski.positions[s.heavy_partition].access(s.walk_id);
        s.walk_id = s.heavy_begin + pos_in_bucket;
        auto const& hlb = ski.heavy_load_buckets;
        auto const* bytes = reinterpret_cast<char const*>(hlb.data().data());
        __builtin_prefetch(bytes + ((s.walk_id * hlb.width()) >> 3));
        s.st = state::heavy_off;
    }

    void step_heavy_off(slot_t& s) {
        auto const& strings = m_dict->m_spss.strings;
        const uint64_t offset = m_dict->m_ssi.ski.heavy_load_buckets.access(s.walk_id);
        s.cache.positions[0] = m_dict->m_spss.strings_offsets.decode(offset);
        const uint64_t abs = s.cache.positions[0].absolute_offset;
        const uint64_t pos_in_kmer = s.info.pos_in_kmer;
        const uint64_t char_offset = abs >= pos_in_kmer ? abs - pos_in_kmer : 0;
        __builtin_prefetch(&strings.data()[(bits_per_char * char_offset) >> 6]);
        __builtin_prefetch(&strings.data()[(bits_per_char * abs) >> 6]);
        s.walk_size = 1;
        s.walk_heavy = true;
        s.st = state::walk_verify;
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
