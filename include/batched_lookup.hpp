#pragma once

#include <cstring>

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
                return util::string_to_uint_kmer<kmer_t>(queries[i], k);
            }
        };
        run(source{queries, m_k}, n, out);
    }

private:
    Dict const* m_dict;
    uint64_t m_k, m_m;
    endpoints_probe m_probe;

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

        while (active != 0) {
            bool advanced = false;
            for (auto& s : m_slots) {
                switch (s.st) {
                    case state::idle:
                    case state::mphf_wait:
                    case state::heavy_wait: continue;
                    case state::walk_pf: step_walk_pf(s); break;
                    case state::walk_code: step_walk_code(s); break;
                    case state::walk_bucket: step_walk_bucket(s); break;
                    case state::walk_verify: step_walk_verify(s); break;
                    case state::walk_hint: step_walk_hint(s); break;
                    case state::walk_locate: step_walk_locate(s); break;
                    case state::heavy_pf: step_heavy_pf(s); break;
                    case state::heavy_pos: step_heavy_pos(s); break;
                    case state::heavy_off: step_heavy_off(s); break;
                }
                advanced = true;
                if (s.st == state::idle) {
                    if (next != n) {
                        load(s, source, next++);
                    } else {
                        active -= 1;
                    }
                }
            }
            if (m_mphf_queue.size() >= group_size or (!advanced and !m_mphf_queue.empty())) {
                start_walks();
            }
            if (m_heavy_queue.size() >= group_size or (!advanced and !m_heavy_queue.empty())) {
                start_heavies();
            }
        }
        assert(m_mphf_queue.empty());
        assert(m_heavy_queue.empty());
    }

    template <typename Source>
    void load(slot_t& s, Source const& source, const uint64_t i) {
        s.out_index = i;
        s.kmer = source.get(i);
        s.kmer_rc = s.kmer;
        s.kmer_rc.reverse_complement_inplace(m_k);
        s.info = util::compute_minimizer<kmer_t>(s.kmer, s.kmer_rc, m_k, m_m, m_dict->m_hasher);
        s.st = state::mphf_wait;
        m_mphf_queue.push_back(&s - m_slots.data());
    }

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
