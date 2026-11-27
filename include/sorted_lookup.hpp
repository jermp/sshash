#pragma once

#include <cstdlib>
#include <cstring>
#include <vector>

#include "include/dictionary.hpp"
#include "include/lookup_kernels.hpp"
#include "include/util.hpp"

namespace sshash {

/*
    Batched lookups of independent kmers for batches that are large
    relative to the index: instead of hiding the latency of random
    accesses (lookup_query_batched), eliminate the randomness.

    The batch is processed in chunks. Within a chunk, every data structure
    of the index is visited in one pass, in the order of its own layout:

      1. compute all minimizers (the shared vector front end), and group
         the queries by minimizer, ordered by minimizer hash -- which is
         also the order in which the minimizers MPHF lays out its buckets,
         so the batched MPHF lookups walk its memory near-sequentially;
      2. resolve the distinct minimizers to their MPHF ids, then visit the
         control codewords in id order, which is their storage order;
      3. visit the mid-load locate sets in bucket-address order;
      4. resolve the kmers of heavy (skew) buckets through the per-
         partition MPHFs, their position and offset arrays visited in
         index order;
      5. verify every candidate text locus in ascending text order -- one
         ordered pass over the packed strings;
      6. resolve the matched loci to their strings in ascending offset
         order -- one ordered pass over the string endpoints.

    Between passes, requests are reordered by LSD radix sort on the next
    pass's address. The per-query cost is therefore a few sorts of small
    records at memory bandwidth plus an amortized share of sequential
    sweeps, instead of a chain of random accesses; the crossover against
    lookup_query_batched is an empirical matter of batch size versus index
    size.

    Results are written in input order and are identical to
    `dictionary::lookup` on every kmer (debug assert).
*/
template <typename Dict>
struct lookup_query_sorted  //
{
    using kmer_t = typename Dict::kmer_type;
    using decoded_t = decltype(std::declval<Dict>().m_spss.strings_offsets.decode(0));

    static constexpr uint64_t bits_per_char = kmer_t::bits_per_char;
    static constexpr uint64_t group_size = 16;
    static constexpr uint64_t default_chunk_size = uint64_t(1) << 23;

    lookup_query_sorted(Dict const* dict)
        : m_dict(dict), m_k(dict->m_k), m_m(dict->m_m)  //
    {
        m_chunk_size = default_chunk_size;
        if (char const* e = std::getenv("SSHASH_SORTED_CHUNK")) {
            const uint64_t v = std::strtoull(e, nullptr, 10);
            if (v >= 1024) m_chunk_size = v;
        }
        m_no_avx512 = std::getenv("SSHASH_NO_AVX512") != nullptr;
        /* key widths for the radix sorts, from the index's own sizes */
        const uint64_t num_chars = dict->m_spss.strings.num_bits() / bits_per_char;
        m_offset_bytes = bytes_for(num_chars);
        m_id_bytes = bytes_for(dict->m_ssi.codewords.control_codewords.size());
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
                return lookup_kernels::pack_string<kmer_t>(queries[i], k);
            }
        };
        run(source{queries, m_k}, n, out);
    }

private:
    Dict const* m_dict;
    uint64_t m_k, m_m;
    uint64_t m_chunk_size;
    bool m_no_avx512 = false;
    int m_offset_bytes = 8, m_id_bytes = 8;

    static int bytes_for(uint64_t max_value) {
        int b = 1;
        while (b < 8 and (max_value >> (8 * b)) != 0) b += 1;
        return b;
    }

    /* ---------------- sorting ---------------- */

    struct rec_t {
        uint64_t key;
        uint64_t payload;
    };

    /* LSD radix sort of (key, payload) records on the low key_bytes bytes
       of the key; stable, two buffers. */
    static void radix_sort(std::vector<rec_t>& v, std::vector<rec_t>& tmp, const int key_bytes) {
        const uint64_t n = v.size();
        if (n <= 1) return;
        tmp.resize(n);
        uint64_t count[256];
        rec_t* src = v.data();
        rec_t* dst = tmp.data();
        for (int b = 0; b != key_bytes; ++b) {
            std::memset(count, 0, sizeof(count));
            const int shift = 8 * b;
            for (uint64_t i = 0; i != n; ++i) count[(src[i].key >> shift) & 255] += 1;
            uint64_t sum = 0;
            for (int d = 0; d != 256; ++d) {
                const uint64_t c = count[d];
                count[d] = sum;
                sum += c;
            }
            for (uint64_t i = 0; i != n; ++i) {
                dst[count[(src[i].key >> shift) & 255]++] = src[i];
            }
            std::swap(src, dst);
        }
        if (src != v.data()) std::memcpy(v.data(), src, n * sizeof(rec_t));
    }

    /* ---------------- per-chunk scratch ---------------- */

    std::vector<kmer_t> m_kmers, m_rcs;
    std::vector<uint32_t> m_pos_in_kmer;
    std::vector<rec_t> m_recs, m_tmp;       // reused by every sort
    std::vector<rec_t> m_cands, m_cands_tmp;
    /* distinct minimizers (SoA) */
    std::vector<uint64_t> m_dmini;
    std::vector<uint32_t> m_dbegin, m_dcount;  // group in m_qids
    std::vector<uint32_t> m_qids;              // query ids grouped by distinct minimizer
    std::vector<uint64_t> m_did;
    std::vector<uint32_t> m_dpos_begin;
    std::vector<uint16_t> m_dpos_count;
    std::vector<uint8_t> m_dheavy;
    std::vector<uint64_t> m_dheavy_code;
    std::vector<decoded_t> m_pool;
    struct match_t {
        uint64_t offset;
        uint32_t qid;
        uint32_t anchor;
        bool backward;
    };
    std::vector<match_t> m_matches;

    template <typename Source>
    void run(Source&& source, const uint64_t n, lookup_result* out) {
        for (uint64_t begin = 0; begin < n; begin += m_chunk_size) {
            const uint64_t cn = std::min(m_chunk_size, n - begin);
            run_chunk(source, begin, cn, out + begin);
        }
    }

    template <typename Source>
    void run_chunk(Source const& source, const uint64_t begin, const uint64_t cn,
                   lookup_result* out)  //
    {
        auto const& strings = m_dict->m_spss.strings;
        auto const& ssi = m_dict->m_ssi;

        /* ---- pass 0: front end -------------------------------------- */
        m_kmers.resize(cn);
        m_rcs.resize(cn);
        m_pos_in_kmer.resize(cn);
        m_minis.resize(cn);
        m_recs.resize(cn);
        for (uint64_t i = 0; i != cn; ++i) {
            kmer_t kmer = source.get(begin + i);
            kmer_t rc = kmer;
            rc.reverse_complement_inplace(m_k);
            auto info = lookup_kernels::compute_minimizer<kmer_t>(kmer, rc, m_k, m_m,
                                                                  m_dict->m_hasher, m_no_avx512);
            m_kmers[i] = kmer;
            m_rcs[i] = rc;
            m_pos_in_kmer[i] = info.pos_in_kmer;
            m_minis[i] = info.minimizer;
            m_recs[i] = {m_dict->m_hasher.hash(info.minimizer), i};
            out[i] = lookup_result();  // negative by default
        }

        /* ---- pass 1: group by minimizer, in hash order --------------- */
        radix_sort(m_recs, m_tmp, 8);

        m_dmini.clear();
        m_dbegin.clear();
        m_dcount.clear();
        m_qids.resize(cn);
        for (uint64_t i = 0; i != cn; ++i) m_qids[i] = uint32_t(m_recs[i].payload);
        {
            uint64_t i = 0;
            while (i != cn) {
                uint64_t j = i + 1;
                while (j != cn and m_recs[j].key == m_recs[i].key) j += 1;
                /* within an equal-hash run, split by minimizer value
                   (collisions of distinct minimizers are ~never) */
                for (uint64_t a = i; a != j;) {
                    const uint64_t mini = m_minis[m_qids[a]];
                    uint64_t b = a;
                    while (b != j and m_minis[m_qids[b]] == mini) b += 1;
                    m_dmini.push_back(mini);
                    m_dbegin.push_back(uint32_t(a));
                    m_dcount.push_back(uint32_t(b - a));
                    a = b;
                }
                i = j;
            }
        }
        const uint64_t nd = m_dmini.size();

        /* ---- pass 2: minimizers MPHF, in hash order ------------------ */
        m_did.resize(nd);
        {
            auto const& mphf = ssi.codewords.mphf;
            std::array<uint64_t, group_size> keys;
            std::array<uint64_t, group_size> ids;
            for (uint64_t d = 0; d < nd; d += group_size) {
                const uint64_t nn = std::min<uint64_t>(group_size, nd - d);
                for (uint64_t i = 0; i != group_size; ++i) {
                    keys[i] = m_dmini[d + (i < nn ? i : 0)];
                }
                mphf.template lookup_batch<group_size>(keys.data(), ids.data());
                for (uint64_t i = 0; i != nn; ++i) m_did[d + i] = ids[i];
            }
        }

        /* ---- pass 3: control codewords, in id order ------------------ */
        m_recs.resize(nd);
        for (uint64_t d = 0; d != nd; ++d) m_recs[d] = {m_did[d], d};
        radix_sort(m_recs, m_tmp, m_id_bytes);

        m_dpos_begin.assign(nd, 0);
        m_dpos_count.assign(nd, 0);
        m_dheavy.assign(nd, 0);
        m_dheavy_code.assign(nd, 0);
        m_pool.clear();
        m_cands.clear();  // reused here for the mid-load requests
        for (uint64_t r = 0; r != nd; ++r) {
            const uint64_t d = m_recs[r].payload;
            const uint64_t code = ssi.codewords.control_codewords.access(m_recs[r].key);
            if ((code & 1) == bucket_t::SINGLETON) {
                m_dpos_begin[d] = uint32_t(m_pool.size());
                m_dpos_count[d] = 1;
                m_pool.push_back(m_dict->m_spss.strings_offsets.decode(code >> 1));
            } else if ((code & 3) == bucket_t::MIDLOAD) {
                constexpr uint64_t mask = (uint64_t(1) << constants::min_l) - 1;
                const uint64_t c = code >> 2;
                const uint64_t bucket_size = (c & mask) + 2;
                const uint64_t bucket_id = c >> constants::min_l;
                const uint64_t bucket_begin =
                    ssi.begin_buckets_of_size[bucket_size] + bucket_id * bucket_size;
                m_cands.push_back({bucket_begin, (uint64_t(bucket_size) << 32) | d});
            } else {
                m_dheavy[d] = 1;
                m_dheavy_code[d] = code >> 2;  // | bucket begin | partition (3 bits) |
            }
        }

        /* ---- pass 3b: mid-load locate sets, in bucket-address order -- */
        radix_sort(m_cands, m_cands_tmp, m_offset_bytes);
        for (auto const& r : m_cands) {
            const uint64_t d = r.payload & 0xffffffff;
            const uint64_t bucket_size = r.payload >> 32;
            const uint64_t bucket_begin = r.key;
            m_dpos_begin[d] = uint32_t(m_pool.size());
            m_dpos_count[d] = uint16_t(bucket_size);
            typename sparse_and_skew_index<kmer_t>::bucket_iterator it{
                &ssi, bucket_begin, bucket_size, bucket_t::MIDLOAD};
            for (uint64_t t = 0; t != bucket_size; ++t, ++it) {
                m_pool.push_back(m_dict->m_spss.strings_offsets.decode(*it));
            }
        }

        /* ---- pass 4: minimizer presence, in text order --------------- */
        m_recs.clear();
        for (uint64_t d = 0; d != nd; ++d) {
            if (m_dheavy[d]) continue;
            m_recs.push_back({m_pool[m_dpos_begin[d]].absolute_offset, d});
        }
        radix_sort(m_recs, m_tmp, m_offset_bytes);
        m_dverdict.assign(nd, 0);
        for (auto const& r : m_recs) {
            const uint64_t d = r.payload;
            const uint64_t read_mmer =
                uint64_t(util::read_kmer_at<kmer_t>(strings, m_m, bits_per_char * r.key));
            if (read_mmer == m_dmini[d] or
                util::canonical_mmer<kmer_t>(read_mmer, m_m) == m_dmini[d]) {
                m_dverdict[d] = 1;
            } else {
                /* minimizer absent: every query of the group is negative */
                for (uint64_t a = 0; a != m_dcount[d]; ++a) {
                    out[m_qids[m_dbegin[d] + a]] = lookup_result(false);
                }
            }
        }

        /* ---- pass 5: heavy (skew) buckets ---------------------------- */
        run_heavies(out);

        /* ---- pass 6: candidate loci ---------------------------------- */
        m_cands.swap(m_heavy_cands);  // candidates contributed by the heavy pass
        for (uint64_t d = 0; d != nd; ++d) {
            if (m_dheavy[d] or !m_dverdict[d]) continue;
            const uint32_t anchor = m_dpos_begin[d];
            for (uint64_t a = 0; a != m_dcount[d]; ++a) {
                const uint32_t q = m_qids[m_dbegin[d] + a];
                const uint64_t l0 = m_pos_in_kmer[q];
                const uint64_t l1 = m_k - m_m - l0;
                for (uint64_t t = 0; t != m_dpos_count[d]; ++t) {
                    const uint64_t position = m_pool[m_dpos_begin[d] + t].absolute_offset;
                    if (position >= l0) {
                        m_cands.push_back({position - l0, (uint64_t(anchor) << 32) | q});
                    }
                    if constexpr (kmer_t::has_reverse_complement) {
                        if (l1 != l0 and position >= l1) {
                            m_cands.push_back({position - l1, (uint64_t(anchor) << 32) | q});
                        }
                    }
                }
            }
        }

        /* ---- pass 7: verification, in text order --------------------- */
        radix_sort(m_cands, m_cands_tmp, m_offset_bytes);
        m_matches.clear();
        for (auto const& r : m_cands) {
            const uint32_t q = uint32_t(r.payload);
            auto read_kmer = util::read_kmer_at<kmer_t>(strings, m_k, bits_per_char * r.key);
            if (read_kmer != m_kmers[q] and read_kmer != m_rcs[q]) continue;
            m_matches.push_back({r.key, q, uint32_t(r.payload >> 32), read_kmer == m_rcs[q]});
        }

        /* ---- pass 8: locate, in text order --------------------------- */
        for (auto const& mt : m_matches) {
            lookup_result res;
            res.kmer_offset = mt.offset;
            res.kmer_orientation = mt.backward ? constants::backward_orientation
                                               : constants::forward_orientation;
            m_dict->m_spss.strings_offsets.offset_to_id(res, m_pool[mt.anchor], m_k);
            if (res.kmer_offset >= res.string_end - m_k + 1) {
                /* the matched locus spans a string boundary: rare; the
                   one-at-a-time lookup is the result by definition */
                res = m_dict->lookup(m_kmers[mt.qid]);
            }
            out[mt.qid] = res;
        }

#ifndef NDEBUG
        for (uint64_t i = 0; i != cn; ++i) {
            assert(equal_lookup_result(m_dict->lookup(m_kmers[i]), out[i]));
        }
#endif
    }

    std::vector<uint64_t> m_minis;
    std::vector<uint8_t> m_dverdict;
    std::vector<rec_t> m_heavy_cands;

    /* Heavy (skew) buckets: each kmer of a heavy group walks the skew
       index on its own -- per-partition kmers MPHF, position in bucket,
       then the kmer's text offset -- each array visited in sorted order.
       The verification itself joins the common pass 7; the minimizer
       presence check is skipped because a kmer that matches at its
       located offset has, by construction, its minimizer at the anchor,
       and a kmer that does not match is negative either way (and
       equal results to `dictionary::lookup` are asserted per query). */
    void run_heavies(lookup_result* out) {
        (void)out;
        auto const& ski = m_dict->m_ssi.ski;
        m_heavy_cands.clear();

        /* collect (partition, query, bucket begin) */
        std::array<std::vector<rec_t>, 8> parts;
        const uint64_t nd = m_dmini.size();
        for (uint64_t d = 0; d != nd; ++d) {
            if (!m_dheavy[d]) continue;
            const uint64_t c = m_dheavy_code[d];
            const uint64_t partition = c & 7;
            const uint64_t bucket_begin = c >> 3;
            for (uint64_t a = 0; a != m_dcount[d]; ++a) {
                const uint32_t q = m_qids[m_dbegin[d] + a];
                parts[partition].push_back({bucket_begin, q});
            }
        }

        std::array<kmer_t, group_size> keys;
        std::array<uint64_t, group_size> ids;
        std::vector<rec_t> recs, tmp;
        for (uint64_t p = 0; p != parts.size(); ++p) {
            auto& part = parts[p];
            if (part.empty()) continue;
            auto const& mphf = ski.mphfs[p];

            /* kmers MPHF in groups; then its position array in id order */
            recs.resize(part.size());
            for (uint64_t i = 0; i < part.size(); i += group_size) {
                const uint64_t nn = std::min<uint64_t>(group_size, part.size() - i);
                for (uint64_t j = 0; j != group_size; ++j) {
                    const uint32_t q = uint32_t(part[i + (j < nn ? j : 0)].payload);
                    if constexpr (kmer_t::has_reverse_complement) {
                        keys[j] = std::min(m_kmers[q], m_rcs[q]);
                    } else {
                        keys[j] = m_kmers[q];
                    }
                }
                mphf.template lookup_batch<group_size>(keys.data(), ids.data());
                for (uint64_t j = 0; j != nn; ++j) recs[i + j] = {ids[j], i + j};
            }
            radix_sort(recs, tmp, 8);

            /* position in bucket, in id order; then the heavy-load offsets
               in address order */
            auto const& positions = ski.positions[p];
            for (auto& r : recs) {
                const uint64_t pos_in_bucket = positions.access(r.key);
                const uint64_t i = r.payload;
                r = {part[i].key + pos_in_bucket, part[i].payload};  // (hlb address, qid)
            }
            radix_sort(recs, tmp, 8);

            auto const& hlb = ski.heavy_load_buckets;
            for (auto const& r : recs) {
                const uint32_t q = uint32_t(r.payload);
                const uint64_t offset = hlb.access(r.key);
                const uint32_t anchor = uint32_t(m_pool.size());
                m_pool.push_back(m_dict->m_spss.strings_offsets.decode(offset));
                const uint64_t position = m_pool[anchor].absolute_offset;
                const uint64_t l0 = m_pos_in_kmer[q];
                const uint64_t l1 = m_k - m_m - l0;
                if (position >= l0) {
                    m_heavy_cands.push_back({position - l0, (uint64_t(anchor) << 32) | q});
                }
                if constexpr (kmer_t::has_reverse_complement) {
                    if (l1 != l0 and position >= l1) {
                        m_heavy_cands.push_back({position - l1, (uint64_t(anchor) << 32) | q});
                    }
                }
            }
        }
    }
};

}  // namespace sshash
