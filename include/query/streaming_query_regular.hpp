#pragma once

#include "include/dictionary.hpp"
#include "include/minimizer_enumerator.hpp"
#include "include/util.hpp"

namespace sshash {

template <class kmer_t>
struct streaming_query_regular {
    streaming_query_regular(dictionary<kmer_t> const* dict)

        : m_dict(dict)

        , m_minimizer_enum(dict->m_k, dict->m_m, dict->m_seed)
        , m_minimizer_enum_rc(dict->m_k, dict->m_m, dict->m_seed)
        , m_minimizer_not_found(false)
        , m_minimizer_rc_not_found(false)

        , m_start(true)

        , m_curr_minimizer(constants::invalid_uint64)
        , m_curr_minimizer_rc(constants::invalid_uint64)
        , m_prev_minimizer(constants::invalid_uint64)
        , m_prev_minimizer_rc(constants::invalid_uint64)

        , m_kmer(constants::invalid_uint64)

        , m_k(dict->m_k)
        , m_m(dict->m_m)
        , m_seed(dict->m_seed)

        , m_it(dict->m_buckets.strings, m_k)
        , m_begin(0)
        , m_end(0)
        , m_pos_in_window(0)
        , m_window_size(0)

        , m_reverse(false)

        , m_num_searches(0)
        , m_num_extensions(0)
        , m_num_invalid(0)
        , m_num_negative(0)

    {
        assert(!m_dict->m_canonical);
    }

    void reset() {
        m_start = true;
        m_minimizer_not_found = false;
        m_minimizer_rc_not_found = false;
        m_curr_minimizer = constants::invalid_uint64;
        m_prev_minimizer = constants::invalid_uint64;
        m_curr_minimizer_rc = constants::invalid_uint64;
        m_prev_minimizer_rc = constants::invalid_uint64;
        m_kmer = constants::invalid_uint64;
        m_begin = 0;
        m_end = 0;
        m_pos_in_window = 0;
        m_window_size = 0;
        m_res = lookup_result();
        m_reverse = false;
    }

    lookup_result lookup_advanced(const char* kmer) {
        /* 1. validation */
        bool is_valid =
            m_start ? util::is_valid<kmer_t>(kmer, m_k) : kmer_t::is_valid(kmer[m_k - 1]);
        if (!is_valid) {
            m_num_invalid += 1;
            reset();
            return m_res;
        }

        /* 2. compute kmer and minimizer */
        if (!m_start) {
            m_kmer.drop_char();
            m_kmer.set(m_k - 1, kmer_t::char_to_uint(kmer[m_k - 1]));
            assert(m_kmer == util::string_to_uint_kmer<kmer_t>(kmer, m_k));
        } else {
            m_kmer = util::string_to_uint_kmer<kmer_t>(kmer, m_k);
        }
        m_curr_minimizer = m_minimizer_enum.template next<false>(m_kmer, m_start);
        assert(m_curr_minimizer == util::compute_minimizer<kmer_t>(m_kmer, m_k, m_m, m_seed));
        m_kmer_rc = m_kmer;
        m_kmer_rc.reverse_complement_inplace(m_k);
        m_curr_minimizer_rc = m_minimizer_enum_rc.template next<true>(m_kmer_rc, m_start);
        assert(m_curr_minimizer_rc == util::compute_minimizer<kmer_t>(m_kmer_rc, m_k, m_m, m_seed));

        bool both_minimizers_not_found = (same_minimizer() and m_minimizer_not_found) and
                                         (same_minimizer_rc() and m_minimizer_rc_not_found);
        if (both_minimizers_not_found) {
            update_state();
            assert(equal_lookup_result(m_dict->lookup_advanced(kmer), m_res));
            m_num_negative += 1;
            return lookup_result();
        }

        /* 3. compute result */
        if (m_reverse) {
            if (same_minimizer_rc()) {
                if (!m_minimizer_rc_not_found) {
                    if (extends_rc()) {
                        extend_rc();
                    } else {
                        lookup_advanced_rc();
                    }
                }
            } else {
                m_res = lookup_result();
            }
        } else {
            if (same_minimizer()) {
                if (!m_minimizer_not_found) {
                    if (extends()) {
                        extend();
                    } else {
                        lookup_advanced();
                    }
                }
            } else {
                m_res = lookup_result();
            }
        }
        if (!found()) {
            m_minimizer_not_found = false;
            m_minimizer_rc_not_found = false;
            if (m_reverse) {
                /* If previous search was successfull for a reverse complement,
                   then try to search for a reverse complement first. */
                search_rc();
                if (!found()) search();
            } else {
                search();
                if (!found()) search_rc();
            }
        }

        /* 4. update state */
        update_state();

        if (!found()) m_num_negative += 1;

        assert(equal_lookup_result(m_dict->lookup_advanced(kmer), m_res));
        return m_res;
    }

    uint64_t num_searches() const { return m_num_searches; }
    uint64_t num_extensions() const { return m_num_extensions; }
    uint64_t num_positive_lookups() const { return num_searches() + num_extensions(); }
    uint64_t num_negative_lookups() const { return m_num_negative; }
    uint64_t num_invalid_lookups() const { return m_num_invalid; }

private:
    dictionary<kmer_t> const* m_dict;

    /* result */
    lookup_result m_res;

    /* (kmer,minimizer) state */
    minimizer_enumerator<kmer_t> m_minimizer_enum;
    minimizer_enumerator<kmer_t> m_minimizer_enum_rc;
    bool m_minimizer_not_found, m_minimizer_rc_not_found;
    bool m_start;
    uint64_t m_curr_minimizer, m_curr_minimizer_rc;
    uint64_t m_prev_minimizer, m_prev_minimizer_rc;
    kmer_t m_kmer, m_kmer_rc;

    /* constants */
    uint64_t m_k, m_m, m_seed;

    /* string state */
    kmer_iterator<kmer_t> m_it;
    uint64_t m_begin, m_end;
    uint64_t m_pos_in_window, m_window_size;
    bool m_reverse;

    /* performance counts */
    uint64_t m_num_searches;
    uint64_t m_num_extensions;
    uint64_t m_num_invalid;
    uint64_t m_num_negative;

    void update_state() {
        m_prev_minimizer = m_curr_minimizer;
        m_prev_minimizer_rc = m_curr_minimizer_rc;
        m_start = false;
    }

    inline bool found() { return m_res.kmer_id != constants::invalid_uint64; }
    inline bool same_minimizer() const { return m_curr_minimizer == m_prev_minimizer; }
    inline bool same_minimizer_rc() const { return m_curr_minimizer_rc == m_prev_minimizer_rc; }

    void locate_bucket() {
        uint64_t bucket_id = (m_dict->m_minimizers).lookup(m_curr_minimizer);
        std::tie(m_begin, m_end) = (m_dict->m_buckets).locate_bucket(bucket_id);
    }
    void locate_bucket_rc() {
        uint64_t bucket_id = (m_dict->m_minimizers).lookup(m_curr_minimizer_rc);
        std::tie(m_begin, m_end) = (m_dict->m_buckets).locate_bucket(bucket_id);
    }

    void lookup_advanced() {
        bool check_minimizer = !same_minimizer();
        if (!m_dict->m_skew_index.empty()) {
            uint64_t num_super_kmers_in_bucket = m_end - m_begin;
            uint64_t log2_bucket_size = bits::util::ceil_log2_uint32(num_super_kmers_in_bucket);
            if (log2_bucket_size > (m_dict->m_skew_index).min_log2) {
                uint64_t p = m_dict->m_skew_index.lookup(m_kmer, log2_bucket_size);
                if (p < num_super_kmers_in_bucket) {
                    lookup_advanced(m_begin + p, m_begin + p + 1, check_minimizer);
                    if (found()) return;
                }
                m_res = lookup_result();
                return;
            }
        }
        lookup_advanced(m_begin, m_end, check_minimizer);
    }

    void lookup_advanced_rc() {
        bool check_minimizer = !same_minimizer_rc();
        if (!m_dict->m_skew_index.empty()) {
            uint64_t num_super_kmers_in_bucket = m_end - m_begin;
            uint64_t log2_bucket_size = bits::util::ceil_log2_uint32(num_super_kmers_in_bucket);
            if (log2_bucket_size > (m_dict->m_skew_index).min_log2) {
                uint64_t p = m_dict->m_skew_index.lookup(m_kmer_rc, log2_bucket_size);
                if (p < num_super_kmers_in_bucket) {
                    lookup_advanced_rc(m_begin + p, m_begin + p + 1, check_minimizer);
                    if (found()) return;
                }
                m_res = lookup_result();
                return;
            }
        }
        lookup_advanced_rc(m_begin, m_end, check_minimizer);
    }

    void lookup_advanced(uint64_t begin, uint64_t end, bool check_minimizer) {
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            uint64_t offset = (m_dict->m_buckets).offsets.access(super_kmer_id);
            m_reverse = false;
            m_it.at(2 * offset);
            auto res = (m_dict->m_buckets).offset_to_id(offset, m_k);
            m_res = res;
            m_pos_in_window = 0;
            m_window_size =
                std::min<uint64_t>(m_k - m_m + 1, res.contig_end(m_k) - offset - m_k + 1);

            while (m_pos_in_window != m_window_size) {
                auto val = m_it.get();
                if (check_minimizer and super_kmer_id == begin and m_pos_in_window == 0) {
                    auto minimizer = util::compute_minimizer(val, m_k, m_m, m_seed);
                    if (minimizer != m_curr_minimizer) {
                        m_minimizer_not_found = true;
                        m_res = lookup_result();
                        return;
                    }
                }

                m_pos_in_window += 1;
                assert(m_pos_in_window <= m_window_size);

                if (m_kmer == val) {
                    m_it.next();
                    m_num_searches += 1;
                    m_res.kmer_orientation = constants::forward_orientation;
                    return;
                }

                m_res.kmer_id += 1;
                m_res.kmer_id_in_contig += 1;
                m_it.next();
            }
        }
        m_res = lookup_result();
    }

    void lookup_advanced_rc(uint64_t begin, uint64_t end, bool check_minimizer) {
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            uint64_t offset = (m_dict->m_buckets).offsets.access(super_kmer_id);
            uint64_t pos_in_string = 2 * offset;
            m_reverse = false;
            m_it.at(pos_in_string);
            auto res = (m_dict->m_buckets).offset_to_id(offset, m_k);
            m_res = res;
            m_pos_in_window = 0;
            m_window_size =
                std::min<uint64_t>(m_k - m_m + 1, res.contig_end(m_k) - offset - m_k + 1);

            while (m_pos_in_window != m_window_size) {
                auto val = m_it.get();
                if (check_minimizer and super_kmer_id == begin and m_pos_in_window == 0) {
                    auto minimizer = util::compute_minimizer(val, m_k, m_m, m_seed);
                    if (minimizer != m_curr_minimizer_rc) {
                        m_minimizer_rc_not_found = true;
                        m_res = lookup_result();
                        return;
                    }
                }

                m_pos_in_window += 1;
                pos_in_string += 2;
                assert(m_pos_in_window <= m_window_size);

                if (m_kmer_rc == val) {
                    m_reverse = true;
                    pos_in_string -= 2;
                    m_num_searches += 1;
                    m_it.at(pos_in_string + 2 * (m_k - 1));
                    m_res.kmer_orientation = constants::backward_orientation;
                    return;
                }

                m_res.kmer_id += 1;
                m_res.kmer_id_in_contig += 1;
                m_it.next();
            }
        }
        m_res = lookup_result();
    }

    inline void extend() {
        m_it.next();
        m_pos_in_window += 1;
        assert(m_pos_in_window <= m_window_size);
        assert(m_res.kmer_orientation == constants::forward_orientation);
        m_res.kmer_id += 1;
        m_res.kmer_id_in_contig += 1;
    }

    inline void extend_rc() {
        assert(m_reverse);
        m_it.next_reverse();
        m_pos_in_window -= 1;
        assert(m_pos_in_window >= 1);
        assert(m_res.kmer_orientation == constants::backward_orientation);
        m_res.kmer_id -= 1;
        m_res.kmer_id_in_contig -= 1;
    }

    inline bool extends() {
        if (m_pos_in_window == m_window_size) return false;
        if (m_kmer == m_it.get()) {
            ++m_num_extensions;
            return true;
        }
        return false;
    }

    inline bool extends_rc() {
        assert(m_reverse);
        if (m_pos_in_window == 1) return false;
        if (m_kmer_rc == m_it.get()) {
            ++m_num_extensions;
            return true;
        }
        return false;
    }

    void search() {
        locate_bucket();
        lookup_advanced();
    }

    void search_rc() {
        locate_bucket_rc();
        lookup_advanced_rc();
    }
};

}  // namespace sshash