#pragma once

#include "../dictionary.hpp"
#include "../minimizer_enumerator.hpp"
#include "../util.hpp"

namespace sshash {

struct streaming_query_regular_parsing {
    streaming_query_regular_parsing(dictionary const* dict)

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

        , m_shift(2 * (dict->m_k - 1))
        , m_k(dict->m_k)
        , m_m(dict->m_m)
        , m_seed(dict->m_seed)

        , m_string_iterator(dict->m_buckets.strings, 0)
        , m_begin(0)
        , m_end(0)
        , m_pos_in_window(0)
        , m_window_size(0)

        , m_reverse(false)

        , m_num_searches(0)
        , m_num_extensions(0)

    {
        assert(!m_dict->m_canonical_parsing);
    }

    inline void start() { m_start = true; }

    lookup_result lookup_advanced(const char* kmer) {
        /* 1. validation */
        bool is_valid = m_start ? util::is_valid(kmer, m_k) : util::is_valid(kmer[m_k - 1]);
        if (!is_valid) {
            m_start = true;
            return lookup_result();
        }

        /* 2. compute kmer and minimizer */
        if (!m_start) {
            m_kmer >>= 2;
            m_kmer += (util::char_to_uint(kmer[m_k - 1])) << m_shift;
            assert(m_kmer == util::string_to_uint_kmer(kmer, m_k));
        } else {
            m_kmer = util::string_to_uint_kmer(kmer, m_k);
        }
        m_curr_minimizer = m_minimizer_enum.next(m_kmer, m_start);
        assert(m_curr_minimizer == util::compute_minimizer(m_kmer, m_k, m_m, m_seed));
        m_kmer_rc = util::compute_reverse_complement(m_kmer, m_k);
        constexpr bool reverse = true;
        m_curr_minimizer_rc = m_minimizer_enum_rc.next<reverse>(m_kmer_rc, m_start);
        assert(m_curr_minimizer_rc == util::compute_minimizer(m_kmer_rc, m_k, m_m, m_seed));

        bool both_minimizers_not_found = (same_minimizer() and m_minimizer_not_found) and
                                         (same_minimizer_rc() and m_minimizer_rc_not_found);
        if (both_minimizers_not_found) {
            update_state();
            assert(equal_lookup_result(m_dict->lookup_advanced(kmer), m_res));
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

        assert(equal_lookup_result(m_dict->lookup_advanced(kmer), m_res));
        return m_res;
    }

    uint64_t num_searches() const { return m_num_searches; }
    uint64_t num_extensions() const { return m_num_extensions; }

private:
    dictionary const* m_dict;

    /* result */
    lookup_result m_res;

    /* (kmer,minimizer) state */
    minimizer_enumerator<> m_minimizer_enum;
    minimizer_enumerator<> m_minimizer_enum_rc;
    bool m_minimizer_not_found, m_minimizer_rc_not_found;
    bool m_start;
    uint64_t m_curr_minimizer, m_curr_minimizer_rc;
    uint64_t m_prev_minimizer, m_prev_minimizer_rc;
    kmer_t m_kmer, m_kmer_rc;

    /* constants */
    uint64_t m_shift, m_k, m_m, m_seed;

    /* string state */
    bit_vector_iterator m_string_iterator;
    uint64_t m_begin, m_end;
    uint64_t m_pos_in_window, m_window_size;
    bool m_reverse;

    /* performance counts */
    uint64_t m_num_searches;
    uint64_t m_num_extensions;

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
            uint64_t log2_bucket_size = util::ceil_log2_uint32(num_super_kmers_in_bucket);
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
            uint64_t log2_bucket_size = util::ceil_log2_uint32(num_super_kmers_in_bucket);
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
            uint64_t pos_in_string = 2 * offset;
            m_reverse = false;
            m_string_iterator.at(pos_in_string);
            auto [res, offset_end] = (m_dict->m_buckets).offset_to_id(offset, m_k);
            m_res = res;
            m_pos_in_window = 0;
            m_window_size = std::min<uint64_t>(m_k - m_m + 1, offset_end - offset - m_k + 1);

            while (m_pos_in_window != m_window_size) {
                kmer_t val = m_string_iterator.read(2 * m_k);

                if (check_minimizer and super_kmer_id == begin and m_pos_in_window == 0) {
                    uint64_t minimizer = util::compute_minimizer(val, m_k, m_m, m_seed);
                    if (minimizer != m_curr_minimizer) {
                        m_minimizer_not_found = true;
                        m_res = lookup_result();
                        return;
                    }
                }

                m_string_iterator.eat(2);
                m_pos_in_window += 1;
                pos_in_string += 2;
                assert(m_pos_in_window <= m_window_size);

                if (m_kmer == val) {
                    m_num_searches += 1;
                    m_res.kmer_orientation = constants::forward_orientation;
                    return;
                }

                m_res.kmer_id += 1;
                m_res.kmer_id_in_contig += 1;
            }
        }

        m_res = lookup_result();
    }

    void lookup_advanced_rc(uint64_t begin, uint64_t end, bool check_minimizer) {
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            uint64_t offset = (m_dict->m_buckets).offsets.access(super_kmer_id);
            uint64_t pos_in_string = 2 * offset;
            m_reverse = false;
            m_string_iterator.at(pos_in_string);
            auto [res, offset_end] = (m_dict->m_buckets).offset_to_id(offset, m_k);
            m_res = res;
            m_pos_in_window = 0;
            m_window_size = std::min<uint64_t>(m_k - m_m + 1, offset_end - offset - m_k + 1);

            while (m_pos_in_window != m_window_size) {
                kmer_t val = m_string_iterator.read(2 * m_k);

                if (check_minimizer and super_kmer_id == begin and m_pos_in_window == 0) {
                    uint64_t minimizer = util::compute_minimizer(val, m_k, m_m, m_seed);
                    if (minimizer != m_curr_minimizer_rc) {
                        m_minimizer_rc_not_found = true;
                        m_res = lookup_result();
                        return;
                    }
                }

                m_string_iterator.eat(2);
                m_pos_in_window += 1;
                pos_in_string += 2;
                assert(m_pos_in_window <= m_window_size);

                if (m_kmer_rc == val) {
                    m_reverse = true;
                    pos_in_string -= 2;
                    m_num_searches += 1;
                    m_string_iterator.at(pos_in_string + 2 * (m_k - 1));
                    m_res.kmer_orientation = constants::backward_orientation;
                    return;
                }

                m_res.kmer_id += 1;
                m_res.kmer_id_in_contig += 1;
            }
        }

        m_res = lookup_result();
    }

    inline void extend() {
        m_string_iterator.eat(2);
        m_pos_in_window += 1;
        assert(m_pos_in_window <= m_window_size);
        assert(m_res.kmer_orientation == constants::forward_orientation);
        m_res.kmer_id += 1;
        m_res.kmer_id_in_contig += 1;
    }

    inline void extend_rc() {
        assert(m_reverse);
        m_string_iterator.eat_reverse(2);
        m_pos_in_window -= 1;
        assert(m_pos_in_window >= 1);
        assert(m_res.kmer_orientation == constants::backward_orientation);
        m_res.kmer_id -= 1;
        m_res.kmer_id_in_contig -= 1;
    }

    inline bool extends() {
        if (m_pos_in_window == m_window_size) return false;
        if (m_kmer == m_string_iterator.read(2 * m_k)) {
            ++m_num_extensions;
            return true;
        }
        return false;
    }

    inline bool extends_rc() {
        assert(m_reverse);
        if (m_pos_in_window == 1) return false;
        if (m_kmer_rc == m_string_iterator.read_reverse(2 * m_k)) {
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