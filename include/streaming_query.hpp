#pragma once

#include "include/dictionary.hpp"
#include "include/minimizer_enumerator.hpp"
#include "include/util.hpp"

namespace sshash {

template <class kmer_t, bool canonical>
struct streaming_query {
    streaming_query(dictionary<kmer_t> const* dict)

        : m_dict(dict)

        , m_start(true)
        , m_kmer(constants::invalid_uint64)
        , m_kmer_rc(constants::invalid_uint64)
        , m_k(dict->m_k)

        , m_minimizer_enum(dict->m_k, dict->m_m, dict->m_hasher)
        , m_minimizer_enum_rc(dict->m_k, dict->m_m, dict->m_hasher)
        , m_curr_minimizer(constants::invalid_uint64)
        , m_prev_minimizer(constants::invalid_uint64)
        , m_curr_minimizer_rc(constants::invalid_uint64)
        , m_prev_minimizer_rc(constants::invalid_uint64)

        , m_it(dict->m_buckets.strings, m_k)
        , m_remaining_contig_bases(0)

        , m_num_searches(0)
        , m_num_extensions(0)
        , m_num_invalid(0)
        , m_num_negative(0)

    {
        if (canonical != m_dict->m_canonical) {
            std::stringstream ss;
            ss << "dict.canonical() = " << (m_dict->canonical() ? "true" : "false")
               << " but required " << (canonical ? "true" : "false");
            throw std::runtime_error(ss.str());
        }
    }

    void reset() {
        m_start = true;
        m_remaining_contig_bases = 0;
        m_res = lookup_result();
    }

    lookup_result lookup_advanced(char const* kmer)  //
    {
        /* 1. validation */
        bool is_valid =
            m_start ? util::is_valid<kmer_t>(kmer, m_k) : kmer_t::is_valid(kmer[m_k - 1]);
        if (!is_valid) {
            m_num_invalid += 1;
            reset();
            return m_res;
        }

        /* 2. compute uint kmer from input char kmer and minimizer */
        if (!m_start) {
            m_kmer.drop_char();
            m_kmer.set(m_k - 1, kmer_t::char_to_uint(kmer[m_k - 1]));
            assert(m_kmer == util::string_to_uint_kmer<kmer_t>(kmer, m_k));
            m_kmer_rc.pad_char();
            m_kmer_rc.set(0, kmer_t::char_to_uint(
                                 kmer_t::canonicalize_basepair_reverse_map[int(kmer[m_k - 1])]));
            m_kmer_rc.take(m_k * kmer_t::bits_per_char);
        } else {
            m_kmer = util::string_to_uint_kmer<kmer_t>(kmer, m_k);
            m_kmer_rc = m_kmer;
            m_kmer_rc.reverse_complement_inplace(m_k);
        }

        uint64_t minimizer = m_minimizer_enum.template next<false>(m_kmer, m_start);
        uint64_t minimizer_rc = m_minimizer_enum_rc.template next<true>(m_kmer_rc, m_start);
        if constexpr (canonical) {
            m_curr_minimizer = std::min(minimizer, minimizer_rc);
        } else {
            m_curr_minimizer = minimizer;
            m_curr_minimizer_rc = minimizer_rc;
        }

        /* 3. compute result */
        if (m_remaining_contig_bases == 0) {
            seed();
        } else {
            auto expected_kmer = (m_res.kmer_orientation == constants::forward_orientation)
                                     ? (m_it.next(), m_it.get())
                                     : (m_it.next_reverse(), m_it.get_reverse());
            if ((expected_kmer == m_kmer) or (expected_kmer == m_kmer_rc)) {
                ++m_num_extensions;
                m_res.kmer_id += m_res.kmer_orientation;
                m_res.kmer_id_in_contig += m_res.kmer_orientation;
                m_remaining_contig_bases -= 1;
            } else {
                seed();
            }
        }

        /* 4. update state */
        m_prev_minimizer = m_curr_minimizer;
        m_prev_minimizer_rc = m_curr_minimizer_rc;
        m_start = false;

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

    /* kmer state */
    bool m_start;
    kmer_t m_kmer, m_kmer_rc;
    uint64_t m_k;

    /* minimizer state */
    minimizer_enumerator<kmer_t> m_minimizer_enum;
    minimizer_enumerator<kmer_t> m_minimizer_enum_rc;
    uint64_t m_curr_minimizer, m_prev_minimizer;
    uint64_t m_curr_minimizer_rc, m_prev_minimizer_rc;

    /* string state */
    kmer_iterator<kmer_t> m_it;
    uint64_t m_remaining_contig_bases;

    /* performance counts */
    uint64_t m_num_searches;
    uint64_t m_num_extensions;
    uint64_t m_num_invalid;
    uint64_t m_num_negative;

    void seed()  //
    {
        m_remaining_contig_bases = 0;

        /* if minimizer does not change and previous minimizer was not found,
           surely any kmer having the same minimizer cannot be found as well */
        if (m_curr_minimizer == m_prev_minimizer and        //
            m_curr_minimizer_rc == m_prev_minimizer_rc and  // always true for canonical = true
            m_res.minimizer_found == false)                 //
        {
            assert(m_res.kmer_id == constants::invalid_uint64);
            m_num_negative += 1;
            return;
        }

        if constexpr (canonical) {
            m_res = m_dict->lookup_uint_canonical(m_kmer, m_kmer_rc, m_curr_minimizer);
            if (m_res.kmer_id == constants::invalid_uint64) {
                m_num_negative += 1;
                return;
            }
        } else {
            m_res = m_dict->lookup_uint_regular(m_kmer, m_curr_minimizer);
            bool minimizer_found = m_res.minimizer_found;
            if (m_res.kmer_id == constants::invalid_uint64) {
                assert(m_res.kmer_orientation == constants::forward_orientation);
                m_res = m_dict->lookup_uint_regular(m_kmer_rc, m_curr_minimizer_rc);
                m_res.kmer_orientation = constants::backward_orientation;
                bool minimizer_rc_found = m_res.minimizer_found;
                m_res.minimizer_found = minimizer_rc_found or minimizer_found;
                if (m_res.kmer_id == constants::invalid_uint64) {
                    m_num_negative += 1;
                    return;
                }
            }
        }

        assert(m_res.minimizer_found == true);
        m_num_searches += 1;
        uint64_t kmer_offset = 2 * (m_res.kmer_id + m_res.contig_id * (m_k - 1));
        m_remaining_contig_bases = (m_res.contig_size - 1) - m_res.kmer_id_in_contig;
        if (m_res.kmer_orientation == constants::backward_orientation) {
            kmer_offset += 2 * m_k;
            m_remaining_contig_bases = m_res.kmer_id_in_contig;
        }
        m_it.at(kmer_offset);
    }
};

}  // namespace sshash