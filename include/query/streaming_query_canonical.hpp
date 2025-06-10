#pragma once

#include "include/dictionary.hpp"
#include "include/util.hpp"

namespace sshash {

template <class kmer_t>
struct streaming_query_canonical {
    streaming_query_canonical(dictionary<kmer_t> const* dict)

        : m_dict(dict)

        , m_kmer(constants::invalid_uint64)
        , m_k(dict->m_k)

        , m_it(dict->m_buckets.strings, m_k)
        , m_remaining_contig_bases(constants::invalid_uint64)

        , m_num_searches(0)
        , m_num_extensions(0)
        , m_num_invalid(0)
        , m_num_negative(0)

    {
        start();
        assert(m_dict->m_canonical);
    }

    void start() { m_start = true; }

    void reset_state() {
        start();
        m_kmer = constants::invalid_uint64;
        m_remaining_contig_bases = constants::invalid_uint64;
        m_res = lookup_result();
    }

    lookup_result lookup_advanced(const char* kmer) {
        /* 1. validation */
        bool is_valid =
            m_start ? util::is_valid<kmer_t>(kmer, m_k) : kmer_t::is_valid(kmer[m_k - 1]);
        if (!is_valid) {
            m_num_invalid += 1;
            m_start = true;
            return lookup_result();
        }

        /* 2. compute uint kmer from input char kmer */
        if (!m_start) {
            m_kmer.drop_char();
            m_kmer.set(m_k - 1, kmer_t::char_to_uint(kmer[m_k - 1]));
            assert(m_kmer == util::string_to_uint_kmer<kmer_t>(kmer, m_k));
        } else {
            m_kmer = util::string_to_uint_kmer<kmer_t>(kmer, m_k);
        }

        if (m_remaining_contig_bases == 0) {
            m_start = true;
            // if (m_res.kmer_orientation == constants::forward_orientation) {
            //     m_it.at(m_it.position() + 2 * (m_k - 1));
            // } else {
            //     m_it.at(m_it.position() - 2 * (m_k - 1));
            // }
        }

        /* 3. compute result */
        if (m_start or m_res.kmer_id == constants::invalid_uint64) {
            /* if at the start of a new query or previous kmer was not found */
            seed();
        } else {
            if (m_res.kmer_orientation == constants::forward_orientation) {
                m_it.next();
            } else {
                m_it.next_reverse();
            }
            auto expected_kmer = m_it.get();
            auto expected_kmer_rc = expected_kmer;
            expected_kmer_rc.reverse_complement_inplace(m_k);
            if ((m_kmer == expected_kmer) or (m_kmer == expected_kmer_rc)) {
                ++m_num_extensions;
                m_res.kmer_id += m_res.kmer_orientation;
                m_res.kmer_id_in_contig += m_res.kmer_orientation;
                m_remaining_contig_bases -= 1;
            } else {
                seed();
            }
        }

        m_start = false;

        // std::cout << m_res << std::endl;

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
    kmer_t m_kmer;
    uint64_t m_k;

    /* string state */
    kmer_iterator<kmer_t> m_it;
    uint64_t m_remaining_contig_bases;

    /* performance counts */
    uint64_t m_num_searches;
    uint64_t m_num_extensions;
    uint64_t m_num_invalid;
    uint64_t m_num_negative;

    void seed() {
        m_res = m_dict->lookup_uint_canonical(m_kmer);
        if (m_res.kmer_id == constants::invalid_uint64) {
            m_num_negative += 1;
            return;
        }
        m_num_searches += 1;
        uint64_t kmer_offset = 2 * (m_res.kmer_id + m_res.contig_id * (m_k - 1));
        m_remaining_contig_bases = (m_res.contig_size - 1) - m_res.kmer_id_in_contig;
        bool reverse = false;
        if (m_res.kmer_orientation == constants::backward_orientation) {
            kmer_offset += 2 * m_k;
            m_remaining_contig_bases = m_res.kmer_id_in_contig;
            reverse = true;
        }
        m_it.at(kmer_offset, reverse);
    }
};

}  // namespace sshash