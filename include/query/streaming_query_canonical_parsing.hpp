#pragma once

#include "include/dictionary.hpp"
#include "include/util.hpp"

namespace sshash {

template <class kmer_t>
struct streaming_query_canonical_parsing {
    streaming_query_canonical_parsing(dictionary<kmer_t> const* dict)

        : m_dict(dict)

        , m_kmer(constants::invalid_uint64)
        , m_kmer_rc(constants::invalid_uint64)

        , m_k(dict->m_k)
        , m_m(dict->m_m)
        , m_seed(dict->m_seed)

        , m_string_iterator(dict->m_buckets.strings, 0)
        , m_remaining_contig_bases(0)
        , m_direction(0)

        , m_num_searches(0)
        , m_num_extensions(0)

    {
        start();
        assert(m_dict->m_canonical_parsing);
    }

    void start() { m_start = true; }

    void reset_state() {
        start();
        m_kmer = constants::invalid_uint64;
        m_kmer_rc = constants::invalid_uint64;
        m_string_iterator.at(0);
        m_remaining_contig_bases = 0;
        m_direction = 0;
        m_res = lookup_result();
    }

    lookup_result lookup_advanced(const char* kmer) {
        // std::cout << "querying '" << std::string(kmer, kmer + m_k) << "'..." << std::endl;

        /* 1. validation */
        bool is_valid =
            m_start ? util::is_valid<kmer_t>(kmer, m_k) : kmer_t::is_valid(kmer[m_k - 1]);
        if (!is_valid) {
            // std::cout << "kmer is invalid" << std::endl;
            m_start = true;
            return lookup_result();
        }

        /* 2. compute kmer */
        if (!m_start) {
            m_kmer.drop_char();
            m_kmer.kth_char_or(m_k - 1, kmer_t::char_to_uint(kmer[m_k - 1]));
            assert(m_kmer == util::string_to_uint_kmer<kmer_t>(kmer, m_k));
        } else {
            m_kmer = util::string_to_uint_kmer<kmer_t>(kmer, m_k);
        }
        m_kmer_rc = m_kmer;
        m_kmer_rc.reverse_complement_inplace(m_k);

        // std::cout << "m_remaining_contig_bases = " << m_remaining_contig_bases << std::endl;

        if (m_remaining_contig_bases == 0) m_start = true;

        /* 3. compute result */
        if (m_start or m_res.kmer_id == constants::invalid_uint64) {
            /* if at the start of a new query or previous kmer was not found */
            seed();
        } else {
            if (m_res.kmer_orientation == constants::backward_orientation) {
                m_string_iterator.eat_reverse(2);
                if (m_kmer_rc == m_string_iterator.read_reverse(2 * m_k)) {
                    ++m_num_extensions;
                    m_res.kmer_id += m_direction;
                    m_res.kmer_id_in_contig += m_direction;
                    m_remaining_contig_bases -= 1;
                } else {
                    seed();
                }
            } else {
                m_string_iterator.eat(2);
                if (m_kmer == m_string_iterator.read(2 * m_k)) {
                    ++m_num_extensions;
                    m_res.kmer_id += m_direction;
                    m_res.kmer_id_in_contig += m_direction;
                    m_remaining_contig_bases -= 1;
                } else {
                    seed();
                }
            }
        }

        /* 4. update state */
        m_start = false;

        // std::cout << m_res << std::endl;

        // auto exp_res = m_dict->lookup_advanced(kmer);
        // if (!equal_lookup_result(exp_res, m_res)) {
        //     // std::cout << "expected:" << std::endl;
        //     exp_res.print();
        //     assert(false);
        // }

        assert(equal_lookup_result(m_dict->lookup_advanced(kmer), m_res));

        return m_res;
    }

    uint64_t num_searches() const { return m_num_searches; }
    uint64_t num_extensions() const { return m_num_extensions; }

private:
    dictionary<kmer_t> const* m_dict;

    /* result */
    lookup_result m_res;

    /* kmer state */
    bool m_start;
    kmer_t m_kmer, m_kmer_rc;

    /* constants */
    uint64_t m_k, m_m, m_seed;

    /* string state */
    bit_vector_iterator<kmer_t> m_string_iterator;
    uint64_t m_remaining_contig_bases;
    int64_t m_direction;

    /* performance counts */
    uint64_t m_num_searches;
    uint64_t m_num_extensions;

    void seed() {
        // std::cout << "seeding" << std::endl;

        // warning: this internally recomputes m_kmer_rc
        m_res = m_dict->lookup_uint_canonical_parsing(m_kmer);
        if (m_res.kmer_id == constants::invalid_uint64) {
            // std::cout << "kmer not found..." << std::endl;
            return;
        }

        // std::cout << "--> kmer FOUND!" << std::endl;

        m_num_searches += 1;
        uint64_t kmer_offset = 2 * (m_res.kmer_id + m_res.contig_id * (m_k - 1));
        m_direction = m_res.kmer_orientation == constants::forward_orientation ? 1 : -1;
        kmer_offset += (m_direction == 1) ? 0 : (2 * m_k);
        m_string_iterator.at(kmer_offset);

        m_remaining_contig_bases = (m_direction == 1)
                                       ? ((m_res.contig_size - 1) - m_res.kmer_id_in_contig)
                                       : (m_res.kmer_id_in_contig);
    }
};

}  // namespace sshash