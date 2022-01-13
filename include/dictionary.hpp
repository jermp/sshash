#pragma once

#include "minimizers.hpp"
#include "buckets.hpp"
#include "skew_index.hpp"

namespace sshash {

struct dictionary {
    dictionary() : m_size(0), m_seed(0), m_k(0), m_m(0), m_canonical_parsing(0) {}

    void build(std::string const& filename, build_configuration const& build_config);

    uint64_t size() const { return m_size; }
    uint64_t seed() const { return m_seed; }
    uint64_t k() const { return m_k; }
    uint64_t m() const { return m_m; }
    bool canonicalized() const { return m_canonical_parsing; }

    uint64_t lookup_uint64(uint64_t uint64_kmer) const {
        uint64_t minimizer = util::compute_minimizer(uint64_kmer, m_k, m_m, m_seed);
        uint64_t bucket_id = m_minimizers.lookup(minimizer);

        if (m_skew_index.empty()) return m_buckets.lookup(bucket_id, uint64_kmer, m_k, m_m);

        auto [begin, end] = m_buckets.locate_bucket(bucket_id);
        uint64_t num_strings_in_bucket = end - begin;
        uint64_t log2_num_strings_in_bucket = util::ceil_log2_uint32(num_strings_in_bucket);
        if (log2_num_strings_in_bucket > m_skew_index.min_log2) {
            uint64_t pos = m_skew_index.lookup(uint64_kmer, log2_num_strings_in_bucket);
            /* It must hold pos < num_strings_in_bucket for the kmer to exist. */
            if (pos < num_strings_in_bucket) {
                return m_buckets.lookup_in_string(begin + pos, uint64_kmer, m_k, m_m);
            }
            return constants::invalid;
        }

        return m_buckets.lookup(begin, end, uint64_kmer, m_k, m_m);
    }

    uint64_t lookup(char const* string_kmer, bool check_reverse_complement_too = true) const {
        uint64_t uint64_kmer = util::string_to_uint64_no_reverse(string_kmer, m_k);

        if (m_canonical_parsing) {  // canonical parsing
            uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, m_k);
            uint64_t minimizer = util::compute_minimizer(uint64_kmer, m_k, m_m, m_seed);
            uint64_t minimizer_rc = util::compute_minimizer(uint64_kmer_rc, m_k, m_m, m_seed);
            uint64_t bucket_id = m_minimizers.lookup(std::min<uint64_t>(minimizer, minimizer_rc));

            if (m_skew_index.empty()) {
                return m_buckets.lookup_canonical(bucket_id, uint64_kmer, uint64_kmer_rc, m_k, m_m);
            }

            auto [begin, end] = m_buckets.locate_bucket(bucket_id);
            uint64_t num_strings_in_bucket = end - begin;
            uint64_t log2_num_strings_in_bucket = util::ceil_log2_uint32(num_strings_in_bucket);
            if (log2_num_strings_in_bucket > m_skew_index.min_log2) {
                uint64_t pos = m_skew_index.lookup(uint64_kmer, log2_num_strings_in_bucket);
                if (pos < num_strings_in_bucket) {
                    uint64_t kmer_id =
                        m_buckets.lookup_in_string(begin + pos, uint64_kmer, m_k, m_m);
                    if (kmer_id != constants::invalid) return kmer_id;
                }
                uint64_t pos_rc = m_skew_index.lookup(uint64_kmer_rc, log2_num_strings_in_bucket);
                if (pos_rc < num_strings_in_bucket) {
                    uint64_t kmer_id =
                        m_buckets.lookup_in_string(begin + pos_rc, uint64_kmer_rc, m_k, m_m);
                    return kmer_id;
                }
                return constants::invalid;
            }
            return m_buckets.lookup_canonical(begin, end, uint64_kmer, uint64_kmer_rc, m_k, m_m);
        }

        // regular parsing
        uint64_t kmer_id = lookup_uint64(uint64_kmer);
        if (check_reverse_complement_too and kmer_id == constants::invalid) {
            uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, m_k);
            kmer_id = lookup_uint64(uint64_kmer_rc);
        }
        return kmer_id;
    }

    void access(uint64_t kmer_id, char* string_kmer) const {
        assert(kmer_id < size());
        m_buckets.access(kmer_id, string_kmer, m_k);
    }

    bool is_member(char const* string_kmer, bool check_reverse_complement_too = true) const {
        return lookup(string_kmer, check_reverse_complement_too) != constants::invalid;
    }

    friend struct membership_query_canonical_parsing;
    friend struct membership_query_regular_parsing;

    struct membership_query_result {
        membership_query_result() : num_kmers(0), num_valid_kmers(0), num_positive_kmers(0) {}
        uint64_t num_kmers;
        uint64_t num_valid_kmers;
        uint64_t num_positive_kmers;
        uint64_t num_searches;
        uint64_t num_extensions;
    };

    membership_query_result membership_query_from_file(std::string const& filename,
                                                       bool multiline) const;

    struct iterator {
        iterator(dictionary const* ptr, uint64_t kmer_id = 0) {
            it = ptr->m_buckets.at(kmer_id, ptr->m_k, ptr->m_size);
        }

        bool has_next() const { return it.has_next(); }
        std::pair<uint64_t, std::string> next() { return it.next(); }

    private:
        typename buckets::iterator it;
    };

    iterator begin() const { return iterator(this); }

    iterator at(uint64_t kmer_id) const {
        assert(kmer_id < size());
        return iterator(this, kmer_id);
    }

    uint64_t num_bits() const {
        return 8 * (sizeof(m_size) + sizeof(m_seed) + sizeof(m_k) + sizeof(m_m) +
                    sizeof(m_canonical_parsing)) +
               m_minimizers.num_bits() + m_buckets.num_bits() + m_skew_index.num_bits();
    }

    void print_info() const;
    void print_space_breakdown() const;

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_size);
        visitor.visit(m_seed);
        visitor.visit(m_k);
        visitor.visit(m_m);
        visitor.visit(m_canonical_parsing);
        visitor.visit(m_minimizers);
        visitor.visit(m_buckets);
        visitor.visit(m_skew_index);
    }

private:
    uint64_t m_size;
    uint64_t m_seed;
    uint16_t m_k;
    uint16_t m_m;
    uint16_t m_canonical_parsing;
    minimizers m_minimizers;
    buckets m_buckets;
    skew_index m_skew_index;
};

}  // namespace sshash