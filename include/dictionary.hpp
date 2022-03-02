#pragma once

#include "util.hpp"
#include "minimizers.hpp"
#include "buckets.hpp"
#include "skew_index.hpp"
#include "abundances.hpp"

namespace sshash {

struct dictionary {
    dictionary() : m_size(0), m_seed(0), m_k(0), m_m(0), m_canonical_parsing(0) {}

    void build(std::string const& filename, build_configuration const& build_config);

    uint64_t size() const { return m_size; }
    uint64_t seed() const { return m_seed; }
    uint64_t k() const { return m_k; }
    uint64_t m() const { return m_m; }
    bool canonicalized() const { return m_canonical_parsing; }
    bool weighted() const { return !m_abundances.empty(); }

    uint64_t lookup(char const* string_kmer, bool check_reverse_complement_too = true) const {
        uint64_t uint64_kmer = util::string_to_uint64_no_reverse(string_kmer, m_k);
        return lookup_uint64(uint64_kmer, check_reverse_complement_too);
    }

    uint64_t lookup_uint64(uint64_t uint64_kmer, bool check_reverse_complement_too = true) const {
        if (m_canonical_parsing) return lookup_uint64_canonical_parsing(uint64_kmer);
        uint64_t kmer_id = lookup_uint64_regular_parsing(uint64_kmer);
        if (check_reverse_complement_too and kmer_id == constants::invalid) {
            uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, m_k);
            kmer_id = lookup_uint64_regular_parsing(uint64_kmer_rc);
        }
        return kmer_id;
    }

    void access(uint64_t kmer_id, char* string_kmer) const {
        assert(kmer_id < size());
        m_buckets.access(kmer_id, string_kmer, m_k);
    }

    uint64_t abundance(uint64_t kmer_id) const {
        assert(kmer_id < size());
        return m_abundances.abundance(kmer_id);
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
               m_minimizers.num_bits() + m_buckets.num_bits() + m_skew_index.num_bits() +
               m_abundances.num_bits();
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
        visitor.visit(m_abundances);
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
    abundances m_abundances;

    uint64_t lookup_uint64_regular_parsing(uint64_t uint64_kmer) const;
    uint64_t lookup_uint64_canonical_parsing(uint64_t uint64_kmer) const;
};

}  // namespace sshash