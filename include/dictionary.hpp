#pragma once

#include "util.hpp"
#include "minimizers.hpp"
#include "buckets.hpp"
#include "skew_index.hpp"
#include "weights.hpp"

namespace sshash {

struct dictionary {
    dictionary() : m_size(0), m_seed(0), m_k(0), m_m(0), m_canonical_parsing(0) {}

    /* Build from input file. */
    void build(std::string const& input_filename, build_configuration const& build_config);

    /* Write super-k-mers to output file in FASTA format. */
    void dump(std::string const& output_filename) const;

    uint64_t size() const { return m_size; }
    uint64_t seed() const { return m_seed; }
    uint64_t k() const { return m_k; }
    uint64_t m() const { return m_m; }
    uint64_t num_contigs() const { return m_buckets.pieces.size() - 1; }
    bool canonicalized() const { return m_canonical_parsing; }
    bool weighted() const { return !m_weights.empty(); }

    /* Lookup queries. Return the kmer_id of the kmer or -1 if it is not found in the dictionary. */
    uint64_t lookup(char const* string_kmer, bool check_reverse_complement = true) const;
    uint64_t lookup_uint(kmer_t uint_kmer, bool check_reverse_complement = true) const;

    /* Advanced lookup queries. Return also contig information. */
    lookup_result lookup_advanced(char const* string_kmer,
                                  bool check_reverse_complement = true) const;
    lookup_result lookup_advanced_uint(kmer_t uint_kmer,
                                       bool check_reverse_complement = true) const;

    /* Return the number of kmers in contig. Since contigs do not have duplicates,
       the length of the contig is always size + k - 1. */
    uint64_t contig_size(uint64_t contig_id) const;

    /* Navigational queries. */
    neighbourhood kmer_forward_neighbours(char const* string_kmer,
                                          bool check_reverse_complement = true) const;
    neighbourhood kmer_forward_neighbours(kmer_t uint_kmer,
                                          bool check_reverse_complement = true) const;
    neighbourhood kmer_backward_neighbours(char const* string_kmer,
                                           bool check_reverse_complement = true) const;
    neighbourhood kmer_backward_neighbours(kmer_t uint_kmer,
                                           bool check_reverse_complement = true) const;

    /* forward and backward */
    neighbourhood kmer_neighbours(char const* string_kmer,
                                  bool check_reverse_complement = true) const;
    neighbourhood kmer_neighbours(kmer_t uint_kmer, bool check_reverse_complement = true) const;
    neighbourhood contig_neighbours(uint64_t contig_id, bool check_reverse_complement = true) const;

    /* Return the weight of the kmer given its id. */
    uint64_t weight(uint64_t kmer_id) const;

    /* Return the string of the kmer whose id is kmer_id. */
    void access(uint64_t kmer_id, char* string_kmer) const;

    /* Membership queries. */
    bool is_member(char const* string_kmer, bool check_reverse_complement = true) const;
    bool is_member_uint(kmer_t uint_kmer, bool check_reverse_complement = true) const;

    /* Streaming queries. */
    friend struct streaming_query_canonical_parsing;
    friend struct streaming_query_regular_parsing;
    streaming_query_report streaming_query_from_file(std::string const& filename,
                                                     bool multiline) const;

    struct iterator {
        iterator(dictionary const* ptr, const uint64_t begin_kmer_id, const uint64_t end_kmer_id) {
            it = ptr->m_buckets.at(begin_kmer_id, end_kmer_id, ptr->m_k);
        }

        bool has_next() const { return it.has_next(); }

        std::pair<uint64_t, std::string>  // (kmer-id, kmer)
        next() {
            return it.next();
        }

    private:
        typename buckets::iterator it;
    };

    iterator begin() const { return iterator(this, 0, size()); }

    iterator at_kmer_id(const uint64_t kmer_id) const {
        assert(kmer_id < size());
        return iterator(this, kmer_id, size());
    }

    std::pair<uint64_t, uint64_t>  // [begin, end)
    contig_offsets(const uint64_t contig_id) const {
        return m_buckets.contig_offsets(contig_id);
    }

    iterator at_contig_id(const uint64_t contig_id) const {
        assert(contig_id < num_contigs());
        auto [begin, end] = contig_offsets(contig_id);
        uint64_t contig_length = end - begin;  // in bases
        assert(contig_length >= m_k);
        uint64_t contig_size = contig_length - m_k + 1;  // in kmers
        uint64_t begin_kmer_id = begin - contig_id * (m_k - 1);
        uint64_t end_kmer_id = begin_kmer_id + contig_size;
        return iterator(this, begin_kmer_id, end_kmer_id);
    }

    pthash::bit_vector const& strings() const { return m_buckets.strings; }

    uint64_t num_bits() const;
    void print_info() const;
    void print_space_breakdown() const;
    void compute_statistics() const;

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_size);
        visitor.visit(t.m_seed);
        visitor.visit(t.m_k);
        visitor.visit(t.m_m);
        visitor.visit(t.m_canonical_parsing);
        visitor.visit(t.m_minimizers);
        visitor.visit(t.m_buckets);
        visitor.visit(t.m_skew_index);
        visitor.visit(t.m_weights);
    }

    uint64_t m_size;
    uint64_t m_seed;
    uint16_t m_k;
    uint16_t m_m;
    uint16_t m_canonical_parsing;
    minimizers m_minimizers;
    buckets m_buckets;
    skew_index m_skew_index;
    weights m_weights;

    lookup_result lookup_uint_regular_parsing(kmer_t uint_kmer) const;
    lookup_result lookup_uint_canonical_parsing(kmer_t uint_kmer) const;
    void forward_neighbours(kmer_t suffix, neighbourhood& res, bool check_reverse_complement) const;
    void backward_neighbours(kmer_t prefix, neighbourhood& res,
                             bool check_reverse_complement) const;
};

}  // namespace sshash
