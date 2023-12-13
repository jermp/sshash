#pragma once

#include "util.hpp"
#include "minimizers.hpp"
#include "buckets.hpp"
#include "skew_index.hpp"
#include "weights.hpp"

namespace sshash {

// Forward declarations of the friend template classes
template <class kmer_t>
struct streaming_query_canonical_parsing;

template <class kmer_t>
struct streaming_query_regular_parsing;

template <class kmer_t = default_kmer_t>
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
    neighbourhood kmer_forward_neighbours(char const* string_kmer) const;
    neighbourhood kmer_forward_neighbours(kmer_t uint_kmer) const;
    neighbourhood kmer_backward_neighbours(char const* string_kmer) const;
    neighbourhood kmer_backward_neighbours(kmer_t uint_kmer) const;

    /* forward and backward */
    neighbourhood kmer_neighbours(char const* string_kmer) const;
    neighbourhood kmer_neighbours(kmer_t uint_kmer) const;
    neighbourhood contig_neighbours(uint64_t contig_id) const;

    /* Return the weight of the kmer given its id. */
    uint64_t weight(uint64_t kmer_id) const;

    /* Return the string of the kmer whose id is kmer_id. */
    void access(uint64_t kmer_id, char* string_kmer) const;

    /* Membership queries. */
    bool is_member(char const* string_kmer, bool check_reverse_complement = true) const;
    bool is_member_uint(kmer_t uint_kmer, bool check_reverse_complement = true) const;

    /* Streaming queries. */
    friend struct streaming_query_canonical_parsing<kmer_t>;
    friend struct streaming_query_regular_parsing<kmer_t>;
    streaming_query_report streaming_query_from_file(std::string const& filename,
                                                     bool multiline) const;

    struct iterator {
        iterator(dictionary const* ptr, uint64_t kmer_id = 0) {
            it = ptr->m_buckets.at(kmer_id, ptr->m_k, ptr->m_size);
        }

        bool has_next() const { return it.has_next(); }
        std::pair<uint64_t, std::string> next() { return it.next(); }

    private:
        typename buckets<kmer_t>::iterator it;
    };

    iterator begin() const { return iterator(this); }

    iterator at(uint64_t kmer_id) const {
        assert(kmer_id < size());
        return iterator(this, kmer_id);
    }

    pthash::bit_vector const& strings() const { return m_buckets.strings; }

    uint64_t num_bits() const;
    void print_info() const;
    void print_space_breakdown() const;
    void compute_statistics() const;

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
        visitor.visit(m_weights);
    }

private:
    uint64_t m_size;
    uint64_t m_seed;
    uint16_t m_k;
    uint16_t m_m;
    uint16_t m_canonical_parsing;
    minimizers m_minimizers;
    buckets<kmer_t> m_buckets;
    skew_index<kmer_t> m_skew_index;
    weights m_weights;

    lookup_result lookup_uint_regular_parsing(kmer_t uint_kmer) const;
    lookup_result lookup_uint_canonical_parsing(kmer_t uint_kmer) const;
    void forward_neighbours(kmer_t suffix, neighbourhood& res) const;
    void backward_neighbours(kmer_t prefix, neighbourhood& res) const;
};

}  // namespace sshash

#include "builder/build.impl"
#include "dictionary.impl"
#include "dump.impl"
#include "info.impl"
#include "statistics.impl"