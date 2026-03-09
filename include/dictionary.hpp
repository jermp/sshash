#pragma once

#include "util.hpp"
#include "spectrum_preserving_string_set.hpp"
#include "sparse_and_skew_index.hpp"
#include "weights.hpp"

namespace sshash {

template <typename Kmer, typename Offsets>
struct dictionary  //
{
    using kmer_type = Kmer;

    template <typename, typename>
    friend struct dictionary_builder;

    dictionary()
        : m_vnum(constants::current_version_number::x,  //
                 constants::current_version_number::y,  //
                 constants::current_version_number::z)
        , m_num_kmers(0)
        , m_num_strings(0)
        , m_k(0)
        , m_m(0)
        , m_canonical(false) {}

    /* Build from input file. */
    void build(std::string const& input_filename, build_configuration const& build_config);

    essentials::version_number vnum() const { return m_vnum; }
    uint64_t num_kmers() const { return m_num_kmers; }
    uint64_t num_strings() const { return m_num_strings; }
    uint64_t k() const { return m_k; }
    uint64_t m() const { return m_m; }
    bool canonical() const { return m_canonical; }
    bool weighted() const { return !m_weights.empty(); }
    hasher_type const& hasher() const { return m_hasher; }

    /* Lookup queries. */
    lookup_result lookup(char const* string_kmer, bool check_reverse_complement = true) const;
    lookup_result lookup(Kmer uint_kmer, bool check_reverse_complement = true) const;

    /* Return the number of kmers in string. Since strings do not have duplicates,
       the length of the string is always size + k - 1. */
    uint64_t string_size(uint64_t string_id) const;

    /* Navigational queries. */
    neighbourhood<Kmer> kmer_forward_neighbours(char const* string_kmer,
                                                bool check_reverse_complement = true) const;
    neighbourhood<Kmer> kmer_forward_neighbours(Kmer uint_kmer,
                                                bool check_reverse_complement = true) const;
    neighbourhood<Kmer> kmer_backward_neighbours(char const* string_kmer,
                                                 bool check_reverse_complement = true) const;
    neighbourhood<Kmer> kmer_backward_neighbours(Kmer uint_kmer,
                                                 bool check_reverse_complement = true) const;

    /* forward and backward */
    neighbourhood<Kmer> kmer_neighbours(char const* string_kmer,
                                        bool check_reverse_complement = true) const;
    neighbourhood<Kmer> kmer_neighbours(Kmer uint_kmer, bool check_reverse_complement = true) const;
    neighbourhood<Kmer> string_neighbours(uint64_t string_id,
                                          bool check_reverse_complement = true) const;

    /* Return the weight of the kmer given its id. */
    uint64_t weight(uint64_t kmer_id) const;

    /* Return the string of the kmer whose id is kmer_id. */
    void access(uint64_t kmer_id, char* string_kmer) const;

    /* Accessor for internal bit vector */
    bits::bit_vector const& strings() const { return m_spss.strings; }

    /* Membership queries. */
    bool is_member(char const* string_kmer, bool check_reverse_complement = true) const;
    bool is_member(Kmer uint_kmer, bool check_reverse_complement = true) const;

    template <typename, bool>
    friend struct streaming_query;

    streaming_query_report  //
    streaming_query_from_file(std::string const& filename, bool multiline) const;

    struct iterator {
        iterator(dictionary const* ptr, const uint64_t begin_kmer_id, const uint64_t end_kmer_id) {
            m_it = ptr->m_spss.at(begin_kmer_id, end_kmer_id);
        }

        bool has_next() const { return m_it.has_next(); }

        /* (kmer-id, encoded kmer) */
        std::pair<uint64_t, Kmer> next() { return m_it.next(); }

    private:
        typename spectrum_preserving_string_set<Kmer, Offsets>::iterator m_it;
    };

    iterator begin() const { return iterator(this, 0, num_kmers()); }

    iterator at_kmer_id(const uint64_t kmer_id) const {
        assert(kmer_id < num_kmers());
        return iterator(this, kmer_id, num_kmers());
    }

    std::pair<uint64_t, uint64_t>  // [begin, end)
    string_offsets(const uint64_t string_id) const {
        return m_spss.string_offsets(string_id);
    }

    /* Accessor for internal offsets structure */
    Offsets const& strings_offsets() const { return m_spss.strings_offsets; }

    iterator at_string_id(const uint64_t string_id) const {
        assert(string_id < num_strings());
        auto [begin, end] = string_offsets(string_id);
        uint64_t string_length = end - begin;  // in bases
        assert(string_length >= m_k);
        uint64_t string_size = string_length - m_k + 1;  // in kmers
        uint64_t begin_kmer_id = begin - string_id * (m_k - 1);
        uint64_t end_kmer_id = begin_kmer_id + string_size;
        return iterator(this, begin_kmer_id, end_kmer_id);
    }

    uint64_t num_bits() const;
    void print_info() const;
    void print_space_breakdown() const;

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
        visitor.visit(t.m_vnum);
        util::check_version_number(t.m_vnum);
        visitor.visit(t.m_num_kmers);
        visitor.visit(t.m_num_strings);
        visitor.visit(t.m_k);
        visitor.visit(t.m_m);
        visitor.visit(t.m_canonical);
        visitor.visit(t.m_hasher);
        visitor.visit(t.m_spss);
        visitor.visit(t.m_ssi);
        visitor.visit(t.m_weights);
    }

    essentials::version_number m_vnum;
    uint64_t m_num_kmers;
    uint64_t m_num_strings;
    uint16_t m_k;
    uint16_t m_m;
    bool m_canonical;
    hasher_type m_hasher;

    spectrum_preserving_string_set<Kmer, Offsets> m_spss;
    sparse_and_skew_index<Kmer> m_ssi;

    weights m_weights;

    lookup_result lookup_regular(Kmer uint_kmer) const;
    lookup_result lookup_regular(Kmer uint_kmer, minimizer_info mini_info) const;

    lookup_result lookup_canonical(Kmer uint_kmer) const;
    lookup_result lookup_canonical(Kmer uint_kmer, Kmer uint_kmer_rc,
                                   minimizer_info mini_info) const;

    void forward_neighbours(Kmer suffix, neighbourhood<Kmer>& res,
                            bool check_reverse_complement) const;
    void backward_neighbours(Kmer prefix, neighbourhood<Kmer>& res,
                             bool check_reverse_complement) const;

    Kmer get_prefix(Kmer kmer) const;
    Kmer get_suffix(Kmer kmer) const;
};

}  // namespace sshash
