#include "include/dictionary.hpp"

namespace sshash {

template <class kmer_t>
lookup_result dictionary<kmer_t>::lookup_uint_regular(kmer_t uint_kmer) const {
    uint64_t minimizer = util::compute_minimizer(uint_kmer, m_k, m_m, m_seed);
    uint64_t bucket_id = m_minimizers.lookup(minimizer);

    if (m_skew_index.empty()) return m_buckets.lookup(bucket_id, uint_kmer, m_k, m_m);

    auto [begin, end] = m_buckets.locate_bucket(bucket_id);
    const uint64_t num_super_kmers_in_bucket = end - begin;
    const uint64_t log2_bucket_size = bits::util::ceil_log2_uint32(num_super_kmers_in_bucket);
    if (log2_bucket_size > m_skew_index.min_log2) {
        uint64_t pos = m_skew_index.lookup(uint_kmer, log2_bucket_size);
        /* It must hold pos < num_super_kmers_in_bucket for the kmer to exist. */
        if (pos < num_super_kmers_in_bucket) {
            return m_buckets.lookup_in_super_kmer(begin + pos, uint_kmer, m_k, m_m);
        }
        return lookup_result();
    }

    return m_buckets.lookup(begin, end, uint_kmer, m_k, m_m);
}

template <class kmer_t>
lookup_result dictionary<kmer_t>::lookup_uint_canonical(kmer_t uint_kmer,
                                                        bool check_minimizer) const {
    kmer_t uint_kmer_rc = uint_kmer;
    uint_kmer_rc.reverse_complement_inplace(m_k);
    uint64_t minimizer = std::min(util::compute_minimizer(uint_kmer, m_k, m_m, m_seed),
                                  util::compute_minimizer(uint_kmer_rc, m_k, m_m, m_seed));
    uint64_t bucket_id = m_minimizers.lookup(minimizer);

    if (m_skew_index.empty()) {
        return m_buckets.lookup_canonical(bucket_id, uint_kmer, uint_kmer_rc, minimizer,  //
                                          m_k, m_m, m_seed, check_minimizer);
    }

    auto [begin, end] = m_buckets.locate_bucket(bucket_id);
    const uint64_t num_super_kmers_in_bucket = end - begin;
    const uint64_t log2_bucket_size = bits::util::ceil_log2_uint32(num_super_kmers_in_bucket);
    if (log2_bucket_size > m_skew_index.min_log2) {
        auto uint_kmer_canon = std::min(uint_kmer, uint_kmer_rc);
        uint64_t pos = m_skew_index.lookup(uint_kmer_canon, log2_bucket_size);
        if (pos < num_super_kmers_in_bucket) {
            auto res = m_buckets.lookup_canonical_in_super_kmer(begin + pos, uint_kmer,
                                                                uint_kmer_rc, m_k, m_m);
            if (res.kmer_id != constants::invalid_uint64) return res;
        }
        return lookup_result();
    }

    return m_buckets.lookup_canonical(begin, end, uint_kmer, uint_kmer_rc, minimizer,  //
                                      m_k, m_m, m_seed, check_minimizer);
}

template <class kmer_t>
uint64_t dictionary<kmer_t>::lookup(char const* string_kmer, bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(string_kmer, m_k);
    return lookup_uint(uint_kmer, check_reverse_complement);
}
template <class kmer_t>
uint64_t dictionary<kmer_t>::lookup_uint(kmer_t uint_kmer, bool check_reverse_complement) const {
    auto res = lookup_advanced_uint(uint_kmer, check_reverse_complement);
    return res.kmer_id;
}

template <class kmer_t>
lookup_result dictionary<kmer_t>::lookup_advanced(char const* string_kmer,
                                                  bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(string_kmer, m_k);
    return lookup_advanced_uint(uint_kmer, check_reverse_complement);
}
template <class kmer_t>
lookup_result dictionary<kmer_t>::lookup_advanced_uint(kmer_t uint_kmer,
                                                       bool check_reverse_complement) const {
    if (m_canonical) return lookup_uint_canonical(uint_kmer);
    auto res = lookup_uint_regular(uint_kmer);
    assert(res.kmer_orientation == constants::forward_orientation);
    if (check_reverse_complement and res.kmer_id == constants::invalid_uint64) {
        kmer_t uint_kmer_rc = uint_kmer;
        uint_kmer_rc.reverse_complement_inplace(m_k);
        res = lookup_uint_regular(uint_kmer_rc);
        res.kmer_orientation = constants::backward_orientation;
    }
    return res;
}

template <class kmer_t>
bool dictionary<kmer_t>::is_member(char const* string_kmer, bool check_reverse_complement) const {
    return lookup(string_kmer, check_reverse_complement) != constants::invalid_uint64;
}
template <class kmer_t>
bool dictionary<kmer_t>::is_member_uint(kmer_t uint_kmer, bool check_reverse_complement) const {
    return lookup_uint(uint_kmer, check_reverse_complement) != constants::invalid_uint64;
}

template <class kmer_t>
void dictionary<kmer_t>::access(uint64_t kmer_id, char* string_kmer) const {
    assert(kmer_id < size());
    m_buckets.access(kmer_id, string_kmer, m_k);
}

template <class kmer_t>
uint64_t dictionary<kmer_t>::weight(uint64_t kmer_id) const {
    assert(kmer_id < size());
    return m_weights.weight(kmer_id);
}

template <class kmer_t>
uint64_t dictionary<kmer_t>::contig_size(uint64_t contig_id) const {
    assert(contig_id < num_contigs());
    auto [begin, end] = m_buckets.contig_offsets(contig_id);
    uint64_t contig_length = end - begin;
    assert(contig_length >= m_k);
    return contig_length - m_k + 1;
}

template <class kmer_t>
void dictionary<kmer_t>::forward_neighbours(kmer_t suffix, neighbourhood<kmer_t>& res,
                                            bool check_reverse_complement) const {
    for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
        kmer_t new_kmer = suffix;
        new_kmer.set(m_k - 1, kmer_t::char_to_uint(kmer_t::alphabet[i]));
        res.forward[i] = lookup_advanced_uint(new_kmer, check_reverse_complement);
    }
}
template <class kmer_t>
void dictionary<kmer_t>::backward_neighbours(kmer_t prefix, neighbourhood<kmer_t>& res,
                                             bool check_reverse_complement) const {
    for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
        kmer_t new_kmer = prefix;
        new_kmer.set(0, kmer_t::char_to_uint(kmer_t::alphabet[i]));
        res.backward[i] = lookup_advanced_uint(new_kmer, check_reverse_complement);
    }
}

template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::kmer_forward_neighbours(
    char const* string_kmer, bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(string_kmer, m_k);
    return kmer_forward_neighbours(uint_kmer, check_reverse_complement);
}

template <class kmer_t>
kmer_t dictionary<kmer_t>::get_suffix(kmer_t kmer) const {
    kmer_t suffix = kmer;
    suffix.drop_char();
    return suffix;
}
template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::kmer_forward_neighbours(
    kmer_t uint_kmer, bool check_reverse_complement) const {
    neighbourhood<kmer_t> res;
    forward_neighbours(get_suffix(uint_kmer), res, check_reverse_complement);
    return res;
}

template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::kmer_backward_neighbours(
    char const* string_kmer, bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(string_kmer, m_k);
    return kmer_backward_neighbours(uint_kmer, check_reverse_complement);
}

template <class kmer_t>
kmer_t dictionary<kmer_t>::get_prefix(kmer_t kmer) const {
    kmer_t prefix = kmer;
    prefix.pad_char();
    prefix.take_chars(m_k);
    return prefix;
}

template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::kmer_backward_neighbours(
    kmer_t uint_kmer, bool check_reverse_complement) const {
    neighbourhood<kmer_t> res;
    backward_neighbours(get_prefix(uint_kmer), res, check_reverse_complement);
    return res;
}

template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::kmer_neighbours(char const* string_kmer,
                                                          bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(string_kmer, m_k);
    return kmer_neighbours(uint_kmer, check_reverse_complement);
}

template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::kmer_neighbours(kmer_t uint_kmer,
                                                          bool check_reverse_complement) const {
    neighbourhood<kmer_t> res;
    forward_neighbours(get_suffix(uint_kmer), res, check_reverse_complement);
    backward_neighbours(get_prefix(uint_kmer), res, check_reverse_complement);
    return res;
}

template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::contig_neighbours(uint64_t contig_id,
                                                            bool check_reverse_complement) const {
    assert(contig_id < num_contigs());
    neighbourhood<kmer_t> res;
    kmer_t suffix = m_buckets.contig_suffix(contig_id, m_k);
    forward_neighbours(suffix, res, check_reverse_complement);
    kmer_t prefix = m_buckets.contig_prefix(contig_id, m_k);
    prefix.pad_char();
    backward_neighbours(prefix, res, check_reverse_complement);
    return res;
}

template <class kmer_t>
uint64_t dictionary<kmer_t>::num_bits() const {
    return 8 * (sizeof(m_vnum) + sizeof(m_size) + sizeof(m_seed) + sizeof(m_k) + sizeof(m_m) +
                sizeof(m_canonical)) +
           m_minimizers.num_bits() + m_buckets.num_bits() + m_skew_index.num_bits() +
           m_weights.num_bits();
}

}  // namespace sshash
