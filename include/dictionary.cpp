#include "dictionary.hpp"

namespace sshash {

template <bool advanced>
contig_query_result dictionary::lookup_uint64_regular_parsing(uint64_t uint64_kmer) const {
    uint64_t minimizer = util::compute_minimizer(uint64_kmer, m_k, m_m, m_seed);
    uint64_t bucket_id = m_minimizers.lookup(minimizer);

    if (m_skew_index.empty()) return m_buckets.lookup<advanced>(bucket_id, uint64_kmer, m_k, m_m);

    auto [begin, end] = m_buckets.locate_bucket(bucket_id);
    uint64_t num_super_kmers_in_bucket = end - begin;
    uint64_t log2_bucket_size = util::ceil_log2_uint32(num_super_kmers_in_bucket);
    if (log2_bucket_size > m_skew_index.min_log2) {
        uint64_t pos = m_skew_index.lookup(uint64_kmer, log2_bucket_size);
        /* It must hold pos < num_super_kmers_in_bucket for the kmer to exist. */
        if (pos < num_super_kmers_in_bucket) {
            return m_buckets.lookup_in_super_kmer<advanced>(begin + pos, uint64_kmer, m_k, m_m);
        }
        return contig_query_result();
    }

    return m_buckets.lookup<advanced>(begin, end, uint64_kmer, m_k, m_m);
}

template <bool advanced>
contig_query_result dictionary::lookup_uint64_canonical_parsing(uint64_t uint64_kmer) const {
    uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, m_k);
    uint64_t minimizer = util::compute_minimizer(uint64_kmer, m_k, m_m, m_seed);
    uint64_t minimizer_rc = util::compute_minimizer(uint64_kmer_rc, m_k, m_m, m_seed);
    uint64_t bucket_id = m_minimizers.lookup(std::min<uint64_t>(minimizer, minimizer_rc));

    if (m_skew_index.empty()) {
        return m_buckets.lookup_canonical<advanced>(bucket_id, uint64_kmer, uint64_kmer_rc, m_k,
                                                    m_m);
    }

    auto [begin, end] = m_buckets.locate_bucket(bucket_id);
    uint64_t num_super_kmers_in_bucket = end - begin;
    uint64_t log2_bucket_size = util::ceil_log2_uint32(num_super_kmers_in_bucket);
    if (log2_bucket_size > m_skew_index.min_log2) {
        uint64_t pos = m_skew_index.lookup(uint64_kmer, log2_bucket_size);
        if (pos < num_super_kmers_in_bucket) {
            auto res = m_buckets.lookup_in_super_kmer<advanced>(begin + pos, uint64_kmer, m_k, m_m);
            assert(res.kmer_orientation == constants::forward_orientation);
            if (res.kmer_id != constants::invalid_uint64) return res;
        }
        uint64_t pos_rc = m_skew_index.lookup(uint64_kmer_rc, log2_bucket_size);
        if (pos_rc < num_super_kmers_in_bucket) {
            auto res =
                m_buckets.lookup_in_super_kmer<advanced>(begin + pos_rc, uint64_kmer_rc, m_k, m_m);
            res.kmer_orientation = constants::backward_orientation;
            return res;
        }
        return contig_query_result();
    }
    return m_buckets.lookup_canonical<advanced>(begin, end, uint64_kmer, uint64_kmer_rc, m_k, m_m);
}

uint64_t dictionary::lookup(char const* string_kmer, bool check_reverse_complement_too) const {
    uint64_t uint64_kmer = util::string_to_uint64_no_reverse(string_kmer, m_k);
    return lookup_uint64(uint64_kmer, check_reverse_complement_too);
}
uint64_t dictionary::lookup_uint64(uint64_t uint64_kmer, bool check_reverse_complement_too) const {
    constexpr bool advanced = false;
    if (m_canonical_parsing) {
        auto res = lookup_uint64_canonical_parsing<advanced>(uint64_kmer);
        return res.kmer_id;
    }
    auto res = lookup_uint64_regular_parsing<advanced>(uint64_kmer);
    if (check_reverse_complement_too and res.kmer_id == constants::invalid_uint64) {
        uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, m_k);
        res = lookup_uint64_regular_parsing<advanced>(uint64_kmer_rc);
    }
    return res.kmer_id;
}

contig_query_result dictionary::lookup_advanced(char const* string_kmer,
                                                bool check_reverse_complement_too) const {
    uint64_t uint64_kmer = util::string_to_uint64_no_reverse(string_kmer, m_k);
    return lookup_advanced_uint64(uint64_kmer, check_reverse_complement_too);
}
contig_query_result dictionary::lookup_advanced_uint64(uint64_t uint64_kmer,
                                                       bool check_reverse_complement_too) const {
    constexpr bool advanced = true;
    if (m_canonical_parsing) return lookup_uint64_canonical_parsing<advanced>(uint64_kmer);
    auto res = lookup_uint64_regular_parsing<advanced>(uint64_kmer);
    assert(res.kmer_orientation == constants::forward_orientation);
    if (check_reverse_complement_too and res.kmer_id == constants::invalid_uint64) {
        uint64_t uint64_kmer_rc = util::compute_reverse_complement(uint64_kmer, m_k);
        res = lookup_uint64_regular_parsing<advanced>(uint64_kmer_rc);
        res.kmer_orientation = constants::backward_orientation;
    }
    return res;
}

bool dictionary::is_member(char const* string_kmer, bool check_reverse_complement_too) const {
    return lookup(string_kmer, check_reverse_complement_too) != constants::invalid_uint64;
}
bool dictionary::is_member_uint64(uint64_t uint64_kmer, bool check_reverse_complement_too) const {
    return lookup_uint64(uint64_kmer, check_reverse_complement_too) != constants::invalid_uint64;
}

void dictionary::access(uint64_t kmer_id, char* string_kmer) const {
    assert(kmer_id < size());
    m_buckets.access(kmer_id, string_kmer, m_k);
}

uint64_t dictionary::weight(uint64_t kmer_id) const {
    assert(kmer_id < size());
    return m_weights.weight(kmer_id);
}

uint64_t dictionary::contig_size(uint64_t contig_id) const {
    uint64_t contig_length = m_buckets.contig_length(contig_id);
    assert(contig_length >= m_k);
    return contig_length - m_k + 1;
}

// std::vector<uint32_t> contig_neighbours(uint64_t contig_id) const;

uint64_t dictionary::num_bits() const {
    return 8 * (sizeof(m_size) + sizeof(m_seed) + sizeof(m_k) + sizeof(m_m) +
                sizeof(m_canonical_parsing)) +
           m_minimizers.num_bits() + m_buckets.num_bits() + m_skew_index.num_bits() +
           m_weights.num_bits();
}

}  // namespace sshash