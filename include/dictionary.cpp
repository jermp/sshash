#include "dictionary.hpp"

namespace sshash {

lookup_result dictionary::lookup_uint_regular_parsing(kmer_t uint_kmer) const {
    uint64_t minimizer = util::compute_minimizer(uint_kmer, m_k, m_m, m_seed);
    uint64_t bucket_id = m_minimizers.lookup(minimizer);

    if (m_skew_index.empty()) return m_buckets.lookup(bucket_id, uint_kmer, m_k, m_m);

    auto [begin, end] = m_buckets.locate_bucket(bucket_id);
    uint64_t num_super_kmers_in_bucket = end - begin;
    uint64_t log2_bucket_size = util::ceil_log2_uint32(num_super_kmers_in_bucket);
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

lookup_result dictionary::lookup_uint_canonical_parsing(kmer_t uint_kmer) const {
    kmer_t uint_kmer_rc = util::compute_reverse_complement(uint_kmer, m_k);
    uint64_t minimizer = util::compute_minimizer(uint_kmer, m_k, m_m, m_seed);
    uint64_t minimizer_rc = util::compute_minimizer(uint_kmer_rc, m_k, m_m, m_seed);
    uint64_t bucket_id = m_minimizers.lookup(std::min<uint64_t>(minimizer, minimizer_rc));

    if (m_skew_index.empty()) {
        return m_buckets.lookup_canonical(bucket_id, uint_kmer, uint_kmer_rc, m_k, m_m);
    }

    auto [begin, end] = m_buckets.locate_bucket(bucket_id);
    uint64_t num_super_kmers_in_bucket = end - begin;
    uint64_t log2_bucket_size = util::ceil_log2_uint32(num_super_kmers_in_bucket);
    if (log2_bucket_size > m_skew_index.min_log2) {
        uint64_t pos = m_skew_index.lookup(uint_kmer, log2_bucket_size);
        if (pos < num_super_kmers_in_bucket) {
            auto res = m_buckets.lookup_in_super_kmer(begin + pos, uint_kmer, m_k, m_m);
            assert(res.kmer_orientation == constants::forward_orientation);
            if (res.kmer_id != constants::invalid_uint64) return res;
        }
        uint64_t pos_rc = m_skew_index.lookup(uint_kmer_rc, log2_bucket_size);
        if (pos_rc < num_super_kmers_in_bucket) {
            auto res = m_buckets.lookup_in_super_kmer(begin + pos_rc, uint_kmer_rc, m_k, m_m);
            res.kmer_orientation = constants::backward_orientation;
            return res;
        }
        return lookup_result();
    }

    return m_buckets.lookup_canonical(begin, end, uint_kmer, uint_kmer_rc, m_k, m_m);
}

uint64_t dictionary::lookup(char const* string_kmer, bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer(string_kmer, m_k);
    return lookup_uint(uint_kmer, check_reverse_complement);
}
uint64_t dictionary::lookup_uint(kmer_t uint_kmer, bool check_reverse_complement) const {
    auto res = lookup_advanced_uint(uint_kmer, check_reverse_complement);
    return res.kmer_id;
}

lookup_result dictionary::lookup_advanced(char const* string_kmer,
                                          bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer(string_kmer, m_k);
    return lookup_advanced_uint(uint_kmer, check_reverse_complement);
}
lookup_result dictionary::lookup_advanced_uint(kmer_t uint_kmer,
                                               bool check_reverse_complement) const {
    if (m_canonical_parsing) return lookup_uint_canonical_parsing(uint_kmer);
    auto res = lookup_uint_regular_parsing(uint_kmer);
    assert(res.kmer_orientation == constants::forward_orientation);
    if (check_reverse_complement and res.kmer_id == constants::invalid_uint64) {
        kmer_t uint_kmer_rc = util::compute_reverse_complement(uint_kmer, m_k);
        res = lookup_uint_regular_parsing(uint_kmer_rc);
        res.kmer_orientation = constants::backward_orientation;
    }
    return res;
}

bool dictionary::is_member(char const* string_kmer, bool check_reverse_complement) const {
    return lookup(string_kmer, check_reverse_complement) != constants::invalid_uint64;
}
bool dictionary::is_member_uint(kmer_t uint_kmer, bool check_reverse_complement) const {
    return lookup_uint(uint_kmer, check_reverse_complement) != constants::invalid_uint64;
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
    assert(contig_id < num_contigs());
    auto [begin, end] = m_buckets.contig_offsets(contig_id);
    uint64_t contig_length = end - begin;
    assert(contig_length >= m_k);
    return contig_length - m_k + 1;
}

void dictionary::forward_neighbours(kmer_t suffix, neighbourhood& res,
                                    bool check_reverse_complement) const {
    res.forward_A = lookup_advanced_uint(suffix + (util::char_to_uint('A') << (2 * (m_k - 1))),
                                         check_reverse_complement);
    res.forward_C = lookup_advanced_uint(suffix + (util::char_to_uint('C') << (2 * (m_k - 1))),
                                         check_reverse_complement);
    res.forward_G = lookup_advanced_uint(suffix + (util::char_to_uint('G') << (2 * (m_k - 1))),
                                         check_reverse_complement);
    res.forward_T = lookup_advanced_uint(suffix + (util::char_to_uint('T') << (2 * (m_k - 1))),
                                         check_reverse_complement);
}
void dictionary::backward_neighbours(kmer_t prefix, neighbourhood& res,
                                     bool check_reverse_complement) const {
    res.backward_A =
        lookup_advanced_uint(prefix + util::char_to_uint('A'), check_reverse_complement);
    res.backward_C =
        lookup_advanced_uint(prefix + util::char_to_uint('C'), check_reverse_complement);
    res.backward_G =
        lookup_advanced_uint(prefix + util::char_to_uint('G'), check_reverse_complement);
    res.backward_T =
        lookup_advanced_uint(prefix + util::char_to_uint('T'), check_reverse_complement);
}

neighbourhood dictionary::kmer_forward_neighbours(char const* string_kmer,
                                                  bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer(string_kmer, m_k);
    return kmer_forward_neighbours(uint_kmer, check_reverse_complement);
}
neighbourhood dictionary::kmer_forward_neighbours(kmer_t uint_kmer,
                                                  bool check_reverse_complement) const {
    neighbourhood res;
    kmer_t suffix = uint_kmer >> 2;
    forward_neighbours(suffix, res, check_reverse_complement);
    return res;
}

neighbourhood dictionary::kmer_backward_neighbours(char const* string_kmer,
                                                   bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer(string_kmer, m_k);
    return kmer_backward_neighbours(uint_kmer, check_reverse_complement);
}
neighbourhood dictionary::kmer_backward_neighbours(kmer_t uint_kmer,
                                                   bool check_reverse_complement) const {
    neighbourhood res;
    kmer_t prefix = (uint_kmer << 2) & ((kmer_t(1) << (2 * m_k)) - 1);
    backward_neighbours(prefix, res, check_reverse_complement);
    return res;
}

neighbourhood dictionary::kmer_neighbours(char const* string_kmer,
                                          bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer(string_kmer, m_k);
    return kmer_neighbours(uint_kmer, check_reverse_complement);
}
neighbourhood dictionary::kmer_neighbours(kmer_t uint_kmer, bool check_reverse_complement) const {
    neighbourhood res;
    kmer_t suffix = uint_kmer >> 2;
    forward_neighbours(suffix, res, check_reverse_complement);
    kmer_t prefix = (uint_kmer << 2) & ((kmer_t(1) << (2 * m_k)) - 1);
    backward_neighbours(prefix, res, check_reverse_complement);
    return res;
}

neighbourhood dictionary::contig_neighbours(uint64_t contig_id,
                                            bool check_reverse_complement) const {
    assert(contig_id < num_contigs());
    neighbourhood res;
    kmer_t suffix = m_buckets.contig_suffix(contig_id, m_k);
    forward_neighbours(suffix, res, check_reverse_complement);
    kmer_t prefix = m_buckets.contig_prefix(contig_id, m_k) << 2;
    backward_neighbours(prefix, res, check_reverse_complement);
    return res;
}

uint64_t dictionary::num_bits() const {
    return 8 * (sizeof(m_size) + sizeof(m_seed) + sizeof(m_k) + sizeof(m_m) +
                sizeof(m_canonical_parsing)) +
           m_minimizers.num_bits() + m_buckets.num_bits() + m_skew_index.num_bits() +
           m_weights.num_bits();
}

}  // namespace sshash