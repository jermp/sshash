#include "include/dictionary.hpp"
#include "include/buckets_statistics.hpp"

namespace sshash {

template <class kmer_t>
void dictionary<kmer_t>::compute_statistics() const {
    uint64_t num_kmers = size();
    uint64_t num_minimizers = m_minimizers.size();
    uint64_t num_super_kmers = m_buckets.offsets.size();

    buckets_statistics buckets_stats(num_minimizers, num_kmers, num_super_kmers);

    std::cout << "computing buckets statistics..." << std::endl;

    for (uint64_t bucket_id = 0; bucket_id != num_minimizers; ++bucket_id) {
        auto [begin, end] = m_buckets.locate_bucket(bucket_id);
        uint64_t num_super_kmers_in_bucket = end - begin;
        buckets_stats.add_num_super_kmers_in_bucket(num_super_kmers_in_bucket);
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            uint64_t offset = m_buckets.offsets.access(super_kmer_id);
            auto res = m_buckets.offset_to_id(offset, m_k);
            kmer_iterator<kmer_t> it(m_buckets.strings, m_k, 2 * offset);
            uint64_t window_size =
                std::min<uint64_t>(m_k - m_m + 1, res.contig_end(m_k) - offset - m_k + 1);
            uint64_t prev_minimizer = constants::invalid_uint64;
            uint64_t w = 0;
            for (; w != window_size; ++w) {
                auto kmer = it.get();
                uint64_t minimizer = util::compute_minimizer(kmer, m_k, m_m, m_seed);
                if (m_canonical) {
                    kmer_t kmer_rc = kmer;
                    kmer_rc.reverse_complement_inplace(m_k);
                    uint64_t minimizer_rc = util::compute_minimizer(kmer_rc, m_k, m_m, m_seed);
                    if (minimizer_rc < minimizer) minimizer = minimizer_rc;
                }
                if (prev_minimizer != constants::invalid_uint64 and minimizer != prev_minimizer) {
                    break;
                }
                prev_minimizer = minimizer;
                it.next();
            }
            buckets_stats.add_num_kmers_in_super_kmer(num_super_kmers_in_bucket, w);
        }
    }
    buckets_stats.print_full();
    std::cout << "DONE" << std::endl;
}

}  // namespace sshash
