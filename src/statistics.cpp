#include "include/dictionary.hpp"
#include "include/buckets_statistics.hpp"
#include "include/minimizer_iterator.hpp"

namespace sshash {

template <class kmer_t>
void dictionary<kmer_t>::compute_statistics() const  //
{
    std::cout << "computing bucket statistics..." << std::endl;

    const uint64_t num_kmers = size();
    const uint64_t num_minimizers = m_minimizers.size();
    const uint64_t num_super_kmers = m_buckets.offsets.size();

    buckets_statistics buckets_stats(num_minimizers, num_kmers, num_super_kmers);
    minimizer_iterator<kmer_t> minimizer_it(m_k, m_m, m_hasher);

    for (uint64_t bucket_id = 0; bucket_id != num_minimizers; ++bucket_id) {
        const auto [begin, end] = m_buckets.locate_bucket(bucket_id);
        const uint64_t bucket_size = end - begin;
        buckets_stats.add_bucket_size(bucket_size);
        for (uint64_t i = begin; i != end; ++i) {
            const uint64_t pos_in_seq = m_buckets.offsets.access(i);
            auto p = m_buckets.pieces.locate(pos_in_seq);
            const uint64_t contig_begin = p.first.val;
            const uint64_t contig_end = p.second.val;
            uint64_t offset = pos_in_seq;
            if (offset <= m_k - m_m) {
                assert(contig_begin == 0);
                offset = 0;
            } else if (offset - (m_k - m_m) < contig_begin) {
                offset = contig_begin;
            } else {
                offset -= m_k - m_m;
            }
            kmer_iterator<kmer_t> it(m_buckets.strings, m_k, kmer_t::bits_per_char * offset);
            minimizer_it.set_position(offset);
            uint64_t num_kmers_in_super_kmer = 0;
            auto kmer = it.get();
            auto mini_info = minimizer_it.next(kmer);
            while (mini_info.pos_in_seq < pos_in_seq) {
                it.next();
                kmer = it.get();
                mini_info = minimizer_it.next(kmer);
            }
            while (mini_info.pos_in_seq == pos_in_seq and
                   (mini_info.pos_in_seq - mini_info.pos_in_kmer + m_k) <= contig_end)  //
            {
                num_kmers_in_super_kmer += 1;
                it.next();
                kmer = it.get();
                mini_info = minimizer_it.next(kmer);
            }
            assert(num_kmers_in_super_kmer > 0);
            buckets_stats.add_num_kmers_in_super_kmer(bucket_size, num_kmers_in_super_kmer);
        }
    }

    buckets_stats.print_full();

    std::cout << "DONE" << std::endl;
}

}  // namespace sshash
