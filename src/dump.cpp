#include "include/dictionary.hpp"

namespace sshash {

namespace util {

template <class kmer_t, typename Hasher = murmurhash2_64>
std::pair<kmer_t, uint64_t> compute_minimizer_pos(kmer_t kmer, uint64_t k, uint64_t m,
                                                  uint64_t seed) {
    assert(m <= kmer_t::max_m);
    assert(m <= k);
    uint64_t min_hash = uint64_t(-1);
    kmer_t minimizer = kmer_t(-1);
    uint64_t pos = 0;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        kmer_t mmer = kmer;
        mmer.take_chars(m);
        uint64_t hash = Hasher::hash(uint64_t(mmer), seed);
        if (hash < min_hash) {
            min_hash = hash;
            minimizer = mmer;
            pos = i;
        }
        kmer.drop_char();
    }
    return {minimizer, pos};
}

}  // namespace util

template <class kmer_t>
void dictionary<kmer_t>::dump(std::string const& filename) const {
    uint64_t num_kmers = size();
    uint64_t num_minimizers = m_minimizers.size();
    uint64_t num_super_kmers = m_buckets.offsets.size();

    std::ofstream out(filename);
    std::cout << "dumping super-k-mers to file '" << filename << "'..." << std::endl;

    /*
        Write header and a dummy empty line "N".
        Header is:
        [k]:[m]:[num_kmers]:[num_minimizers]:[num_super_kmers]
    */
    out << '>' << m_k << ':' << m_m << ':' << num_kmers << ':' << num_minimizers << ':'
        << num_super_kmers << "\nN\n";

    for (uint64_t bucket_id = 0; bucket_id != num_minimizers; ++bucket_id) {
        auto [begin, end] = m_buckets.locate_bucket(bucket_id);
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            uint64_t offset = m_buckets.offsets.access(super_kmer_id);
            auto res = m_buckets.offset_to_id(offset, m_k);
            kmer_iterator<kmer_t> it(m_buckets.strings, m_k, 2 * offset);
            uint64_t window_size =
                std::min<uint64_t>(m_k - m_m + 1, res.contig_end(m_k) - offset - m_k + 1);
            kmer_t prev_minimizer = constants::invalid_uint64;
            bool super_kmer_header_written = false;
            for (uint64_t w = 0; w != window_size; ++w) {
                auto kmer = it.get();
                auto [minimizer, pos] = util::compute_minimizer_pos(kmer, m_k, m_m, m_seed);
                if (m_canonical) {
                    kmer_t kmer_rc = kmer;
                    kmer_rc.reverse_complement_inplace(m_k);
                    auto [minimizer_rc, pos_rc] =
                        util::compute_minimizer_pos(kmer_rc, m_k, m_m, m_seed);
                    if (minimizer_rc < minimizer) {
                        minimizer = minimizer_rc;
                        pos = pos_rc;
                    }
                }
                if (!super_kmer_header_written) {
                    /*
                        Write super-kmer header:
                        [minimizer_id]:[super_kmer_id]:[minimizer_string]:[position_of_minimizer_in_super_kmer]
                    */
                    out << '>' << bucket_id << ':' << super_kmer_id - begin << ':'
                        << util::uint_kmer_to_string(minimizer, m_m) << ':' << pos << '\n';
                    out << util::uint_kmer_to_string(kmer, m_k);
                    super_kmer_header_written = true;
                } else {
                    if (minimizer != prev_minimizer) {
                        break;
                    } else {
                        kmer_t last_char = kmer;
                        last_char.drop_chars(m_k - 1);
                        out << kmer_t::uint64_to_char(last_char.pop_char());
                    }
                }
                prev_minimizer = minimizer;
                it.next();
            }
            out << '\n';
        }
    }

    out.close();
    std::cout << "DONE" << std::endl;
}

}  // namespace sshash