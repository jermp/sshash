#include "dictionary.hpp"

#include <unordered_set>

namespace sshash {

void dictionary::dump(std::string const& filename) const {
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

    std::vector<std::pair<std::string, uint32_t>> super_kmers;  // (super-kmer,minimizer-position)
    std::string super_kmer;
    uint32_t p;
    uint64_t resolved_kmers_in_large_buckets = 0;  // "large" means > 1 super-kmer
    uint64_t total_kmers_in_large_buckets = 0;

    for (uint64_t bucket_id = 0; bucket_id != num_minimizers; ++bucket_id) {
        auto [begin, end] = m_buckets.locate_bucket(bucket_id);
        uint64_t num_super_kmers = end - begin;
        if (num_super_kmers == 1) continue;
        assert(num_super_kmers > 1);
        super_kmers.clear();
        super_kmers.reserve(num_super_kmers);
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            uint64_t offset = m_buckets.offsets.access(super_kmer_id);
            auto [_, contig_end] = m_buckets.offset_to_id(offset, m_k);
            (void)_;
            bit_vector_iterator bv_it(m_buckets.strings, 2 * offset);
            uint64_t window_size = std::min<uint64_t>(m_k - m_m + 1, contig_end - offset - m_k + 1);
            uint64_t prev_minimizer = constants::invalid_uint64;
            bool super_kmer_header_written = false;
            for (uint64_t w = 0; w != window_size; ++w) {
                uint64_t kmer = bv_it.read_and_advance_by_two(2 * m_k);
                auto [minimizer, pos] = util::compute_minimizer_pos(kmer, m_k, m_m, m_seed);
                if (m_canonical_parsing) {
                    uint64_t kmer_rc = util::compute_reverse_complement(kmer, m_k);
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
                    super_kmer = util::uint_kmer_to_string_no_reverse(kmer, m_k);
                    p = pos;
                    out << '>' << bucket_id << ':' << super_kmer_id - begin << ':'
                        << util::uint_kmer_to_string_no_reverse(minimizer, m_m) << ':' << pos
                        << '\n';
                    out << super_kmer;
                    super_kmer_header_written = true;
                } else {
                    if (minimizer != prev_minimizer) {
                        break;
                    } else {
                        char c = util::uint64_to_char(kmer >> (2 * (m_k - 1)));
                        super_kmer.push_back(c);
                        out << c;
                    }
                }
                prev_minimizer = minimizer;
            }
            super_kmers.push_back({super_kmer, p});
            out << '\n';
        }

        // std::cout << "super_kmers are:\n";
        std::sort(super_kmers.begin(), super_kmers.end(),
                  [](auto const& x, auto const& y) { return x.second < y.second; });
        assert(super_kmers.size() == num_super_kmers);
        // uint32_t max_pos = super_kmers.back().second;
        for (auto const& p : super_kmers) {
            auto const& s = p.first;
            total_kmers_in_large_buckets += s.length() - m_k + 1;
            // uint32_t shift = max_pos - p.second;
            // std::cout << std::string(shift, ' ');
            // for (uint32_t i = 0; i != s.length(); ++i) {
            //     if (i == p.second or i == p.second + m_m) std::cout << '-';
            //     std::cout << s[i];
            // }
            // std::cout << '\n';
        }

        if (num_super_kmers > 16) {
            /* cannot disambiguate more than 16 super-kmers with a sketch */
            continue;
        }

        /****/
        std::unordered_set<uint16_t> sketches;
        sketches.reserve(num_super_kmers);

        /*
            Example with 2 super-kmers.

                    TACGTGTT-TGGCACTGACCTG-TAGGTCTGATAAGACGCG
            TTGAGCAGGGTATCGA-TGGCACTGACCTG-AAAGCCGG
                      ...321 0           0 123...
                    ----------***********----------
                     negative             positive

            In this case, the result is i = -0 because
            'TT' from first super-kmer and 'AT' from the second, can distinguish them.

            Example with 4 super-kmers.

                              -CTGGCGCTGGTGA-GCCGTAAGAAAAAAGAAG
            GATGCGCTGCCAGAAGCG-CTGGCGCTGGTGA
            CTGATGAATGGGCGTCCG-CTGGCGCTGGTGA
            GTGGTCTTCGGTTCGCTG-CTGGCGCTGGTGA-TTTATCTCTTCTGGCTCT

            In this case, the result is i = -2, with
            ??, GC, CC, CT. The sketch ?? indicates a "don't-care" value that ca be chosen
            as we want but different from {GC, CC, CT}.
        */

        /* search backward */
        uint32_t i = 0;
        uint32_t num_dont_cares = 0;
        // uint32_t N = 1;  // for tri-grams
        uint32_t N = 0;  // for bi-grams
        while (true) {
            for (auto const& p : super_kmers) {
                if (num_dont_cares == num_super_kmers) break;
                if (p.second <= i + N) {
                    num_dont_cares += 1;
                    continue;
                }
                uint32_t pos = p.second - i;
                auto const& s = p.first;

                /* bigram */
                uint16_t sketch =
                    (util::char_to_uint(s[pos - 1]) << 2) + (util::char_to_uint(s[pos]));
                // /* trigram */
                // uint16_t sketch = (util::char_to_uint(s[pos - 2]) << 4) +
                //                   (util::char_to_uint(s[pos - 1]) << 2) +
                //                   (util::char_to_uint(s[pos]));

                // std::cout << s[pos - 1] << s[pos] << "(" << sketch << ") ";
                if (auto it = sketches.find(sketch); it == sketches.cend()) {
                    sketches.insert(sketch);
                } else {  // collision
                    num_dont_cares = 0;
                    break;
                }
            }
            if (num_dont_cares == num_super_kmers) break;
            if (sketches.size() + num_dont_cares == num_super_kmers) break;
            sketches.clear();
            i += 1;
        }

        if (num_dont_cares == num_super_kmers) {
            /* search forward */
            sketches.clear();
            i = 0;
            num_dont_cares = 0;
            while (true) {
                for (auto const& p : super_kmers) {
                    uint32_t pos = p.second + (m_m - 1) + i;
                    auto const& s = p.first;
                    if (pos >= s.length() - 1 + N) {
                        num_dont_cares += 1;
                        continue;
                    }
                    /* bigram */
                    uint16_t sketch =
                        (util::char_to_uint(s[pos]) << 2) + util::char_to_uint(s[pos + 1]);
                    // /* trigram */
                    // uint16_t sketch = (util::char_to_uint(s[pos]) << 4) +
                    //                   (util::char_to_uint(s[pos + 1]) << 2) +
                    //                   (util::char_to_uint(s[pos + 2]));

                    // std::cout << s[pos] << s[pos + 1] << "(" << sketch << ") ";
                    if (auto it = sketches.find(sketch); it == sketches.cend()) {
                        sketches.insert(sketch);
                    } else {  // collision
                        num_dont_cares = 0;
                        break;
                    }
                }
                if (num_dont_cares == num_super_kmers) break;
                if (sketches.size() + num_dont_cares == num_super_kmers) break;
                sketches.clear();
                i += 1;
            }

            if (num_dont_cares == num_super_kmers) {
                std::cout << " == sketch NOT FOUND" << std::endl;
                uint32_t max_pos = super_kmers.back().second;
                for (auto const& p : super_kmers) {
                    auto const& s = p.first;
                    uint32_t shift = max_pos - p.second;
                    std::cout << std::string(shift, ' ');
                    for (uint32_t i = 0; i != s.length(); ++i) {
                        if (i == p.second or i == p.second + m_m) std::cout << '-';
                        std::cout << s[i];
                    }
                    std::cout << '\n';
                }
                std::cout << std::endl;
            } else {
                // if (num_dont_cares > 0) {
                //     std::cout << "i = +" << i << " (with DON'T-CARE) " << std::endl << std::endl;
                // } else {
                //     std::cout << "i = +" << i << std::endl << std::endl;
                // }
                for (auto const& p : super_kmers) {
                    resolved_kmers_in_large_buckets += p.first.length() - m_k + 1;
                }
            }

        } else {
            // if (num_dont_cares > 0) {
            //     std::cout << "i = -" << i << " (with DON'T-CARE) " << std::endl << std::endl;
            // } else {
            //     std::cout << "i = -" << i << std::endl << std::endl;
            // }
            for (auto const& p : super_kmers) {
                resolved_kmers_in_large_buckets += p.first.length() - m_k + 1;
            }
        }
    }

    out.close();
    std::cout << "DONE" << std::endl;
    std::cerr << "kmers resolved by sketching with don't-cares in buckets of size > 1: "
              << (resolved_kmers_in_large_buckets * 100.0) / total_kmers_in_large_buckets << "%"
              << std::endl;
    std::cerr << "kmers in fall-back mphf: "
              << (total_kmers_in_large_buckets - resolved_kmers_in_large_buckets) * 100.0 /
                     num_kmers
              << "%" << std::endl;
    std::cerr << "fall-back mphf: "
              << (total_kmers_in_large_buckets - resolved_kmers_in_large_buckets) * 3.0 / num_kmers
              << " bits/kmer" << std::endl;
}

}  // namespace sshash
