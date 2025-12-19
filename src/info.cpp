#include "include/dictionary.hpp"

namespace sshash {

template <typename Kmer, typename Offsets>
void dictionary<Kmer, Offsets>::print_space_breakdown() const {
    const uint64_t num_bytes = (num_bits() + 7) / 8;

    auto perc = [](uint64_t amount, uint64_t total) -> double { return (amount * 100.0) / total; };

    std::cout << "total index size: " << num_bytes << " [B] -- "
              << essentials::convert(num_bytes, essentials::MB) << " [MB]" << '\n';
    std::cout << "SPACE BREAKDOWN:\n";
    std::cout << "  mphf: " << static_cast<double>(m_ssi.codewords.mphf.num_bits()) / num_kmers()
              << " [bits/kmer] ("
              << static_cast<double>(m_ssi.codewords.mphf.num_bits()) /
                     m_ssi.codewords.mphf.num_keys()
              << " [bits/key]) -- " << perc(m_ssi.codewords.mphf.num_bits(), num_bits()) << "%\n";
    std::cout << "  strings_offsets: " << (8.0 * m_spss.strings_offsets.num_bytes()) / num_kmers()
              << " [bits/kmer] -- " << perc(m_spss.strings_offsets.num_bytes() * 8, num_bits())
              << "%\n";

    std::cout << "  control_codewords: "
              << (8.0 * m_ssi.codewords.control_codewords.num_bytes()) / num_kmers()
              << " [bits/kmer] -- "
              << perc(8 * m_ssi.codewords.control_codewords.num_bytes(), num_bits()) << "%\n";
    std::cout << "  mid_load_buckets: " << (8.0 * m_ssi.mid_load_buckets.num_bytes()) / num_kmers()
              << " [bits/kmer] -- " << perc(8 * m_ssi.mid_load_buckets.num_bytes(), num_bits())
              << "%\n";
    std::cout << "  begin_buckets_of_size: "
              << (8.0 * essentials::vec_bytes(m_ssi.begin_buckets_of_size)) / num_kmers()
              << " [bits/kmer] -- "
              << perc(8 * essentials::vec_bytes(m_ssi.begin_buckets_of_size), num_bits()) << "%\n";

    std::cout << "  strings: " << (8.0 * m_spss.strings.num_bytes()) / num_kmers()
              << " [bits/kmer] -- " << perc(8 * m_spss.strings.num_bytes(), num_bits()) << "%\n";
    std::cout << "  skew_index: " << static_cast<double>(m_ssi.ski.num_bits()) / num_kmers()
              << " [bits/kmer] -- " << perc(m_ssi.ski.num_bits(), num_bits()) << "%\n";
    std::cout << "  weights: " << static_cast<double>(m_weights.num_bits()) / num_kmers()
              << " [bits/kmer] -- " << perc(m_weights.num_bits(), num_bits()) << "%\n";

    if (weighted()) m_weights.print_space_breakdown(num_kmers());

    std::cout << "  --------------\n";
    std::cout << "  total: " << static_cast<double>(num_bits()) / num_kmers() << " [bits/kmer]"
              << std::endl;
}

template <typename Kmer, typename Offsets>
void dictionary<Kmer, Offsets>::print_info() const {
    std::cout << "=== dictionary info:\n";
    std::cout << "version number = " << m_vnum.to_string() << '\n';
    std::cout << "num_kmers = " << num_kmers() << '\n';
    std::cout << "num_strings = " << num_strings() << '\n';
    std::cout << "k = " << k() << '\n';
    std::cout << "num_minimizers = " << m_ssi.codewords.size() << std::endl;
    std::cout << "m = " << m() << '\n';
    std::cout << "canonical = " << (canonical() ? "true" : "false") << '\n';
    std::cout << "weighted = " << (weighted() ? "true" : "false") << '\n';
    print_space_breakdown();
}

}  // namespace sshash
