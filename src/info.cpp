#include "include/dictionary.hpp"

namespace sshash {

double perc(uint64_t amount, uint64_t total) { return (amount * 100.0) / total; }

template <class kmer_t>
void dictionary<kmer_t>::print_space_breakdown() const {
    const uint64_t num_bytes = (num_bits() + 7) / 8;
    std::cout << "total index size: " << num_bytes << " [B] -- "
              << essentials::convert(num_bytes, essentials::MB) << " [MB]" << '\n';
    std::cout << "SPACE BREAKDOWN:\n";
    std::cout << "  minimizers: " << static_cast<double>(m_minimizers.num_bits()) / size()
              << " [bits/kmer] ("
              << static_cast<double>(m_minimizers.num_bits()) / m_minimizers.size()
              << " [bits/key]) -- " << perc(m_minimizers.num_bits(), num_bits()) << "%\n";
    std::cout << "  strings_endpoints: " << (8.0 * m_buckets.strings_endpoints.num_bytes()) / size()
              << " [bits/kmer] -- " << perc(m_buckets.strings_endpoints.num_bytes() * 8, num_bits())
              << "%\n";

    std::cout << "  offsets: " << (8.0 * m_buckets.offsets.num_bytes()) / size()
              << " [bits/kmer] -- " << perc(8 * m_buckets.offsets.num_bytes(), num_bits()) << "%\n";
    std::cout << "  offsets2: " << (8.0 * m_buckets.offsets2.num_bytes()) / size()
              << " [bits/kmer] -- " << perc(8 * m_buckets.offsets2.num_bytes(), num_bits())
              << "%\n";
    std::cout << "  offsets3: " << (8.0 * m_buckets.offsets3.num_bytes()) / size()
              << " [bits/kmer] -- " << perc(8 * m_buckets.offsets3.num_bytes(), num_bits())
              << "%\n";
    std::cout << "  start_lists_of_size: "
              << (8.0 * essentials::vec_bytes(m_buckets.start_lists_of_size)) / size()
              << " [bits/kmer] -- "
              << perc(8 * essentials::vec_bytes(m_buckets.start_lists_of_size), num_bits())
              << "%\n";

    std::cout << "  strings: " << (8.0 * m_buckets.strings.num_bytes()) / size()
              << " [bits/kmer] -- " << perc(8 * m_buckets.strings.num_bytes(), num_bits()) << "%\n";
    std::cout << "  skew_index: " << static_cast<double>(m_skew_index.num_bits()) / size()
              << " [bits/kmer] -- " << perc(m_skew_index.num_bits(), num_bits()) << "%\n";
    std::cout << "  weights: " << static_cast<double>(m_weights.num_bits()) / size()
              << " [bits/kmer] -- " << perc(m_weights.num_bits(), num_bits()) << "%\n";
    if (weighted()) m_weights.print_space_breakdown(size());
    std::cout << "  --------------\n";
    std::cout << "  total: " << static_cast<double>(num_bits()) / size() << " [bits/kmer]"
              << std::endl;
}

template <class kmer_t>
void dictionary<kmer_t>::print_info() const {
    std::cout << "=== dictionary info:\n";
    std::cout << "version number = " << m_vnum.to_string() << '\n';
    std::cout << "num_kmers = " << size() << '\n';
    std::cout << "k = " << k() << '\n';
    std::cout << "num_minimizers = " << m_minimizers.size() << std::endl;
    std::cout << "m = " << m() << '\n';
    std::cout << "canonical = " << (canonical() ? "true" : "false") << '\n';
    std::cout << "weighted = " << (weighted() ? "true" : "false") << '\n';
    print_space_breakdown();
}

}  // namespace sshash
