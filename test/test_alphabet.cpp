#include "../include/util.hpp"

using namespace sshash;

std::ostream& operator<<(std::ostream& os, __uint128_t x) {
    os << *(reinterpret_cast<uint64_t*>(&x) + 0);
    os << *(reinterpret_cast<uint64_t*>(&x) + 1);
    return os;
}

template <typename T>
void expect(T got, T expected) {
    if (got != expected) {
        std::cerr << "got '" << got << "' but expected '" << expected << "'" << std::endl;
    }
}

using kmer_t = default_kmer_t;

int main(int /*argc*/, char** argv) {
    uint64_t k = std::stoull(argv[1]);
    std::string read(
        "ttgttagcaaatgaagtcttagaacttttaagtcaagaaatttcaaaaaatgagatggaaaactatatatctcaaatcaaattcaatgaa"
        "aaactatcaaataacgaaactgctatttttacagcaccaaacgaacttatggctaaatttatacaaactagatatgcttctaagatcgct"
        "catctttttgagataaaaacaggaaataaaccaaatataagcatcactactcaaaaaaataaactatctatcaaaacaaaagacgtagat"
        "gtaaaacagatcagaactcaaagttcgcttttaaatccaagctatacttttgaaagcttcgtcgtaggcgactcaaatcaattcgcatat"
        "attagttcaaaacaagtagcagcaaatccaggccttgtttataatccactttttatatatggctcaactggacttggcaaaactcacctt"
        "ttacaatccatcggaaattactgtttagaacacggaaaaacagttatatgtgtaactagcgaacaatttatgagcgattttatgagaaa"
        "c");

    for (uint64_t i = 0; i != read.length() - k + 1; ++i) {
        bool is_valid = util::is_valid(read.data() + i, k);
        if (!is_valid) {
            std::cerr << "ERROR: '" << std::string(read.data() + i, k) << "' is valid!"
                      << std::endl;
        }
        std::cout << "read: '" << std::string(read.data() + i, k) << "'; ";
        kmer_t x = util::string_to_uint_kmer<kmer_t>(read.data() + i, k);
        std::string kmer = util::uint_kmer_to_string<kmer_t>(x, k);
        std::cout << "capitalized: '" << kmer << "'" << std::endl;
    }

    /****/

    expect(kmer_t::char_to_uint('A'), uint64_t(0));
    expect(kmer_t::char_to_uint('a'), uint64_t(0));

    expect(kmer_t::char_to_uint('C'), uint64_t(1));
    expect(kmer_t::char_to_uint('c'), uint64_t(1));

    expect(kmer_t::char_to_uint('T'), uint64_t(2));
    expect(kmer_t::char_to_uint('t'), uint64_t(2));

    expect(kmer_t::char_to_uint('G'), uint64_t(3));
    expect(kmer_t::char_to_uint('g'), uint64_t(3));

    /****/

    expect(util::canonicalize_basepair_forward_map[int('A')], 'A');
    expect(util::canonicalize_basepair_forward_map[int('a')], 'a');

    expect(util::canonicalize_basepair_forward_map[int('C')], 'C');
    expect(util::canonicalize_basepair_forward_map[int('c')], 'c');

    expect(util::canonicalize_basepair_forward_map[int('T')], 'T');
    expect(util::canonicalize_basepair_forward_map[int('t')], 't');

    expect(util::canonicalize_basepair_forward_map[int('G')], 'G');
    expect(util::canonicalize_basepair_forward_map[int('g')], 'g');

    /****/

    expect(util::canonicalize_basepair_reverse_map[int('A')], 'T');
    expect(util::canonicalize_basepair_reverse_map[int('a')], 't');

    expect(util::canonicalize_basepair_reverse_map[int('C')], 'G');
    expect(util::canonicalize_basepair_reverse_map[int('c')], 'g');

    expect(util::canonicalize_basepair_reverse_map[int('T')], 'A');
    expect(util::canonicalize_basepair_reverse_map[int('t')], 'a');

    expect(util::canonicalize_basepair_reverse_map[int('G')], 'C');
    expect(util::canonicalize_basepair_reverse_map[int('g')], 'c');

    return 0;
}