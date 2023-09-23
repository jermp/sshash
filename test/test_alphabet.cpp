#include "../include/util.hpp"

using namespace sshash;

template <typename T>
void expect(T got, T expected) {
    if (got != expected) {
        std::cerr << "got '" << got << "' but expected '" << expected << "'" << std::endl;
    }
}

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
        kmer_t x = util::string_to_uint_kmer(read.data() + i, k);
        std::string kmer = util::uint_kmer_to_string(x, k);
        std::cout << "capitalized: '" << kmer << "'" << std::endl;
    }

    /****/

    expect(util::char_to_uint('A'), kmer_t(0));
    expect(util::char_to_uint('a'), kmer_t(0));

    expect(util::char_to_uint('C'), kmer_t(1));
    expect(util::char_to_uint('c'), kmer_t(1));

    expect(util::char_to_uint('T'), kmer_t(2));
    expect(util::char_to_uint('t'), kmer_t(2));

    expect(util::char_to_uint('G'), kmer_t(3));
    expect(util::char_to_uint('g'), kmer_t(3));

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