#include "include/util.hpp"
#include "src/common.hpp"  // for random_kmer
#include <cctype>

using namespace sshash;

std::ostream& operator<<(std::ostream& os, __uint128_t x) {
    os << static_cast<uint64_t>(x) << static_cast<uint64_t>(x >> 64);
    return os;
}

template <typename T>
void expect(T got, T expected) {
    if (got != expected) {
        std::cerr << "got '" << got << "' but expected '" << expected << "'" << std::endl;
    }
}

using kmer_t = default_kmer_t;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Usage " << argv[0] << " k" << std::endl;
        return 1;
    }
    uint64_t k = std::stoull(argv[1]);
    if (k > kmer_t::max_k) {
        std::cerr << "k must be less <= " << kmer_t::max_k << " but got k = " << k << '\n';
        return 1;
    }
    std::string read(
        "ttgttagcaaatgaagtcttagaacttttaagtcaagaaatttcaaaaaatgagatggaaaactatatatctcaaatcaaattcaatgaa"
        "aaactatcaaataacgaaactgctatttttacagcaccaaacgaacttatggctaaatttatacaaactagatatgcttctaagatcgct"
        "catctttttgagataaaaacaggaaataaaccaaatataagcatcactactcaaaaaaataaactatctatcaaaacaaaagacgtagat"
        "gtaaaacagatcagaactcaaagttcgcttttaaatccaagctatacttttgaaagcttcgtcgtaggcgactcaaatcaattcgcatat"
        "attagttcaaaacaagtagcagcaaatccaggccttgtttataatccactttttatatatggctcaactggacttggcaaaactcacctt"
        "ttacaatccatcggaaattactgtttagaacacggaaaaacagttatatgtgtaactagcgaacaatttatgagcgattttatgagaaa"
        "c");

    for (uint64_t i = 0; i != read.length() - k + 1; ++i) {
        bool is_valid = util::is_valid<kmer_t>(read.data() + i, k);
        if (!is_valid) {
            std::cerr << "ERROR: '" << std::string(read.data() + i, k) << "' is NOT valid!"
                      << std::endl;
        }
        kmer_t x = util::string_to_uint_kmer<kmer_t>(read.data() + i, k);
        std::string kmer = util::uint_kmer_to_string<kmer_t>(x, k);
        std::transform(kmer.begin(), kmer.end(), kmer.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        expect(std::string(read.data() + i, k), kmer);
    }

    /****/

#ifdef SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING
    expect(kmer_t::char_to_uint('A'), 0);
    expect(kmer_t::char_to_uint('a'), 0);
    expect(kmer_t::char_to_uint('C'), 1);
    expect(kmer_t::char_to_uint('c'), 1);
    expect(kmer_t::char_to_uint('G'), 2);
    expect(kmer_t::char_to_uint('g'), 2);
    expect(kmer_t::char_to_uint('T'), 3);
    expect(kmer_t::char_to_uint('t'), 3);
#else
    expect(kmer_t::char_to_uint('A'), 0);
    expect(kmer_t::char_to_uint('a'), 0);
    expect(kmer_t::char_to_uint('C'), 1);
    expect(kmer_t::char_to_uint('c'), 1);
    expect(kmer_t::char_to_uint('T'), 2);
    expect(kmer_t::char_to_uint('t'), 2);
    expect(kmer_t::char_to_uint('G'), 3);
    expect(kmer_t::char_to_uint('g'), 3);
#endif

    for (uint64_t kmer_len = 1; kmer_len <= k; ++kmer_len) {
        std::string kmer, rc;
        kmer.resize(kmer_len);
        rc.resize(kmer_len);
        for (uint64_t i = 0; i != 1000; ++i) {
            // generate a random kmer of length kmer_len
            random_kmer(kmer.data(), kmer_len);
            kmer_t::compute_reverse_complement(kmer.data(), rc.data(), kmer_len);
            kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(kmer.data(), kmer_len);
            uint_kmer.reverse_complement_inplace(kmer_len);
            expect(util::uint_kmer_to_string(uint_kmer, kmer_len), rc);
        }
    }

    /****/

    expect(kmer_t::canonicalize_basepair_forward_map[int('A')], 'A');
    expect(kmer_t::canonicalize_basepair_forward_map[int('a')], 'a');

    expect(kmer_t::canonicalize_basepair_forward_map[int('C')], 'C');
    expect(kmer_t::canonicalize_basepair_forward_map[int('c')], 'c');

    expect(kmer_t::canonicalize_basepair_forward_map[int('T')], 'T');
    expect(kmer_t::canonicalize_basepair_forward_map[int('t')], 't');

    expect(kmer_t::canonicalize_basepair_forward_map[int('G')], 'G');
    expect(kmer_t::canonicalize_basepair_forward_map[int('g')], 'g');

    /****/

    expect(kmer_t::canonicalize_basepair_reverse_map[int('A')], 'T');
    expect(kmer_t::canonicalize_basepair_reverse_map[int('a')], 't');

    expect(kmer_t::canonicalize_basepair_reverse_map[int('C')], 'G');
    expect(kmer_t::canonicalize_basepair_reverse_map[int('c')], 'g');

    expect(kmer_t::canonicalize_basepair_reverse_map[int('T')], 'A');
    expect(kmer_t::canonicalize_basepair_reverse_map[int('t')], 'a');

    expect(kmer_t::canonicalize_basepair_reverse_map[int('G')], 'C');
    expect(kmer_t::canonicalize_basepair_reverse_map[int('g')], 'c');

    return 0;
}
