#include "include/kmer.hpp"
#include "include/util.hpp"
#include "tools/common.hpp"  // for random_kmer

#include <cctype>

using namespace sshash;

std::ostream& operator<<(std::ostream& os, __uint128_t x) {
    os << static_cast<uint64_t>(x) << static_cast<uint64_t>(x >> 64);
    return os;
}

static uint64_t num_failed = 0;

template <typename T>
void expect(T got, T expected) {
    if (got != expected) {
        ++num_failed;
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

    auto check_take = [](auto all_bits) {
        uint_kmer_t<decltype(all_bits), 4> value;
        value.bits = all_bits;
        value.take_chars(value.max_k);
        if (value.bits != all_bits) return false;
        value.take(value.uint_kmer_bits - 1);
        if (value.bits != (all_bits >> 1)) return false;
        value.take(0);
        return value.bits == 0;
    };
    if (!check_take(~uint64_t(0)) || !check_take(~__uint128_t(0))) {
        std::cerr << "incorrect k-mer mask\n";
        return 1;
    }

    using dna = dna_uint_kmer_t<uint64_t>;
    using protein = aa_uint_kmer_t<uint64_t>;
    for (unsigned byte = 128; byte < 256; ++byte) {
        char c = static_cast<char>(byte);
        char reverse = 1;
        dna::compute_reverse_complement(&c, &reverse, 1);
        if (dna::is_valid(c) || protein::is_valid(c) ||
            protein::char_to_uint(c) != constants::invalid_uint64 || reverse != 0) {
            std::cerr << "incorrect alphabet mapping for byte " << byte << std::endl;
            return 1;
        }
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
            ++num_failed;
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
    expect(kmer_t::char_to_uint('A'), uint64_t(0));
    expect(kmer_t::char_to_uint('a'), uint64_t(0));
    expect(kmer_t::char_to_uint('C'), uint64_t(1));
    expect(kmer_t::char_to_uint('c'), uint64_t(1));
    expect(kmer_t::char_to_uint('G'), uint64_t(2));
    expect(kmer_t::char_to_uint('g'), uint64_t(2));
    expect(kmer_t::char_to_uint('T'), uint64_t(3));
    expect(kmer_t::char_to_uint('t'), uint64_t(3));
#else
    expect(kmer_t::char_to_uint('A'), uint64_t(0));
    expect(kmer_t::char_to_uint('a'), uint64_t(0));
    expect(kmer_t::char_to_uint('C'), uint64_t(1));
    expect(kmer_t::char_to_uint('c'), uint64_t(1));
    expect(kmer_t::char_to_uint('T'), uint64_t(2));
    expect(kmer_t::char_to_uint('t'), uint64_t(2));
    expect(kmer_t::char_to_uint('G'), uint64_t(3));
    expect(kmer_t::char_to_uint('g'), uint64_t(3));
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

    return num_failed != 0;
}
