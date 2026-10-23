#include <chrono>
#include <filesystem>

#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/info.cpp"
#include "include/streaming_query.hpp"

using namespace sshash;

struct test_directory {
    std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("sshash-protein-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    test_directory() { std::filesystem::create_directory(path); }
    ~test_directory() {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};

template <typename Kmer, typename Offsets>
void check(std::filesystem::path const& directory, std::string const& sequence,
           uint64_t k, uint64_t m) {
    auto input = directory / "input.fa";
    {
        std::ofstream out(input);
        out.exceptions(std::ios::failbit | std::ios::badbit);
        out << ">test\n" << sequence << '\n';
    }
    build_configuration config;
    config.k = k;
    config.m = m;
    config.verbose = false;
    config.tmp_dirname = directory.string();
    dictionary<Kmer, Offsets> dict;
    dict.build(input.string(), config);
    streaming_query<decltype(dict)> stream(&dict);
    for (uint64_t i = 0; i + k <= sequence.size(); ++i) {
        auto result = dict.lookup(sequence.data() + i, false);
        if (result.kmer_id != i || result.kmer_orientation != constants::forward_orientation) {
            throw std::runtime_error("incorrect forward lookup for " + sequence.substr(i, k));
        }
        if (!equal_lookup_result(result, stream.lookup(sequence.data() + i))) {
            throw std::runtime_error("streaming lookup differs from direct lookup");
        }
    }
}

template <typename Offsets>
void check_alphabets(std::filesystem::path const& directory) {
    check<dna_uint_kmer_t<uint64_t>, Offsets>(directory, "AAACG", 3, 1);
    check<aa_uint_kmer_t<uint64_t>, Offsets>(directory, "MKWVTFISRLDPQEAHNYCG", 3, 1);
    check<aa_uint_kmer_t<__uint128_t>, Offsets>(directory,
        "MKWVTFISRLDPQEAHNYCGMKWVTFISRLDPQEAHNYCG", 21, 3);
}

int main() {
    try {
        test_directory directory;
        check_alphabets<decoded_offsets>(directory.path);
        check_alphabets<encoded_offsets>(directory.path);
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
