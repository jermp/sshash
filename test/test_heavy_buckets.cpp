#include <chrono>
#include <filesystem>

#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/info.cpp"

using namespace sshash;

struct test_directory {
    std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("sshash-heavy-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    test_directory() { std::filesystem::create_directory(path); }
    ~test_directory() {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};

template <typename Offsets>
void check(std::filesystem::path const& directory, std::vector<std::string> const& words) {
    auto input = directory / "input.fa";
    {
        std::ofstream out(input);
        out.exceptions(std::ios::failbit | std::ios::badbit);
        for (auto const& word : words) out << ">test\n" << word << '\n';
    }
    build_configuration config;
    config.k = 3;
    config.m = 1;
    config.verbose = false;
    config.tmp_dirname = directory.string();
    dictionary<aa_uint_kmer_t<uint64_t>, Offsets> dict;
    dict.build(input.string(), config);
    for (uint64_t i = 0; i < words.size(); ++i) {
        if (dict.lookup(words[i].data()).kmer_id != i) {
            throw std::runtime_error("incorrect heavy-bucket lookup for " + words[i]);
        }
    }
}

int main() {
    try {
        test_directory directory;
        // One record per distinct protein 3-mer produces dense minimizer buckets.
        std::vector<std::string> words;
        const std::string alphabet = "ACDEFGHIKLMNPQRSTVWY";
        for (char a : alphabet) for (char b : alphabet) for (char c : alphabet) {
            words.push_back(std::string{a, b, c});
        }
        check<decoded_offsets>(directory.path, words);
        check<encoded_offsets>(directory.path, words);
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
