#include <chrono>
#include <filesystem>
#include <map>
#include <set>
#include <sstream>
#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/info.cpp"
using namespace sshash;
template <typename Offsets>
void check(std::filesystem::path const& directory, std::vector<std::string> const& words, uint64_t k) {
    using Kmer = dna_uint_kmer_t<uint64_t>;
    dictionary<Kmer, Offsets> dict;
    dict.call_minimizer_positions(0, [](uint64_t) { throw std::runtime_error("empty dictionary emitted a position"); });
    auto input = directory / "input.fa";
    {
        std::ofstream out(input);
        out.exceptions(std::ios::failbit | std::ios::badbit);
        for (auto const& word : words) out << ">test\n" << word << '\n';
    }
    build_configuration c;
    c.k = k;
    c.m = 1;
    c.verbose = false;
    c.tmp_dirname = directory.string();
    dict.build(input.string(), c);
    std::map<uint64_t, std::set<uint64_t>> expected;
    uint64_t base = 0;
    for (auto const& word : words) {
        for (uint64_t i = 0; i + k <= word.size(); ++i) {
            auto encoded = util::string_to_uint_kmer<Kmer>(word.data() + i, k);
            auto mini = util::compute_minimizer(encoded, k, 1, dict.hasher());
            expected[mini.minimizer].insert(base + i + mini.pos_in_kmer);
        }
        base += word.size();
    }
    auto verify = [&] {
        for (uint64_t mini = 0; mini < 4; ++mini) {
            if (util::canonical_mmer<Kmer>(mini, 1) != mini) continue;
            std::set<uint64_t> observed;
            dict.call_minimizer_positions(mini, [&](uint64_t offset) {
                if (!observed.insert(offset).second)
                    throw std::runtime_error("minimizer enumeration repeated a position");
            });
            if (observed != expected[mini]) throw std::runtime_error("incorrect minimizer positions");
        }
    };
    verify();
    std::stringstream bytes;
    essentials::generic_saver saver(bytes);
    saver.visit(dict);
    essentials::generic_loader loader(bytes);
    loader.visit(dict);
    verify();
}
int main() {
    auto dir = std::filesystem::temp_directory_path() /
        ("sshash-positions-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(dir);
    try {
        std::vector<std::string> words;
        using Kmer = dna_uint_kmer_t<uint64_t>;
        for (uint64_t i = 0; i < 1024; ++i) {
            std::string word(5, 'A');
            for (uint64_t j = 0; j < 5; ++j) word[j] = "ACGT"[(i >> (2*j)) & 3];
            auto encoded = util::string_to_uint_kmer<Kmer>(word.data(), 5);
            auto rc = encoded; rc.reverse_complement_inplace(5);
            if (encoded < rc) words.push_back(word);
        }
        check<decoded_offsets>(dir, {"AAACG"}, 3);
        check<encoded_offsets>(dir, {"AAACG"}, 3);
        check<decoded_offsets>(dir, words, 5);
        check<encoded_offsets>(dir, words, 5);
        std::filesystem::remove_all(dir);
    } catch (std::exception const& e) {
        std::cerr << e.what() << '\n';
        std::filesystem::remove_all(dir);
        return 1;
    }
}
