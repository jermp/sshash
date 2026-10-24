#include <chrono>
#include <filesystem>
#include <map>
#include <set>
#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/info.cpp"
using namespace sshash;
template <typename Offsets>
void check(std::filesystem::path const& directory, std::vector<std::string> const& words, uint64_t k,
           uint64_t m, bool require_heavy_and_absent = false) {
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
    c.m = m;
    c.verbose = false;
    c.tmp_dirname = directory.string();
    dict.build(input.string(), c);
    std::map<uint64_t, std::set<uint64_t>> expected;
    uint64_t base = 0, kmer_id = 0;
    for (auto const& word : words) {
        for (uint64_t i = 0; i + k <= word.size(); ++i) {
            auto encoded = util::string_to_uint_kmer<Kmer>(word.data() + i, k);
            auto mini = util::compute_minimizer(encoded, k, m, dict.hasher());
            expected[mini.minimizer].insert(base + i + mini.pos_in_kmer);
            std::string accessed(k, ' ');
            dict.access(kmer_id, accessed.data());
            if (accessed != word.substr(i, k))
                throw std::runtime_error("incorrect k-mer access");
            ++kmer_id;
        }
        base += word.size();
    }
    uint64_t heavy = 0, absent = 0;
    for (uint64_t mini = 0; mini < (uint64_t(1) << (Kmer::bits_per_char * m)); ++mini) {
        if (util::canonical_mmer<Kmer>(mini, m) != mini) continue;
        heavy += expected[mini].size() > (uint64_t(1) << constants::min_l);
        absent += expected[mini].empty();
        std::set<uint64_t> observed;
        dict.call_minimizer_positions(mini, [&](uint64_t offset) {
            if (!observed.insert(offset).second)
                throw std::runtime_error("minimizer enumeration repeated a position");
        });
        if (observed != expected[mini]) throw std::runtime_error("incorrect minimizer positions");
    }
    if (require_heavy_and_absent && (heavy < 2 || absent == 0))
        throw std::runtime_error("fixture must cover adjacent heavy buckets and absent minimizers");
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
        check<decoded_offsets>(dir, {"AAACG"}, 3, 1);
        check<encoded_offsets>(dir, {"AAACG"}, 3, 1);
        // This subset has two heavy buckets (65 and 146 positions) and two absent minimizers.
        words.resize(384);
        check<decoded_offsets>(dir, words, 5, 2, true);
        check<encoded_offsets>(dir, words, 5, 2, true);
        std::filesystem::remove_all(dir);
    } catch (std::exception const& e) {
        std::cerr << e.what() << '\n';
        std::filesystem::remove_all(dir);
        return 1;
    }
}
