#include <chrono>
#include <filesystem>
#include <random>
#include <unordered_map>

#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/info.cpp"

using namespace sshash;

template <typename Kmer, typename Offsets>
void check_forward_heavy(std::filesystem::path const& directory, uint64_t k) {
    auto reverse = [](std::string s) {
        std::reverse(s.begin(), s.end());
        for (char& c : s) c = c == 'A' ? 'T' : c == 'T' ? 'A' : c == 'C' ? 'G' : 'C';
        return s;
    };
    const std::string anchor = "TGACGATAAACCAACACGC";
    std::mt19937_64 random(20261104);
    auto flank = [&] {
        std::string s;
        for (uint64_t i = 0; i < (k - anchor.size()) / 2; ++i) s += "ACGT"[random() % 4];
        return s;
    };
    build_configuration config;
    config.k = k;
    config.m = anchor.size();
    config.preserve_orientation = true;
    config.verbose = false;
    config.tmp_dirname = directory.string();
    hasher_type hasher(config.seed);
    auto anchor_kmer = util::string_to_uint_kmer<Kmer>(anchor.data(), anchor.size());
    uint64_t minimizer = util::canonical_mmer<Kmer>(uint64_t(anchor_kmer), anchor.size());
    std::vector<std::string> words;
    std::unordered_map<std::string, uint64_t> ids;
    for (uint64_t family = 0, attempts = 0; family < 512; ++attempts) {
        if (attempts > 10000) throw std::runtime_error("could not populate heavy bucket");
        auto word = flank() + anchor + flank();
        auto reversed = reverse(word);
        auto value = util::string_to_uint_kmer<Kmer>(word.data(), k);
        auto value_rc = value;
        value_rc.reverse_complement_inplace(k);
        if (ids.count(word) || ids.count(reversed) ||
            util::compute_minimizer(value, value_rc, k, config.m, hasher).minimizer != minimizer)
            continue;
        ids.emplace(word, words.size());
        words.push_back(word);
        if (family % 4 < 2) {
            ids.emplace(reversed, words.size());
            words.push_back(reversed);
        }
        ++family;
    }
    auto input = directory / "forward.fa";
    {
        std::ofstream out(input);
        for (auto const& word : words) out << ">test\n" << word << '\n';
    }
    using dict_type = dictionary<Kmer, Offsets>;
    dict_type dict;
    dict.build(input.string(), config);
    auto verify = [&](dict_type const& index, std::string const& word) {
        for (bool with_rc : {false, true}) {
            auto result = index.lookup(word.data(), with_rc);
            bool present = ids.count(word) || (with_rc && ids.count(reverse(word)));
            if ((result.kmer_id != constants::invalid_uint64) != present)
                throw std::runtime_error("incorrect forward-only membership");
            if (!present) continue;
            std::string found(k, ' ');
            index.access(result.kmer_id, found.data());
            if (result.kmer_orientation == constants::backward_orientation) {
                if (!with_rc)
                    throw std::runtime_error("forward-only lookup returned reverse strand");
                found = reverse(found);
            }
            if (found != word || (!with_rc && result.kmer_id != ids.at(word)))
                throw std::runtime_error("incorrect forward-only node or orientation");
        }
    };
    auto path = (directory / "forward.sshash").string();
    essentials::save(dict, path.c_str());
    dict_type loaded, mapped;
    essentials::load(loaded, path.c_str());
    essentials::mmap(mapped, path.c_str());
    for (auto word : words) {
        for (auto const* index : {&dict, &loaded, &mapped}) {
            verify(*index, word);
            verify(*index, reverse(word));
        }
        word[k / 2] = word[k / 2] == 'A' ? 'C' : 'A';
        verify(dict, word);
        verify(mapped, word);
    }
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    std::filesystem::path dir(argv[1]);
    std::filesystem::create_directory(dir);
    check_forward_heavy<dna_uint_kmer_t<uint64_t>, decoded_offsets>(dir, 31);
    check_forward_heavy<dna_uint_kmer_t<uint64_t>, encoded_offsets>(dir, 31);
    check_forward_heavy<dna_uint_kmer_t<__uint128_t>, decoded_offsets>(dir, 41);
    check_forward_heavy<dna_uint_kmer_t<__uint128_t>, encoded_offsets>(dir, 41);
}
