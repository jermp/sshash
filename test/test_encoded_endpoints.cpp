#include <chrono>
#include <filesystem>
#include <sstream>

#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/info.cpp"

using namespace sshash;

struct test_directory {
    std::filesystem::path path =
        std::filesystem::temp_directory_path() /
        ("sshash-endpoints-" +
         std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));

    test_directory() { std::filesystem::create_directory(path); }
    ~test_directory() {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};

using dict_type = dictionary<dna_uint_kmer_t<uint64_t>, encoded_offsets>;

void verify(dict_type const& dict, std::vector<std::string> const& strings) {
    uint64_t base = 0, kmer_id = 0, string_id = 0;
    for (auto const& string : strings) {
        // Check the raw endpoints before lookup can assert on corrupted bounds.
        auto begin = dict.strings_offsets().access(string_id);
        auto end = dict.strings_offsets().access(++string_id);
        if (begin != base || end != base + string.size()) {
            throw std::runtime_error("incorrect encoded string endpoint");
        }
        for (uint64_t i = 0; i + dict.k() <= string.size(); ++i, ++kmer_id) {
            if (dict.lookup(string.data() + i).kmer_id != kmer_id) {
                throw std::runtime_error("incorrect k-mer ID near encoded endpoint");
            }
        }
        base += string.size();
    }
}

int main() {
    try {
        test_directory directory;
        build_configuration config;
        config.k = 3;
        config.m = 1;
        config.verbose = false;
        config.tmp_dirname = directory.path.string();
        auto input = directory.path / "input.fa";

        // Cover powers of two, their neighboring sizes, and multiple strings.
        for (auto const& strings : std::vector<std::vector<std::string>>{{"AAA"},
                                                                         {"AAAC"},
                                                                         {"AAACC"},
                                                                         {"AAACCGT"},
                                                                         {"AAACCGTA"},
                                                                         {"AAACCGTAT"},
                                                                         {"AAAC", "CCGT"}}) {
            {
                std::ofstream out(input);
                out.exceptions(std::ios::failbit | std::ios::badbit);
                for (auto const& string : strings) out << ">test\n" << string << '\n';
            }
            dict_type dict;
            dict.build(input.string(), config);
            verify(dict, strings);

            std::stringstream bytes;
            essentials::generic_saver saver(bytes);
            saver.visit(dict);
            dict_type loaded;
            essentials::generic_loader loader(bytes);
            loader.visit(loaded);
            verify(loaded, strings);
        }
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
