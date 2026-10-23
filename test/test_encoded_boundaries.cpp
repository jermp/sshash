#include <chrono>
#include <filesystem>

#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/info.cpp"

using namespace sshash;

struct test_directory {
    std::filesystem::path path =
        std::filesystem::temp_directory_path() /
        ("sshash-boundaries-" +
         std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));

    test_directory() { std::filesystem::create_directory(path); }
    ~test_directory() {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};

using dict_type = dictionary<dna_uint_kmer_t<uint64_t>, encoded_offsets>;

void verify(dict_type const& dict, std::vector<std::string> const& strings) {
    // The reverse complement TCGTC spans the second/third string boundary.
    if (dict.lookup("GACGA").kmer_id != constants::invalid_uint64) {
        throw std::runtime_error("absent k-mer GACGA reported present");
    }
    uint64_t kmer_id = 0;
    for (auto const& string : strings) {
        for (uint64_t i = 0; i + dict.k() <= string.size(); ++i, ++kmer_id) {
            if (dict.lookup(string.data() + i).kmer_id != kmer_id) {
                throw std::runtime_error("stored k-mer lost at string boundary");
            }
        }
    }
}

int main() {
    try {
        test_directory directory;
        std::vector<std::string> strings = {"CGGGGGC", "GTGTATT", "CGTCAGA", "CCACGAA"};
        auto input = directory.path / "input.fa";
        {
            std::ofstream out(input);
            out.exceptions(std::ios::failbit | std::ios::badbit);
            for (auto const& string : strings) out << ">test\n" << string << '\n';
        }
        build_configuration config;
        config.k = 5;
        config.m = 1;
        config.verbose = false;
        config.tmp_dirname = directory.path.string();
        dict_type dict;
        dict.build(input.string(), config);
        verify(dict, strings);
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
