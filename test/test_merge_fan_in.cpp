#include <chrono>
#include <filesystem>
#include <random>

#include <sys/resource.h>

#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/info.cpp"

using namespace sshash;

struct test_directory {
    std::filesystem::path path =
        std::filesystem::temp_directory_path() /
        ("sshash-fan-in-" +
         std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));

    test_directory() { std::filesystem::create_directory(path); }
    ~test_directory() {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};

static void require(bool condition, std::string const& message) {
    if (!condition) throw std::runtime_error(message);
}

int main() {
    /* the fan-in follows the limits, so typical limits merge thousands of runs in one pass */
    require(merge_fan_in(1024, 65530) == 512, "fan-in for 1024 open files");
    require(merge_fan_in(1048576, 65530) == 32765, "fan-in for the default vm.max_map_count");
    require(merge_fan_in(1048576, 262144) == 131072, "fan-in for 1048576 open files");
    require(merge_fan_in(3, 65530) == 2, "fan-in is at least 2");

    /* a low soft limit is raised to the hard limit before merging */
    rlimit limit;
    getrlimit(RLIMIT_NOFILE, &limit);
    if (limit.rlim_max != RLIM_INFINITY and limit.rlim_max > 32) {
        limit.rlim_cur = 32;
        require(setrlimit(RLIMIT_NOFILE, &limit) == 0, "cannot lower the soft RLIMIT_NOFILE");
        uint64_t fan_in = max_merge_fan_in();
        rlimit raised;
        getrlimit(RLIMIT_NOFILE, &raised);
        require(raised.rlim_cur == limit.rlim_max, "soft RLIMIT_NOFILE not raised to the hard limit");
        require(fan_in > 16, "fan-in ignores the raised soft limit");
    }

    /*
        Every build thread flushes at least one sorted run of minimizer tuples,
        so 49 threads write 49 runs. With at most 32 open files (soft and hard),
        the runs used to be merged all at once and the build failed with "cannot
        open file"; they must now be merged 16 at a time (49 -> 4 -> final merge).
    */
    limit.rlim_cur = limit.rlim_max = 32;
    require(setrlimit(RLIMIT_NOFILE, &limit) == 0, "cannot lower RLIMIT_NOFILE");

    const uint64_t k = 31;
    const uint64_t num_strings = 500;
    const uint64_t string_length = 200;
    std::mt19937_64 rng(42);
    std::vector<std::string> strings(num_strings);
    for (auto& string : strings) {
        for (uint64_t i = 0; i != string_length; ++i) string += "ACGT"[rng() % 4];
    }

    try {
        test_directory directory;
        auto input = directory.path / "input.fa";
        {
            std::ofstream out(input);
            for (auto const& string : strings) out << ">test\n" << string << '\n';
        }
        build_configuration config;
        config.k = k;
        config.m = 15;
        config.verbose = false;
        config.num_threads = 49;
        config.ram_limit_in_GiB = 1;
        config.tmp_dirname = directory.path.string();
        dictionary<dna_uint_kmer_t<uint64_t>, decoded_offsets> dict;
        dict.build(input.string(), config);

        const uint64_t num_kmers = num_strings * (string_length - k + 1);
        require(dict.num_kmers() == num_kmers, "expected " + std::to_string(num_kmers) +
                                                   " k-mers but got " +
                                                   std::to_string(dict.num_kmers()));
        uint64_t expected_id = 0;
        for (auto const& string : strings) {
            for (uint64_t i = 0; i + k <= string.size(); ++i, ++expected_id) {
                require(dict.lookup(string.data() + i, false).kmer_id == expected_id,
                        "wrong id for k-mer " + string.substr(i, k));
            }
        }
        for (auto const& entry : std::filesystem::directory_iterator(directory.path)) {
            require(entry.path() == input, "temporary file left behind: " + entry.path().string());
        }
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
