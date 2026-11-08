#include <chrono>
#include <filesystem>

#include <sys/resource.h>

#include "src/builder/build.cpp"
#include "src/dictionary.cpp"
#include "src/info.cpp"

using namespace sshash;

struct test_directory {
    std::filesystem::path path =
        std::filesystem::temp_directory_path() /
        ("sshash-malformed-" +
         std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));

    test_directory() { std::filesystem::create_directory(path); }
    ~test_directory() {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};

using dict_type = dictionary<dna_uint_kmer_t<uint64_t>, decoded_offsets>;

void expect_rejected(test_directory const& directory, std::string const& contents,
                     std::string const& reason) {
    auto input = directory.path / "input.fa";
    {
        std::ofstream out(input);
        out.exceptions(std::ios::failbit | std::ios::badbit);
        out << contents;
    }
    build_configuration config;
    config.k = 21;
    config.m = 11;
    config.verbose = false;
    config.tmp_dirname = directory.path.string();
    try {
        dict_type dict;
        dict.build(input.string(), config);
    } catch (std::runtime_error const& error) {
        if (std::string(error.what()).find(reason) == std::string::npos) {
            throw std::runtime_error("expected error containing '" + reason + "' but got '" +
                                     error.what() + "'");
        }
        return;
    }
    throw std::runtime_error("malformed input accepted, expected: " + reason);
}

void expect_num_kmers(test_directory const& directory, std::string const& contents,
                      uint64_t num_kmers) {
    auto input = directory.path / "input.fa";
    {
        std::ofstream out(input);
        out.exceptions(std::ios::failbit | std::ios::badbit);
        out << contents;
    }
    build_configuration config;
    config.k = 21;
    config.m = 11;
    config.verbose = false;
    config.tmp_dirname = directory.path.string();
    dict_type dict;
    dict.build(input.string(), config);
    if (dict.num_kmers() != num_kmers) {
        throw std::runtime_error("expected " + std::to_string(num_kmers) + " k-mers but got " +
                                 std::to_string(dict.num_kmers()));
    }
}

int main() {
    /*
        A sequence shorter than k used to underflow its k-mer count, so the build
        never terminated and wrote temp files until the disk was full. Cap the size
        of written files so that a regression fails with SIGXFSZ instead.
    */
    rlimit limit;
    getrlimit(RLIMIT_FSIZE, &limit);
    if (limit.rlim_cur == RLIM_INFINITY or limit.rlim_cur > (64 << 20)) {
        limit.rlim_cur = 64 << 20;
        setrlimit(RLIMIT_FSIZE, &limit);
    }

    const std::string seq = "ACGTTGCAAGGCTTACCATGGTACCAGTTA";
    try {
        test_directory directory;
        expect_rejected(directory, ">a\n" + seq + "\n>b\nACGTA\n",
                        "sequence 2 has length 5, shorter than k = 21");
        // a wrapped record is read as one sequence, even with a tail shorter than k
        expect_num_kmers(directory,
                         ">a\n" + seq + "\nTTGACCATGCAGTCAGGTACCATGACTTGA\nAC\n>b\n"
                         "GATCCAGTTGCAAGCTTAGGCATCGATAGC\n",
                         (62 - 21 + 1) + (30 - 21 + 1));
        expect_rejected(directory, "@a\n" + seq + "\n+\n" + std::string(30, 'I') + "\n",
                        "expected a FASTA header starting with '>' before sequence 1");
        expect_rejected(directory, ">a\n" + seq, "sequence 1 is not terminated by a newline");
        expect_rejected(directory, ">a\n" + seq + "\n>b\n", "record 2 has a header but no sequence");
        expect_rejected(directory, ">a\n>b\n" + seq + "\n", "record 1 has a header but no sequence");
        expect_rejected(directory, ">a\n" + seq + "\n" + seq, "sequence 1 is not terminated by a newline");
        expect_rejected(directory, ">a\nACGTTGCAAGGCNTACCATGGTACCAGTTA\n",
                        "sequence 1 has character 'N' at position 12, outside the alphabet");
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
