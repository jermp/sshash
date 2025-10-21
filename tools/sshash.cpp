#include <iostream>

#include "common.hpp"
#include "perf.hpp"

#include "test/check.hpp"
#include "test/check_from_file.hpp"

#include "src/build.cpp"
#include "src/dictionary.cpp"
#include "src/query.cpp"
#include "src/info.cpp"

#include "build.cpp"
#include "query.cpp"
#include "permute.cpp"

using namespace sshash;

int check(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with the tool 'build'.", "-i", true);
    parser.add("verbose", "Verbose output.", "--verbose", false, true);
    if (!parser.parse()) return 0;
    auto index_filename = parser.get<std::string>("index_filename");
    bool verbose = parser.get<bool>("verbose");
    dictionary_type dict;
    load_dictionary(dict, index_filename, verbose);
    check_dictionary(dict);
    check_correctness_navigational_string_query(dict);
    check_correctness_kmer_iterator(dict);
    check_correctness_string_iterator(dict);
    return 0;
}

int bench(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with the tool 'build'.", "-i", true);
    parser.add("verbose", "Verbose output.", "--verbose", false, true);
    if (!parser.parse()) return 0;
    auto index_filename = parser.get<std::string>("index_filename");
    bool verbose = parser.get<bool>("verbose");
    dictionary_type dict;
    load_dictionary(dict, index_filename, verbose);

    // perf_test_lookup_by_list_size(dict);

    perf_test_lookup_access(dict);
    if (dict.weighted()) perf_test_lookup_weight(dict);
    perf_test_iterator(dict);

    return 0;
}

int help(char* arg0) {
    std::cout << "== SSHash: (S)parse and (S)kew (Hash)ing of k-mers ========================="
              << std::endl
              << std::endl;
    std::cout << "Usage: " << arg0 << " <tool> ...\n\n"
              << "Available tools:\n"
              << "  build              \t build a dictionary \n"
              << "  query              \t query a dictionary \n"
              << "  check              \t check correctness of a dictionary \n"
              << "  bench              \t run performance tests for a dictionary \n"
              << "  permute            \t permute a weighted input file \n"
              << std::endl;

    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) return help(argv[0]);
    print_cmd(argc, argv);
    auto tool = std::string(argv[1]);
    if (tool == "build") {
        return build(argc - 1, argv + 1);
    } else if (tool == "query") {
        return query(argc - 1, argv + 1);
    } else if (tool == "check") {
        return check(argc - 1, argv + 1);
    } else if (tool == "bench") {
        return bench(argc - 1, argv + 1);
    } else if (tool == "permute") {
        return permute(argc - 1, argv + 1);
    }
    std::cout << "Unsupported tool '" << tool << "'.\n" << std::endl;
    return help(argv[0]);
}
