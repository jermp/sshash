#include <iostream>

#include "common.hpp"
#include "bench_utils.hpp"
#include "check_utils.hpp"
#include "build.cpp"
#include "query.cpp"
#include "permute.cpp"

using namespace sshash;

int check(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with the tool 'build'.", "-i", true);
    parser.add("verbose", "Verbose output.", "--verbose", false, true);
    if (!parser.parse()) return 1;
    auto index_filename = parser.get<std::string>("index_filename");
    bool verbose = parser.get<bool>("verbose");
    dictionary dict;
    load_dictionary(dict, index_filename, verbose);
    check_dictionary(dict);
    check_correctness_navigational_contig_query(dict);
    return 0;
}

int bench(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with the tool 'build'.", "-i", true);
    parser.add("verbose", "Verbose output.", "--verbose", false, true);
    if (!parser.parse()) return 1;
    auto index_filename = parser.get<std::string>("index_filename");
    bool verbose = parser.get<bool>("verbose");
    dictionary dict;
    load_dictionary(dict, index_filename, verbose);
    perf_test_lookup_access(dict);
    if (dict.weighted()) perf_test_lookup_weight(dict);
    perf_test_iterator(dict);
    return 0;
}

int dump(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with the tool 'build'.", "-i", true);
    parser.add("output_filename", "A FASTA file where the output will be saved.", "-o", true);
    parser.add("verbose", "Verbose output.", "--verbose", false, true);
    if (!parser.parse()) return 1;
    auto index_filename = parser.get<std::string>("index_filename");
    auto output_filename = parser.get<std::string>("output_filename");
    bool verbose = parser.get<bool>("verbose");
    dictionary dict;
    load_dictionary(dict, index_filename, verbose);
    dict.dump(output_filename);
    return 0;
}

int compute_statistics(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with the tool 'build'.", "-i", true);
    parser.add("verbose", "Verbose output.", "--verbose", false, true);
    if (!parser.parse()) return 1;
    auto index_filename = parser.get<std::string>("index_filename");
    bool verbose = parser.get<bool>("verbose");
    dictionary dict;
    load_dictionary(dict, index_filename, verbose);
    dict.compute_statistics();
    return 0;
}

int help(char* arg0) {
    std::cout << "(S)parse and (S)kew (Hash)ing of k-mers" << std::endl << std::endl;
    std::cout << "Usage: " << arg0 << " <tool> ...\n\n"
              << "Available tools:\n"
              << "  build              \t build a dictionary \n"
              << "  query              \t query a dictionary \n"
              << "  check              \t check correctness of a dictionary \n"
              << "  bench              \t run performance tests for a dictionary \n"
              << "  dump               \t write super-k-mers of a dictionary to a fasta file \n"
              << "  permute            \t permute a weighted input file \n"
              << "  compute-statistics \t compute index statistics " << std::endl;
    return 1;
}

int main(int argc, char** argv) {
    if (argc < 2) return help(argv[0]);
    auto tool = std::string(argv[1]);
    if (tool == "build") {
        return build(argc - 1, argv + 1);
    } else if (tool == "query") {
        return query(argc - 1, argv + 1);
    } else if (tool == "check") {
        return check(argc - 1, argv + 1);
    } else if (tool == "bench") {
        return bench(argc - 1, argv + 1);
    } else if (tool == "dump") {
        return dump(argc - 1, argv + 1);
    } else if (tool == "permute") {
        return permute(argc - 1, argv + 1);
    } else if (tool == "compute-statistics") {
        return compute_statistics(argc - 1, argv + 1);
    }
    std::cout << "Unsupported tool '" << tool << "'.\n" << std::endl;
    return help(argv[0]);
}
