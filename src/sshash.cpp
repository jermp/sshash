#include <iostream>

#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/dictionary.hpp"
#include "../include/query/streaming_query.hpp"
#include "permute.hpp"
#include "bench_utils.hpp"
#include "check_utils.hpp"

using namespace sshash;

int build(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);

    /* Required arguments. */
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:\n"
               "\t- without duplicate nor invalid kmers\n"
               "\t- one DNA sequence per line.\n"
               "\tFor example, it could be the de Bruijn graph topology output by BCALM.",
               "-i", true);
    parser.add("k", "K-mer length (must be <= " + std::to_string(constants::max_k) + ").", "-k",
               true);
    parser.add("m", "Minimizer length (must be < k).", "-m", true);

    /* Optional arguments. */
    parser.add("seed",
               "Seed for construction (default is " + std::to_string(constants::seed) + ").", "-s",
               false);
    parser.add("l",
               "A (integer) constant that controls the space/time trade-off of the dictionary. "
               "A reasonable values lies between 2 and 12 (default is " +
                   std::to_string(constants::min_l) + ").",
               "-l", false);
    parser.add("c",
               "A (floating point) constant that trades construction speed for space effectiveness "
               "of minimal perfect hashing. "
               "A reasonable value lies between 3.0 and 10.0 (default is " +
                   std::to_string(constants::c) + ").",
               "-c", false);
    parser.add("output_filename", "Output file name where the data structure will be serialized.",
               "-o", false);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory. Default is directory '" +
            constants::default_tmp_dirname + "'.",
        "-d", false);
    parser.add("canonical_parsing",
               "Canonical parsing of k-mers. This option changes the parsing and results in a "
               "trade-off between index space and lookup time.",
               "--canonical-parsing", false, true);
    parser.add("weighted", "Also store the weights in compressed format.", "--weighted", false,
               true);
    parser.add("check", "Check correctness after construction.", "--check", false, true);
    parser.add("bench", "Run benchmark after construction.", "--bench", false, true);
    parser.add("verbose", "Verbose output during construction.", "--verbose", false, true);

    if (!parser.parse()) return 1;

    auto input_filename = parser.get<std::string>("input_filename");
    auto k = parser.get<uint64_t>("k");
    auto m = parser.get<uint64_t>("m");

    dictionary dict;

    build_configuration build_config;
    build_config.k = k;
    build_config.m = m;

    if (parser.parsed("seed")) build_config.seed = parser.get<uint64_t>("seed");
    if (parser.parsed("l")) build_config.l = parser.get<double>("l");
    if (parser.parsed("c")) build_config.c = parser.get<double>("c");
    build_config.canonical_parsing = parser.get<bool>("canonical_parsing");
    build_config.weighted = parser.get<bool>("weighted");
    build_config.verbose = parser.get<bool>("verbose");
    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    build_config.print();

    dict.build(input_filename, build_config);
    assert(dict.k() == k);

    bool check = parser.get<bool>("check");
    if (check) {
        check_correctness_lookup_access(dict, input_filename);
        check_correctness_navigational_kmer_query(dict, input_filename);
        check_correctness_navigational_contig_query(dict);
        if (build_config.weighted) check_correctness_weights(dict, input_filename);
        check_correctness_iterator(dict);
    }
    bool bench = parser.get<bool>("bench");
    if (bench) {
        perf_test_lookup_access(dict);
        if (dict.weighted()) perf_test_lookup_weight(dict);
        perf_test_iterator(dict);
    }
    if (parser.parsed("output_filename")) {
        auto output_filename = parser.get<std::string>("output_filename");
        essentials::logger("saving data structure to disk...");
        essentials::save(dict, output_filename.c_str());
        essentials::logger("DONE");
    }

    return 0;
}

int query(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with src/build.cpp.", "-i", true);
    parser.add("query_filename",
               "Must be a FASTA/FASTQ file (.fa/fasta or .fq/fastq extension) compressed with gzip "
               "or not.",
               "-q", true);
    parser.add("multiline",
               "Use this option if more the one DNA line must be parsed after each header."
               " Only valid for FASTA files (not FASTQ).",
               "--multiline", false, true);
    parser.add("print_index_info", "Print index information.", "--print-index-info", false, true);
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");
    auto query_filename = parser.get<std::string>("query_filename");

    dictionary dict;
    essentials::logger("loading index from file '" + index_filename + "'...");
    uint64_t num_bytes_read = essentials::load(dict, index_filename.c_str());
    std::cout << "index size: " << essentials::convert(num_bytes_read, essentials::MB) << " [MB] ("
              << (num_bytes_read * 8.0) / dict.size() << " [bits/kmer])" << std::endl;
    if (parser.get<bool>("print_index_info")) dict.print_info();

    bool multiline = parser.get<bool>("multiline");

    essentials::logger("performing queries from file '" + query_filename + "'...");
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> t;
    t.start();
    auto report = dict.streaming_query_from_file(query_filename, multiline);
    t.stop();
    essentials::logger("DONE");

    std::cout << "==== query report:\n";
    std::cout << "num_kmers = " << report.num_kmers << std::endl;
    std::cout << "num_positive_kmers = " << report.num_positive_kmers << " ("
              << (report.num_positive_kmers * 100.0) / report.num_kmers << "%)" << std::endl;
    std::cout << "num_searches = " << report.num_searches << "/" << report.num_positive_kmers
              << " (" << (report.num_searches * 100.0) / report.num_positive_kmers << "%)"
              << std::endl;
    std::cout << "num_extensions = " << report.num_extensions << "/" << report.num_positive_kmers
              << " (" << (report.num_extensions * 100.0) / report.num_positive_kmers << "%)"
              << std::endl;
    std::cout << "elapsed = " << t.elapsed() / 1000 << " millisec / ";
    std::cout << t.elapsed() / 1000000 << " sec / ";
    std::cout << t.elapsed() / 1000000 / 60 << " min / ";
    std::cout << (t.elapsed() * 1000) / report.num_kmers << " ns/kmer" << std::endl;

    return 0;
}

int check(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with src/build.cpp.", "-i", true);
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");

    dictionary dict;
    uint64_t num_bytes_read = essentials::load(dict, index_filename.c_str());
    std::cout << "index size: " << essentials::convert(num_bytes_read, essentials::MB) << " [MB] ("
              << (num_bytes_read * 8.0) / dict.size() << " [bits/kmer])" << std::endl;
    dict.print_info();

    check_dictionary(dict);
    check_correctness_navigational_contig_query(dict);

    return 0;
}

int bench(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with src/build.cpp.", "-i", true);
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");

    dictionary dict;
    uint64_t num_bytes_read = essentials::load(dict, index_filename.c_str());
    std::cout << "index size: " << essentials::convert(num_bytes_read, essentials::MB) << " [MB] ("
              << (num_bytes_read * 8.0) / dict.size() << " [bits/kmer])" << std::endl;
    dict.print_info();

    perf_test_lookup_access(dict);
    if (dict.weighted()) perf_test_lookup_weight(dict);
    perf_test_iterator(dict);

    return 0;
}

int dump(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with src/build.cpp.", "-i", true);
    parser.add("output_filename", "A FASTA file where the output will be saved.", "-o", true);
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");
    auto output_filename = parser.get<std::string>("output_filename");

    dictionary dict;
    uint64_t num_bytes_read = essentials::load(dict, index_filename.c_str());
    std::cout << "index size: " << essentials::convert(num_bytes_read, essentials::MB) << " [MB] ("
              << (num_bytes_read * 8.0) / dict.size() << " [bits/kmer])" << std::endl;
    dict.print_info();
    dict.dump(output_filename);

    return 0;
}

int help(char* arg0) {
    std::cout << "(S)parse and (S)kew (Hash)ing of k-mers" << std::endl << std::endl;
    std::cout << "Usage: " << arg0 << " <tool> ...\n\n"
              << "Available tools:\n"
              << "  build   \t build a dictionary \n"
              << "  query   \t query a dictionary \n"
              << "  check   \t check correctness of a dictionary \n"
              << "  bench   \t run performance tests for a dictionary \n"
              << "  dump    \t write super-k-mers of a dictionary to a fasta file \n"
              << "  permute \t permute an weighted input file" << std::endl;
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
    }
    std::cout << "Unsupported tool '" << tool << "'.\n" << std::endl;
    return help(argv[0]);
}