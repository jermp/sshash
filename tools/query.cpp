#include "common.hpp"
#include "include/streaming_query.hpp"

using namespace sshash;

template <class kmer_t>
int query(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with the tool 'build'.", "-i", true);
    parser.add("query_filename",
               "Must be a FASTA/FASTQ file (.fa/fasta or .fq/fastq extension) compressed with gzip "
               "or not.",
               "-q", true);
    parser.add("multiline",
               "Use this option if more the one DNA line must be parsed after each header."
               " Only valid for FASTA files (not FASTQ).",
               "--multiline", false, true);
    parser.add("verbose", "Verbose output.", "--verbose", false, true);
    if (!parser.parse()) return 1;

    auto index_filename = parser.get<std::string>("index_filename");
    auto query_filename = parser.get<std::string>("query_filename");
    bool verbose = parser.get<bool>("verbose");
    bool multiline = parser.get<bool>("multiline");

    dictionary<kmer_t> dict;
    load_dictionary(dict, index_filename, verbose);

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
    std::cout << "num_negative_kmers = " << report.num_negative_kmers << " ("
              << (report.num_negative_kmers * 100.0) / report.num_kmers << "%)" << std::endl;
    std::cout << "num_invalid_kmers = " << report.num_invalid_kmers << " ("
              << (report.num_invalid_kmers * 100.0) / report.num_kmers << "%)" << std::endl;
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