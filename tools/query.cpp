#include "common.hpp"
#include "include/streaming_query.hpp"

using namespace sshash;

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
    if (!parser.parse()) return 0;

    auto index_filename = parser.get<std::string>("index_filename");
    auto query_filename = parser.get<std::string>("query_filename");
    bool verbose = parser.get<bool>("verbose");
    bool multiline = parser.get<bool>("multiline");

    dictionary_type dict;
    load_dictionary(dict, index_filename, verbose);

    essentials::logger("performing queries from file '" + query_filename + "'...");
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::milliseconds> t;
    t.start();
    auto report = dict.streaming_query_from_file(query_filename, multiline);
    t.stop();
    essentials::logger("DONE");

    essentials::json_lines query_stats;
    query_stats.add("index_filename", index_filename.c_str());
    query_stats.add("query_filename", query_filename.c_str());
    query_stats.add("num_kmers", report.num_kmers);
    query_stats.add("num_positive_kmers", report.num_positive_kmers);
    query_stats.add("num_negative_kmers", report.num_negative_kmers);
    query_stats.add("num_invalid_kmers", report.num_invalid_kmers);
    query_stats.add("num_searches", report.num_searches);
    query_stats.add("num_extensions", report.num_extensions);
    query_stats.add("elapsed_millisec", uint64(t.elapsed()));

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
    std::cout << "elapsed = " << t.elapsed() / 1000 << " sec / ";
    std::cout << t.elapsed() / 1000 / 60 << " min / ";
    std::cout << (t.elapsed() * 1e6) / report.num_kmers << " ns/kmer" << std::endl;

    query_stats.print();

    return 0;
}