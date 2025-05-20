#include <iostream>
#include <unordered_set>

#include "include/util.hpp"
#include "external/gz/zip_stream.hpp"
#include "external/pthash/external/cmd_line_parser/include/parser.hpp"

using namespace sshash;

std::unordered_set<uint64_t> parse_file(std::istream& is, const uint64_t k) {
    std::unordered_set<uint64_t> kmers;

    std::string sequence;
    uint64_t num_kmers = 0;
    uint64_t num_bases = 0;
    uint64_t num_sequences = 0;

    while (!is.eof())  //
    {
        std::getline(is, sequence);  // header sequence
        std::getline(is, sequence);  // DNA sequence

        if (sequence.length() < k) continue;

        if (++num_sequences % 100000 == 0) {
            std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
                      << num_kmers << " kmers" << std::endl;
        }

        num_bases += sequence.length();

        for (uint64_t end = 0; end != sequence.length() - k + 1; ++end) {
            char const* kmer = sequence.data() + end;
            assert(util::is_valid<default_kmer_t>(kmer, k));
            default_kmer_t uint_kmer = util::string_to_uint_kmer<default_kmer_t>(kmer, k);

            kmers.insert(uint_kmer.kmer);

            ++num_kmers;
        }
    }

    std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, " << num_kmers
              << " kmers" << std::endl;

    return kmers;
}

std::unordered_set<uint64_t> parse_file(std::string const& filename, const uint64_t k) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    zip_istream zis(is);
    auto kmers = parse_file(zis, k);
    is.close();
    return kmers;
}

void query_from_fastq_file(std::string const& query_filename,
                           std::unordered_set<uint64_t> const& kmers, const uint64_t k) {
    std::ifstream is(query_filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + query_filename + "'");
    zip_istream zis(is);
    std::string line;
    uint64_t num_positive_kmers = 0;
    uint64_t num_kmers = 0;
    while (!zis.eof())  //
    {
        /* We assume the file is well-formed, i.e., there are exactly 4 lines per read. */
        std::getline(zis, line);  // skip first header line
        std::getline(zis, line);
        if (line.size() >= k) {
            for (uint64_t i = 0; i != line.size() - k + 1; ++i) {
                char const* kmer = line.data() + i;
                if (util::is_valid<default_kmer_t>(kmer, k))  //
                {
                    default_kmer_t uint_kmer = util::string_to_uint_kmer<default_kmer_t>(kmer, k);
                    if (auto it = kmers.find(uint_kmer.kmer); it != kmers.end()) {
                        num_positive_kmers += 1;
                    }
                    default_kmer_t uint_kmer_rc = uint_kmer;
                    uint_kmer_rc.reverse_complement_inplace(k);
                    if (auto it = kmers.find(uint_kmer_rc.kmer); it != kmers.end()) {
                        num_positive_kmers += 1;
                    }
                }
                num_kmers += 1;
            }
        }
        std::getline(zis, line);  // skip '+'
        std::getline(zis, line);  // skip score
    }
    is.close();
    std::cout << "num_kmers " << num_kmers << std::endl;
    std::cout << "num_positive_kmers " << num_positive_kmers << std::endl;
}

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("input_filename", "Must be in the same input format as that used by 'sshash build'.",
               "-i", true);
    parser.add("query_filename", "Must be a FASTQ file (.fq/fastq extension) compressed with gzip.",
               "-q", true);
    parser.add("k", "K-mer length (must be <= 31).", "-k", true);
    if (!parser.parse()) return 1;

    auto k = parser.get<uint64_t>("k");
    if (k > 31) {
        std::cerr << "K-mer length (must be <= 31)." << std::endl;
        return 1;
    }
    auto input_filename = parser.get<std::string>("input_filename");
    auto query_filename = parser.get<std::string>("query_filename");

    essentials::logger("reading input kmers...");
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> t;
    t.start();
    auto kmers = parse_file(input_filename, k);
    t.stop();
    essentials::logger("DONE");

    essentials::logger("reading input kmers...");
    t.reset();
    t.start();
    query_from_fastq_file(query_filename, kmers, k);
    t.stop();
    essentials::logger("DONE");

    return 0;
}