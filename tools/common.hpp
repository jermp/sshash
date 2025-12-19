#pragma once

#include "external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "include/dictionary_types.hpp"

#include <vector>

namespace sshash {

void print_cmd(int argc, char** argv) {
    for (int i = 0; i != argc; ++i) std::cout << argv[i] << ' ';
    std::cout << std::endl;
}

void random_kmer(char* kmer, uint64_t k) {
    for (uint64_t i = 0; i != k; ++i) kmer[i] = "ACGT"[rand() % 4];
}

template <typename Dict>
void load_dictionary(Dict& dict, std::string const& index_filename, bool verbose) {
    const uint64_t num_bytes_read = essentials::load(dict, index_filename.c_str());
    if (verbose) {
        std::cout << "total index size: " << num_bytes_read << " [B] -- "
                  << essentials::convert(num_bytes_read, essentials::MB) << " [MB] ("
                  << (num_bytes_read * 8.0) / dict.num_kmers() << " [bits/kmer])" << std::endl;
        dict.print_info();
    }
}

}  // namespace sshash