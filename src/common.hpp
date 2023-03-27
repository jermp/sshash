#pragma once

#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/dictionary.hpp"

#include <vector>

namespace sshash {

void random_kmer(char* kmer, uint64_t k) {
    for (uint64_t i = 0; i != k; ++i) kmer[i] = "ACGT"[rand() % 4];
}

void load_dictionary(dictionary& dict, std::string const& index_filename, bool verbose) {
    uint64_t num_bytes_read = essentials::load(dict, index_filename.c_str());
    if (verbose) {
        std::cout << "index size: " << essentials::convert(num_bytes_read, essentials::MB)
                  << " [MB] (" << (num_bytes_read * 8.0) / dict.size() << " [bits/kmer])"
                  << std::endl;
        dict.print_info();
    }
}

}  // namespace sshash