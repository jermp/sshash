#include <iostream>

#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/dictionary.hpp"
#include "check_utils.hpp"

using namespace sshash;

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with src/build.cpp");
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
