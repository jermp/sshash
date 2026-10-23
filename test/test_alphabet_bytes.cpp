#include <iostream>
#include <stdexcept>

#include "include/kmer.hpp"

using namespace sshash;

void require(bool condition, char const* message) {
    if (!condition) throw std::runtime_error(message);
}

int main() {
    try {
        using dna = dna_uint_kmer_t<uint64_t>;
        using protein = aa_uint_kmer_t<uint64_t>;
        for (unsigned byte = 0; byte < 256; ++byte) {
            char c = static_cast<char>(byte);
            bool is_dna = std::string("ACGTacgt").find(c) != std::string::npos;
            bool is_protein = (byte >= 'A' && byte <= 'Z') || (byte >= 'a' && byte <= 'z');
            require(dna::is_valid(c) == is_dna, "incorrect DNA byte validation");
            require(protein::is_valid(c) == is_protein, "incorrect protein byte validation");

            char reverse = 1;
            dna::compute_reverse_complement(&c, &reverse, 1);
            if (!is_dna) require(reverse == 0, "invalid DNA byte has a complement");
        }

        std::string input = "ACGTacgt";
        std::string reverse(input.size(), '\0');
        dna::compute_reverse_complement(input.data(), reverse.data(), input.size());
        require(reverse == "acgtACGT", "ASCII reverse complement changed");
        for (uint64_t i = 0; i < sizeof(amino_acids) - 1; ++i) {
            char upper = amino_acids[i];
            char lower = upper + ('a' - 'A');
            require(protein::char_to_uint(upper) == i, "uppercase protein mapping changed");
            require(protein::char_to_uint(lower) == i, "lowercase protein mapping changed");
        }
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
