#pragma once

#include <vector>

namespace sshash {

void random_kmer(char* kmer, uint64_t k) {
    for (uint64_t i = 0; i != k; ++i) kmer[i] = "ACGT"[rand() % 4];
}

}  // namespace sshash