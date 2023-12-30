#pragma once

#include "kmer.hpp"

namespace sshash::constants {

constexpr uint64_t uint_kmer_bits = sizeof(kmer_t) * 8;

/* max *odd* size that can be packed into sizeof(kmer_t)*8 bits */
constexpr uint64_t max_k = sizeof(kmer_t) * 4 - 1;

/* max *odd* size that can be packed into 64 bits */
constexpr uint64_t max_m = 31;

constexpr uint64_t invalid_uint64 = uint64_t(-1);

constexpr uint64_t seed = 1;
constexpr double c = 3.0;  // for PTHash

constexpr uint64_t min_l = 6;
constexpr uint64_t max_l = 12;
static_assert(min_l < max_l);

static const std::string default_tmp_dirname(".");
constexpr bool forward_orientation = 0;
constexpr bool backward_orientation = 1;

}  // namespace sshash::constants