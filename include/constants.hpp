#pragma once

#include "kmer.hpp"

namespace sshash::constants {
constexpr uint64_t invalid_uint64 = uint64_t(-1);

constexpr uint64_t seed = 1;

/* for PTHash */
constexpr double lambda = 7.0;
constexpr uint64_t avg_partition_size = 3000000;

constexpr uint64_t min_l = 6;
constexpr uint64_t max_l = 12;
static_assert(min_l < max_l);

static const std::string default_tmp_dirname(".");
constexpr bool forward_orientation = 0;
constexpr bool backward_orientation = 1;

}  // namespace sshash::constants