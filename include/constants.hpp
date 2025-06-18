#pragma once

#include "kmer.hpp"

namespace sshash::constants {

constexpr uint64_t invalid_uint64 = uint64_t(-1);
constexpr uint64_t default_ram_limit_in_GiB = 8;
constexpr uint64_t seed = 1;

/* for PTHash */
constexpr double lambda = 7.0;
constexpr uint64_t avg_partition_size = 3000000;

constexpr uint64_t min_l = 6;
constexpr uint64_t max_l = 12;
static_assert(min_l < max_l);

static const std::string default_tmp_dirname(".");
constexpr int forward_orientation = 1;
constexpr int backward_orientation = -1;

namespace current_version_number {
constexpr uint8_t x = 4;
constexpr uint8_t y = 0;
constexpr uint8_t z = 0;
}  // namespace current_version_number

}  // namespace sshash::constants