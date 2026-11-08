#pragma once

namespace sshash::constants {

constexpr uint64_t invalid_uint64 = uint64_t(-1);
constexpr uint64_t default_ram_limit_in_GiB = 8;
constexpr uint64_t seed = 1;

/* for PTHash */
constexpr double lambda = 5.0;
constexpr uint64_t avg_partition_size = 3000000;

constexpr uint64_t min_l = 6;
constexpr uint64_t max_l = 13;
static_assert(min_l < max_l);
static_assert(max_l - min_l + 1 <= 8);

static const std::string default_tmp_dirname(".");
constexpr int forward_orientation = 1;
constexpr int backward_orientation = -1;

namespace current_version_number {
constexpr uint8_t x = 6;
constexpr uint8_t y = 0;
constexpr uint8_t z = 0;
}  // namespace current_version_number

/* Indexes built with preserve_orientation use this major version. Readers that
   only know the current version reject them instead of treating them as
   canonical. Their layout is the current one plus a trailing orientation flag. */
constexpr uint8_t oriented_version_major = 7;

}  // namespace sshash::constants
