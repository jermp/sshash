#include "src/builder/build.cpp"

using namespace sshash;

void check(uint64_t count, uint64_t expected) {
    // Exercise the builder's calculation without allocating enormous strings.
    auto nb = util::offset_widths(count, count, 1, count);
    if (nb.per_absolute_offset != expected || nb.per_relative_offset != expected ||
        nb.per_string_id != expected) {
        throw std::runtime_error("incorrect offset width for count " + std::to_string(count));
    }
    // Relative positions range from 0 through max_len - m, inclusive.
    nb = util::offset_widths(100, 20, 5, 3);
    if (nb.per_absolute_offset != 7 || nb.per_relative_offset != 4 || nb.per_string_id != 2) {
        throw std::runtime_error("incorrect offset ranges");
    }
}

int main() {
    try {
        check(1, 0);
        check(2, 1);
        // In particular, log2(2^49 + 1) can round down to exactly 49.
        for (unsigned exponent = 2; exponent < 64; ++exponent) {
            uint64_t power = uint64_t(1) << exponent;
            check(power - 1, exponent);
            check(power, exponent);
            check(power + 1, exponent + 1);
        }
        check(UINT64_MAX, 64);
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
