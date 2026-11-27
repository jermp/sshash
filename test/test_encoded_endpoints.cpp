#include "include/util.hpp"
#include "include/offsets.hpp"

using namespace sshash;

int main() {
    // Cover small and large endpoints without allocating the sequences.
    for (unsigned exponent : {2, 3, 16, 32, 40}) {
        for (int delta : {-1, 0, 1}) {
            uint64_t end = (uint64_t(1) << exponent) + delta;
            num_bits nb;
            nb.per_absolute_offset = bits::util::ceil_log2_uint64(end);
            encoded_offsets::builder builder;
            builder.set_num_bits(nb);
            builder.push_back(0);
            builder.push_back(end / 2);
            builder.push_back(end);
            encoded_offsets offsets;
            builder.build(offsets);
            if (offsets.access(0) != 0 || offsets.access(1) != end / 2 ||
                offsets.access(2) != end) {
                std::cerr << "incorrect encoded endpoint " << end << std::endl;
                return 1;
            }
        }
    }
}
