#include "include/builder/dictionary_builder.hpp"

#include "encode_strings.cpp"
#include "compute_minimizer_tuples.cpp"
#include "build_sparse_and_skew_index.cpp"

namespace sshash {

namespace {

inline void validate_and_normalize_build_config(build_configuration& bc, uint64_t max_k,
                                                uint64_t max_m) {
    if (bc.k == 0) throw std::runtime_error("k must be > 0");
    if (bc.k > max_k) {
        throw std::runtime_error("k must be less <= " + std::to_string(max_k) +
                                 " but got k = " + std::to_string(bc.k));
    }
    if (bc.m == 0) throw std::runtime_error("m must be > 0");
    if (bc.m > max_m) {
        throw std::runtime_error("m must be less <= " + std::to_string(max_m) +
                                 " but got m = " + std::to_string(bc.m));
    }
    if (bc.m > bc.k) throw std::runtime_error("m must be <= k");

    /* Clamp --ram-limit to the floor. Below this, the streaming buffers
       plus pthash's internal working memory can't usefully be made to
       fit; rather than try to squeeze further we treat the floor as the
       effective budget. */
    if (bc.ram_limit_in_GiB < constants::min_ram_limit_in_GiB) {
        if (bc.verbose) {
            std::cout << "  --> NOTE: --ram-limit raised from " << bc.ram_limit_in_GiB
                      << " GiB to the floor of " << constants::min_ram_limit_in_GiB << " GiB"
                      << std::endl;
        }
        bc.ram_limit_in_GiB = constants::min_ram_limit_in_GiB;
    }
}

}  // namespace

template <typename Kmer, typename Offsets>
void dictionary<Kmer, Offsets>::build(std::string const& filename,
                                      build_configuration const& build_config)  //
{
    build_configuration bc = build_config;
    validate_and_normalize_build_config(bc, Kmer::max_k, Kmer::max_m);
    dictionary_builder<Kmer, Offsets> builder(bc);
    builder.build(*this, filename);
}

template <typename Kmer, typename Offsets>
void dictionary<Kmer, Offsets>::build_streaming_save(std::string const& input_filename,
                                                     build_configuration const& build_config,
                                                     std::string const& output_filename)  //
{
    build_configuration bc = build_config;
    validate_and_normalize_build_config(bc, Kmer::max_k, Kmer::max_m);
    dictionary_builder<Kmer, Offsets> builder(bc);
    builder.build_streaming_save(*this, input_filename, output_filename);
}

}  // namespace sshash
