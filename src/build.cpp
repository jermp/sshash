#include "include/builder/dictionary_builder.hpp"

#include "include/builder/parse_file.cpp"
#include "include/builder/build_sparse_and_skew_index.cpp"
#include "include/builder/compute_minimizer_tuples.cpp"

namespace sshash {

template <typename Kmer, typename Offsets>
void dictionary<Kmer, Offsets>::build(std::string const& filename,
                                      build_configuration const& build_config)  //
{
    /* Validate the build configuration. */
    if (build_config.k == 0) throw std::runtime_error("k must be > 0");
    if (build_config.k > Kmer::max_k) {
        throw std::runtime_error("k must be less <= " + std::to_string(Kmer::max_k) +
                                 " but got k = " + std::to_string(build_config.k));
    }
    if (build_config.m == 0) throw std::runtime_error("m must be > 0");
    if (build_config.m > Kmer::max_m) {
        throw std::runtime_error("m must be less <= " + std::to_string(Kmer::max_m) +
                                 " but got m = " + std::to_string(build_config.m));
    }
    if (build_config.m > build_config.k) throw std::runtime_error("m must be <= k");

    dictionary_builder<Kmer, Offsets> builder(build_config);
    builder.build(*this, filename);
}

}  // namespace sshash
