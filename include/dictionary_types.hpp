#pragma once

#include "dictionary.hpp"
#include "offsets.hpp"
#include "kmer.hpp"

namespace sshash {

using dictionary_type = dictionary<default_kmer_t, decoded_offsets>;
// using dictionary_type = dictionary<default_kmer_t, encoded_offsets>;

}  // namespace sshash