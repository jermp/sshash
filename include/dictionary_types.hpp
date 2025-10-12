#pragma once

#include "dictionary.hpp"
#include "endpoints.hpp"
#include "kmer.hpp"

namespace sshash {

// using dictionary_type = dictionary<default_kmer_t, endpoints>;
using dictionary_type = dictionary<default_kmer_t, _endpoints>;

}  // namespace sshash