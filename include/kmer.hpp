#pragma once

namespace sshash {

#ifdef SSHASH_USE_MAX_KMER_LENGTH_63
typedef __uint128_t kmer_t;
#else
typedef uint64_t kmer_t;
#endif

}  // namespace sshash