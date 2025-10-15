#pragma once

#include "util.hpp"

namespace sshash {

template <typename Kmer>
struct skew_index  //
{
    skew_index() {
        mphfs.resize(0);
        positions.resize(0);
    }

    /* Returns the number of kmers in the index. */
    uint64_t print_info() const {
        uint64_t num_partitions = mphfs.size();
        uint64_t lower = 1ULL << constants::min_l;
        uint64_t upper = 2 * lower;
        uint64_t num_kmers_in_skew_index = 0;
        for (uint64_t partition_id = 0; partition_id != num_partitions; ++partition_id) {
            uint64_t n = mphfs[partition_id].num_keys();
            assert(n == positions[partition_id].size());
            std::cout << "num_kmers belonging to buckets of size > " << lower << " and <= " << upper
                      << ": " << n << "; ";
            std::cout << "bits/kmer = " << static_cast<double>(mphfs[partition_id].num_bits()) / n
                      << " (mphf) + " << (positions[partition_id].num_bytes() * 8.0) / n
                      << " (positions)\n";
            num_kmers_in_skew_index += n;
            lower = upper;
            upper = 2 * lower;
        }
        return num_kmers_in_skew_index;
    }

    uint64_t lookup(const Kmer uint_kmer, uint64_t code) const {
        code >>= 2;
        uint64_t partition_id = code & 7;
        uint64_t begin = code >> 3;
        assert(partition_id < mphfs.size());
        auto const& f = mphfs[partition_id];
        auto const& p = positions[partition_id];
        uint64_t pos_in_bucket = p.access(f(uint_kmer));
        uint64_t offset = heavy_load_buckets.access(begin + pos_in_bucket);
        return offset;
    }

    uint64_t num_bits() const {
        uint64_t n = (2 * sizeof(size_t)) * 8; /* for std::vector::size */
        for (uint64_t partition_id = 0; partition_id != mphfs.size(); ++partition_id) {
            auto const& f = mphfs[partition_id];
            auto const& p = positions[partition_id];
            n += f.num_bits() + p.num_bytes() * 8;
        }
        return n + 8 * heavy_load_buckets.num_bytes();
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

    std::vector<kmers_pthash_type<Kmer>> mphfs;
    std::vector<bits::compact_vector> positions;
    bits::compact_vector heavy_load_buckets;

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.mphfs);
        visitor.visit(t.positions);
        visitor.visit(t.heavy_load_buckets);
    }
};

}  // namespace sshash
