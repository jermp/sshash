#pragma once

#include "util.hpp"

namespace sshash {

template <class kmer_t>
struct skew_index {
    skew_index() : min_log2(constants::min_l), max_log2(constants::max_l), log2_max_bucket_size(0) {
        mphfs.resize(0);
        positions.resize(0);
    }

    /* Returns the number of kmers in the index. */
    uint64_t print_info() const {
        uint64_t num_partitions = mphfs.size();
        uint64_t lower = 1ULL << min_log2;
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

    bool empty() const { return mphfs.empty(); }

    uint64_t lookup(kmer_t uint_kmer, uint64_t log2_bucket_size) const {
        assert(log2_bucket_size >= uint64_t(min_log2 + 1));
        assert(log2_bucket_size <= log2_max_bucket_size);
        uint64_t partition_id = log2_bucket_size - (min_log2 + 1);
        if (log2_bucket_size == log2_max_bucket_size or log2_bucket_size > max_log2) {
            partition_id = mphfs.size() - 1;
        }
        assert(partition_id < mphfs.size());
        auto const& f = mphfs[partition_id];
        auto const& p = positions[partition_id];
        uint64_t position = p.access(f(uint_kmer));
        return position;
    }

    uint64_t num_bits() const {
        uint64_t n = (sizeof(min_log2) + sizeof(max_log2) + sizeof(log2_max_bucket_size) +
                      2 * sizeof(size_t) /* for std::vector::size */) *
                     8;
        for (uint64_t partition_id = 0; partition_id != mphfs.size(); ++partition_id) {
            auto const& f = mphfs[partition_id];
            auto const& p = positions[partition_id];
            n += f.num_bits() + p.num_bytes() * 8;
        }
        return n;
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

    uint16_t min_log2;
    uint16_t max_log2;
    uint32_t log2_max_bucket_size;
    std::vector<kmers_pthash_type<kmer_t>> mphfs;
    std::vector<bits::compact_vector> positions;

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.min_log2);
        visitor.visit(t.max_log2);
        visitor.visit(t.log2_max_bucket_size);
        visitor.visit(t.mphfs);
        visitor.visit(t.positions);
    }
};

}  // namespace sshash
