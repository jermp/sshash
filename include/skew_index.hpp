#pragma once

#include "util.hpp"

namespace sshash {

struct skew_index {
    skew_index()
        : min_log2(constants::min_l)
        , max_log2(constants::max_l)
        , log2_max_num_super_kmers_in_bucket(0) {
        mphfs.resize(0);
        positions.resize(0);
    }

    /* Returns the number of k-mers in the index. */
    uint64_t print_info() const;

    bool empty() const { return mphfs.empty(); }

    uint64_t lookup(kmer_t uint_kmer, uint64_t log2_bucket_size) const {
        assert(log2_bucket_size >= uint64_t(min_log2 + 1));
        assert(log2_bucket_size <= log2_max_num_super_kmers_in_bucket);
        uint64_t partition_id = log2_bucket_size - (min_log2 + 1);
        if (log2_bucket_size == log2_max_num_super_kmers_in_bucket or log2_bucket_size > max_log2) {
            partition_id = mphfs.size() - 1;
        }
        assert(partition_id < mphfs.size());
        auto const& mphf = mphfs[partition_id];
        auto const& P = positions[partition_id];
        uint64_t position = P.access(mphf(uint_kmer));
        return position;
    }

    uint64_t num_bits() const {
        uint64_t n =
            (sizeof(min_log2) + sizeof(max_log2) + sizeof(log2_max_num_super_kmers_in_bucket) +
             2 * sizeof(size_t) /* for std::vector::size */) *
            8;
        for (uint64_t partition_id = 0; partition_id != mphfs.size(); ++partition_id) {
            auto const& mphf = mphfs[partition_id];
            auto const& P = positions[partition_id];
            n += mphf.num_bits() + P.bytes() * 8;
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
    uint32_t log2_max_num_super_kmers_in_bucket;
    std::vector<kmers_pthash_type> mphfs;
    std::vector<pthash::compact_vector> positions;

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.min_log2);
        visitor.visit(t.max_log2);
        visitor.visit(t.log2_max_num_super_kmers_in_bucket);
        visitor.visit(t.mphfs);
        visitor.visit(t.positions);
    }
};

}  // namespace sshash