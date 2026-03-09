#pragma once

#include "util.hpp"
#include "minimizers_control_map.hpp"

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

template <typename Kmer>
struct sparse_and_skew_index  //
{
    struct bucket_iterator {
        bucket_iterator(sparse_and_skew_index const* ssi, uint64_t pos, uint64_t size,
                        bucket_t bucket_type)
            : m_size(size)
            , m_offset(constants::invalid_uint64)
            , m_bucket_type(bucket_type)  //
        {
            assert(size > 0);
            if (size == 1) {
                m_offset = pos;
            } else {
                m_it = ssi->mid_load_buckets.get_iterator_at(pos);
            }
        }

        uint64_t operator*() const { return m_size == 1 ? m_offset : *m_it; }
        void operator++() {
            if (m_size != 1) ++m_it;
        }

        uint64_t size() const { return m_size; }
        bucket_t bucket_type() const { return m_bucket_type; }

    private:
        uint64_t m_size;
        uint64_t m_offset;
        bucket_t m_bucket_type;
        bits::compact_vector::iterator m_it;
    };

    bucket_iterator lookup(const Kmer uint_kmer, const minimizer_info mini_info) const  //
    {
        uint64_t code = codewords.lookup(mini_info.minimizer);

        uint64_t status = code & 1;
        if (status == bucket_t::SINGLETON) {  // minimizer occurs once
            uint64_t offset = code >> 1;
            return {this, offset, 1, bucket_t::SINGLETON};
        }

        status = code & 3;
        if (status == bucket_t::MIDLOAD) {  // minimizer occurs more than once, but is not part of
                                            // the skew index
            constexpr uint64_t mask = (uint64_t(1) << constants::min_l) - 1;
            code >>= 2;
            uint64_t bucket_size = (code & mask) + 2;
            uint64_t bucket_id = code >> constants::min_l;
            assert(bucket_size < begin_buckets_of_size.size());
            uint64_t begin = begin_buckets_of_size[bucket_size] + bucket_id * bucket_size;
            return {this, begin, bucket_size, bucket_t::MIDLOAD};
        }

        assert(status == bucket_t::HEAVYLOAD);  // minimizer is part of the skew index
        uint64_t offset = ski.lookup(uint_kmer, code);
        return {this, offset, 1, bucket_t::HEAVYLOAD};
    }

    uint64_t num_bits() const {
        return codewords.num_bits() +
               8 * (essentials::vec_bytes(begin_buckets_of_size) + mid_load_buckets.num_bytes()) +
               ski.num_bits();
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

    minimizers_control_map codewords;
    std::vector<uint32_t> begin_buckets_of_size;
    bits::compact_vector mid_load_buckets;
    skew_index<Kmer> ski;

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.codewords);
        visitor.visit(t.begin_buckets_of_size);
        visitor.visit(t.mid_load_buckets);
        visitor.visit(t.ski);
    }
};

}  // namespace sshash