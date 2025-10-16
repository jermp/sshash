#pragma once

#include "minimizers_control_map.hpp"
#include "skew_index.hpp"

namespace sshash {

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
        uint64_t offset = skew_index.lookup(uint_kmer, code);
        return {this, offset, 1, bucket_t::HEAVYLOAD};
    }

    uint64_t num_bits() const {
        return codewords.num_bits() +
               8 * (essentials::vec_bytes(begin_buckets_of_size) + mid_load_buckets.num_bytes()) +
               skew_index.num_bits();
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
    skew_index<Kmer> skew_index;

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.codewords);
        visitor.visit(t.begin_buckets_of_size);
        visitor.visit(t.mid_load_buckets);
        visitor.visit(t.skew_index);
    }
};

}  // namespace sshash