#pragma once

#include "util.hpp"
#include "bit_vector_iterator.hpp"
#include "ef_sequence.hpp"

namespace sshash {

struct buckets {
    std::pair<lookup_result, uint64_t> offset_to_id(uint64_t offset, uint64_t k) const {
        auto [pos, contig_begin, contig_end] = pieces.locate(offset);

        /* The following two facts hold. */
        assert(pieces.access(pos) == contig_end);
        assert(contig_end >= offset);

        bool shift = contig_end > offset;
        assert(pos >= shift);
        uint64_t contig_id = pos - shift;

        if (contig_end == offset) {
            assert(shift == false);
            assert(pos + 1 < pieces.size());
            contig_begin = contig_end;
            contig_end = pieces.access(pos + 1);
        }

        /* Now, the following facts hold. */
        assert(offset >= contig_id * (k - 1));
        assert(contig_begin <= offset);
        assert(offset < contig_end);

        uint64_t absolute_kmer_id = offset - contig_id * (k - 1);
        uint64_t relative_kmer_id = offset - contig_begin;
        uint64_t contig_length = contig_end - contig_begin;
        assert(contig_length >= k);
        uint64_t contig_size = contig_length - k + 1;

        lookup_result res;
        res.kmer_id = absolute_kmer_id;
        res.kmer_id_in_contig = relative_kmer_id;
        res.contig_id = contig_id;
        res.contig_size = contig_size;

        return {res, contig_end};
    }

    /* Return where the contig begins and ends in strings. */
    std::pair<uint64_t, uint64_t>  // [begin, end)
    contig_offsets(const uint64_t contig_id) const {
        uint64_t begin = pieces.access(contig_id);
        uint64_t end = pieces.access(contig_id + 1);
        assert(end > begin);
        return {begin, end};
    }

    kmer_t contig_prefix(uint64_t contig_id, uint64_t k) const {
        uint64_t contig_begin = pieces.access(contig_id);
        bit_vector_iterator bv_it(strings, 2 * contig_begin);
        return bv_it.read(2 * (k - 1));
    }

    kmer_t contig_suffix(uint64_t contig_id, uint64_t k) const {
        uint64_t contig_end = pieces.access(contig_id + 1);
        bit_vector_iterator bv_it(strings, 2 * (contig_end - k + 1));
        return bv_it.read(2 * (k - 1));
    }

    std::pair<uint64_t, uint64_t> locate_bucket(uint64_t bucket_id) const {
        uint64_t begin = num_super_kmers_before_bucket.access(bucket_id) + bucket_id;
        uint64_t end = num_super_kmers_before_bucket.access(bucket_id + 1) + bucket_id + 1;
        assert(begin < end);
        return {begin, end};
    }

    lookup_result lookup(uint64_t bucket_id, kmer_t target_kmer, uint64_t k, uint64_t m) const {
        auto [begin, end] = locate_bucket(bucket_id);
        return lookup(begin, end, target_kmer, k, m);
    }

    lookup_result lookup(uint64_t begin, uint64_t end, kmer_t target_kmer, uint64_t k,
                         uint64_t m) const {
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            auto res = lookup_in_super_kmer(super_kmer_id, target_kmer, k, m);
            if (res.kmer_id != constants::invalid_uint64) {
                assert(is_valid(res));
                return res;
            }
        }
        return lookup_result();
    }

    lookup_result lookup_in_super_kmer(uint64_t super_kmer_id, kmer_t target_kmer, uint64_t k,
                                       uint64_t m) const {
        uint64_t offset = offsets.access(super_kmer_id);
        auto [res, contig_end] = offset_to_id(offset, k);
        bit_vector_iterator bv_it(strings, 2 * offset);
        uint64_t window_size = std::min<uint64_t>(k - m + 1, contig_end - offset - k + 1);
        for (uint64_t w = 0; w != window_size; ++w) {
            kmer_t read_kmer = bv_it.read_and_advance_by_two(2 * k);
            if (read_kmer == target_kmer) {
                res.kmer_id += w;
                res.kmer_id_in_contig += w;
                assert(is_valid(res));
                return res;
            }
        }
        return lookup_result();
    }

    lookup_result lookup_canonical(uint64_t bucket_id, kmer_t target_kmer, kmer_t target_kmer_rc,
                                   uint64_t k, uint64_t m) const {
        auto [begin, end] = locate_bucket(bucket_id);
        return lookup_canonical(begin, end, target_kmer, target_kmer_rc, k, m);
    }

    lookup_result lookup_canonical(uint64_t begin, uint64_t end, kmer_t target_kmer,
                                   kmer_t target_kmer_rc, uint64_t k, uint64_t m) const {
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            uint64_t offset = offsets.access(super_kmer_id);
            auto [res, contig_end] = offset_to_id(offset, k);
            bit_vector_iterator bv_it(strings, 2 * offset);
            uint64_t window_size = std::min<uint64_t>(k - m + 1, contig_end - offset - k + 1);
            for (uint64_t w = 0; w != window_size; ++w) {
                kmer_t read_kmer = bv_it.read_and_advance_by_two(2 * k);
                if (read_kmer == target_kmer) {
                    res.kmer_id += w;
                    res.kmer_id_in_contig += w;
                    res.kmer_orientation = constants::forward_orientation;
                    assert(is_valid(res));
                    return res;
                }
                if (read_kmer == target_kmer_rc) {
                    res.kmer_id += w;
                    res.kmer_id_in_contig += w;
                    res.kmer_orientation = constants::backward_orientation;
                    assert(is_valid(res));
                    return res;
                }
            }
        }
        return lookup_result();
    }

    uint64_t id_to_offset(uint64_t id, uint64_t k) const {
        constexpr uint64_t linear_scan_threshold = 8;
        uint64_t lo = 0;
        uint64_t hi = pieces.size() - 1;
        assert(pieces.access(0) == 0);
        while (lo < hi) {
            if (hi - lo <= linear_scan_threshold) {
                for (; lo < hi; ++lo) {
                    uint64_t val = pieces.access(lo) - lo * (k - 1);
                    if (val > id) break;
                }
                break;
            }
            uint64_t mid = lo + (hi - lo) / 2;
            uint64_t val = pieces.access(mid);
            assert(val >= mid * (k - 1));
            if (id <= val - mid * (k - 1)) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        if (lo < pieces.size() and pieces.access(lo) - lo * (k - 1) > id) --lo;
        return id + lo * (k - 1);
    }

    void access(uint64_t kmer_id, char* string_kmer, uint64_t k) const {
        uint64_t offset = id_to_offset(kmer_id, k);
        bit_vector_iterator bv_it(strings, 2 * offset);
        kmer_t read_kmer = bv_it.read(2 * k);
        util::uint_kmer_to_string(read_kmer, string_kmer, k);
    }

    struct iterator {
        iterator() {}

        iterator(buckets const* ptr,                                        //
                 const uint64_t begin_kmer_id, const uint64_t end_kmer_id,  // [begin,end)
                 const uint64_t k)
            : m_buckets(ptr)
            , m_begin_kmer_id(begin_kmer_id)
            , m_end_kmer_id(end_kmer_id)
            , m_k(k)  //
        {
            m_bv_it = bit_vector_iterator(m_buckets->strings, -1);
            m_offset = m_buckets->id_to_offset(m_begin_kmer_id, k);
            auto [pos, piece_end] = m_buckets->pieces.next_geq(m_offset);
            if (piece_end == m_offset) pos += 1;
            m_pieces_it = m_buckets->pieces.at(pos);
            next_piece();
            m_ret.second.resize(m_k, 0);
        }

        bool has_next() const { return m_begin_kmer_id != m_end_kmer_id; }

        std::pair<uint64_t, std::string> next() {
            if (m_offset == m_next_offset - m_k + 1) {
                m_offset = m_next_offset;
                next_piece();
            }

            while (m_offset != m_next_offset - m_k + 1) {
                m_ret.first = m_begin_kmer_id;
                if (m_clear) {
                    util::uint_kmer_to_string(m_read_kmer, m_ret.second.data(), m_k);
                } else {
                    memmove(m_ret.second.data(), m_ret.second.data() + 1, m_k - 1);
                    m_ret.second[m_k - 1] = util::uint64_to_char(m_last_two_bits);
                }
                m_clear = false;
                m_read_kmer >>= 2;
                m_last_two_bits = m_bv_it.get_next_two_bits();
                m_read_kmer += m_last_two_bits << (2 * (m_k - 1));
                ++m_begin_kmer_id;
                ++m_offset;
                return m_ret;
            }

            return next();
        }

    private:
        std::pair<uint64_t, std::string> m_ret;
        buckets const* m_buckets;
        uint64_t m_begin_kmer_id, m_end_kmer_id;
        uint64_t m_k;
        uint64_t m_offset;
        uint64_t m_next_offset;
        bit_vector_iterator m_bv_it;
        ef_sequence<true>::iterator m_pieces_it;

        kmer_t m_read_kmer;
        uint64_t m_last_two_bits;
        bool m_clear;

        void next_piece() {
            m_bv_it.at(2 * m_offset);
            m_next_offset = m_pieces_it.next();
            assert(m_next_offset > m_offset);
            m_read_kmer = m_bv_it.take(2 * m_k);
            m_clear = true;
        }
    };

    iterator at(const uint64_t begin_kmer_id, const uint64_t end_kmer_id, const uint64_t k) const {
        return iterator(this, begin_kmer_id, end_kmer_id, k);
    }

    uint64_t num_bits() const {
        return pieces.num_bits() + num_super_kmers_before_bucket.num_bits() +
               8 * (offsets.bytes() + strings.bytes());
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(pieces);
        visitor.visit(num_super_kmers_before_bucket);
        visitor.visit(offsets);
        visitor.visit(strings);
    }

    ef_sequence<true> pieces;
    ef_sequence<false> num_super_kmers_before_bucket;
    pthash::compact_vector offsets;
    pthash::bit_vector strings;

private:
    bool is_valid(lookup_result res) const {
        return (res.contig_size != constants::invalid_uint64 and
                res.kmer_id_in_contig < res.contig_size) and
               (res.contig_id != constants::invalid_uint64 or res.contig_id < pieces.size());
    }
};

}  // namespace sshash