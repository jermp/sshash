#pragma once

#include "util.hpp"
#include "bit_vector_iterator.hpp"
#include "ef_sequence.hpp"

namespace sshash {

struct buckets {
    std::pair<uint64_t, uint64_t> offset_to_id(uint64_t offset, uint64_t k) const {
        auto [pos, piece_end] = pieces.next_geq(offset);
        uint64_t p = pos - (piece_end > offset);
        if (piece_end == offset) {
            assert(pos + 1 < pieces.size());
            piece_end = pieces.access(pos + 1);
        }
        assert(offset >= p * (k - 1));
        assert(piece_end > offset);
        return {offset - p * (k - 1), piece_end};
    }

    uint64_t contig_length(uint64_t contig_id) const {
        uint64_t length = pieces.access(contig_id + 1) - pieces.access(contig_id);
        return length;
    }

    std::pair<uint64_t, uint64_t> locate_bucket(uint64_t bucket_id) const {
        uint64_t begin = num_super_kmers_before_bucket.access(bucket_id) + bucket_id;
        uint64_t end = num_super_kmers_before_bucket.access(bucket_id + 1) + bucket_id + 1;
        assert(begin < end);
        return {begin, end};
    }

    uint64_t lookup(uint64_t bucket_id, uint64_t target_kmer, uint64_t k, uint64_t m) const {
        auto [begin, end] = locate_bucket(bucket_id);
        return lookup(begin, end, target_kmer, k, m);
    }

    uint64_t lookup(uint64_t begin, uint64_t end, uint64_t target_kmer, uint64_t k,
                    uint64_t m) const {
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            uint64_t kmer_id = lookup_in_super_kmer(super_kmer_id, target_kmer, k, m);
            if (kmer_id != constants::invalid) return kmer_id;
        }
        return constants::invalid;
    }

    uint64_t lookup_in_super_kmer(uint64_t super_kmer_id, uint64_t target_kmer, uint64_t k,
                                  uint64_t m) const {
        uint64_t offset = offsets.access(super_kmer_id);
        auto [kmer_id, offset_end] = offset_to_id(offset, k);
        bit_vector_iterator bv_it(strings, 2 * offset);
        uint64_t window_size = std::min<uint64_t>(k - m + 1, offset_end - offset - k + 1);
        for (uint64_t w = 0; w != window_size; ++w) {
            uint64_t read_kmer = bv_it.read_and_advance_by_two(2 * k);
            if (read_kmer == target_kmer) return kmer_id + w;
        }
        return constants::invalid;
    }

    /* used by canonical parsing */
    uint64_t lookup_canonical(uint64_t bucket_id, uint64_t target_kmer, uint64_t target_kmer_rc,
                              uint64_t k, uint64_t m) const {
        auto [begin, end] = locate_bucket(bucket_id);
        return lookup_canonical(begin, end, target_kmer, target_kmer_rc, k, m);
    }
    uint64_t lookup_canonical(uint64_t begin, uint64_t end, uint64_t target_kmer,
                              uint64_t target_kmer_rc, uint64_t k, uint64_t m) const {
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            uint64_t offset = offsets.access(super_kmer_id);
            auto [kmer_id, offset_end] = offset_to_id(offset, k);
            bit_vector_iterator bv_it(strings, 2 * offset);
            uint64_t window_size = std::min<uint64_t>(k - m + 1, offset_end - offset - k + 1);
            for (uint64_t w = 0; w != window_size; ++w) {
                uint64_t read_kmer = bv_it.read_and_advance_by_two(2 * k);
                if (read_kmer == target_kmer or read_kmer == target_kmer_rc) return kmer_id + w;
            }
        }
        return constants::invalid;
    }
    /****************************/

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
        uint64_t read_kmer = bv_it.read(2 * k);
        util::uint64_to_string_no_reverse(read_kmer, string_kmer, k);
    }

    struct iterator {
        iterator() {}

        iterator(buckets const* ptr, uint64_t kmer_id, uint64_t k, uint64_t num_kmers)
            : m_buckets(ptr), m_kmer_id(kmer_id), m_k(k), m_num_kmers(num_kmers) {
            bv_it = bit_vector_iterator(m_buckets->strings, -1);
            offset = m_buckets->id_to_offset(m_kmer_id, k);
            auto [pos, piece_end] = m_buckets->pieces.next_geq(offset);
            if (piece_end == offset) pos += 1;
            pieces_it = m_buckets->pieces.at(pos);
            next_piece();
            ret.second.resize(k, 0);
        }

        bool has_next() const { return m_kmer_id != m_num_kmers; }

        std::pair<uint64_t, std::string> next() {
            if (offset == next_offset - m_k + 1) {
                offset = next_offset;
                next_piece();
            }

            while (offset != next_offset - m_k + 1) {
                ret.first = m_kmer_id;
                if (clear) {
                    util::uint64_to_string_no_reverse(read_kmer, ret.second.data(), m_k);
                } else {
                    memmove(ret.second.data(), ret.second.data() + 1, m_k - 1);
                    ret.second[m_k - 1] = util::uint64_to_char(last_two_bits);
                }
                clear = false;
                read_kmer >>= 2;
                last_two_bits = bv_it.get_next_two_bits();
                read_kmer += last_two_bits << (2 * (m_k - 1));
                ++m_kmer_id;
                ++offset;
                return ret;
            }

            return next();
        }

    private:
        std::pair<uint64_t, std::string> ret;
        buckets const* m_buckets;
        uint64_t m_kmer_id, m_k, m_num_kmers;
        uint64_t offset;
        uint64_t next_offset;
        bit_vector_iterator bv_it;
        ef_sequence<true>::iterator pieces_it;

        uint64_t read_kmer;
        uint64_t last_two_bits;
        bool clear;

        void next_piece() {
            bv_it.at(2 * offset);
            next_offset = pieces_it.next();
            assert(next_offset > offset);
            read_kmer = bv_it.take(2 * m_k);
            clear = true;
        }
    };

    iterator at(uint64_t kmer_id, uint64_t k, uint64_t size) const {
        return iterator(this, kmer_id, k, size);
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
};

}  // namespace sshash