#pragma once

#include "external/pthash/external/bits/include/compact_vector.hpp"
#include "external/pthash/external/bits/include/endpoints_sequence.hpp"

namespace sshash {

struct num_bits {
    num_bits() : per_absolute_offset(0), per_relative_offset(0), per_string_id(0) {}
    uint64_t per_absolute_offset;
    uint64_t per_relative_offset;
    uint64_t per_string_id;
};

template <typename Seq>
struct offsets  //
{
    struct builder {
        builder() {}

        void reserve(uint64_t n) { m_v.reserve(n); }
        void push_back(uint64_t val) { m_v.push_back(val); }

        uint64_t operator[](uint64_t i) {
            assert(i < m_v.size());
            return m_v[i];
        }

        uint64_t front() const { return m_v.front(); }
        uint64_t back() const { return m_v.back(); }
        uint64_t size() const { return m_v.size(); }

        void set_num_bits(num_bits nb) { m_nb = nb; }

    protected:
        num_bits m_nb;
        std::vector<uint64_t> m_v;
    };

    std::pair<uint64_t, uint64_t> id_to_offset(const uint64_t kmer_id, const uint64_t k) const  //
    {
        constexpr uint64_t linear_scan_threshold = 32;
        uint64_t lo = 0;
        uint64_t hi = m_seq.size() - 1;
        assert(m_seq.access(0) == 0);
        while (hi - lo > linear_scan_threshold) {
            uint64_t mid = lo + (hi - lo) / 2;
            uint64_t val = m_seq.access(mid);
            assert(val >= mid * (k - 1));
            if (kmer_id <= val - mid * (k - 1)) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        assert(lo < hi);
        assert(hi < m_seq.size());
        for (auto it = m_seq.get_iterator_at(lo); lo <= hi; ++lo, ++it) {
            uint64_t val = *it - lo * (k - 1);
            if (val > kmer_id) break;
        }
        assert(lo > 0);
        return {lo, kmer_id + (lo - 1) * (k - 1)};
    }

    uint64_t access(uint64_t i) const {
        assert(i < size());
        return m_seq.access(i);
    }

    uint64_t size() const { return m_seq.size(); }

    uint64_t num_bytes() const {
        return sizeof(m_num_bits_per_relative_offset) + m_seq.num_bytes();
    }

    struct iterator {
        iterator() {}
        iterator(offsets const* e, uint64_t pos) { m_it = e->m_seq.get_iterator_at(pos); }

        uint64_t value() const { return *m_it; }
        void next() { ++m_it; }

    private:
        typename Seq::iterator m_it;
    };

    iterator get_iterator_at(uint64_t pos) const {
        assert(pos < size());
        return {this, pos};
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

protected:
    Seq m_seq;
    uint64_t m_num_bits_per_relative_offset;

    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_num_bits_per_relative_offset);
        visitor.visit(t.m_seq);
    }
};

struct decoded_offsets
    : public offsets<bits::endpoints_sequence<>>  //
{
    struct builder
        : public offsets<bits::endpoints_sequence<>>::builder  //
    {
        uint64_t num_bits_per_offset() const { return m_nb.per_absolute_offset; }

        uint64_t encode(uint64_t offset, uint64_t, uint64_t) { return offset; }

        void build(decoded_offsets& e) {
            assert(std::is_sorted(m_v.begin(), m_v.end()));
            e.m_seq.encode(m_v.begin(), m_v.size(), m_v.back());
            std::vector<uint64_t>().swap(m_v);
        }
    };

    struct decoded_offset {
        uint64_t absolute_offset;
    };

    decoded_offset decode(const uint64_t encoded_offset) const { return {encoded_offset}; }

    void offset_to_id(lookup_result& res, decoded_offset /* p */, const uint64_t k) const  //
    {
        assert(res.kmer_offset != constants::invalid_uint64);

        auto q = m_seq.locate(res.kmer_offset);
        res.string_id = q.first.pos;
        res.string_begin = q.first.val;
        res.string_end = q.second.val;
        res.kmer_id = res.kmer_offset - res.string_id * (k - 1);     // absolute kmer id
        res.kmer_id_in_string = res.kmer_offset - res.string_begin;  // relative kmer id

        assert(res.string_id < m_seq.size());
        assert(res.string_begin < res.string_end);
        assert(res.kmer_offset >= res.string_id * (k - 1));
        assert(res.string_begin <= res.kmer_offset);
        assert(res.kmer_offset < res.string_end);
    }
};

struct encoded_offsets
    : public offsets<bits::compact_vector>  //
{
    struct builder
        : public offsets<bits::compact_vector>::builder  //
    {
        uint64_t num_bits_per_offset() const {
            return m_nb.per_string_id + m_nb.per_relative_offset;
        }

        uint64_t encode(uint64_t offset, uint64_t begin, uint64_t string_id) {
            /* encode offset as | string-id | relative offset | */
            assert(string_id < m_v.size());
            assert(offset >= begin);
            assert((offset - begin) < (1ULL << m_nb.per_relative_offset));
            uint64_t relative_offset = offset - begin;
            return (string_id << m_nb.per_relative_offset) + relative_offset;
        }

        void build(encoded_offsets& e) {
            assert(std::is_sorted(m_v.begin(), m_v.end()));
            e.m_seq.build(m_v.begin(), m_v.size(), m_nb.per_absolute_offset);
            e.m_num_bits_per_relative_offset = m_nb.per_relative_offset;
            std::vector<uint64_t>().swap(m_v);
        }
    };

    struct decoded_offset {
        uint64_t absolute_offset;
        uint64_t relative_offset;
        uint64_t string_id;
        uint64_t string_begin;
        uint64_t string_end;
    };

    decoded_offset decode(const uint64_t encoded_offset) const {
        uint64_t relative_offset = encoded_offset & ((1ULL << m_num_bits_per_relative_offset) - 1);
        uint64_t string_id = encoded_offset >> m_num_bits_per_relative_offset;
        assert(string_id + 1 < m_seq.size());
        uint64_t begin = m_seq.access(string_id);
        uint64_t end = m_seq.access(string_id + 1);
        return {begin + relative_offset, relative_offset, string_id, begin, end};
    }

    void offset_to_id(lookup_result& res, decoded_offset p, const uint64_t k) const  //
    {
        assert(res.kmer_offset != constants::invalid_uint64);

        res.string_id = p.string_id;
        res.string_begin = p.string_begin;
        res.string_end = p.string_end;
        res.kmer_id = res.kmer_offset - res.string_id * (k - 1);     // absolute kmer id
        res.kmer_id_in_string = res.kmer_offset - res.string_begin;  // relative kmer id

        assert(res.string_id < m_seq.size());
        assert(res.string_begin < res.string_end);
    }
};

}  // namespace sshash