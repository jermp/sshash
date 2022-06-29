#pragma once

#include "../external/pthash/include/encoders/bit_vector.hpp"
#include "../external/pthash/include/encoders/compact_vector.hpp"
#include "../external/pthash/include/encoders/darray.hpp"

namespace sshash {

template <bool index_zeros = false>
struct ef_sequence {
    ef_sequence() : m_universe(0) {}

    template <typename ForwardIterator>
    void encode(ForwardIterator begin, uint64_t n, uint64_t u) {
        if (n == 0) return;
        m_universe = u;

        uint64_t l = uint64_t((n && u / n) ? pthash::util::msb(u / n) : 0);
        pthash::bit_vector_builder bvb_high_bits(n + (u >> l) + 1);
        pthash::compact_vector::builder cv_builder_low_bits(n, l);

        uint64_t low_mask = (uint64_t(1) << l) - 1;
        uint64_t last = 0;
        for (size_t i = 0; i != n; ++i, ++begin) {
            auto v = *begin;
            if (i and v < last) {  // check the order
                std::cerr << "error at " << i << "/" << n << ":\n";
                std::cerr << "last " << last << "\n";
                std::cerr << "current " << v << "\n";
                throw std::runtime_error("ef_sequence is not sorted");
            }
            if (l) cv_builder_low_bits.push_back(v & low_mask);
            bvb_high_bits.set((v >> l) + i, 1);
            last = v;
        }

        pthash::bit_vector(&bvb_high_bits).swap(m_high_bits);
        cv_builder_low_bits.build(m_low_bits);
        m_high_bits_d1.build(m_high_bits);
        if constexpr (index_zeros) m_high_bits_d0.build(m_high_bits);
    }

    struct iterator {
        iterator() : m_ef(nullptr) {}

        iterator(ef_sequence const* ef, uint64_t pos = 0)
            : m_ef(ef), m_pos(pos), m_l(ef->m_low_bits.width()) {
            assert(m_pos < m_ef->size());
            assert(m_l < 64);
            if (m_ef->m_high_bits_d1.num_positions() == 0) return;
            uint64_t begin = m_ef->m_high_bits_d1.select(m_ef->m_high_bits, m_pos);
            m_high_enum = pthash::bit_vector::unary_iterator(m_ef->m_high_bits, begin);
            m_low_enum = m_ef->m_low_bits.at(m_pos);
        }

        bool good() const { return m_ef != nullptr; }
        bool has_next() const { return m_pos < m_ef->size(); }

        uint64_t next() {
            assert(good() and has_next());
            uint64_t high = m_high_enum.next();
            assert(high == m_ef->m_high_bits_d1.select(m_ef->m_high_bits, m_pos));
            uint64_t low = m_low_enum.value();
            uint64_t val = (((high - m_pos) << m_l) | low);
            ++m_pos;
            return val;
        }

    private:
        ef_sequence const* m_ef;
        uint64_t m_pos;
        uint64_t m_l;
        pthash::bit_vector::unary_iterator m_high_enum;
        pthash::compact_vector::iterator m_low_enum;
    };

    iterator at(uint64_t pos) const {
        assert(pos < size());
        return iterator(this, pos);
    }

    inline uint64_t access(uint64_t i) const {
        assert(i < size());
        return ((m_high_bits_d1.select(m_high_bits, i) - i) << m_low_bits.width()) |
               m_low_bits.access(i);
    }

    // inline uint64_t diff(uint64_t i) const {
    //     assert(i < size() && encode_prefix_sum);
    //     uint64_t low1 = m_low_bits.access(i);
    //     uint64_t low2 = m_low_bits.access(i + 1);
    //     uint64_t l = m_low_bits.width();
    //     uint64_t pos = m_high_bits_d1.select(m_high_bits, i);
    //     uint64_t h1 = pos - i;
    //     uint64_t h2 = bit_vector::unary_iterator(m_high_bits, pos + 1).next() - i - 1;
    //     uint64_t val1 = (h1 << l) | low1;
    //     uint64_t val2 = (h2 << l) | low2;
    //     return val2 - val1;
    // }

    // Return [position,value] of the rightmost smallest element >= x.
    // Return [size(),back()] if x > back() (largest element).
    inline std::pair<uint64_t, uint64_t> next_geq(uint64_t x) const {
        static_assert(index_zeros == true, "must build index on zeros");
        assert(m_high_bits_d0.num_positions());

        if (x >= back()) return {size() - (x == back()), back()};

        uint64_t h_x = x >> m_low_bits.width();
        uint64_t begin = h_x ? m_high_bits_d0.select(m_high_bits, h_x - 1) - h_x + 1 : 0;
        assert(begin < size());

        // uint64_t end = m_high_bits_d0.select(m_high_bits, h_x) - h_x;
        // assert(end <= size());
        // assert(begin <= end);
        // return binary search for x in [begin, end)

        auto it = at(begin);
        uint64_t pos = begin;
        uint64_t val = it.next();
        while (val < x) {
            ++pos;
            val = it.next();
        }
        assert(val >= x);

        // now pos is the position of the first element (the leftmost one)
        // that is >= x

        if (val == x) {
            // keep scanning to pick the rightmost one
            while (val == x) {
                ++pos;
                val = it.next();
            }
            assert(val > x);
            return {pos - 1, x};
        }

        return {pos, val};
    }

    // Return [pos_next,prev,next] such that prev < x <= next.
    // where pos_next is the position of next.
    // Assumptions:
    // - all elements of the sequence are distinct;
    // - x <= back();
    // - first element of the sequence is 0.
    // Because of these assumption, this is specific for our use-case, NOT for generic use.
    inline std::tuple<uint64_t, uint64_t, uint64_t> locate(uint64_t x) const {
        static_assert(index_zeros == true, "must build index on zeros");
        assert(m_high_bits_d0.num_positions());

        if (x == 0) return {0, 0, 0};

        uint64_t h_x = x >> m_low_bits.width();
        uint64_t begin = h_x ? m_high_bits_d0.select(m_high_bits, h_x - 1) - h_x + 1 : 0;
        assert(begin < size());

        auto it = at(begin);
        uint64_t pos_next = begin;
        uint64_t next = it.next();
        uint64_t prev = next;
        while (next < x) {
            ++pos_next;
            prev = next;
            next = it.next();
        }
        assert(next >= x);
        assert(pos_next > 0);
        return {pos_next, pos_next != begin ? prev : access(pos_next - 1), next};
    }

    // Return the position of the rightmost largest element <= x.
    // Return size() if x > back() (largest element).
    inline uint64_t prev_leq(uint64_t x) const {
        auto [pos, val] = next_geq(x);
        return pos - (val > x);
    }

    inline uint64_t back() const { return m_universe; }
    inline uint64_t size() const { return m_low_bits.size(); }

    uint64_t num_bits() const {
        return 8 * (sizeof(m_universe) + m_high_bits.bytes() + m_high_bits_d1.bytes() +
                    m_high_bits_d0.bytes() + m_low_bits.bytes());
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_universe);
        visitor.visit(m_high_bits);
        visitor.visit(m_high_bits_d1);
        visitor.visit(m_high_bits_d0);
        visitor.visit(m_low_bits);
    }

private:
    uint64_t m_universe;
    pthash::bit_vector m_high_bits;
    pthash::darray1 m_high_bits_d1;
    pthash::darray0 m_high_bits_d0;
    pthash::compact_vector m_low_bits;
};

}  // namespace sshash