#pragma once

#include "external/pthash/include/encoders/bit_vector.hpp"
#include "util.hpp"

namespace sshash {

template <class kmer_t>
struct bit_vector_iterator {
    bit_vector_iterator() : m_bv(nullptr) {}

    bit_vector_iterator(pthash::bit_vector const& bv, uint64_t pos) : m_bv(&bv) { at(pos); }

    void at(uint64_t pos) {
        m_pos = pos;
        m_avail = 0;
        m_buf = 0;
    }

    inline kmer_t read(uint64_t l) {
        assert(l <= kmer_t::uint_kmer_bits);
        if (m_avail < l) fill_buf();
        kmer_t val = m_buf;
        val.take(l);
        return val;
    }

    inline kmer_t read_reverse(uint64_t l) {
        assert(l <= kmer_t::uint_kmer_bits);
        if (m_avail < l) fill_buf_reverse();
        kmer_t val = m_buf;
        val.drop(kmer_t::uint_kmer_bits - l);
        return val;
    }

    inline void eat(uint64_t l) {
        assert(l <= kmer_t::uint_kmer_bits);
        if (m_avail < l) fill_buf();
        m_buf.drop(l);
        m_avail -= l;
        m_pos += l;
    }

    inline void eat_reverse(uint64_t l) {
        assert(l <= kmer_t::uint_kmer_bits);
        if (m_avail < l) fill_buf_reverse();
        m_buf.pad(l);
        m_avail -= l;
        m_pos -= l;
    }

    inline kmer_t read_and_advance_by_char(uint64_t l) {
        assert(l <= kmer_t::uint_kmer_bits);
        if (m_avail < l) fill_buf();
        kmer_t val = m_buf;
        val.take(l);
        m_buf.drop_char();
        m_avail -= kmer_t::bits_per_char;
        m_pos += kmer_t::bits_per_char;
        return val;
    }

    inline uint64_t get_next_char() {
        if (m_avail < kmer_t::bits_per_char) fill_buf();
        m_avail -= kmer_t::bits_per_char;
        m_pos += kmer_t::bits_per_char;
        return m_buf.pop_char();
    }

    inline kmer_t take(uint64_t l) {
        assert(l <= kmer_t::uint_kmer_bits);
        if (m_avail < l) fill_buf();
        kmer_t val = m_buf;
        val.take(l);
        m_buf.drop(l);
        m_avail -= l;
        m_pos += l;
        return val;
    }

    inline uint64_t position() const { return m_pos; }

private:
    inline void fill_buf() {
        static_assert(kmer_t::uint_kmer_bits % 64 == 0);
        for (int i = kmer_t::uint_kmer_bits - 64; i >= 0; i -= 64) {
            if (m_pos + i < m_bv->size()) { m_buf.append64(m_bv->get_word64(m_pos + i)); }
        }
        m_avail = kmer_t::uint_kmer_bits;
    }

    inline void fill_buf_reverse() {
        static_assert(kmer_t::uint_kmer_bits % 64 == 0);
        for (int i = kmer_t::uint_kmer_bits; i > 0; i -= 64) {
            m_buf.append64(m_bv->get_word64(std::max<uint64_t>(m_pos, kmer_t::uint_kmer_bits) - i));
        }
        m_avail = std::min<uint64_t>(m_pos, kmer_t::uint_kmer_bits);
        m_buf.pad(kmer_t::uint_kmer_bits - m_avail);
    }

    pthash::bit_vector const* m_bv;
    uint64_t m_pos;
    uint64_t m_avail;
    kmer_t m_buf;
};

}  // namespace sshash