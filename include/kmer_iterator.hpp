#pragma once

#include "util.hpp"

namespace sshash {

template <typename Kmer, typename BitVector>
struct kmer_iterator  //
{
    kmer_iterator() {}

    kmer_iterator(BitVector const& bv, const uint64_t k)
        : m_bv(&bv), m_uint_kmer_bits(Kmer::bits_per_char * k), m_pos(0), m_avail(0), m_buff(0) {}

    kmer_iterator(BitVector const& bv, const uint64_t k, const uint64_t pos)
        : kmer_iterator(bv, k)  //
    {
        at(pos);
    }

    void at(uint64_t pos) {
        m_pos = pos;
        m_avail = 0;
        m_buff = 0;
    }

    Kmer get() {
        if (m_avail < m_uint_kmer_bits) fill_buff();
        auto kmer = m_buff;
        kmer.take(m_uint_kmer_bits);
        return kmer;
    }

    Kmer get_reverse() {
        if (m_avail < m_uint_kmer_bits) fill_buff_reverse();
        auto kmer = m_buff;
        kmer.drop(Kmer::uint_kmer_bits - m_uint_kmer_bits);
        return kmer;
    }

    void next() {
        if (m_avail < Kmer::bits_per_char) fill_buff();
        m_buff.drop_char();
        m_avail -= Kmer::bits_per_char;
        m_pos += Kmer::bits_per_char;
    }

    void next_reverse() {
        if (m_avail < Kmer::bits_per_char) fill_buff_reverse();
        m_buff.pad_char();
        m_avail -= Kmer::bits_per_char;
        m_pos -= Kmer::bits_per_char;
    }

    inline uint64_t get_next_char() {
        if (m_avail < Kmer::bits_per_char) fill_buff();
        m_avail -= Kmer::bits_per_char;
        m_pos += Kmer::bits_per_char;
        return m_buff.pop_char();
    }

    inline uint64_t position() const { return m_pos; }

private:
    inline void fill_buff() {
        static_assert(Kmer::uint_kmer_bits % 64 == 0);
        for (int i = Kmer::uint_kmer_bits - 64; i >= 0; i -= 64) {
            m_buff.append64(m_bv->get_word64(m_pos + i));
        }
        m_avail = Kmer::uint_kmer_bits;
    }

    inline void fill_buff_reverse() {
        static_assert(Kmer::uint_kmer_bits % 64 == 0);
        for (int i = Kmer::uint_kmer_bits; i > 0; i -= 64) {
            m_buff.append64(m_bv->get_word64(std::max<uint64_t>(m_pos, Kmer::uint_kmer_bits) - i));
        }
        m_avail = std::min<uint64_t>(m_pos, Kmer::uint_kmer_bits);
        m_buff.pad(Kmer::uint_kmer_bits - m_avail);
    }

    BitVector const* m_bv;
    uint64_t m_uint_kmer_bits;
    uint64_t m_pos, m_avail;
    Kmer m_buff;
};

}  // namespace sshash