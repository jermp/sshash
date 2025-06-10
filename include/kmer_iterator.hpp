#pragma once

#include "external/pthash/external/bits/include/bit_vector.hpp"
#include "util.hpp"

namespace sshash {

template <class kmer_t>
struct kmer_iterator  //
{
    kmer_iterator() {}

    kmer_iterator(bits::bit_vector const& bv, const uint64_t k)
        : m_bv(&bv)
        , m_k(k)
        , m_uint_kmer_bits(kmer_t::bits_per_char * m_k)
        , m_pos(0)
        , m_avail(0)
        , m_kmer(0)
        , m_buff(0) {}

    kmer_iterator(bits::bit_vector const& bv, const uint64_t k,  //
                  const uint64_t pos, bool reverse = false)
        : kmer_iterator(bv, k)  //
    {
        at(pos, reverse);
    }

    void at(uint64_t pos, bool reverse = false) {
        assert(m_k > 0);
        m_pos = pos;
        m_avail = 0;
        m_buff = 0;

        if (!reverse) {
            fill_buff();
            m_kmer = m_buff;
            m_kmer.take(m_uint_kmer_bits);
            m_buff.drop(m_uint_kmer_bits);
            m_pos += m_uint_kmer_bits;
        } else {
            fill_buff_reverse();
            m_kmer = m_buff;
            m_kmer.drop(kmer_t::uint_kmer_bits - m_uint_kmer_bits);
            m_buff.pad(m_uint_kmer_bits);
            m_pos -= m_uint_kmer_bits;
        }

        m_avail -= m_uint_kmer_bits;
    }

    inline kmer_t get() const { return m_kmer; }

    void next() {
        uint64_t c = get_next_char();
        m_kmer.drop_char();
        m_kmer.set(m_k - 1, c);
    }

    void next_reverse() {
        if (m_avail == 0) fill_buff_reverse();
        uint64_t c = m_buff.at(kmer_t::uint_kmer_bits / kmer_t::bits_per_char - 1);
        m_buff.pad_char();
        m_avail -= kmer_t::bits_per_char;
        m_pos -= kmer_t::bits_per_char;
        m_kmer.pad_char();
        m_kmer.set(0, c);
        m_kmer.take(m_uint_kmer_bits);
    }

    inline uint64_t get_next_char() {
        if (m_avail == 0) fill_buff();
        m_avail -= kmer_t::bits_per_char;
        m_pos += kmer_t::bits_per_char;
        return m_buff.pop_char();
    }

    inline uint64_t position() const { return m_pos; }

private:
    inline void fill_buff() {
        static_assert(kmer_t::uint_kmer_bits % 64 == 0);
        for (int i = kmer_t::uint_kmer_bits - 64; i >= 0; i -= 64) {
            m_buff.append64(m_bv->get_word64(m_pos + i));
        }
        m_avail = kmer_t::uint_kmer_bits;
    }

    inline void fill_buff_reverse() {
        static_assert(kmer_t::uint_kmer_bits % 64 == 0);
        for (int i = kmer_t::uint_kmer_bits; i > 0; i -= 64) {
            m_buff.append64(
                m_bv->get_word64(std::max<uint64_t>(m_pos, kmer_t::uint_kmer_bits) - i));
        }
        m_avail = std::min<uint64_t>(m_pos, kmer_t::uint_kmer_bits);
        m_buff.pad(kmer_t::uint_kmer_bits - m_avail);
    }

    bits::bit_vector const* m_bv;
    uint64_t m_k, m_uint_kmer_bits;
    uint64_t m_pos, m_avail;
    kmer_t m_kmer, m_buff;
};

}  // namespace sshash