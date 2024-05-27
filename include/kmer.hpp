#pragma once

#include "bitpack.hpp"
#include <bitset>
#include <string>

template <size_t N>
bool operator<(std::bitset<N> const& a, std::bitset<N> const& b) {
    return a.to_string() < b.to_string();
}

namespace sshash {
template <typename Kmer, uint8_t BitsPerChar>
struct uint_kmer_t {
    using uint_t = Kmer;
    Kmer kmer;

    uint_kmer_t() {}
    uint_kmer_t(uint64_t kmer) : kmer(kmer) {}

    explicit operator uint64_t() const {
        if constexpr (std::is_constructible_v<uint64_t, Kmer>) {
            return static_cast<uint64_t>(kmer);
        } else {  // std::bitset?
            return (kmer & Kmer(uint64_t(-1))).to_ulong();
        }
    }

    // TODO: change to <=> when switching to C++20
    bool operator==(uint_kmer_t const& t) const { return kmer == t.kmer; }
    bool operator!=(uint_kmer_t const& t) const { return kmer != t.kmer; }
    bool operator<(uint_kmer_t const& t) const { return kmer < t.kmer; }

    void pad(uint16_t b) {
        assert(b < uint_kmer_bits);
        kmer <<= b;
    }
    void pad_char() { pad(bits_per_char); }

    void drop(uint16_t b) {
        assert(b < uint_kmer_bits);
        kmer >>= b;
    }
    void drop64() {
        if constexpr (uint_kmer_bits == 64) {
            kmer = 0;
        } else {
            drop(64);
        }
    }
    void drop_char() { drop(bits_per_char); }
    void drop_chars(uint16_t k) { drop(k * bits_per_char); }

    void take(uint16_t b) { kmer &= ~(~Kmer(0) << b); }
    void take64() { kmer &= Kmer(uint64_t(-1)); }
    void take_char() { take(bits_per_char); }
    void take_chars(uint16_t k) { take(k * bits_per_char); }

    uint64_t pop64() {
        uint64_t res(*this);
        drop64();
        return res;
    }
    uint64_t pop_char() {
        uint64_t res(*this);
        res &= (uint64_t(1) << bits_per_char) - 1;
        drop_char();
        return res;
    }

    void append(uint16_t b, uint64_t n) {
        assert(b < uint_kmer_bits);
        kmer = (kmer << b) | Kmer(n);
    }
    void append64(uint64_t n) {
        if constexpr (uint_kmer_bits == 64) {
            kmer = n;
        } else {
            append(64, n);
        }
    }
    void append_char(uint64_t c) { append(bits_per_char, c); }

    // assigns a character at k-th position
    // assuming that the position is empty
    void kth_char_or(uint16_t k, uint64_t c) { kmer |= Kmer(c) << (k * bits_per_char); }

    static constexpr uint16_t uint_kmer_bits = 8 * sizeof(Kmer);
    static constexpr uint8_t bits_per_char = BitsPerChar;
    // max *odd* size that can be packed into uint_kmer_bits bits
    static constexpr uint16_t max_k = []() {
        uint16_t max_k_any = uint_kmer_bits / bits_per_char;
        return max_k_any % 2 == 0 ? max_k_any - 1 : max_k_any;
    }();

    static_assert(uint_kmer_bits % 64 == 0, "Kmer must use 64*k bits");
    static_assert(bits_per_char < 64, "Less than 64 bits per character");
};

template <typename Kmer, uint8_t BitsPerChar, char const* Alphabet>
struct alpha_kmer_t : uint_kmer_t<Kmer, BitsPerChar> {
    using base = uint_kmer_t<Kmer, BitsPerChar>;
    using base::base;
    static constexpr char const* alphabet = Alphabet;
    static constexpr uint8_t alphabet_size = std::char_traits<char>::length(Alphabet);

    static uint64_t char_to_uint(char c);
    static char uint64_to_char(uint64_t x) { return alphabet[x]; }
};

#ifdef SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING
inline constexpr char nucleotides[] = "ACGT";
#else
inline constexpr char nucleotides[] = "ACTG";
#endif

template <typename Kmer>
struct dna_uint_kmer_t : alpha_kmer_t<Kmer, 2, nucleotides> {
    using base = alpha_kmer_t<Kmer, 2, nucleotides>;
    using base::uint_kmer_bits;
    using base::bits_per_char;
    using base::max_k;
    using base::base;
    /*
        This works with the map:
        A -> 00; C -> 01; G -> 11; T -> 10.

        Example.
        reverse_complement("ACTCACG") = CGTGAGT, in binary:
        reverse_complement("00.01.10.01.00.01.11") = 01.11.10.11.00.11.10.
    */
    [[maybe_unused]] static uint64_t crc64(uint64_t x) {
        /* complement */
#ifdef SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING
        uint64_t c = ~x;
#else
        uint64_t c = x ^ 0xaaaaaaaaaaaaaaaa;  // ...1010.1010.1010.1010
#endif
        /* swap byte order */
        uint64_t res = __builtin_bswap64(c);

        /* Swap nuc order in bytes */
        const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
        const uint64_t c2 = 0x3333333333333333;
        res = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4);  // swap 2-nuc order in bytes
        res = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2);  // swap nuc order in 2-nuc
        return res;
    }

    [[maybe_unused]] dna_uint_kmer_t reverse_complement(uint64_t k) {
        dna_uint_kmer_t x(*this);
        assert(k <= max_k);
        dna_uint_kmer_t res(0);
        for (uint16_t i = 0; i < uint_kmer_bits; i += 64) { res.append64(crc64(x.pop64())); }
        // res is full reverse-complement to x
        res.drop(uint_kmer_bits - k * bits_per_char);
        return res;
    }

#ifdef SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING
    /*
    char decimal  binary
    A     65     01000001 -> 00
    C     67     01000011 -> 01
    G     71     01000111 -> 10
    T     84     01010100 -> 11

    a     97     01100001 -> 00
    c     99     01100011 -> 01
    g    103     01100111 -> 10
    t    116     01110100 -> 11
    */
    static uint64_t char_to_uint(char c) { return (((c >> 1) ^ (c >> 2)) & 3); }
#else
    /*
    char decimal  binary
    A     65     01000.00.1 -> 00
    C     67     01000.01.1 -> 01
    G     71     01000.11.1 -> 11
    T     84     01010.10.0 -> 10

    a     97     01100.00.1 -> 00
    c     99     01100.01.1 -> 01
    g    103     01100.11.1 -> 11
    t    116     01110.10.0 -> 10
    */
    static uint64_t char_to_uint(char c) { return (c >> 1) & 3; }
#endif
};

// also supports bitpack<__uint128_t, 1>, std::bitset<256>, etc
#ifdef SSHASH_USE_MAX_KMER_LENGTH_63
using default_kmer_t = dna_uint_kmer_t<__uint128_t>;
#else
using default_kmer_t = dna_uint_kmer_t<uint64_t>;
#endif

}  // namespace sshash
