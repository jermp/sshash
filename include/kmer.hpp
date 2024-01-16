#pragma once

namespace sshash {

template <typename Int, uint8_t BitsPerChar>
struct uint_kmer_t {
    using uint_t = Int;
    Int kmer;

    uint_kmer_t() {}
    uint_kmer_t(Int kmer) : kmer(kmer) {}

    // Only if you're sure the result is uint64_t
    // e.g. with minimizers
    explicit operator uint64_t() const {
        assert(uint64_t(kmer) == kmer);
        return kmer;
    }

    // TODO: change to <=> when switching to C++20
    bool operator==(uint_kmer_t const& t) const { return kmer == t.kmer; }
    bool operator!=(uint_kmer_t const& t) const { return kmer != t.kmer; }
    bool operator<(uint_kmer_t const& t) const { return kmer < t.kmer; }

    void pad(uint16_t b) { kmer <<= b; }
    void pad_char() { kmer <<= bits_per_char; }

    void drop(uint16_t b) { kmer >>= b; }
    void drop64() { drop(64); }
    void drop_char() { drop(bits_per_char); }
    void drop_chars(uint16_t k) { drop(k * bits_per_char); }

    void take(uint16_t b) { kmer &= (Int(1) << b) - 1; }
    void take64() { kmer = uint64_t(kmer); }
    void take_char() { take(bits_per_char); }
    void take_chars(uint16_t k) { take(k * bits_per_char); }

    uint64_t pop64() {
        uint64_t res(kmer);
        drop64();
        return res;
    }
    uint64_t pop_char() {
        uint64_t res(kmer);
        res &= (uint64_t(1) << bits_per_char) - 1;
        drop_char();
        return res;
    }

    void append_bits(uint64_t n, uint16_t b) {
        pad(b);
        kmer |= n;
    }
    void append64(uint64_t n) { append_bits(n, 64); }
    void append_char(uint64_t c) { append_bits(c, bits_per_char); }

    // assigns a character at k-th position
    // assuming that the position is empty
    // use append_char instead, if possible
    void add_kth_char(uint16_t k, uint64_t c) { kmer |= Int(c) << (k * bits_per_char); }

    static constexpr uint16_t uint_kmer_bits = sizeof(Int) * 8;
    static constexpr uint8_t bits_per_char = BitsPerChar;
    // max *odd* size that can be packed into uint_kmer_bits bits
    static constexpr uint16_t max_k = []() {
        uint16_t max_k_any = uint_kmer_bits / bits_per_char;
        return max_k_any % 2 == 0 ? max_k_any - 1 : max_k_any;
    }();

    static_assert(uint_kmer_bits % 64 == 0, "Int must use 64*k bits");
    static_assert(bits_per_char <= 64, "At most 64 bits per character");
};

template <typename Int, uint8_t BitsPerChar, char const* Alphabet>
struct alpha_kmer_t : uint_kmer_t<Int, BitsPerChar> {
    using base = uint_kmer_t<Int, BitsPerChar>;
    using base::base;
    static constexpr char const* alphabet = Alphabet;
    static constexpr uint8_t alphabet_size = strlen(Alphabet);

    static uint64_t char_to_uint(char c);

    static char uint64_to_char(uint64_t x) { return alphabet[x]; }
};

constexpr char dna_alphabet[] = "ACTG";

template <typename Int>
struct dna_uint_kmer_t : alpha_kmer_t<Int, 2, dna_alphabet> {
    using base = alpha_kmer_t<Int, 2, dna_alphabet>;
    using base::uint_kmer_bits;
    using base::base;
    /*
        This works with the map:
        A -> 00; C -> 01; G -> 11; T -> 10.

        Example.
        reverse_complement("ACTCACG") = CGTGAGT, in binary:
        reverse_complement("00.01.10.01.00.01.11") = 01.11.10.11.00.11.10.
    */
    [[maybe_unused]] static uint64_t crc64(uint64_t x) {
        /* Complement, swap byte order */
        uint64_t res = __builtin_bswap64(x ^ 0xaaaaaaaaaaaaaaaa);

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
        res.drop(uint_kmer_bits - k);
        return res;
    }

    /*
    char decimal  binary
    A     65     01000-00-1 -> 00
    C     67     01000-01-1 -> 01
    G     71     01000-11-1 -> 11
    T     84     01010-10-0 -> 10

    a     97     01100-00-1 -> 00
    c     99     01100-01-1 -> 01
    g    103     01100-11-1 -> 11
    t    116     01110-10-0 -> 10
    */
    static uint64_t char_to_uint(char c) { return (c >> 1) & 3; }
};

using default_kmer_t = dna_uint_kmer_t<__uint128_t>;

}  // namespace sshash
