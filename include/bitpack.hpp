#pragma once

namespace sshash {
// full binary tree of given height
// with Int type in its leafs
template <typename Int, uint16_t height>
struct bitpack {
    static_assert(height > 0);
    using halfpack = std::conditional_t<height == 1, Int, bitpack<Int, height - 1>>;
    static constexpr uint16_t hsize = 8 * sizeof(halfpack);
    halfpack a, b;

    bitpack() {}
    bitpack(uint64_t x) : a(x), b(0) {}
    bitpack(halfpack a, halfpack b) : a(a), b(b) {}
    explicit operator uint64_t() const { return (uint64_t)a; }

    bool operator==(bitpack const& t) const { return std::pair{a, b} == std::pair{t.a, t.b}; }
    bool operator!=(bitpack const& t) const { return std::pair{a, b} != std::pair{t.a, t.b}; }
    bool operator<(bitpack const& t) const { return std::pair{a, b} < std::pair{t.a, t.b}; }

    // shift in [0, size)
    bitpack& operator>>=(uint16_t shift) {
        if (shift < hsize) {
            a = (a >> shift) | (b << (hsize - shift));
            b >>= shift;
        } else {
            a = b >> (shift - hsize);
            b = 0;
        }
        return *this;
    }
    bitpack& operator<<=(uint16_t shift) {
        if (shift < hsize) {
            b = (b << shift) | (a >> (hsize - shift));
            a <<= shift;
        } else {
            b = a << (shift - hsize);
            a = 0;
        }
        return *this;
    }
    bitpack operator<<(uint16_t shift) const { return bitpack(*this) <<= shift; }
    bitpack operator>>(uint16_t shift) const { return bitpack(*this) >>= shift; }

    bitpack& operator|=(bitpack const& t) {
        a |= t.a;
        b |= t.b;
        return *this;
    }
    bitpack& operator&=(bitpack const& t) {
        a &= t.a;
        b &= t.b;
        return *this;
    }
    bitpack operator|(bitpack const& t) const { return bitpack(*this) |= t; }
    bitpack operator&(bitpack const& t) const { return bitpack(*this) &= t; }

    bitpack operator~() const { return {~a, ~b}; }
};
}  // namespace sshash
