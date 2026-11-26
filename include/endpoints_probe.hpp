#pragma once

#include <type_traits>

namespace sshash {

/*
    Addresses inside the string-endpoints sequence, captured through its
    own public `visit` traversal, so that the dependent reads of
    `offset_to_id` -- the hint, then the entry lines of the search --
    can be prefetched one sweep before they are read. These are pure
    prefetch hints: if the traversal shape ever changes and the capture
    comes back incomplete, the engine runs identically without them.
*/
struct endpoints_probe {
    bool valid = false;
    char const* hints_bytes = nullptr;
    uint64_t hints_width = 0;
    uint64_t hints_size = 0;
    uint64_t const* high_words = nullptr;
    uint8_t const* low_bytes = nullptr;

    template <typename T, typename V, typename = void>
    struct has_visit : std::false_type {};
    template <typename T, typename V>
    struct has_visit<T, V,
                     std::void_t<decltype(std::declval<T const&>().visit(std::declval<V&>()))>>
        : std::true_type {};

    template <typename T>
    void visit(T const& t) {
        if constexpr (std::is_same_v<T, bits::compact_vector>) {
            if (hints_bytes == nullptr) {
                hints_bytes = reinterpret_cast<char const*>(t.data().data());
                hints_width = t.width();
                hints_size = t.size();
            }
        } else if constexpr (std::is_same_v<T, bits::bit_vector>) {
            if (high_words == nullptr) high_words = t.data().data();
        } else if constexpr (std::is_same_v<T, essentials::owning_span<uint8_t>>) {
            if (low_bytes == nullptr) low_bytes = t.data();
        } else if constexpr (std::is_arithmetic_v<T> or std::is_enum_v<T>) {
            /* nothing */
        } else if constexpr (has_visit<T, endpoints_probe>::value) {
            t.visit(*this);
        }
    }

    template <typename Offsets>
    void capture(Offsets const& offsets) {
        offsets.visit(*this);
        valid = hints_bytes != nullptr and high_words != nullptr and
                low_bytes != nullptr and hints_width > 0 and hints_width < 57;
    }
};

}  // namespace sshash
