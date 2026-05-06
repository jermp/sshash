#pragma once

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <vector>

#include "essentials.hpp"
#include "external/pthash/external/bits/include/bit_vector.hpp"

#include "include/builder/disk_backed_strings.hpp"

namespace sshash {

/*
    A typed substitution: the saver replaces the visit byte content of an
    object at `address` with the bytes of `filename` only if the visited
    type T satisfies `std::type_index(typeid(T)) == type`.

    Type discrimination is necessary because in C++ a struct's address
    coincides with the address of its first member when the struct has
    standard layout. Without the type check, a substitution registered
    for an inner field would also fire (incorrectly) on every enclosing
    parent that shares its address.
*/
struct typed_address_sub {
    std::string filename;
    std::type_index type;
};

/*
    A saver that mirrors `essentials::generic_saver`, except for two
    interception mechanisms used during streaming save:

    1. The `bits::bit_vector` instance whose address matches `strings_addr`
       has its bytes streamed from `strings_storage` (which writes the
       same on-disk format `bits::bit_vector::visit_impl` produces).

    2. Any object whose address appears in `address_subs` has its visit
       byte content replaced by a copy of the corresponding tmp file.
       This is type-agnostic — it works for `bits::compact_vector`, for
       pthash MPHFs, or anything else whose serialized form has been
       saved to a file via `essentials::save`.

    The substitution check is performed at the start of every visit<T>
    call (whatever T is); if no match, the call falls through to the
    regular `essentials::generic_saver` logic (POD via save_pod, or
    recursion via val.visit(*this)).
*/
struct streaming_strings_saver {
    streaming_strings_saver(std::ostream& os,                            //
                            bits::bit_vector const* strings_addr,        //
                            disk_backed_strings const* strings_storage,  //
                            std::unordered_map<void const*, typed_address_sub> address_subs)
        : m_os(os)
        , m_strings_addr(strings_addr)
        , m_strings_storage(strings_storage)
        , m_address_subs(std::move(address_subs)) {
        if (m_strings_addr == nullptr || m_strings_storage == nullptr) {
            throw std::runtime_error("streaming_strings_saver requires non-null arguments");
        }
    }

    template <typename T>
    void visit(T const& val) {
        /* Type+address substitution (compact_vectors, MPHFs, etc.).
           Both must match: address alone is ambiguous when a struct
           shares its address with its first member. */
        void const* addr = static_cast<void const*>(&val);
        auto it = m_address_subs.find(addr);
        if (it != m_address_subs.end() && it->second.type == std::type_index(typeid(T))) {
            stream_file_into_os(it->second.filename);
            return;
        }
        /* Strings: dedicated callback because the on-disk strings file
           holds raw words (not the bits::bit_vector serialized form);
           `disk_backed_strings::save_to(os)` writes the visit_impl format
           on the fly. */
        if constexpr (std::is_same_v<T, bits::bit_vector>) {
            if (&val == m_strings_addr) {
                m_strings_storage->save_to(m_os);
                return;
            }
        }
        if constexpr (essentials::is_pod<T>::value) {
            essentials::save_pod(m_os, val);
        } else {
            val.visit(*this);
        }
    }

    template <typename T, typename Allocator>
    void visit(std::vector<T, Allocator> const& vec) {
        visit_seq(vec);
    }

    template <typename T>
    void visit(essentials::owning_span<T> const& vec) {
        visit_seq(vec);
    }

    std::size_t bytes() { return static_cast<std::size_t>(m_os.tellp()); }

private:
    std::ostream& m_os;
    bits::bit_vector const* m_strings_addr;
    disk_backed_strings const* m_strings_storage;
    std::unordered_map<void const*, typed_address_sub> m_address_subs;

    template <typename Vec>
    void visit_seq(Vec const& vec) {
        using T = typename Vec::value_type;
        const std::size_t n = vec.size();
        visit(n);
        if constexpr (essentials::is_pod<T>::value) {
            m_os.write(reinterpret_cast<char const*>(vec.data()),
                       static_cast<std::streamsize>(sizeof(T) * n));
        } else {
            for (auto const& v : vec) visit(v);
        }
    }

    void stream_file_into_os(std::string const& filename) {
        std::ifstream in(filename, std::ifstream::binary);
        if (!in.is_open()) {
            throw std::runtime_error("cannot open spilled component file '" + filename + "'");
        }
        char buf[64 * 1024];
        while (in.good()) {
            in.read(buf, sizeof(buf));
            const std::streamsize got = in.gcount();
            if (got > 0) m_os.write(buf, got);
        }
        in.close();
    }
};

/*
    Save `t` to `filename`. Any embedded bits::bit_vector matching
    `strings_addr` is streamed from `strings_storage`; any embedded
    object whose address appears in `address_subs` has its bytes copied
    from the corresponding tmp file. Other fields are saved via the
    standard essentials path.
*/
template <typename T>
void save_streaming(T const& t, char const* filename,            //
                    bits::bit_vector const* strings_addr,        //
                    disk_backed_strings const& strings_storage,  //
                    std::unordered_map<void const*, typed_address_sub> address_subs = {}) {
    std::ofstream out(filename, std::ios::binary);
    if (!out.good()) {
        throw std::runtime_error(std::string("error opening file '") + filename + "' for writing");
    }
    streaming_strings_saver saver(out, strings_addr, &strings_storage, std::move(address_subs));
    saver.visit(t);
    out.close();
}

/* Helper: register a typed substitution at the address of `addr`. */
template <typename T>
inline void register_sub(std::unordered_map<void const*, typed_address_sub>& subs, T const* addr,
                         std::string filename) {
    subs.insert_or_assign(static_cast<void const*>(addr),
                          typed_address_sub{std::move(filename), std::type_index(typeid(T))});
}

}  // namespace sshash
