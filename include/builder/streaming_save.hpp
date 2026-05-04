#pragma once

#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "essentials.hpp"
#include "external/pthash/external/bits/include/bit_vector.hpp"

#include "include/builder/disk_backed_strings.hpp"

namespace sshash {

/*
    A saver that mirrors `essentials::generic_saver`, except that any visit
    to a specific `bits::bit_vector` instance (identified by address) is
    redirected to `disk_backed_strings::save_to`, which streams the strings
    bytes from the on-disk tmp file. All other visits go through the regular
    `essentials` path.

    Using address-based identification means we don't need to add any
    intermediate type or marker to `bits::bit_vector` itself.
*/
struct streaming_strings_saver {
    streaming_strings_saver(std::ostream& os,                                 //
                            bits::bit_vector const* strings_addr,             //
                            disk_backed_strings const* strings_storage)       //
        : m_os(os), m_strings_addr(strings_addr), m_strings_storage(strings_storage) {
        if (m_strings_addr == nullptr || m_strings_storage == nullptr) {
            throw std::runtime_error("streaming_strings_saver requires non-null arguments");
        }
    }

    template <typename T>
    void visit(T const& val) {
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
};

/*
    Save `t` to `filename`, streaming any embedded `bits::bit_vector` whose
    address matches `strings_addr` from `strings_storage` instead of from
    RAM. Other fields are saved using the standard `essentials` path.
*/
template <typename T>
void save_streaming(T const& t, char const* filename,                    //
                    bits::bit_vector const* strings_addr,                //
                    disk_backed_strings const& strings_storage)          //
{
    std::ofstream out(filename, std::ios::binary);
    if (!out.good()) {
        throw std::runtime_error(std::string("error opening file '") + filename + "' for writing");
    }
    streaming_strings_saver saver(out, strings_addr, &strings_storage);
    saver.visit(t);
    out.close();
}

}  // namespace sshash
