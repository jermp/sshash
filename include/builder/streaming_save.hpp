#pragma once

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "essentials.hpp"
#include "external/pthash/external/bits/include/bit_vector.hpp"
#include "external/pthash/external/bits/include/compact_vector.hpp"

#include "include/builder/disk_backed_strings.hpp"

namespace sshash {

/*
    A saver that mirrors `essentials::generic_saver`, except that any visit
    to a specific `bits::bit_vector` instance (identified by address) is
    redirected to `disk_backed_strings::save_to`, which streams the strings
    bytes from the on-disk tmp file. Likewise, visits to `bits::compact_vector`
    instances whose addresses appear in `compact_vector_subs` are replaced
    with byte-for-byte streaming from the corresponding tmp file (which is
    expected to be in `bits::compact_vector::visit_impl`'s on-disk format).

    Address-based identification means we don't need to add any intermediate
    type or marker to bits::bit_vector / bits::compact_vector themselves.
*/
struct streaming_strings_saver {
    streaming_strings_saver(
        std::ostream& os,                                                                //
        bits::bit_vector const* strings_addr,                                            //
        disk_backed_strings const* strings_storage,                                      //
        std::unordered_map<bits::compact_vector const*, std::string> compact_vector_subs //
        )
        : m_os(os)
        , m_strings_addr(strings_addr)
        , m_strings_storage(strings_storage)
        , m_compact_vector_subs(std::move(compact_vector_subs)) {
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
        if constexpr (std::is_same_v<T, bits::compact_vector>) {
            auto it = m_compact_vector_subs.find(&val);
            if (it != m_compact_vector_subs.end()) {
                stream_file_into_os(it->second);
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
    std::unordered_map<bits::compact_vector const*, std::string> m_compact_vector_subs;

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
    bits::compact_vector whose address appears in `compact_vector_subs`
    has its bytes copied from the corresponding tmp file. Other fields are
    saved via the standard essentials path.
*/
template <typename T>
void save_streaming(T const& t, char const* filename,                                  //
                    bits::bit_vector const* strings_addr,                              //
                    disk_backed_strings const& strings_storage,                        //
                    std::unordered_map<bits::compact_vector const*, std::string>       //
                        compact_vector_subs = {})                                      //
{
    std::ofstream out(filename, std::ios::binary);
    if (!out.good()) {
        throw std::runtime_error(std::string("error opening file '") + filename + "' for writing");
    }
    streaming_strings_saver saver(out, strings_addr, &strings_storage,
                                  std::move(compact_vector_subs));
    saver.visit(t);
    out.close();
}

}  // namespace sshash
