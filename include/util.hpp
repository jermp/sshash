#pragma once

#include <array>
#include <cassert>
#include <fstream>
#include <thread>
#include <cmath>  // for std::ceil on linux

#include "hash_util.hpp"

namespace sshash {

enum bucket_t : int {
    SINGLETON = 0,  // minimizer appears only once
    MIDLOAD = 1,    // minimizer appears > 1 but < 2^l times
    HEAVYLOAD = 3   // minimizer appears >= 2^l times
};

enum input_file_t { fasta, cf_seg };

struct streaming_query_report {
    streaming_query_report()
        : num_kmers(0)
        , num_positive_kmers(0)
        , num_negative_kmers(0)
        , num_invalid_kmers(0)
        , num_searches(0)
        , num_extensions(0) {}

    uint64_t num_kmers;
    uint64_t num_positive_kmers;
    uint64_t num_negative_kmers;
    uint64_t num_invalid_kmers;
    uint64_t num_searches;
    uint64_t num_extensions;
};

struct lookup_result {
    lookup_result(bool mf = true)
        : kmer_id(constants::invalid_uint64)
        , kmer_id_in_string(constants::invalid_uint64)
        , kmer_offset(constants::invalid_uint64)
        , kmer_orientation(constants::forward_orientation)

        , string_id(constants::invalid_uint64)
        , string_begin(constants::invalid_uint64)
        , string_end(constants::invalid_uint64)

        , minimizer_found(mf) {}

    uint64_t kmer_id;            // "absolute" kmer-id
    uint64_t kmer_id_in_string;  // "relative" kmer-id: 0 <= kmer_id_in_string < string_size,
                                 // where string_size = string_end - string_begin - k + 1
    uint64_t kmer_offset;
    int64_t kmer_orientation;

    uint64_t string_id;
    uint64_t string_begin;
    uint64_t string_end;

    bool minimizer_found;
};

inline std::ostream& operator<<(std::ostream& os, lookup_result const& res) {
    os << "  == kmer_id = " << res.kmer_id << '\n';
    os << "  == kmer_id_in_string = " << res.kmer_id_in_string << '\n';
    os << "  == kmer_offset = " << res.kmer_offset << '\n';
    os << "  == kmer_orientation = " << res.kmer_orientation << '\n';
    os << "  == string_id = " << res.string_id << '\n';
    os << "  == string_begin = " << res.string_begin << '\n';
    os << "  == string_end = " << res.string_end << '\n';
    os << "  == string_length = " << (res.string_end - res.string_begin) << '\n';
    os << "  == minimizer_found = " << (res.minimizer_found ? "true" : "false") << '\n';
    return os;
}

template <class kmer_t>
struct neighbourhood {
    std::array<lookup_result, kmer_t::alphabet_size> forward;
    std::array<lookup_result, kmer_t::alphabet_size> backward;
};

struct minimizer_info {
    minimizer_info()
        : minimizer(constants::invalid_uint64)
        , pos_in_seq(constants::invalid_uint64)
        , pos_in_kmer(constants::invalid_uint64) {}

    minimizer_info(uint64_t mm, uint64_t pk)
        : minimizer(mm), pos_in_seq(constants::invalid_uint64), pos_in_kmer(pk) {}

    minimizer_info(uint64_t mm, uint64_t ps, uint64_t pk)
        : minimizer(mm), pos_in_seq(ps), pos_in_kmer(pk) {}

    uint64_t minimizer;
    uint64_t pos_in_seq;
    uint64_t pos_in_kmer;

    bool operator==(minimizer_info rhs) const {
        return minimizer == rhs.minimizer and    //
               pos_in_seq == rhs.pos_in_seq and  //
               pos_in_kmer == rhs.pos_in_kmer;   //
    }
    bool operator!=(minimizer_info rhs) const { return !(*this == rhs); }
};

[[maybe_unused]] static bool equal_lookup_result(lookup_result expected, lookup_result got) {
    bool good = true;
    if (expected.kmer_id != got.kmer_id) {
        std::cout << "expected kmer_id " << expected.kmer_id << " but got " << got.kmer_id
                  << std::endl;
        good = false;
    }
    if (expected.kmer_id_in_string != got.kmer_id_in_string) {
        std::cout << "expected kmer_id_in_string " << expected.kmer_id_in_string << " but got "
                  << got.kmer_id_in_string << std::endl;
        good = false;
    }
    if (got.kmer_id != constants::invalid_uint64 and
        expected.kmer_orientation != got.kmer_orientation) {
        std::cout << "expected kmer_orientation " << expected.kmer_orientation << " but got "
                  << got.kmer_orientation << std::endl;
        good = false;
    }
    if (expected.string_id != got.string_id) {
        std::cout << "expected string_id " << expected.string_id << " but got " << got.string_id
                  << std::endl;
        good = false;
    }
    if (expected.string_begin != got.string_begin) {
        std::cout << "expected string_begin " << expected.string_begin << " but got "
                  << got.string_begin << std::endl;
        good = false;
    }
    if (expected.string_end != got.string_end) {
        std::cout << "expected string_end " << expected.string_end << " but got " << got.string_end
                  << std::endl;
        good = false;
    }
    return good;
}

struct build_configuration {
    build_configuration()
        : k(31)
        , m(20)
        , seed(constants::seed)
        , num_threads(1)
        , ram_limit_in_GiB(constants::default_ram_limit_in_GiB)

        , lambda(constants::lambda)

        , canonical(false)
        , weighted(false)
        , verbose(true)

        , tmp_dirname(constants::default_tmp_dirname)

    {}

    uint64_t k;  // kmer length
    uint64_t m;  // minimizer length
    uint64_t seed;
    uint64_t num_threads;
    uint64_t ram_limit_in_GiB;

    double lambda;  // drive PTHash trade-off

    bool canonical;
    bool weighted;
    bool verbose;

    std::string tmp_dirname;

    void print() const {
        std::cout << "k = " << k                                              //
                  << ", m = " << m                                            //
                  << ", seed = " << seed                                      //
                  << ", num_threads = " << num_threads                        //
                  << ", ram_limit_in_GiB = " << ram_limit_in_GiB              //
                  << ", lambda = " << lambda                                  //
                  << ", canonical = " << (canonical ? "true" : "false")       //
                  << ", weighted = " << (weighted ? "true" : "false")         //
                  << ", verbose = " << (verbose ? "true" : "false")           //
                  << ", tmp_dirname = '" << tmp_dirname << "'" << std::endl;  //
    }
};

namespace util {

static inline void check_version_number(essentials::version_number const& vnum) {
    if (vnum.x != constants::current_version_number::x) {
        throw std::runtime_error("MAJOR index version mismatch: SSHash index needs rebuilding");
    }
}

static inline uint64_t get_seed_for_hash_function(build_configuration const& build_config) {
    static const uint64_t my_favourite_seed = 1234567890;
    return build_config.seed != my_favourite_seed ? my_favourite_seed : ~my_favourite_seed;
}

[[maybe_unused]] static bool ends_with(std::string const& str, std::string const& pattern) {
    if (pattern.size() > str.size()) return false;
    return std::equal(pattern.begin(), pattern.end(), str.end() - pattern.size());
}

template <class kmer_t>
[[maybe_unused]] static kmer_t string_to_uint_kmer(char const* str, uint64_t k) {
    assert(k <= kmer_t::max_k);
    kmer_t x = 0;
    for (uint64_t i = 0; i != k; ++i) x.set(i, kmer_t::char_to_uint(str[i]));
    return x;
}

template <class kmer_t>
static void uint_kmer_to_string(kmer_t x, char* str, uint64_t k) {
    assert(k <= kmer_t::max_k);
    for (uint64_t i = 0; i != k; ++i) str[i] = kmer_t::uint64_to_char(x.pop_char());
}

template <class kmer_t>
[[maybe_unused]] static std::string uint_kmer_to_string(kmer_t x, uint64_t k) {
    assert(k <= kmer_t::max_k);
    std::string str;
    str.resize(k);
    uint_kmer_to_string(x, str.data(), k);
    return str;
}

template <class kmer_t>
[[maybe_unused]] static std::string uint_minimizer_to_string(uint64_t minimizer, uint64_t m) {
    assert(m <= kmer_t::max_m);
    std::string str;
    str.resize(m);
    kmer_t x = minimizer;
    uint_kmer_to_string(x, str.data(), m);
    return str;
}

template <class kmer_t>
[[maybe_unused]] static bool is_valid(char const* str, uint64_t size) {
    for (uint64_t i = 0; i != size; ++i) {
        if (!kmer_t::is_valid(str[i])) return false;
    }
    return true;
}

template <class kmer_t>
static kmer_t read_kmer_at(bits::bit_vector const& bv, const uint64_t k, const uint64_t pos) {
    static_assert(kmer_t::uint_kmer_bits % 64 == 0);
    kmer_t kmer = 0;
    for (int i = kmer_t::uint_kmer_bits - 64; i >= 0; i -= 64) {
        if (pos + i < bv.num_bits()) kmer.append64(bv.get_word64(pos + i));
    }
    kmer.take(kmer_t::bits_per_char * k);
    return kmer;
}

/*
    This implements the random minimizer.
*/
template <class kmer_t>
minimizer_info compute_minimizer(kmer_t kmer, const uint64_t k, const uint64_t m,
                                 hasher_type const& hasher)  //
{
    assert(m <= kmer_t::max_m);
    assert(m <= k);
    uint64_t min_hash = constants::invalid_uint64;
    kmer_t minimizer = kmer_t(-1);
    uint64_t pos = 0;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        kmer_t mmer = kmer;
        mmer.take_chars(m);
        uint64_t hash = hasher.hash(uint64_t(mmer));
        if (hash < min_hash) {
            min_hash = hash;
            minimizer = mmer;
            pos = i;
        }
        kmer.drop_char();
    }
    return {uint64_t(minimizer), pos};
}

}  // namespace util

struct buffered_lines_iterator {
    static const uint64_t BUFFER_SIZE = 1024;

    buffered_lines_iterator(std::istream& is, uint64_t buffer_size = BUFFER_SIZE)
        : m_is(is), m_buffer_size(buffer_size), m_read_chars(0) {}

    bool fill_buffer(std::string& buffer) {
        if (buffer.size() >= m_buffer_size) return false;

        bool empty_line_was_read = false;
        uint64_t size = buffer.size();
        buffer.resize(m_buffer_size);

        char* ptr = buffer.data() + size;
        while (size < m_buffer_size) {
            // read until '\n' or rest is filled
            uint64_t rest = m_buffer_size - size;
            m_is.getline(ptr, rest, '\n');
            uint64_t read_chars = m_is.gcount();
            m_read_chars += read_chars;

            if (!m_is) {
                if (read_chars + 1 == rest) {  // '\n' not found
                    m_is.clear();
                    size += read_chars;
                    break;
                }
            } else if (!eof()) {
                assert(read_chars > 0);
                --read_chars;  // discard the delimiter
            }

            if (read_chars == 0) {  // empty line was read
                empty_line_was_read = true;
                break;
            }

            size += read_chars;
            ptr += read_chars;
        }

        buffer.resize(size);
        return empty_line_was_read;
    }

    bool eof() const { return m_is.eof(); }

    uint64_t read_chars() const { return m_read_chars; }

private:
    std::istream& m_is;
    uint64_t m_buffer_size;
    uint64_t m_read_chars;
};

}  // namespace sshash
