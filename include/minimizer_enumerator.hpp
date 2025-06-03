#pragma once

#include <vector>

#include "kmer.hpp"

namespace sshash {

// template <typename T>
// struct fixed_size_deque {
//     fixed_size_deque(uint64_t size) : m_begin(0), m_end(0), m_size(0), m_buffer(size) {
//         assert(size > 0);
//     }

//     bool empty() const { return m_size == 0; }

//     void push_back(T const& val) {
//         assert(m_size < m_buffer.size());
//         m_buffer[m_end] = val;
//         if (++m_end == m_buffer.size()) m_end = 0;
//         ++m_size;
//     }

//     void pop_front() {
//         assert(!empty());
//         if (++m_begin == m_buffer.size()) m_begin = 0;
//         --m_size;
//     }

//     void pop_back() {
//         assert(!empty());
//         if (m_end == 0) m_end = m_buffer.size();
//         --m_end;
//         --m_size;
//     }

//     T const& front() const { return m_buffer[m_begin]; }

//     T const& back() const {
//         if (m_end == 0) return m_buffer.back();
//         return m_buffer[m_end - 1];
//     }

// private:
//     uint64_t m_begin;
//     uint64_t m_end;
//     uint64_t m_size;
//     std::vector<T> m_buffer;
// };

// /*
//     "Monotone queue" method.
// */
// template <class kmer_t, typename Hasher = murmurhash2_64>
// struct minimizer_enumerator {
//     minimizer_enumerator() {}

//     minimizer_enumerator(uint64_t k, uint64_t m, uint64_t seed)
//         : m_k(k)
//         , m_m(m)
//         , m_seed(seed)
//         , m_position(0)
//         , m_q(k - m + 1) /* deque cannot contain more than k - m + 1 elements  */
//     {}

//     template <bool reverse = false>
//     uint64_t next(kmer_t kmer, bool clear) {
//         if (clear) {
//             for (uint64_t i = 0; i != m_k - m_m + 1; ++i) {
//                 kmer_t mmer = kmer;
//                 if constexpr (reverse) {
//                     mmer.drop_chars(m_k - m_m - i);
//                 } else {
//                     kmer.drop_char();
//                 }
//                 mmer.take_chars(m_m);
//                 eat(uint64_t(mmer));
//             }
//         } else {
//             kmer_t mmer = kmer;
//             if constexpr (reverse) {
//                 mmer.take_chars(m_m);
//             } else {
//                 mmer.drop_chars(m_k - m_m);
//             }
//             eat(uint64_t(mmer));
//         }
//         return m_q.front().value;
//     }

// private:
//     uint64_t m_k;
//     uint64_t m_m;
//     uint64_t m_seed;
//     uint64_t m_position;

//     struct mmer_t {
//         mmer_t() {}
//         mmer_t(uint64_t hash, uint64_t position, uint64_t value)
//             : hash(hash), position(position), value(value) {}
//         uint64_t hash, position;
//         uint64_t value;
//     };

//     /* NOTE: we could use a std::deque<mmer_t> here,
//              but std::deque is terribly space-inefficient. */
//     fixed_size_deque<mmer_t> m_q;

//     void eat(uint64_t mmer) {
//         uint64_t hash = Hasher::hash(mmer, m_seed);

//         /* Removes from front elements which are no longer in the window */
//         while (!m_q.empty() and m_position + m_m - 1 >= m_k and
//                m_q.front().position <= (m_position + m_m - 1) - m_k) {
//             m_q.pop_front();
//         }
//         /* Removes from back elements which are no longer useful */
//         while (!m_q.empty() and hash < m_q.back().hash) m_q.pop_back();

//         m_q.push_back({hash, m_position, mmer});
//         m_position += 1;
//     }
// };

/*
    "Re-scan" method.
*/
template <class kmer_t, typename Hasher = murmurhash2_64>
struct minimizer_enumerator {
    minimizer_enumerator() {}

    minimizer_enumerator(uint64_t k, uint64_t m, uint64_t seed)
        : m_k(k)
        , m_m(m)
        , m_seed(seed)
        , m_position(0)
        , m_min_value(constants::invalid_uint64)
        , m_min_position(0)
        , m_min_hash(constants::invalid_uint64) {}

    template <bool reverse = false>
    uint64_t next(kmer_t kmer, bool clear) {
        if (clear) {
            rescan<reverse>(kmer);
        } else {
            m_position += 1;
            if (m_min_position <= (m_position + m_m - 1) - m_k) {
                /* min leaves the window: re-scan to compute the new min */
                rescan<reverse>(kmer);
            } else {
                kmer_t mmer = kmer;
                if constexpr (reverse) {
                    mmer.take_chars(m_m);
                } else {
                    mmer.drop_chars(m_k - m_m);
                }
                uint64_t hash = Hasher::hash(uint64_t(mmer), m_seed);
                if (hash < m_min_hash) {
                    m_min_hash = hash;
                    m_min_value = uint64_t(mmer);
                    m_min_position = m_position;
                }
            }
        }
        return m_min_value;
    }

private:
    uint64_t m_k;
    uint64_t m_m;
    uint64_t m_seed;
    uint64_t m_position;
    uint64_t m_min_value, m_min_position, m_min_hash;

    template <bool reverse = false>
    void rescan(kmer_t kmer) {
        m_min_hash = constants::invalid_uint64;
        m_min_value = constants::invalid_uint64;
        for (uint64_t i = 0; i != m_k - m_m + 1; ++i) {
            kmer_t mmer = kmer;
            if constexpr (reverse) {
                mmer.drop_chars(m_k - m_m - i);
            } else {
                kmer.drop_char();
            }
            mmer.take_chars(m_m);
            uint64_t hash = Hasher::hash(uint64_t(mmer), m_seed);
            if (hash < m_min_hash) {
                m_min_hash = hash;
                m_min_value = uint64_t(mmer);
                m_min_position = m_position;
            }
            m_position += 1;
        }
    }
};

}  // namespace sshash