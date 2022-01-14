#pragma once

#include <vector>

namespace sshash {

template <typename T>
struct fixed_size_deque {
    fixed_size_deque(uint64_t size) : m_begin(0), m_end(0), m_size(0), m_buffer(size) {
        assert(size > 0);
    }

    bool empty() const { return m_size == 0; }

    void push_back(T const& val) {
        assert(m_size < m_buffer.size());
        m_buffer[m_end] = val;
        if (++m_end == m_buffer.size()) m_end = 0;
        ++m_size;
    }

    void pop_front() {
        assert(!empty());
        if (++m_begin == m_buffer.size()) m_begin = 0;
        --m_size;
    }

    void pop_back() {
        assert(!empty());
        if (m_end == 0) m_end = m_buffer.size();
        --m_end;
        --m_size;
    }

    T const& front() const { return m_buffer[m_begin]; }

    T const& back() const {
        if (m_end == 0) return m_buffer.back();
        return m_buffer[m_end - 1];
    }

private:
    uint64_t m_begin;
    uint64_t m_end;
    uint64_t m_size;
    std::vector<T> m_buffer;
};

template <typename Hasher = util::murmurhash2_64>
struct minimizer_enumerator {
    minimizer_enumerator() {}

    minimizer_enumerator(uint64_t k, uint64_t m, uint64_t seed)
        : m_k(k)
        , m_m(m)
        , m_seed(seed)
        , m_position(0)
        , m_mask((1ULL << (2 * m_m)) - 1)
        , m_q(k - m + 1) /* deque cannot contain more than k - m + 1 elements  */
    {}

    template <bool reverse = false>
    uint64_t next(uint64_t kmer, bool clear) {
        if (clear) {
            if constexpr (reverse) {
                for (uint64_t i = 0; i != m_k - m_m + 1; ++i) {
                    uint64_t minimizer = (kmer >> (2 * (m_k - m_m - i))) & m_mask;
                    eat(minimizer);
                }
            } else {
                for (uint64_t i = 0; i != m_k - m_m + 1; ++i) {
                    uint64_t minimizer = kmer & m_mask;
                    kmer >>= 2;
                    eat(minimizer);
                }
            }
        } else {
            if constexpr (reverse) {
                uint64_t minimizer = kmer & m_mask;
                eat(minimizer);
            } else {
                uint64_t minimizer = kmer >> (2 * (m_k - m_m));
                eat(minimizer);
            }
        }
        return m_q.front().minimizer;
    }

private:
    uint64_t m_k;
    uint64_t m_m;
    uint64_t m_seed;
    uint64_t m_position;
    uint64_t m_mask;

    struct minimizer_t {
        minimizer_t() {}
        minimizer_t(uint64_t h, uint64_t p, uint64_t m) : hash(h), position(p), minimizer(m) {}
        uint64_t hash, position, minimizer;
    };

    /* NOTE: we could use a std::deque<minimizer_t> here,
             but std::deque is terribly space-inefficient. */
    fixed_size_deque<minimizer_t> m_q;

    void eat(uint64_t minimizer) {
        uint64_t hash = Hasher::hash(minimizer, m_seed);

        /* Removes from front elements which are no longer in the window */
        while (!m_q.empty() and m_position + m_m - 1 >= m_k and
               m_q.front().position <= (m_position + m_m - 1) - m_k) {
            m_q.pop_front();
        }
        /* Removes from back elements which are no longer useful */
        while (!m_q.empty() and hash < m_q.back().hash) m_q.pop_back();

        m_q.push_back({hash, m_position, minimizer});
        m_position += 1;
    }
};

}  // namespace sshash