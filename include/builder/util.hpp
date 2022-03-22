#pragma once

namespace sshash {

struct empty_bucket_runtime_error : public std::runtime_error {
    empty_bucket_runtime_error()
        : std::runtime_error("try a different choice of l or change seed") {}
};

struct parse_runtime_error : public std::runtime_error {
    parse_runtime_error()
        : std::runtime_error("did you provide an input file with abundance counts?") {}
};

void expect(char got, char expected) {
    if (got != expected) {
        std::cout << "got '" << got << "' but expected '" << expected << "'" << std::endl;
        throw parse_runtime_error();
    }
}

struct compact_string_pool {
    compact_string_pool() {}

    struct builder {
        builder(uint64_t k) : k(k), offset(0), num_strings(0) {}

        void build(compact_string_pool& pool) {
            pool.k = k;
            pool.num_strings = num_strings;
            pool.pieces.swap(pieces);
            pool.strings.build(&bvb_strings);
        }

        void append(char const* string, uint64_t size, bool glue) {
            assert(size >= k);
            uint64_t prefix = 0;
            if (glue) {
                prefix = k - 1;
            } else {
                /* otherwise, start a new piece */
                pieces.push_back(bvb_strings.size() / 2);
            }
            for (uint64_t i = prefix; i != size; ++i) {
                bvb_strings.append_bits(util::char_to_uint64(string[i]), 2);
            }
            num_strings += 1;
            offset = bvb_strings.size() / 2;
        }

        void finalize() {
            /* So pieces will be of size p+1, where p is the number of DNA strings
               in the input file. */
            pieces.push_back(bvb_strings.size() / 2);
            assert(pieces.front() == 0);

            /* Push a final sentinel (dummy) symbol to avoid bounds' checking. */
            bvb_strings.append_bits(0, 2);
        }

        uint64_t k;
        uint64_t offset;
        uint64_t num_strings;
        std::vector<uint64_t> pieces;
        pthash::bit_vector_builder bvb_strings;
    };

    uint64_t num_bits() const { return strings.size(); }
    uint64_t size() const { return num_strings; }

    uint64_t k;
    uint64_t num_strings;
    std::vector<uint64_t> pieces;
    pthash::bit_vector strings;
};

typedef uint8_t num_kmers_in_string_uint_type;

#pragma pack(push, 1)
struct minimizer_tuple {
    minimizer_tuple(uint64_t minimizer, uint64_t offset, uint64_t num_kmers_in_string)
        : minimizer(minimizer), offset(offset), num_kmers_in_string(num_kmers_in_string) {}
    uint64_t minimizer;
    uint64_t offset;
    num_kmers_in_string_uint_type num_kmers_in_string;
};
#pragma pack(pop)

struct list_type {
    list_type(std::vector<minimizer_tuple>::iterator begin,
              std::vector<minimizer_tuple>::iterator end)
        : m_begin(begin), m_size(std::distance(begin, end)) {}

    struct iterator {
        iterator(std::vector<minimizer_tuple>::iterator begin) : m_begin(begin) {}

        inline std::pair<uint64_t, uint64_t> operator*() const {
            return {(*m_begin).offset, (*m_begin).num_kmers_in_string};
        }

        inline void operator++() { ++m_begin; }
        bool operator==(iterator const& other) const { return m_begin == other.m_begin; }
        bool operator!=(iterator const& other) const { return !(*this == other); }

    private:
        std::vector<minimizer_tuple>::iterator m_begin;
    };

    iterator begin() const { return iterator(m_begin); }
    iterator end() const { return iterator(m_begin + m_size); }
    uint64_t size() const { return m_size; }

private:
    std::vector<minimizer_tuple>::iterator m_begin;
    uint64_t m_size;
};

struct minimizers_tuples {
    minimizers_tuples() {}

    // void reserve(uint64_t n) {
    //     tuples.reserve(n);
    // }

    void emplace_back(uint64_t minimizer, uint64_t offset, uint64_t num_kmers_in_string) {
        tuples.emplace_back(minimizer, offset, num_kmers_in_string);
    }

    minimizer_tuple& back() { return tuples.back(); }

    void sort() {
        std::sort(tuples.begin(), tuples.end(),
                  [](minimizer_tuple const& x, minimizer_tuple const& y) {
                      return (x.minimizer < y.minimizer) or
                             (x.minimizer == y.minimizer and x.offset < y.offset);
                  });
    }

    struct iterator {
        iterator(std::vector<minimizer_tuple>::iterator b, std::vector<minimizer_tuple>::iterator e)
            : begin(b), end(e) {}

        inline uint64_t minimizer() const { return (*begin).minimizer; }
        inline uint64_t operator*() const { return minimizer(); }
        list_type list() const { return list_type(begin, next_begin()); }
        bool has_next() const { return begin != end; }
        inline void next() { begin = next_begin(); }
        inline void operator++() { next(); }

        inline iterator operator+(uint64_t /*offset*/) const {
            assert(false);
            return iterator(end, end);
        }

    private:
        std::vector<minimizer_tuple>::iterator begin, end;

        std::vector<minimizer_tuple>::iterator next_begin() const {
            auto i = begin;
            uint64_t prev = (*i).minimizer;
            while (true) {
                ++i;
                if (i == end) break;
                uint64_t curr = (*i).minimizer;
                if (curr != prev) break;
            }
            assert(i > begin);
            return i;
        }
    };

    iterator begin() { return iterator(tuples.begin(), tuples.end()); }

private:
    std::vector<minimizer_tuple> tuples;
};

}  // namespace sshash