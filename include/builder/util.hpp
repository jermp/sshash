#pragma once

namespace sshash {

void print_time(double time, uint64_t num_kmers, std::string const& message) {
    std::cout << "=== " << message << " " << time / 1000000 << " [sec] ("
              << (time * 1000) / num_kmers << " [ns/kmer])" << std::endl;
}

struct empty_bucket_runtime_error : public std::runtime_error {
    empty_bucket_runtime_error()
        : std::runtime_error("try a different choice of l or change seed") {}
};

struct parse_runtime_error : public std::runtime_error {
    parse_runtime_error() : std::runtime_error("did you provide an input file with weights?") {}
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
        builder(uint64_t k) : k(k), offset(0), num_super_kmers(0) {}

        void build(compact_string_pool& pool) {
            pool.m_num_super_kmers = num_super_kmers;
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
            num_super_kmers += 1;
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
        uint64_t num_super_kmers;
        std::vector<uint64_t> pieces;
        pthash::bit_vector_builder bvb_strings;
    };

    uint64_t num_bits() const { return strings.size(); }
    uint64_t num_super_kmers() const { return m_num_super_kmers; }

    std::vector<uint64_t> pieces;
    pthash::bit_vector strings;

private:
    uint64_t m_num_super_kmers;
};

typedef uint8_t num_kmers_in_super_kmer_uint_type;

#pragma pack(push, 1)
struct minimizer_tuple {
    minimizer_tuple(uint64_t minimizer, uint64_t offset, uint64_t num_kmers_in_super_kmer)
        : minimizer(minimizer), offset(offset), num_kmers_in_super_kmer(num_kmers_in_super_kmer) {}
    uint64_t minimizer;
    uint64_t offset;
    num_kmers_in_super_kmer_uint_type num_kmers_in_super_kmer;
};
#pragma pack(pop)

struct list_type {
    list_type(std::vector<minimizer_tuple>::iterator begin,
              std::vector<minimizer_tuple>::iterator end)
        : m_begin(begin), m_size(std::distance(begin, end)) {}

    struct iterator {
        iterator(std::vector<minimizer_tuple>::iterator begin) : m_begin(begin) {}

        inline std::pair<uint64_t, uint64_t> operator*() const {
            return {(*m_begin).offset, (*m_begin).num_kmers_in_super_kmer};
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

/* suitable for external-memory iterator */
struct minimizers_tuples_iterator : std::forward_iterator_tag {
    typedef minimizer_tuple value_type;

    minimizers_tuples_iterator(uint8_t const* begin, uint8_t const* end)
        : m_begin(begin), m_end(end) {}

    inline uint64_t operator*() const { return minimizer(); }
    inline void operator++() { next_begin(); }

private:
    uint8_t const* m_begin;
    uint8_t const* m_end;

    inline uint64_t minimizer() const { return *reinterpret_cast<uint64_t const*>(m_begin); }

    void next_begin() {
        uint64_t prev = minimizer();
        while (true) {
            m_begin += sizeof(minimizer_tuple);
            if (m_begin == m_end) break;
            uint64_t curr = minimizer();
            if (curr != prev) break;
        }
    }
};

struct minimizers_tuples {
    minimizers_tuples() : m_run_identifier(pthash::clock_type::now().time_since_epoch().count()) {}

    // void reserve(uint64_t n) {
    //     m_tuples.reserve(n);
    // }

    void emplace_back(uint64_t minimizer, uint64_t offset, uint64_t num_kmers_in_super_kmer) {
        m_tuples.emplace_back(minimizer, offset, num_kmers_in_super_kmer);
    }

    minimizer_tuple& back() { return m_tuples.back(); }

    void sort() {
        std::sort(m_tuples.begin(), m_tuples.end(),
                  [](minimizer_tuple const& x, minimizer_tuple const& y) {
                      return (x.minimizer < y.minimizer) or
                             (x.minimizer == y.minimizer and x.offset < y.offset);
                  });
    }

    /* internal-memory iterator */
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

    iterator begin() { return iterator(m_tuples.begin(), m_tuples.end()); }

    std::string get_minimizers_filename() const {
        std::stringstream filename;
        // TODO: use a tmp_dir
        filename << "sshash.tmp.run_" << m_run_identifier << ".minimizers.bin";
        return filename.str();
    }

    void flush() {
        std::ofstream out(get_minimizers_filename().c_str(), std::ofstream::binary);
        if (!out.is_open()) throw std::runtime_error("cannot open file");
        out.write(reinterpret_cast<char const*>(m_tuples.data()),
                  m_tuples.size() * sizeof(minimizer_tuple));
        out.close();
    }

private:
    uint64_t m_run_identifier;
    std::vector<minimizer_tuple> m_tuples;
};

}  // namespace sshash