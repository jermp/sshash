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
    list_type(minimizer_tuple const* begin, minimizer_tuple const* end)
        : m_begin(begin), m_end(end), m_size(std::distance(begin, end)) {}

    struct iterator {
        iterator(minimizer_tuple const* begin) : m_begin(begin) {}

        inline std::pair<uint64_t, uint64_t> operator*() const {
            return {(*m_begin).offset, (*m_begin).num_kmers_in_super_kmer};
        }

        inline void operator++() { ++m_begin; }
        bool operator==(iterator const& other) const { return m_begin == other.m_begin; }
        bool operator!=(iterator const& other) const { return !(*this == other); }

    private:
        minimizer_tuple const* m_begin;
    };

    iterator begin() const { return iterator(m_begin); }
    iterator end() const { return iterator(m_end); }
    uint64_t size() const { return m_size; }

    minimizer_tuple const* begin_ptr() const { return m_begin; }
    minimizer_tuple const* end_ptr() const { return m_end; }

private:
    minimizer_tuple const* m_begin;
    minimizer_tuple const* m_end;
    uint64_t m_size;
};

struct minimizers_tuples_iterator : std::forward_iterator_tag {
    typedef minimizer_tuple value_type;

    minimizers_tuples_iterator(minimizer_tuple const* begin, minimizer_tuple const* end)
        : m_list_begin(begin), m_list_end(begin), m_end(end) {
        m_list_end = next_begin();
    }

    inline uint64_t minimizer() const { return (*m_list_begin).minimizer; }
    inline uint64_t operator*() const { return minimizer(); }
    inline void next() {
        m_list_begin = m_list_end;
        m_list_end = next_begin();
    }
    inline void operator++() { next(); }
    bool has_next() const { return m_list_begin != m_end; }
    list_type list() const { return list_type(m_list_begin, m_list_end); }

private:
    minimizer_tuple const* m_list_begin;
    minimizer_tuple const* m_list_end;
    minimizer_tuple const* m_end;

    minimizer_tuple const* next_begin() {
        minimizer_tuple const* begin = m_list_begin;
        uint64_t prev_minimizer = (*begin).minimizer;
        while (begin != m_end) {
            ++begin;
            uint64_t curr_minimizer = (*begin).minimizer;
            if (curr_minimizer != prev_minimizer) break;
        }
        return begin;
    }
};

struct minimizers_tuples {
    static constexpr uint64_t ram_limit = 0.5 * essentials::GB;

    minimizers_tuples(std::string tmp_dirname = constants::default_tmp_dirname)
        : m_buffer_size(0)
        , m_num_files_to_merge(0)
        , m_run_identifier(pthash::clock_type::now().time_since_epoch().count())
        , m_tmp_dirname(tmp_dirname) {
        m_buffer_size = ram_limit / sizeof(minimizer_tuple);
        std::cout << "m_buffer_size " << m_buffer_size << std::endl;
    }

    void emplace_back(uint64_t minimizer, uint64_t offset, uint64_t num_kmers_in_super_kmer) {
        if (m_buffer.size() == m_buffer_size) sort_and_flush();
        m_buffer.emplace_back(minimizer, offset, num_kmers_in_super_kmer);
    }

    minimizer_tuple& back() { return m_buffer.back(); }

    void sort_and_flush() {
        std::cout << "sorting buffer..." << std::endl;
        std::sort(m_buffer.begin(), m_buffer.end(),
                  [](minimizer_tuple const& x, minimizer_tuple const& y) {
                      return (x.minimizer < y.minimizer) or
                             (x.minimizer == y.minimizer and x.offset < y.offset);
                  });
        auto tmp_output_filename = get_tmp_output_filename(m_num_files_to_merge);
        std::cout << "saving to file '" << tmp_output_filename << "'..." << std::endl;
        std::ofstream out(tmp_output_filename.c_str(), std::ofstream::binary);
        if (!out.is_open()) throw std::runtime_error("cannot open file");
        out.write(reinterpret_cast<char const*>(m_buffer.data()),
                  m_buffer.size() * sizeof(minimizer_tuple));
        out.close();

        m_buffer.clear();
        ++m_num_files_to_merge;
    }

    void finalize() {
        if (!m_buffer.empty()) sort_and_flush();
    }

    std::string get_minimizers_filename() const {
        std::stringstream filename;
        filename << m_tmp_dirname << "/sshash.tmp.run_" << m_run_identifier << ".minimizers.bin";
        return filename.str();
    }

    void merge() {
        if (m_num_files_to_merge == 0) return;

        assert(m_num_files_to_merge > 0);
        std::cout << "files to merge = " << m_num_files_to_merge << std::endl;

        struct iterator_type {
            iterator_type(minimizer_tuple const* b, minimizer_tuple const* e) : begin(b), end(e) {}
            minimizer_tuple const* begin;
            minimizer_tuple const* end;
        };
        std::vector<iterator_type> iterators;
        std::vector<uint32_t> idx_heap;
        iterators.reserve(m_num_files_to_merge);
        idx_heap.reserve(m_num_files_to_merge);
        std::vector<mm::file_source<minimizer_tuple>> mm_files(m_num_files_to_merge);

        auto heap_idx_comparator = [&](uint32_t i, uint32_t j) {
            minimizer_tuple const* begin_i = iterators[i].begin;
            minimizer_tuple const* begin_j = iterators[j].begin;
            if ((*begin_i).minimizer != (*begin_j).minimizer) {
                return (*begin_i).minimizer > (*begin_j).minimizer;
            }
            return (*begin_i).offset > (*begin_j).offset;
        };

        auto advance_heap_head = [&]() {
            uint32_t idx = idx_heap.front();
            iterators[idx].begin += 1;
            if (iterators[idx].begin != iterators[idx].end) {  // percolate down the head
                uint64_t pos = 0;
                uint64_t size = idx_heap.size();
                while (2 * pos + 1 < size) {
                    uint64_t i = 2 * pos + 1;
                    if (i + 1 < size and heap_idx_comparator(idx_heap[i], idx_heap[i + 1])) ++i;
                    if (heap_idx_comparator(idx_heap[i], idx_heap[pos])) break;
                    std::swap(idx_heap[pos], idx_heap[i]);
                    pos = i;
                }
            } else {
                std::pop_heap(idx_heap.begin(), idx_heap.end(), heap_idx_comparator);
                idx_heap.pop_back();
            }
        };

        /* create the input iterators and make the heap */
        for (uint64_t i = 0; i != m_num_files_to_merge; ++i) {
            auto tmp_output_filename = get_tmp_output_filename(i);
            mm_files[i].open(tmp_output_filename, mm::advice::sequential);
            iterators.emplace_back(mm_files[i].data(), mm_files[i].data() + mm_files[i].size());
            idx_heap.push_back(i);
        }
        std::make_heap(idx_heap.begin(), idx_heap.end(), heap_idx_comparator);

        std::ofstream out(get_minimizers_filename().c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open file");

        uint64_t num_written_tuples = 0;
        while (!idx_heap.empty()) {
            minimizer_tuple const* begin = iterators[idx_heap.front()].begin;
            out.write(reinterpret_cast<char const*>(begin), sizeof(minimizer_tuple));
            num_written_tuples += 1;
            if (num_written_tuples % 100000 == 0) {
                std::cout << "num_written_tuples = " << num_written_tuples << std::endl;
            }
            advance_heap_head();
        }
        std::cout << "num_written_tuples = " << num_written_tuples << std::endl;
        out.close();

        /* remove tmp files */
        for (uint64_t i = 0; i != m_num_files_to_merge; ++i) {
            mm_files[i].close();
            auto tmp_output_filename = get_tmp_output_filename(i);
            std::remove(tmp_output_filename.c_str());
        }

        std::vector<minimizer_tuple>().swap(m_buffer);
        m_num_files_to_merge = 0;  // any other call to merge() will do nothing
    }

    void remove_tmp_file() { std::remove(get_minimizers_filename().c_str()); }

private:
    uint64_t m_buffer_size;
    uint64_t m_num_files_to_merge;
    uint64_t m_run_identifier;
    std::string m_tmp_dirname;
    std::vector<minimizer_tuple> m_buffer;

    std::string get_tmp_output_filename(uint64_t id) {
        std::stringstream filename;
        filename << m_tmp_dirname << "/sshash.tmp.run_" << m_run_identifier << ".minimizers." << id
                 << ".bin";
        return filename.str();
    };
};

}  // namespace sshash