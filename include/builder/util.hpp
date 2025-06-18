#pragma once

#include "file_merging_iterator.hpp"
#include "parallel_sort.hpp"

namespace sshash {

[[maybe_unused]] static void print_time(double time, uint64_t num_kmers,
                                        std::string const& message) {
    std::cout << "=== " << message << " " << time / 1000000 << " [sec] ("
              << (time * 1000) / num_kmers << " [ns/kmer])" << std::endl;
}

struct parse_runtime_error : public std::runtime_error {
    parse_runtime_error() : std::runtime_error("did you provide an input file with weights?") {}
};

[[maybe_unused]] static void expect(char got, char expected) {
    if (got != expected) {
        std::cout << "got '" << got << "' but expected '" << expected << "'" << std::endl;
        throw parse_runtime_error();
    }
}

template <class kmer_t>
struct compact_string_pool {
    compact_string_pool() {}

    struct builder {
        builder(uint64_t k) : k(k), offset(0), num_super_kmers(0) {
            const uint64_t num_bits = 8 * 8 * essentials::GB;  // 8 GB of memory
            bvb_strings.reserve(num_bits);
        }

        void build(compact_string_pool& pool) {
            pool.m_num_super_kmers = num_super_kmers;
            pool.pieces.swap(pieces);
            bvb_strings.build(pool.strings);
        }

        void append(char const* string, uint64_t size, bool glue) {
            assert(size >= k);
            uint64_t prefix = 0;
            if (glue) {
                prefix = k - 1;
            } else { /* otherwise, start a new piece */
                assert(bvb_strings.num_bits() % kmer_t::bits_per_char == 0);
                pieces.push_back(bvb_strings.num_bits() / kmer_t::bits_per_char);
            }
            for (uint64_t i = prefix; i != size; ++i) {
                bvb_strings.append_bits(kmer_t::char_to_uint(string[i]), kmer_t::bits_per_char);
            }
            num_super_kmers += 1;
            offset = bvb_strings.num_bits() / kmer_t::bits_per_char;
        }

        void finalize() {
            /* So pieces will be of size p+1, where p is the number of DNA sequences
               in the input file. */
            pieces.push_back(bvb_strings.num_bits() / kmer_t::bits_per_char);
            assert(pieces.front() == 0);

            /* Push a final sentinel (dummy) value to avoid bounds' checking in
               kmer_iterator::fill_buff(). */
            bvb_strings.append_bits(0, kmer_t::uint_kmer_bits);
        }

        uint64_t k;
        uint64_t offset;
        uint64_t num_super_kmers;
        std::vector<uint64_t> pieces;
        bits::bit_vector::builder bvb_strings;
    };

    uint64_t num_bits() const { return strings.num_bits(); }
    uint64_t num_super_kmers() const { return m_num_super_kmers; }

    std::vector<uint64_t> pieces;
    bits::bit_vector strings;

private:
    uint64_t m_num_super_kmers;
};

typedef uint8_t num_kmers_in_super_kmer_uint_type;

#pragma pack(push, 1)
struct minimizer_tuple {
    minimizer_tuple() {}
    minimizer_tuple(uint64_t minimizer, uint64_t offset, uint64_t num_kmers_in_super_kmer)
        : minimizer(minimizer), offset(offset), num_kmers_in_super_kmer(num_kmers_in_super_kmer) {}
    uint64_t minimizer;
    uint64_t offset;
    num_kmers_in_super_kmer_uint_type num_kmers_in_super_kmer;

    bool operator>(minimizer_tuple other) const {
        if (minimizer != other.minimizer) return minimizer > other.minimizer;
        return offset > other.offset;
    }
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

template <typename ValueType>
struct bytes_iterator {
    bytes_iterator(uint8_t const* begin, uint8_t const* end) : m_begin(begin), m_end(end) {}

    void next() { m_begin += sizeof(ValueType); }
    bool has_next() const { return m_begin != m_end; }
    ValueType operator*() const { return *reinterpret_cast<ValueType const*>(m_begin); }

private:
    uint8_t const* m_begin;
    uint8_t const* m_end;
};

struct minimizers_tuples {
    minimizers_tuples() {}
    minimizers_tuples(build_configuration const& build_config)
        /* allocate half memory because __gnu::parallel_sort is not in-place */
        : m_buffer_size((build_config.ram_limit_in_GiB * essentials::GiB) /
                        (2 * sizeof(minimizer_tuple)))
        , m_num_minimizers(0)
        , m_run_identifier(pthash::clock_type::now().time_since_epoch().count())
        , m_num_threads(build_config.num_threads)
        , m_tmp_dirname(build_config.tmp_dirname)  //
    {
        init();
        m_buffer.reserve(m_buffer_size);
    }

    void init() { m_num_files_to_merge = 0; }

    void emplace_back(uint64_t minimizer, uint64_t offset, uint64_t num_kmers_in_super_kmer) {
        if (m_buffer.size() == m_buffer_size) sort_and_flush();
        m_buffer.emplace_back(minimizer, offset, num_kmers_in_super_kmer);
    }

    minimizer_tuple& back() { return m_buffer.back(); }

    void sort_and_flush() {
        // std::cout << "sorting buffer..." << std::endl;
        m_temp_buffer.resize(m_buffer.size());
        parallel_sort(m_buffer, m_temp_buffer, m_num_threads,
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
        // std::cout << "DONE" << std::endl;
        m_buffer.clear();
        ++m_num_files_to_merge;
    }

    void finalize() {
        if (!m_buffer.empty()) sort_and_flush();
    }

    std::string get_minimizers_filename() const {
        assert(m_num_files_to_merge > 0);
        if (m_num_files_to_merge == 1) return get_tmp_output_filename(0);
        std::stringstream filename;
        filename << m_tmp_dirname << "/sshash.tmp.run_" << m_run_identifier << ".minimizers.bin";
        return filename.str();
    }

    struct files_name_iterator {
        files_name_iterator(minimizers_tuples const* ptr) : m_id(0), m_ptr(ptr) {}

        std::string operator*() { return m_ptr->get_tmp_output_filename(m_id); }
        void operator++() { ++m_id; }

    private:
        uint64_t m_id;
        minimizers_tuples const* m_ptr;
    };

    files_name_iterator files_name_iterator_begin() { return files_name_iterator(this); }

    void merge() {
        assert(m_buffer.empty());
        std::vector<minimizer_tuple>().swap(m_buffer);
        std::vector<minimizer_tuple>().swap(m_temp_buffer);

        if (m_num_files_to_merge == 0) return;

        if (m_num_files_to_merge == 1) {
            if (m_num_minimizers != 0) return;
            /* just count num. distinct minimizers and do not write twice on disk */
            mm::file_source<minimizer_tuple> input(get_minimizers_filename(),
                                                   mm::advice::sequential);
            for (minimizers_tuples_iterator it(input.data(), input.data() + input.size());
                 it.has_next(); it.next()) {
                ++m_num_minimizers;
            }
            input.close();
            return;
        }

        std::cout << " == files to merge = " << m_num_files_to_merge << std::endl;

        assert(m_num_files_to_merge > 1);
        typedef bytes_iterator<minimizer_tuple> bytes_iterator_type;
        file_merging_iterator<bytes_iterator_type> fm_iterator(files_name_iterator_begin(),
                                                               m_num_files_to_merge);

        std::ofstream out(get_minimizers_filename().c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open file");

        m_num_minimizers = 0;
        uint64_t num_written_tuples = 0;
        uint64_t prev_minimizer = constants::invalid_uint64;
        while (fm_iterator.has_next()) {
            auto file_it = *fm_iterator;
            minimizer_tuple mt = *file_it;
            out.write(reinterpret_cast<char const*>(&mt), sizeof(minimizer_tuple));
            num_written_tuples += 1;
            if (mt.minimizer != prev_minimizer) {
                prev_minimizer = mt.minimizer;
                ++m_num_minimizers;
            }
            if (num_written_tuples % 50000000 == 0) {
                std::cout << "num_written_tuples = " << num_written_tuples << std::endl;
            }
            fm_iterator.next();
        }
        std::cout << "num_written_tuples = " << num_written_tuples << std::endl;

        out.close();
        fm_iterator.close();

        /* remove tmp files */
        for (uint64_t i = 0; i != m_num_files_to_merge; ++i) {
            auto tmp_output_filename = get_tmp_output_filename(i);
            std::remove(tmp_output_filename.c_str());
        }
    }

    uint64_t num_minimizers() const { return m_num_minimizers; }
    uint64_t num_files_to_merge() const { return m_num_files_to_merge; }
    uint64_t buffer_size() const { return m_buffer_size; }
    std::vector<minimizer_tuple>& buffer() { return m_buffer; }

    void remove_tmp_file() { std::remove(get_minimizers_filename().c_str()); }

private:
    uint64_t m_buffer_size;
    uint64_t m_num_files_to_merge;
    uint64_t m_num_minimizers;
    uint64_t m_run_identifier;
    uint64_t m_num_threads;
    std::string m_tmp_dirname;
    std::vector<minimizer_tuple> m_buffer;
    std::vector<minimizer_tuple> m_temp_buffer;

    std::string get_tmp_output_filename(uint64_t id) const {
        std::stringstream filename;
        filename << m_tmp_dirname << "/sshash.tmp.run_" << m_run_identifier << ".minimizers." << id
                 << ".bin";
        return filename.str();
    }
};

}  // namespace sshash