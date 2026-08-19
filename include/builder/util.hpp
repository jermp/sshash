#pragma once

#include <vector>
#include <atomic>

#include "file_merging_iterator.hpp"
#include "parallel_sort.hpp"

namespace sshash {

struct parse_runtime_error : public std::runtime_error {
    parse_runtime_error() : std::runtime_error("did you provide an input file with weights?") {}
};

[[maybe_unused]] static void expect(char got, char expected) {
    if (got != expected) {
        std::cout << "got '" << got << "' but expected '" << expected << "'" << std::endl;
        throw parse_runtime_error();
    }
}

typedef uint8_t num_kmers_in_super_kmer_uint_type;

#pragma pack(push, 2)
struct minimizer_tuple {
    minimizer_tuple() {}
    minimizer_tuple(minimizer_info mini_info, uint64_t num_kmers_in_super_kmer)
        : minimizer(mini_info.minimizer)
        , pos_in_seq(mini_info.pos_in_seq)
        , pos_in_kmer(mini_info.pos_in_kmer)
        , num_kmers_in_super_kmer(num_kmers_in_super_kmer) {}

    bool operator>(minimizer_tuple other) const {
        if (minimizer != other.minimizer) return minimizer > other.minimizer;
        return pos_in_seq > other.pos_in_seq;
    }
    bool operator<(minimizer_tuple other) const {
        if (minimizer != other.minimizer) return minimizer < other.minimizer;
        return pos_in_seq < other.pos_in_seq;
    }

    static minimizer_tuple max() {
        return {{uint64_t(-1), uint64_t(-1), uint64_t(-1)}, uint64_t(-1)};
    }

    uint64_t minimizer;
    uint64_t pos_in_seq;
    num_kmers_in_super_kmer_uint_type pos_in_kmer;
    num_kmers_in_super_kmer_uint_type num_kmers_in_super_kmer;
};
#pragma pack(pop)

inline std::ostream& operator<<(std::ostream& os, minimizer_tuple const& mt) {
    os << "minimizer = " << mt.minimizer << std::endl;
    os << "pos_in_seq = " << mt.pos_in_seq << std::endl;
    os << "pos_in_kmer = " << int(mt.pos_in_kmer) << std::endl;
    os << "num_kmers_in_super_kmer = " << int(mt.num_kmers_in_super_kmer) << std::endl;
    return os;
}

/*
    The bucket of a minimizer: its tuples, sorted by pos_in_seq. The minimizer
    is forward (see `util::compute_minimizer`), so a locus is never abandoned
    and later re-selected: each tuple carries a distinct pos_in_seq, and the
    number of tuples (= super-kmers) equals the number of minimizer positions,
    which is what size() returns. `minimizers_tuples::merge` checks this.
*/
struct bucket_type {
    bucket_type(minimizer_tuple const* begin, minimizer_tuple const* end)
        : m_begin(begin)
        , m_end(end)  //
    {
        assert([&] { /* one tuple per minimizer position: the scheme is forward */
                     for (auto it = begin; it + 1 < end; ++it) {
                         if (it->pos_in_seq == (it + 1)->pos_in_seq) return false;
                     }
                     return true;
        }());
    }

    struct iterator {
        iterator(minimizer_tuple const* begin) : m_begin(begin) {}

        inline minimizer_tuple operator*() const { return *m_begin; }
        inline void operator++() { ++m_begin; }
        bool operator==(iterator const& other) const { return m_begin == other.m_begin; }
        bool operator!=(iterator const& other) const { return !(*this == other); }

    private:
        minimizer_tuple const* m_begin;
    };

    iterator begin() const { return iterator(m_begin); }
    iterator end() const { return iterator(m_end); }

    uint64_t size() const { return std::distance(m_begin, m_end); }

    minimizer_tuple const* begin_ptr() const { return m_begin; }
    minimizer_tuple const* end_ptr() const { return m_end; }

private:
    minimizer_tuple const* m_begin;
    minimizer_tuple const* m_end;
};

/*
    Iterate over the "bucket" of a minimizer, i.e.,
    the sorted list of minimizer tuples
    (minimizer, pos_in_seq, pos_in_kmer, num_kmers_in_superkmer).
*/
struct minimizers_tuples_iterator {
    typedef minimizer_tuple value_type;
    using iterator_category = std::forward_iterator_tag;

    minimizers_tuples_iterator(minimizer_tuple const* begin, minimizer_tuple const* end)
        : m_bucket_begin(begin), m_bucket_end(begin), m_end(end) {
        m_bucket_end = next_begin();
    }

    inline uint64_t minimizer() const { return (*m_bucket_begin).minimizer; }
    inline uint64_t operator*() const { return minimizer(); }
    inline void next() {
        m_bucket_begin = m_bucket_end;
        m_bucket_end = next_begin();
    }
    inline void operator++() { next(); }
    bool has_next() const { return m_bucket_begin != m_end; }
    bucket_type bucket() const { return bucket_type(m_bucket_begin, m_bucket_end); }

private:
    minimizer_tuple const* m_bucket_begin;
    minimizer_tuple const* m_bucket_end;
    minimizer_tuple const* m_end;

    minimizer_tuple const* next_begin() {
        if (m_bucket_begin == m_end) return m_end;
        minimizer_tuple const* begin = m_bucket_begin;
        uint64_t prev_minimizer = begin->minimizer;
        while (++begin != m_end and begin->minimizer == prev_minimizer) {}
        return begin;
    }
};

struct minimizers_tuples {
    minimizers_tuples() {}
    minimizers_tuples(build_configuration const& build_config)
        : m_num_minimizers(0)
        , m_num_minimizer_positions(0)
        , m_run_identifier(pthash::clock_type::now().time_since_epoch().count())
        , m_build_config(build_config)  //
    {
        init();
    }

    void init() { m_num_files_to_merge = 0; }

    void sort_and_flush(std::vector<minimizer_tuple>& buffer) {
        parallel_sort(buffer, m_build_config.num_threads,
                      [](minimizer_tuple const& x, minimizer_tuple const& y) {
                          return (x.minimizer < y.minimizer) or
                                 (x.minimizer == y.minimizer and x.pos_in_seq < y.pos_in_seq);
                      });
        uint64_t id = m_num_files_to_merge.fetch_add(1);
        auto tmp_output_filename = get_tmp_output_filename(id);
        if (m_build_config.verbose) {
            std::cout << "saving to file '" << tmp_output_filename << "'..." << std::endl;
        }
        std::ofstream out(tmp_output_filename.c_str(), std::ofstream::binary);
        if (!out.is_open()) throw std::runtime_error("cannot open file");
        out.write(reinterpret_cast<char const*>(buffer.data()),
                  buffer.size() * sizeof(minimizer_tuple));
        out.close();
        buffer.clear();
    }

    std::string get_minimizers_filename() const {
        assert(m_num_files_to_merge > 0);
        std::stringstream filename;
        filename << m_build_config.tmp_dirname << "/sshash.tmp.run_" << m_run_identifier
                 << ".minimizers.bin";
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
        if (m_num_files_to_merge == 0) return;

        if (m_num_files_to_merge == 1) {
            std::rename(get_tmp_output_filename(0).c_str(), get_minimizers_filename().c_str());
            if (m_num_minimizers != 0) return;

            assert(m_num_minimizers == 0);
            assert(m_num_minimizer_positions == 0);
            mm::file_source<minimizer_tuple> input(get_minimizers_filename(),
                                                   mm::advice::sequential);
            count_and_check_forward(input.data(), input.data() + input.size());
            input.close();
            return;
        }

        assert(m_num_files_to_merge > 1);
        file_merging_iterator<minimizer_tuple> fm_iterator(files_name_iterator_begin(),
                                                           m_num_files_to_merge);

        std::ofstream out(get_minimizers_filename().c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open file");

        if (m_build_config.verbose) {
            std::cout << "saving to file '" << get_minimizers_filename() << "'" << std::endl;
        }

        m_num_minimizers = 0;
        m_num_minimizer_positions = 0;
        uint64_t prev_minimizer = constants::invalid_uint64;
        uint64_t prev_pos_in_seq = constants::invalid_uint64;
        while (fm_iterator.has_next()) {
            minimizer_tuple mt = *fm_iterator;
            count_and_check_forward_one(mt, prev_minimizer, prev_pos_in_seq);
            out.write(reinterpret_cast<char const*>(&mt), sizeof(minimizer_tuple));
            if (m_build_config.verbose and m_num_minimizer_positions % 100'000'000 == 0) {
                std::cout << "processed " << m_num_minimizer_positions << " minimizer tuples"
                          << std::endl;
            }
            fm_iterator.next();
        }

        out.close();
        fm_iterator.close();

        /* remove tmp files */
        for (uint64_t i = 0; i != m_num_files_to_merge; ++i) {
            auto tmp_output_filename = get_tmp_output_filename(i);
            std::remove(tmp_output_filename.c_str());
        }
    }

    uint64_t num_files_to_merge() const { return m_num_files_to_merge; }
    uint64_t num_minimizers() const { return m_num_minimizers; }

    /* One tuple per minimizer position, the scheme being forward (checked by
       `merge`), so this is also the number of super-kmers and the number of
       records in the merged tuples file. */
    uint64_t num_minimizer_positions() const { return m_num_minimizer_positions; }

    void remove_tmp_file() { std::remove(get_minimizers_filename().c_str()); }

private:
    std::atomic<uint64_t> m_num_files_to_merge;
    uint64_t m_num_minimizers;
    uint64_t m_num_minimizer_positions;
    uint64_t m_run_identifier;
    build_configuration m_build_config;

    /*
        Count minimizers and minimizer positions over tuples sorted by
        (minimizer, pos_in_seq), checking the forwardness requirement: the
        scheme never re-selects an abandoned locus, so no (minimizer,
        pos_in_seq) pair may appear twice. Everything downstream sizes the
        index by the tuple count, so a violation must stop the build.
    */
    void count_and_check_forward_one(minimizer_tuple const& mt, uint64_t& prev_minimizer,
                                     uint64_t& prev_pos_in_seq)  //
    {
        if (mt.minimizer != prev_minimizer) {
            prev_minimizer = mt.minimizer;
            ++m_num_minimizers;
        } else if (mt.pos_in_seq == prev_pos_in_seq) {
            throw std::runtime_error(
                "the minimizer scheme is not forward: "
                "a (minimizer, position) pair occurs in more than one super-kmer");
        }
        prev_pos_in_seq = mt.pos_in_seq;
        ++m_num_minimizer_positions;
    }

    void count_and_check_forward(minimizer_tuple const* begin, minimizer_tuple const* end) {
        uint64_t prev_minimizer = constants::invalid_uint64;
        uint64_t prev_pos_in_seq = constants::invalid_uint64;
        for (; begin != end; ++begin) {
            count_and_check_forward_one(*begin, prev_minimizer, prev_pos_in_seq);
        }
    }

    std::string get_tmp_output_filename(uint64_t id) const {
        std::stringstream filename;
        filename << m_build_config.tmp_dirname << "/sshash.tmp.run_" << m_run_identifier
                 << ".minimizers." << id << ".bin";
        return filename.str();
    }
};

}  // namespace sshash
