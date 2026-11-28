#pragma once

#include <vector>
#include <atomic>
#include <numeric>

#include <sys/resource.h>

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

/*
    Every run being merged is an open file and a memory mapping: merge at most
    half as many runs as the process may have of either, leaving the rest to the
    process, and at least two.
*/
static inline uint64_t merge_fan_in(uint64_t max_open_files, uint64_t max_map_count) {
    return std::max<uint64_t>(std::min(max_open_files, max_map_count) / 2, 2);
}

/*
    Raise the soft limit on open files to the hard limit, then return the merge
    fan-in for the limits in effect (vm.max_map_count is read where available).
*/
static inline uint64_t max_merge_fan_in() {
    uint64_t max_open_files = uint64_t(-1);
    rlimit limit;
    if (getrlimit(RLIMIT_NOFILE, &limit) == 0) {
        if (limit.rlim_cur < limit.rlim_max) {
            rlimit raised = limit;
            raised.rlim_cur = limit.rlim_max;
            if (setrlimit(RLIMIT_NOFILE, &raised) == 0) limit = raised;
        }
        if (limit.rlim_cur != RLIM_INFINITY) max_open_files = limit.rlim_cur;
    }
    uint64_t max_map_count = uint64_t(-1);
    std::ifstream in("/proc/sys/vm/max_map_count");
    if (uint64_t value; in >> value) max_map_count = value;
    return merge_fan_in(max_open_files, max_map_count);
}

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
        files_name_iterator(minimizers_tuples const* ptr, uint64_t const* ids)
            : m_ids(ids), m_ptr(ptr) {}

        std::string operator*() { return m_ptr->get_tmp_output_filename(*m_ids); }
        void operator++() { ++m_ids; }

    private:
        uint64_t const* m_ids;
        minimizers_tuples const* m_ptr;
    };

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
        std::vector<uint64_t> ids(m_num_files_to_merge);
        std::iota(ids.begin(), ids.end(), 0);

        /*
            Merge at most `fan_in` runs at a time into longer runs until the rest
            fit in one final merge. The new runs take fresh ids.
        */
        const uint64_t fan_in = max_merge_fan_in();
        if (m_build_config.verbose) {
            std::cout << "merging " << ids.size() << " runs, at most " << fan_in << " at a time"
                      << std::endl;
        }
        while (ids.size() > fan_in) {
            std::vector<uint64_t> merged_ids;
            for (uint64_t i = 0; i < ids.size(); i += fan_in) {
                const uint64_t n = std::min<uint64_t>(fan_in, ids.size() - i);
                if (n == 1) {
                    merged_ids.push_back(ids[i]);
                    continue;
                }
                const uint64_t id = m_num_files_to_merge.fetch_add(1);
                std::ofstream out(get_tmp_output_filename(id).c_str(), std::ofstream::binary);
                if (!out.is_open()) throw std::runtime_error("cannot open file");
                merge_runs(ids.data() + i, n, [&](minimizer_tuple const& mt) {
                    out.write(reinterpret_cast<char const*>(&mt), sizeof(minimizer_tuple));
                });
                out.close();
                merged_ids.push_back(id);
            }
            if (m_build_config.verbose) {
                std::cout << "merged " << ids.size() << " runs into " << merged_ids.size()
                          << " runs" << std::endl;
            }
            ids.swap(merged_ids);
        }

        std::ofstream out(get_minimizers_filename().c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open file");

        if (m_build_config.verbose) {
            std::cout << "saving to file '" << get_minimizers_filename() << "'" << std::endl;
        }

        m_num_minimizers = 0;
        m_num_minimizer_positions = 0;
        uint64_t prev_minimizer = constants::invalid_uint64;
        uint64_t prev_pos_in_seq = constants::invalid_uint64;
        merge_runs(ids.data(), ids.size(), [&](minimizer_tuple const& mt) {
            count_and_check_forward_one(mt, prev_minimizer, prev_pos_in_seq);
            out.write(reinterpret_cast<char const*>(&mt), sizeof(minimizer_tuple));
            if (m_build_config.verbose and m_num_minimizer_positions % 100'000'000 == 0) {
                std::cout << "processed " << m_num_minimizer_positions << " minimizer tuples"
                          << std::endl;
            }
        });

        out.close();
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

    /*
        Merge the sorted runs with the given ids, passing each tuple in order
        to `append`, then remove the runs.
    */
    template <typename Append>
    void merge_runs(uint64_t const* ids, uint64_t num_runs, Append append) {
        file_merging_iterator<minimizer_tuple> fm_iterator(files_name_iterator(this, ids),
                                                           num_runs);
        while (fm_iterator.has_next()) {
            append(*fm_iterator);
            fm_iterator.next();
        }
        fm_iterator.close();
        for (uint64_t i = 0; i != num_runs; ++i) {
            std::remove(get_tmp_output_filename(ids[i]).c_str());
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
