#pragma once

#include "file_merging_iterator.hpp"
#include "../buckets_statistics.hpp"

namespace sshash {

#pragma pack(push, 4)
struct bucket_pair {
    bucket_pair(uint64_t id, uint32_t size) : id(id), size(size) {}
    uint64_t id;
    uint32_t size;

    bool operator>(bucket_pair other) const { return id > other.id; }
};
#pragma pack(pop)

struct bucket_pairs_iterator : std::forward_iterator_tag {
    bucket_pairs_iterator(bucket_pair const* begin, bucket_pair const* end)
        : m_begin(begin)
        , m_end(end)
        , m_i(0)
        , m_val(0)  // first returned val is always 0
    {}

    uint64_t operator*() const { return m_val; }
    void operator++() {
        ++m_i;
        if (m_begin != nullptr and m_end != nullptr and m_i == (*m_begin).id) {
            m_val += (*m_begin).size;
            ++m_begin;
            assert(m_begin <= m_end);
        }
    }

private:
    bucket_pair const* m_begin;
    bucket_pair const* m_end;
    uint64_t m_i;
    uint64_t m_val;
};

struct bucket_pairs {
    static constexpr uint64_t ram_limit = 0.25 * essentials::GB;

    bucket_pairs(std::string const& tmp_dirname)
        : m_buffer_size(0)
        , m_num_files_to_merge(0)
        , m_run_identifier(pthash::clock_type::now().time_since_epoch().count())
        , m_tmp_dirname(tmp_dirname) {
        m_buffer_size = ram_limit / sizeof(bucket_pair);
        std::cout << "m_buffer_size " << m_buffer_size << std::endl;
    }

    void emplace_back(uint64_t id, uint32_t size) {
        if (m_buffer.size() == m_buffer_size) sort_and_flush();
        m_buffer.emplace_back(id, size);
    }

    void sort_and_flush() {
        std::cout << "sorting buffer..." << std::endl;
        std::sort(m_buffer.begin(), m_buffer.end(),
                  [](bucket_pair const& x, bucket_pair const& y) { return x.id < y.id; });

        auto tmp_output_filename = get_tmp_output_filename(m_num_files_to_merge);
        std::cout << "saving to file '" << tmp_output_filename << "'..." << std::endl;
        std::ofstream out(tmp_output_filename.c_str(), std::ofstream::binary);
        if (!out.is_open()) throw std::runtime_error("cannot open file");
        out.write(reinterpret_cast<char const*>(m_buffer.data()),
                  m_buffer.size() * sizeof(bucket_pair));
        out.close();

        m_buffer.clear();
        ++m_num_files_to_merge;
    }

    void finalize() {
        if (!m_buffer.empty()) sort_and_flush();
    }

    std::string get_bucket_pairs_filename() const {
        assert(m_num_files_to_merge > 0);
        if (m_num_files_to_merge == 1) return get_tmp_output_filename(0);
        std::stringstream filename;
        filename << m_tmp_dirname << "/sshash.tmp.run_" << m_run_identifier << ".bucket_pairs.bin";
        return filename.str();
    }

    struct files_name_iterator {
        files_name_iterator(bucket_pairs const* ptr) : m_id(0), m_ptr(ptr) {}

        std::string operator*() { return m_ptr->get_tmp_output_filename(m_id); }
        void operator++() { ++m_id; }

    private:
        uint64_t m_id;
        bucket_pairs const* m_ptr;
    };

    files_name_iterator files_name_iterator_begin() { return files_name_iterator(this); }

    void merge() {
        if (m_num_files_to_merge <= 1) return;
        assert(m_num_files_to_merge > 1);

        std::cout << " == files to merge = " << m_num_files_to_merge << std::endl;

        typedef bytes_iterator<bucket_pair> bytes_iterator_type;
        file_merging_iterator<bytes_iterator_type> fm_iterator(files_name_iterator_begin(),
                                                               m_num_files_to_merge);

        std::ofstream out(get_bucket_pairs_filename().c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open file");

        uint64_t num_written_pairs = 0;
        while (fm_iterator.has_next()) {
            auto file_it = *fm_iterator;
            bucket_pair val = *file_it;
            out.write(reinterpret_cast<char const*>(&val), sizeof(bucket_pair));
            num_written_pairs += 1;
            if (num_written_pairs % 50000000 == 0) {
                std::cout << "num_written_pairs = " << num_written_pairs << std::endl;
            }
            fm_iterator.next();
        }
        std::cout << "num_written_pairs = " << num_written_pairs << std::endl;

        out.close();
        fm_iterator.close();

        /* remove tmp files */
        for (uint64_t i = 0; i != m_num_files_to_merge; ++i) {
            auto tmp_output_filename = get_tmp_output_filename(i);
            std::remove(tmp_output_filename.c_str());
        }

        std::vector<bucket_pair>().swap(m_buffer);
    }

    uint64_t num_files_to_merge() const { return m_num_files_to_merge; }
    void remove_tmp_file() { std::remove(get_bucket_pairs_filename().c_str()); }

private:
    uint64_t m_buffer_size;
    uint64_t m_num_files_to_merge;
    uint64_t m_run_identifier;
    std::string m_tmp_dirname;
    std::vector<bucket_pair> m_buffer;

    std::string get_tmp_output_filename(uint64_t id) const {
        std::stringstream filename;
        filename << m_tmp_dirname << "/sshash.tmp.run_" << m_run_identifier << ".bucket_pairs."
                 << id << ".bin";
        return filename.str();
    }
};

buckets_statistics build_index(parse_data& data, minimizers const& m_minimizers, buckets& m_buckets,
                               build_configuration const& build_config) {
    uint64_t num_buckets = m_minimizers.size();
    uint64_t num_kmers = data.num_kmers;
    uint64_t num_super_kmers = data.strings.num_super_kmers();

    pthash::compact_vector::builder offsets;
    offsets.resize(num_super_kmers, std::ceil(std::log2(data.strings.num_bits() / 2)));

    std::cout << "bits_per_offset = ceil(log2(" << data.strings.num_bits() / 2
              << ")) = " << std::ceil(std::log2(data.strings.num_bits() / 2)) << std::endl;

    mm::file_source<minimizer_tuple> input(data.minimizers.get_minimizers_filename(),
                                           mm::advice::sequential);

    bucket_pairs bucket_pairs_manager(build_config.tmp_dirname);
    uint64_t num_singletons = 0;
    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size()); it.has_next();
         it.next()) {
        uint32_t list_size = it.list().size();
        assert(list_size > 0);
        if (list_size != 1) {
            uint64_t bucket_id = m_minimizers.lookup(it.minimizer());
            bucket_pairs_manager.emplace_back(bucket_id + 1, list_size - 1);
        } else {
            ++num_singletons;
        }
    }
    bucket_pairs_manager.finalize();

    std::cout << "num_singletons " << num_singletons << "/" << num_buckets << " ("
              << (num_singletons * 100.0) / num_buckets << "%)" << std::endl;

    if (bucket_pairs_manager.num_files_to_merge() > 0) {
        bucket_pairs_manager.merge();
        mm::file_source<bucket_pair> bucket_pairs_file(
            bucket_pairs_manager.get_bucket_pairs_filename(), mm::advice::sequential);
        bucket_pairs_iterator iterator(bucket_pairs_file.data(),
                                       bucket_pairs_file.data() + bucket_pairs_file.size());
        m_buckets.num_super_kmers_before_bucket.encode(iterator, num_buckets + 1,
                                                       num_super_kmers - num_buckets);
        bucket_pairs_file.close();
        bucket_pairs_manager.remove_tmp_file();
    } else {
        /* all buckets are singletons, thus pass an empty iterator that always returns 0 */
        bucket_pairs_iterator iterator(nullptr, nullptr);
        m_buckets.num_super_kmers_before_bucket.encode(iterator, num_buckets + 1,
                                                       num_super_kmers - num_buckets);
    }

    buckets_statistics buckets_stats(num_buckets, num_kmers, num_super_kmers);

    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size()); it.has_next();
         it.next()) {
        uint64_t bucket_id = m_minimizers.lookup(it.minimizer());
        uint64_t base = m_buckets.num_super_kmers_before_bucket.access(bucket_id) + bucket_id;
        uint64_t num_super_kmers_in_bucket =
            (m_buckets.num_super_kmers_before_bucket.access(bucket_id + 1) + bucket_id + 1) - base;
        assert(num_super_kmers_in_bucket > 0);
        buckets_stats.add_num_super_kmers_in_bucket(num_super_kmers_in_bucket);
        uint64_t offset_pos = 0;
        auto list = it.list();
        for (auto [offset, num_kmers_in_super_kmer] : list) {
            offsets.set(base + offset_pos++, offset);
            buckets_stats.add_num_kmers_in_super_kmer(num_super_kmers_in_bucket,
                                                      num_kmers_in_super_kmer);
        }
        assert(offset_pos == num_super_kmers_in_bucket);
    }

    m_buckets.pieces.encode(data.strings.pieces.begin(), data.strings.pieces.size(),
                            data.strings.pieces.back());
    offsets.build(m_buckets.offsets);
    m_buckets.strings.swap(data.strings.strings);

    input.close();

    return buckets_stats;
}

}  // namespace sshash