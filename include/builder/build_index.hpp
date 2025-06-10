#pragma once

#include "file_merging_iterator.hpp"
#include "include/buckets_statistics.hpp"

namespace sshash {

#pragma pack(push, 4)
struct bucket_pair {
    bucket_pair(uint64_t id, uint32_t size) : id(id), size(size) {}
    uint64_t id;
    uint32_t size;

    bool operator>(bucket_pair other) const { return id > other.id; }
};
#pragma pack(pop)

struct bucket_pairs_iterator {
    using iterator_category = std::forward_iterator_tag;

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
        : m_buffer_size(ram_limit / sizeof(bucket_pair))
        , m_num_files_to_merge(0)
        , m_run_identifier(pthash::clock_type::now().time_since_epoch().count())
        , m_tmp_dirname(tmp_dirname) {}

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

template <class kmer_t>
buckets_statistics build_index(parse_data<kmer_t>& data, buckets<kmer_t>& m_buckets,
                               const uint64_t num_buckets,
                               build_configuration const& build_config)  //
{
    const uint64_t num_kmers = data.num_kmers;
    const uint64_t num_super_kmers = data.strings.num_super_kmers();
    const uint64_t num_threads = build_config.num_threads;

    bits::compact_vector::builder offsets_builder;
    offsets_builder.resize(num_super_kmers,
                           std::ceil(std::log2(data.strings.num_bits() / kmer_t::bits_per_char)));

    std::cout << "bits_per_offset = ceil(log2(" << data.strings.num_bits() / kmer_t::bits_per_char
              << ")) = " << std::ceil(std::log2(data.strings.num_bits() / kmer_t::bits_per_char))
              << std::endl;

    std::cout << "reading from '" << data.minimizers.get_minimizers_filename() << "'..."
              << std::endl;
    mm::file_source<minimizer_tuple> input(data.minimizers.get_minimizers_filename(),
                                           mm::advice::sequential);
    minimizer_tuple const* begin = input.data();
    minimizer_tuple const* end = input.data() + input.size();
    assert(input.size() == num_super_kmers);

    bucket_pairs bucket_pairs_manager(build_config.tmp_dirname);
    uint64_t num_singletons = 0;
    for (minimizers_tuples_iterator it(begin, end); it.has_next(); it.next()) {
        uint32_t list_size = it.list().size();
        assert(list_size > 0);
        if (list_size != 1) {
            uint64_t bucket_id = it.minimizer();
            bucket_pairs_manager.emplace_back(bucket_id + 1, list_size - 1);
        } else {
            ++num_singletons;
        }
    }
    bucket_pairs_manager.finalize();

    std::cout << "num_singletons " << num_singletons << "/" << num_buckets << " ("
              << (num_singletons * 100.0) / num_buckets << "%)" << std::endl;

    essentials::timer_type timer;
    timer.start();
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
    timer.stop();
    std::cout << "building: " << timer.elapsed() / 1000000 << " [sec]" << std::endl;

    timer.reset();

    buckets_statistics buckets_stats(num_buckets, num_kmers, num_super_kmers);

    timer.start();
    uint64_t block_size = (num_super_kmers + num_threads - 1) / num_threads;
    std::vector<uint64_t> offsets;
    offsets.reserve(num_threads + 1);
    for (uint64_t offset = -1; offset != num_super_kmers;) {
        offsets.push_back(offset + 1);
        offset = std::min<uint64_t>((offset + 1) + block_size, num_super_kmers);
        minimizer_tuple const* b = begin + offset;
        uint64_t curr_minimizer = (*b).minimizer;
        while (b + 1 < end) {  // adjust offset
            uint64_t next_minimizer = (*(b + 1)).minimizer;
            if (curr_minimizer != next_minimizer) break;
            b += 1;
            offset += 1;
        }
    }
    offsets.push_back(num_super_kmers);
    assert(offsets.size() == num_threads + 1);

    std::vector<buckets_statistics> threads_buckets_stats(num_threads);

    auto exe = [&](const uint64_t thread_id) {
        assert(thread_id + 1 < offsets.size());
        const uint64_t offset_begin = offsets[thread_id];
        const uint64_t offset_end = offsets[thread_id + 1];
        auto& tbs = threads_buckets_stats[thread_id];
        for (minimizers_tuples_iterator it(begin + offset_begin, begin + offset_end);  //
             it.has_next();                                                            //
             it.next())                                                                //
        {
            uint64_t bucket_id = it.minimizer();
            uint64_t base = m_buckets.num_super_kmers_before_bucket.access(bucket_id) + bucket_id;
            uint64_t num_super_kmers_in_bucket =
                (m_buckets.num_super_kmers_before_bucket.access(bucket_id + 1) + bucket_id + 1) -
                base;
            assert(num_super_kmers_in_bucket > 0);
            tbs.add_num_super_kmers_in_bucket(num_super_kmers_in_bucket);
            uint64_t pos = 0;
            auto list = it.list();
            for (auto [offset, num_kmers_in_super_kmer] : list) {
                offsets_builder.set(base + pos++, offset);
                tbs.add_num_kmers_in_super_kmer(num_super_kmers_in_bucket, num_kmers_in_super_kmer);
            }
            assert(pos == num_super_kmers_in_bucket);
        }
    };

    std::vector<std::thread> threads(num_threads);
    for (uint64_t thread_id = 0; thread_id != num_threads; ++thread_id) {
        threads_buckets_stats[thread_id] =
            buckets_statistics(num_buckets, num_kmers, num_super_kmers);
        threads[thread_id] = std::thread(exe, thread_id);
    }
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }
    for (auto const& tbs : threads_buckets_stats) buckets_stats += tbs;

    input.close();
    timer.stop();
    std::cout << "computing minimizers offsets: " << timer.elapsed() / 1000000 << " [sec]"
              << std::endl;

    timer.reset();
    timer.start();
    m_buckets.pieces.encode(data.strings.pieces.begin(),  //
                            data.strings.pieces.size(),   //
                            data.strings.pieces.back());  //
    offsets_builder.build(m_buckets.offsets);
    m_buckets.strings.swap(data.strings.strings);
    timer.stop();
    std::cout << "encoding: " << timer.elapsed() / 1000000 << " [sec]" << std::endl;

    return buckets_stats;
}

}  // namespace sshash