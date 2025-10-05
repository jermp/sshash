#pragma once

#include "include/buckets_statistics.hpp"

namespace sshash {

struct bucket_size_iterator {
    using iterator_category = std::forward_iterator_tag;

    bucket_size_iterator(minimizer_tuple const* begin, minimizer_tuple const* end)
        : m_val(0)  // first returned value is always 0
        , m_it(begin, end) {}

    uint64_t operator*() const { return m_val; }

    void operator++() {
        if (!m_it.has_next()) return;
        uint64_t size = m_it.bucket().size();
        assert(size > 0);
        m_val += size - 1;  // directly compute the cumulative sum
        m_it.next();
    }

private:
    uint64_t m_val;
    minimizers_tuples_iterator m_it;
};

template <class kmer_t>
buckets_statistics build_sparse_index(parse_data<kmer_t>& data, buckets<kmer_t>& m_buckets,
                                      build_configuration const& build_config)  //
{
    const uint64_t num_kmers = data.num_kmers;
    const uint64_t num_minimizer_positions = data.minimizers.num_minimizer_positions();
    const uint64_t num_super_kmers = data.minimizers.num_super_kmers();
    const uint64_t num_buckets = data.minimizers.num_minimizers();
    const uint64_t num_threads = build_config.num_threads;

    bits::compact_vector::builder offsets_builder;
    offsets_builder.resize(num_minimizer_positions,
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

    essentials::timer_type timer;
    timer.start();
    {
        bucket_size_iterator iterator(begin, end);
        m_buckets.bucket_sizes.encode(iterator, num_buckets + 1,
                                      num_minimizer_positions - num_buckets);
    }
    timer.stop();
    std::cout << "encoding bucket sizes: " << timer.elapsed() / 1000000 << " [sec]" << std::endl;

    timer.reset();

    buckets_statistics buckets_stats(num_buckets, num_kmers, num_minimizer_positions);

    timer.start();
    const uint64_t block_size = (num_super_kmers + num_threads - 1) / num_threads;
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

    std::vector<buckets_statistics> threads_buckets_stats;
    threads_buckets_stats.reserve(num_threads);

    auto exe = [&](const uint64_t thread_id) {
        assert(thread_id + 1 < offsets.size());
        const uint64_t offset_begin = offsets[thread_id];
        const uint64_t offset_end = offsets[thread_id + 1];
        auto& tbs = threads_buckets_stats[thread_id];
        for (minimizers_tuples_iterator it(begin + offset_begin, begin + offset_end);  //
             it.has_next();                                                            //
             it.next())                                                                //
        {
            const uint64_t bucket_id = it.minimizer();
            const auto [begin, end] = m_buckets.locate_bucket(bucket_id);
            assert(end > begin);
            const uint64_t bucket_size = end - begin;
            assert(bucket_size == it.bucket().size());
            tbs.add_bucket_size(bucket_size);
            uint64_t pos = 0;
            auto bucket = it.bucket();
            uint64_t prev_pos_in_seq = constants::invalid_uint64;
            for (auto mt : bucket) {
                if (mt.pos_in_seq != prev_pos_in_seq) {
                    offsets_builder.set(begin + pos++, mt.pos_in_seq);
                    prev_pos_in_seq = mt.pos_in_seq;
                }
                tbs.add_num_kmers_in_super_kmer(bucket_size, mt.num_kmers_in_super_kmer);
            }
            assert(pos == bucket_size);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    assert(offsets.size() <= num_threads + 1);
    for (uint64_t thread_id = 0; thread_id + 1 < size(offsets); ++thread_id) {
        threads_buckets_stats.emplace_back(num_buckets, num_kmers, num_minimizer_positions);
        threads.emplace_back(exe, thread_id);
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
    m_buckets.pieces.encode(data.pieces.begin(), data.pieces.size(), data.pieces.back());
    offsets_builder.build(m_buckets.offsets);
    m_buckets.strings.swap(data.strings);
    timer.stop();
    std::cout << "encoding string boundaries and building offsets: " << timer.elapsed() / 1000000
              << " [sec]" << std::endl;

    return buckets_stats;
}

}  // namespace sshash