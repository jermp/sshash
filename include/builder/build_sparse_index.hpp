#pragma once

#include <mutex>
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
                                      build_configuration const& /*build_config*/)  //
{
    const uint64_t num_kmers = data.num_kmers;
    const uint64_t num_minimizer_positions = data.minimizers.num_minimizer_positions();
    const uint64_t num_buckets = data.minimizers.num_minimizers();


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

    bits::compact_vector::builder offsets_builder;
    offsets_builder.resize(num_minimizer_positions,
                           std::ceil(std::log2(data.strings.num_bits() / kmer_t::bits_per_char)));
    uint64_t prev_minimizer = constants::invalid_uint64, prev_pos_in_seq = constants::invalid_uint64, bucket_size = 0;
    for (auto mt : input) {
        if (mt.minimizer != prev_minimizer) {
            auto [bucket_begin, bucket_end] = m_buckets.locate_bucket(mt.minimizer);
            bucket_size = bucket_end - bucket_begin;
            buckets_stats.add_bucket_size(bucket_size);
            prev_minimizer = mt.minimizer;
            prev_pos_in_seq = constants::invalid_uint64;
        }
        buckets_stats.add_num_kmers_in_super_kmer(bucket_size, mt.num_kmers_in_super_kmer);
        if (mt.pos_in_seq != prev_pos_in_seq) {
            offsets_builder.push_back(mt.pos_in_seq);
            prev_pos_in_seq = mt.pos_in_seq;
        }
    }

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
