#pragma once

#include "external/pthash/include/pthash.hpp"

namespace sshash {

template <class kmer_t>
void build_skew_index(skew_index<kmer_t>& m_skew_index,         //
                      parse_data<kmer_t>& data,                 //
                      buckets<kmer_t> const& m_buckets,         //
                      build_configuration const& build_config,  //
                      buckets_statistics const& buckets_stats)  //
{
    const uint64_t min_log2_size = m_skew_index.min_log2;
    const uint64_t max_log2_size = m_skew_index.max_log2;
    const uint64_t min_size = 1ULL << min_log2_size;
    const uint64_t k = build_config.k;
    assert(build_config.k > 0 and build_config.m <= build_config.k);

    m_skew_index.log2_max_bucket_size = std::ceil(std::log2(buckets_stats.max_bucket_size()));

    std::cout << "max_bucket_size " << buckets_stats.max_bucket_size() << std::endl;
    std::cout << "log2_max_bucket_size " << m_skew_index.log2_max_bucket_size << std::endl;

    mm::file_source<minimizer_tuple> input(data.minimizers.get_minimizers_filename(),
                                           mm::advice::sequential);

    uint64_t num_buckets_in_skew_index = 0;
    uint64_t num_super_kmers_in_skew_index = 0;
    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size()); it.has_next();
         it.next())  //
    {
        auto bucket = it.bucket();
        if (bucket.size() > min_size) {
            num_super_kmers_in_skew_index += bucket.num_super_kmers();
            ++num_buckets_in_skew_index;
        }
    }
    std::cout << "num_buckets_in_skew_index " << num_buckets_in_skew_index << "/"
              << buckets_stats.num_buckets() << " ("
              << (num_buckets_in_skew_index * 100.0) / buckets_stats.num_buckets() << "%)"
              << std::endl;

    if (num_buckets_in_skew_index == 0) {
        input.close();
        return;
    }

    std::vector<bucket_type> buckets;
    buckets.reserve(num_buckets_in_skew_index);
    std::vector<minimizer_tuple> tuples;  // backed memory
    tuples.reserve(num_super_kmers_in_skew_index);
    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size()); it.has_next();
         it.next())  //
    {
        auto bucket = it.bucket();
        if (bucket.size() > min_size) {
            minimizer_tuple const* begin = tuples.data() + tuples.size();
            std::copy(bucket.begin_ptr(), bucket.end_ptr(), std::back_inserter(tuples));
            minimizer_tuple const* end = tuples.data() + tuples.size();
            buckets.push_back(bucket_type(begin, end));
        }
    }
    assert(buckets.size() == num_buckets_in_skew_index);
    input.close();

    std::sort(buckets.begin(), buckets.end(),
              [](bucket_type const& x, bucket_type const& y) { return x.size() < y.size(); });

    uint64_t num_partitions = max_log2_size - min_log2_size + 1;
    if (buckets_stats.max_bucket_size() < (1ULL << max_log2_size)) {
        num_partitions = m_skew_index.log2_max_bucket_size - min_log2_size;
    }
    std::cout << "num_partitions " << num_partitions << std::endl;

    std::vector<uint64_t> num_kmers(num_partitions, 0);
    m_skew_index.mphfs.resize(num_partitions);
    m_skew_index.positions.resize(num_partitions);

    {
        std::cout << "computing sizes of partitions..." << std::endl;

        uint64_t partition_id = 0;
        uint64_t lower = min_size;
        uint64_t upper = 2 * lower;
        uint64_t num_kmers_in_skew_index = 0;
        for (uint64_t i = 0; i <= buckets.size(); ++i) {
            while (i == buckets.size() or buckets[i].size() > upper) {
                std::cout << "  partition_id = " << partition_id
                          << ": num_kmers belonging to buckets of size > " << lower
                          << " and <= " << upper << ": " << num_kmers[partition_id] << std::endl;
                num_kmers_in_skew_index += num_kmers[partition_id];

                if (i == buckets.size()) break;

                lower = upper;
                upper = 2 * lower;
                partition_id += 1;
                if (partition_id == num_partitions - 1) { upper = buckets_stats.max_bucket_size(); }
            }

            if (i == buckets.size()) break;

            assert(buckets[i].size() > lower and buckets[i].size() <= upper);
            for (auto mt : buckets[i]) { num_kmers[partition_id] += mt.num_kmers_in_super_kmer; }
        }
        assert(partition_id == num_partitions - 1);
        std::cout << "num_kmers_in_skew_index " << num_kmers_in_skew_index << " ("
                  << (num_kmers_in_skew_index * 100.0) / buckets_stats.num_kmers() << "%)"
                  << std::endl;
        assert(num_kmers_in_skew_index ==
               std::accumulate(num_kmers.begin(), num_kmers.end(), uint64_t(0)));
    }

    {
        std::cout << "building partitions..." << std::endl;

        pthash::build_configuration mphf_build_config;
        mphf_build_config.lambda = build_config.lambda;
        mphf_build_config.alpha = 0.94;
        mphf_build_config.seed = util::get_seed_for_hash_function(build_config);
        mphf_build_config.verbose = false;
        mphf_build_config.num_threads = build_config.num_threads;
        mphf_build_config.avg_partition_size = constants::avg_partition_size;

        uint64_t partition_id = 0;
        uint64_t lower = min_size;
        uint64_t upper = 2 * lower;
        uint64_t num_bits_per_pos = min_log2_size + 1;

        /* Temporary storage for kmers and positions within a partition. */
        std::vector<kmer_t> kmers;
        std::vector<uint32_t> positions_in_bucket;
        bits::compact_vector::builder cvb_positions;
        kmers.reserve(num_kmers[partition_id]);
        positions_in_bucket.reserve(num_kmers[partition_id]);
        cvb_positions.resize(num_kmers[partition_id], num_bits_per_pos);

        for (uint64_t i = 0; i <= buckets.size(); ++i) {
            while (i == buckets.size() or buckets[i].size() > upper) {
                std::cout << "  lower = " << lower << "; upper = " << upper
                          << "; num_bits_per_pos = " << num_bits_per_pos
                          << "; num_kmers_in_partition = " << kmers.size() << std::endl;
                assert(num_kmers[partition_id] == kmers.size());
                assert(num_kmers[partition_id] == positions_in_bucket.size());

                if (num_kmers[partition_id] > 0)  //
                {
                    if (build_config.verbose) {
                        const uint64_t avg_partition_size =
                            pthash::compute_avg_partition_size(kmers.size(), mphf_build_config);
                        const uint64_t num_partitions =
                            pthash::compute_num_partitions(kmers.size(), avg_partition_size);
                        assert(num_partitions > 0);
                        std::cout << "    building MPHF with " << mphf_build_config.num_threads
                                  << " threads and " << num_partitions
                                  << " partitions (avg. partition size = " << avg_partition_size
                                  << ")..." << std::endl;
                    }

                    auto& mphf = m_skew_index.mphfs[partition_id];
                    mphf.build_in_internal_memory(kmers.begin(), kmers.size(), mphf_build_config);

                    std::cout << "    built mphs[" << partition_id << "] for " << kmers.size()
                              << " kmers; bits/key = "
                              << static_cast<double>(mphf.num_bits()) / mphf.num_keys()
                              << std::endl;

                    for (uint64_t i = 0; i != kmers.size(); ++i) {
                        kmer_t kmer = kmers[i];
                        uint64_t pos = mphf(kmer);
                        uint32_t pos_in_bucket = positions_in_bucket[i];
                        cvb_positions.set(pos, pos_in_bucket);
                    }
                    auto& positions = m_skew_index.positions[partition_id];
                    cvb_positions.build(positions);

                    std::cout << "    built positions[" << partition_id << "] for "
                              << positions.size() << " kmers; bits/key = "
                              << (positions.num_bytes() * 8.0) / positions.size() << std::endl;
                }

                if (i == buckets.size()) break;

                lower = upper;
                upper = 2 * lower;
                num_bits_per_pos += 1;
                partition_id += 1;
                if (partition_id == num_partitions - 1) {
                    upper = buckets_stats.max_bucket_size();
                    num_bits_per_pos = m_skew_index.log2_max_bucket_size;
                }

                kmers.clear();
                positions_in_bucket.clear();
                kmers.reserve(num_kmers[partition_id]);
                positions_in_bucket.reserve(num_kmers[partition_id]);
                cvb_positions.resize(num_kmers[partition_id], num_bits_per_pos);
            }

            if (i == buckets.size()) break;

            assert(buckets[i].size() > lower and buckets[i].size() <= upper);
            uint64_t pos_in_bucket = -1;
            uint64_t prev_pos_in_seq = constants::invalid_uint64;
            for (auto mt : buckets[i])  //
            {
                if (mt.pos_in_seq != prev_pos_in_seq) {
                    prev_pos_in_seq = mt.pos_in_seq;
                    ++pos_in_bucket;
                }
                assert(mt.pos_in_seq >= mt.pos_in_kmer);
                const uint64_t starting_pos_of_super_kmer = mt.pos_in_seq - mt.pos_in_kmer;
                kmer_iterator<kmer_t> it(m_buckets.strings, k,
                                         kmer_t::bits_per_char * starting_pos_of_super_kmer);
                for (uint64_t i = 0; i != mt.num_kmers_in_super_kmer; ++i) {
                    auto kmer = it.get();
                    if (build_config.canonical) { /* take the canonical kmer */
                        auto kmer_rc = kmer;
                        kmer_rc.reverse_complement_inplace(k);
                        kmer = std::min(kmer, kmer_rc);
                    }
                    kmers.push_back(kmer);
                    positions_in_bucket.push_back(pos_in_bucket);
                    it.next();
                }
                assert(pos_in_bucket < (1ULL << cvb_positions.width()));
            }
        }
        assert(partition_id == num_partitions - 1);
    }

    std::cout << "num_bits_for_skew_index " << m_skew_index.num_bits() << "("
              << static_cast<double>(m_skew_index.num_bits()) / buckets_stats.num_kmers()
              << " [bits/kmer])" << std::endl;
}

}  // namespace sshash