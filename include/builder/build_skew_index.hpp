#pragma once

#include "external/pthash/include/pthash.hpp"

namespace sshash {

template <class kmer_t>
void build_skew_index(skew_index<kmer_t>& m_skew_index, parse_data<kmer_t>& data,
                      buckets<kmer_t> const& m_buckets, build_configuration const& build_config,
                      buckets_statistics const& buckets_stats) {
    const uint64_t min_log2_size = m_skew_index.min_log2;
    const uint64_t max_log2_size = m_skew_index.max_log2;
    const uint64_t min_size = 1ULL << min_log2_size;

    m_skew_index.log2_max_num_super_kmers_in_bucket =
        std::ceil(std::log2(buckets_stats.max_num_super_kmers_in_bucket()));

    std::cout << "max_num_super_kmers_in_bucket " << buckets_stats.max_num_super_kmers_in_bucket()
              << std::endl;
    std::cout << "log2_max_num_super_kmers_in_bucket "
              << m_skew_index.log2_max_num_super_kmers_in_bucket << std::endl;

    mm::file_source<minimizer_tuple> input(data.minimizers.get_minimizers_filename(),
                                           mm::advice::sequential);

    uint64_t num_buckets_in_skew_index = 0;
    uint64_t num_super_kmers_in_skew_index = 0;
    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size()); it.has_next();
         it.next()) {
        uint64_t list_size = it.list().size();
        if (list_size > min_size) {
            num_super_kmers_in_skew_index += list_size;
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

    std::vector<list_type> lists;
    lists.reserve(num_buckets_in_skew_index);
    std::vector<minimizer_tuple> lists_tuples;  // backed memory
    lists_tuples.reserve(num_super_kmers_in_skew_index);
    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size()); it.has_next();
         it.next()) {
        auto list = it.list();
        if (list.size() > min_size) {
            minimizer_tuple const* begin = lists_tuples.data() + lists_tuples.size();
            std::copy(list.begin_ptr(), list.end_ptr(), std::back_inserter(lists_tuples));
            minimizer_tuple const* end = lists_tuples.data() + lists_tuples.size();
            lists.push_back(list_type(begin, end));
        }
    }
    assert(lists.size() == num_buckets_in_skew_index);
    input.close();

    std::sort(lists.begin(), lists.end(),
              [](list_type const& x, list_type const& y) { return x.size() < y.size(); });

    uint64_t num_partitions = max_log2_size - min_log2_size + 1;
    if (buckets_stats.max_num_super_kmers_in_bucket() < (1ULL << max_log2_size)) {
        num_partitions = m_skew_index.log2_max_num_super_kmers_in_bucket - min_log2_size;
    }
    std::cout << "num_partitions " << num_partitions << std::endl;

    std::vector<uint64_t> num_kmers_in_partition(num_partitions, 0);
    m_skew_index.mphfs.resize(num_partitions);
    m_skew_index.positions.resize(num_partitions);

    {
        std::cout << "computing sizes of partitions..." << std::endl;

        uint64_t partition_id = 0;
        uint64_t lower = min_size;
        uint64_t upper = 2 * lower;
        uint64_t num_kmers_in_skew_index = 0;
        for (uint64_t i = 0; i <= lists.size(); ++i) {
            while (i == lists.size() or lists[i].size() > upper) {
                std::cout << "  partition_id = " << partition_id
                          << ": num_kmers belonging to buckets of size > " << lower
                          << " and <= " << upper << ": " << num_kmers_in_partition[partition_id]
                          << std::endl;
                num_kmers_in_skew_index += num_kmers_in_partition[partition_id];

                if (i == lists.size()) break;

                lower = upper;
                upper = 2 * lower;
                partition_id += 1;
                if (partition_id == num_partitions - 1) {
                    upper = buckets_stats.max_num_super_kmers_in_bucket();
                }
            }

            if (i == lists.size()) break;

            assert(lists[i].size() > lower and lists[i].size() <= upper);
            for (auto [offset, num_kmers_in_super_kmer] : lists[i]) {
                (void)offset;  // unused
                num_kmers_in_partition[partition_id] += num_kmers_in_super_kmer;
            }
        }
        assert(partition_id == num_partitions - 1);
        std::cout << "num_kmers_in_skew_index " << num_kmers_in_skew_index << " ("
                  << (num_kmers_in_skew_index * 100.0) / buckets_stats.num_kmers() << "%)"
                  << std::endl;
        assert(num_kmers_in_skew_index == std::accumulate(num_kmers_in_partition.begin(),
                                                          num_kmers_in_partition.end(),
                                                          uint64_t(0)));
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

        /* tmp storage for keys and super_kmer_ids ******/
        std::vector<kmer_t> keys_in_partition;
        std::vector<uint32_t> super_kmer_ids_in_partition;
        keys_in_partition.reserve(num_kmers_in_partition[partition_id]);
        super_kmer_ids_in_partition.reserve(num_kmers_in_partition[partition_id]);
        bits::compact_vector::builder cvb_positions;
        cvb_positions.resize(num_kmers_in_partition[partition_id], num_bits_per_pos);
        /*******/

        for (uint64_t i = 0; i <= lists.size(); ++i) {
            while (i == lists.size() or lists[i].size() > upper) {
                std::cout << "  lower " << lower << "; upper " << upper << "; num_bits_per_pos "
                          << num_bits_per_pos << "; keys_in_partition.size() "
                          << keys_in_partition.size() << std::endl;
                assert(num_kmers_in_partition[partition_id] == keys_in_partition.size());
                assert(num_kmers_in_partition[partition_id] == super_kmer_ids_in_partition.size());

                if (num_kmers_in_partition[partition_id] > 0)  //
                {
                    if (build_config.verbose) {
                        const uint64_t avg_partition_size = pthash::compute_avg_partition_size(
                            keys_in_partition.size(), mphf_build_config);
                        const uint64_t num_partitions = pthash::compute_num_partitions(
                            keys_in_partition.size(), avg_partition_size);
                        assert(num_partitions > 0);
                        std::cout << "    building MPHF with " << mphf_build_config.num_threads
                                  << " threads and " << num_partitions
                                  << " partitions (avg. partition size = " << avg_partition_size
                                  << ")..." << std::endl;
                    }

                    auto& mphf = m_skew_index.mphfs[partition_id];
                    mphf.build_in_internal_memory(keys_in_partition.begin(),
                                                  keys_in_partition.size(), mphf_build_config);

                    std::cout << "    built mphs[" << partition_id << "] for "
                              << keys_in_partition.size() << " keys; bits/key = "
                              << static_cast<double>(mphf.num_bits()) / mphf.num_keys()
                              << std::endl;

                    for (uint64_t i = 0; i != keys_in_partition.size(); ++i) {
                        kmer_t kmer = keys_in_partition[i];
                        uint64_t pos = mphf(kmer);
                        uint32_t super_kmer_id = super_kmer_ids_in_partition[i];
                        cvb_positions.set(pos, super_kmer_id);
                    }
                    auto& positions = m_skew_index.positions[partition_id];
                    cvb_positions.build(positions);

                    std::cout << "    built positions[" << partition_id << "] for "
                              << positions.size() << " keys; bits/key = "
                              << (positions.num_bytes() * 8.0) / positions.size() << std::endl;
                }

                if (i == lists.size()) break;

                lower = upper;
                upper = 2 * lower;
                num_bits_per_pos += 1;
                partition_id += 1;
                if (partition_id == num_partitions - 1) {
                    upper = buckets_stats.max_num_super_kmers_in_bucket();
                    num_bits_per_pos = m_skew_index.log2_max_num_super_kmers_in_bucket;
                }

                keys_in_partition.clear();
                super_kmer_ids_in_partition.clear();
                keys_in_partition.reserve(num_kmers_in_partition[partition_id]);
                super_kmer_ids_in_partition.reserve(num_kmers_in_partition[partition_id]);
                cvb_positions.resize(num_kmers_in_partition[partition_id], num_bits_per_pos);
            }

            if (i == lists.size()) break;

            assert(lists[i].size() > lower and lists[i].size() <= upper);
            uint64_t super_kmer_id = 0;
            for (auto [offset, num_kmers_in_super_kmer] : lists[i]) {
                kmer_iterator<kmer_t> it(m_buckets.strings, build_config.k,
                                         kmer_t::bits_per_char * offset);
                for (uint64_t i = 0; i != num_kmers_in_super_kmer; ++i) {
                    auto kmer = it.get();
                    if (build_config.canonical) { /* take the canonical kmer */
                        auto kmer_rc = kmer;
                        kmer_rc.reverse_complement_inplace(build_config.k);
                        kmer = std::min(kmer, kmer_rc);
                    }
                    keys_in_partition.push_back(kmer);
                    super_kmer_ids_in_partition.push_back(super_kmer_id);
                    it.next();
                }
                assert(super_kmer_id < (1ULL << cvb_positions.width()));
                ++super_kmer_id;
            }
        }
        assert(partition_id == num_partitions - 1);
    }

    std::cout << "num_bits_for_skew_index " << m_skew_index.num_bits() << "("
              << static_cast<double>(m_skew_index.num_bits()) / buckets_stats.num_kmers()
              << " [bits/kmer])" << std::endl;
}

}  // namespace sshash