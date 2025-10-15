#include "dictionary_builder.hpp"
#include "include/buckets_statistics.hpp"

namespace sshash {

template <typename Kmer, typename Offsets>
buckets_statistics                                                                            //
dictionary_builder<Kmer, Offsets>::build_sparse_and_skew_index(dictionary<Kmer, Offsets>& d)  //
{
    const uint64_t num_minimizer_positions = minimizers.num_minimizer_positions();
    const uint64_t num_minimizers = minimizers.num_minimizers();
    const uint64_t min_size = 1ULL << constants::min_l;
    const uint64_t num_bits_per_offset = strings_offsets_builder.num_bits_per_offset();

    std::cout << "num_bits_per_offset = " << num_bits_per_offset << std::endl;

    bits::compact_vector::builder control_codewords_builder;
    control_codewords_builder.resize(num_minimizers, num_bits_per_offset + 1);

    mm::file_source<minimizer_tuple> input(minimizers.get_minimizers_filename(),
                                           mm::advice::sequential);

    essentials::timer_type timer;
    timer.start();

    buckets_statistics buckets_stats(num_minimizers, num_kmers, num_minimizer_positions);

    uint64_t num_buckets_larger_than_1_not_in_skew_index = 0;
    uint64_t num_buckets_in_skew_index = 0;
    uint64_t num_super_kmers_in_buckets_larger_than_1 = 0;
    uint64_t num_minimizer_positions_of_buckets_larger_than_1 = 0;
    uint64_t num_minimizer_positions_of_buckets_in_skew_index = 0;

    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size());  //
         it.has_next(); it.next())                                                  //
    {
        const uint64_t bucket_id = it.minimizer();
        assert(bucket_id < num_minimizers);
        auto bucket = it.bucket();
        const uint64_t bucket_size = bucket.size();
        buckets_stats.add_bucket_size(bucket_size);

        if (bucket_size > 1) {
            if (bucket_size <= min_size) {
                ++num_buckets_larger_than_1_not_in_skew_index;
                num_minimizer_positions_of_buckets_larger_than_1 += bucket_size;
            } else {
                ++num_buckets_in_skew_index;
                num_minimizer_positions_of_buckets_in_skew_index += bucket_size;
            }
            num_super_kmers_in_buckets_larger_than_1 += bucket.num_super_kmers();
        }

        uint64_t prev_pos_in_seq = constants::invalid_uint64;
        for (auto mt : bucket) {
            if (bucket_size == 1 and mt.pos_in_seq != prev_pos_in_seq) {
                /*
                    For minimizers occurring once, store a (log(N)+1)-bit
                    code, as follows: |offset|0|, i.e., the LSB is 0.
                */
                uint64_t code = mt.pos_in_seq << 1;  // first LS bit encodes status code: 0
                assert(code < (uint64_t(1) << (num_bits_per_offset + 1)));
                control_codewords_builder.set(bucket_id, code);
                prev_pos_in_seq = mt.pos_in_seq;
            }
            buckets_stats.add_num_kmers_in_super_kmer(bucket_size, mt.num_kmers_in_super_kmer);
        }
    }

    assert(buckets_stats.num_buckets() == num_minimizers);
    std::cout << "num_buckets_larger_than_1_not_in_skew_index "
              << num_buckets_larger_than_1_not_in_skew_index << "/" << buckets_stats.num_buckets()
              << " ("
              << (num_buckets_larger_than_1_not_in_skew_index * 100.0) / buckets_stats.num_buckets()
              << "%)" << std::endl;
    std::cout << "num_buckets_in_skew_index " << num_buckets_in_skew_index << "/"
              << buckets_stats.num_buckets() << " ("
              << (num_buckets_in_skew_index * 100.0) / buckets_stats.num_buckets() << "%)"
              << std::endl;

    strings_offsets_builder.build(d.m_spss.strings_offsets);
    strings_builder.build(d.m_spss.strings);

    /* compute offsets2 and offsets3 */
    assert(buckets_stats.num_buckets() == num_minimizers);

    const uint64_t max_bucket_size = buckets_stats.max_bucket_size();
    const uint64_t log2_max_bucket_size = std::ceil(std::log2(max_bucket_size));
    std::cout << "max_bucket_size " << max_bucket_size << std::endl;
    std::cout << "log2_max_bucket_size " << log2_max_bucket_size << std::endl;

    std::vector<bucket_type> buckets;
    buckets.reserve(num_buckets_larger_than_1_not_in_skew_index + num_buckets_in_skew_index);
    std::vector<minimizer_tuple> tuples;  // backed memory
    tuples.reserve(num_super_kmers_in_buckets_larger_than_1);

    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size());  //
         it.has_next(); it.next())                                                  //
    {
        auto bucket = it.bucket();
        if (bucket.size() > 1) {
            minimizer_tuple const* begin = tuples.data() + tuples.size();
            std::copy(bucket.begin_ptr(), bucket.end_ptr(), std::back_inserter(tuples));
            minimizer_tuple const* end = tuples.data() + tuples.size();
            buckets.push_back(bucket_type(begin, end));
        }
    }
    assert(buckets.size() ==
           num_buckets_larger_than_1_not_in_skew_index + num_buckets_in_skew_index);

    input.close();

    std::sort(buckets.begin(), buckets.end(),
              [](bucket_type const& x, bucket_type const& y) { return x.size() < y.size(); });

    uint64_t num_partitions = constants::max_l - constants::min_l + 1;
    if (max_bucket_size < min_size) {
        num_partitions = 0;
    } else if (max_bucket_size < (1ULL << constants::max_l)) {
        num_partitions = log2_max_bucket_size - constants::min_l;
    }
    std::cout << "skew index num_partitions " << num_partitions << std::endl;
    assert(num_partitions <= 8);  // so that we need 3 bits to encode a partition_id

    std::cout << "num_minimizer_positions_of_buckets_larger_than_1 "
              << num_minimizer_positions_of_buckets_larger_than_1 << "/" << num_minimizer_positions
              << " ("
              << (num_minimizer_positions_of_buckets_larger_than_1 * 100.0) /
                     num_minimizer_positions
              << "%)" << std::endl;
    std::cout << "num_minimizer_positions_of_buckets_in_skew_index "
              << num_minimizer_positions_of_buckets_in_skew_index << "/" << num_minimizer_positions
              << " ("
              << (num_minimizer_positions_of_buckets_in_skew_index * 100.0) /
                     num_minimizer_positions
              << "%)" << std::endl;

    bits::compact_vector::builder mid_load_buckets_builder;
    bits::compact_vector::builder heavy_load_buckets_builder;
    mid_load_buckets_builder.resize(num_minimizer_positions_of_buckets_larger_than_1,
                                    num_bits_per_offset);
    heavy_load_buckets_builder.resize(num_minimizer_positions_of_buckets_in_skew_index,
                                      num_bits_per_offset);

    uint64_t curr_bucket_size = 2;
    uint64_t list_id = 0;

    d.m_ssi.begin_buckets_of_size.resize(min_size + 1, 0);
    uint64_t offsets2_curr_size = 0;
    uint64_t offsets3_curr_size = 0;

    uint64_t partition_id = 0;
    uint64_t lower = min_size;
    uint64_t upper = 2 * lower;

    for (auto bucket : buckets) {
        const uint64_t bucket_size = bucket.size();
        assert(bucket_size >= 2);

        if (bucket_size > curr_bucket_size) {
            while (bucket_size > curr_bucket_size) ++curr_bucket_size;
            if (curr_bucket_size <= min_size) {
                d.m_ssi.begin_buckets_of_size[curr_bucket_size] = offsets2_curr_size;
            } else {
                while (curr_bucket_size > upper) {
                    lower = upper;
                    upper = 2 * lower;
                    partition_id += 1;
                    if (partition_id == num_partitions - 1) upper = max_bucket_size;
                }
            }
            list_id = 0;
        }

        if (curr_bucket_size <= min_size) {
            uint64_t prev_pos_in_seq = constants::invalid_uint64;
            for (auto mt : bucket) {
                if (prev_pos_in_seq == constants::invalid_uint64) {  // only once
                    uint64_t p = (list_id << constants::min_l) | (curr_bucket_size - 2);
                    uint64_t code = (p << 2) | 1;  // first two LS bits encode status code: 01
                    assert(code < (uint64_t(1) << (num_bits_per_offset + 1)));
                    control_codewords_builder.set(mt.minimizer, code);
                }
                if (mt.pos_in_seq != prev_pos_in_seq) {
                    mid_load_buckets_builder.push_back(mt.pos_in_seq);
                    prev_pos_in_seq = mt.pos_in_seq;
                    offsets2_curr_size += 1;
                }
            }
            ++list_id;
        } else {
            uint64_t prev_pos_in_seq = constants::invalid_uint64;
            for (auto mt : bucket) {
                if (prev_pos_in_seq == constants::invalid_uint64) {  // only once
                    assert(partition_id < 8);
                    uint64_t p = (offsets3_curr_size << 3) | partition_id;
                    uint64_t code = (p << 2) | 3;  // first two LS bits encode status code: 11
                    assert(code < (uint64_t(1) << (num_bits_per_offset + 1)));
                    control_codewords_builder.set(mt.minimizer, code);
                }
                if (mt.pos_in_seq != prev_pos_in_seq) {
                    heavy_load_buckets_builder.push_back(mt.pos_in_seq);
                    prev_pos_in_seq = mt.pos_in_seq;
                    offsets3_curr_size += 1;
                }
            }
        }
    }

    control_codewords_builder.build(d.m_ssi.codewords.control_codewords);
    mid_load_buckets_builder.build(d.m_ssi.mid_load_buckets);
    heavy_load_buckets_builder.build(d.m_ssi.m_skew_index.heavy_load_buckets);

    timer.stop();
    std::cout << "computing minimizers offsets: " << timer.elapsed() / 1000000 << " [sec]"
              << std::endl;

    timer.reset();

    if (num_buckets_in_skew_index == 0) return buckets_stats;

    /* build skew index */
    timer.start();
    std::vector<uint64_t> num_kmers_in_partition(num_partitions, 0);
    d.m_ssi.m_skew_index.mphfs.resize(num_partitions);
    d.m_ssi.m_skew_index.positions.resize(num_partitions);

    {
        std::cout << "computing sizes of partitions..." << std::endl;

        uint64_t partition_id = 0;
        uint64_t lower = min_size;
        uint64_t upper = 2 * lower;
        uint64_t num_kmers_in_skew_index = 0;
        for (uint64_t i = buckets.size() - num_buckets_in_skew_index; i <= buckets.size(); ++i)  //
        {
            while (i == buckets.size() or buckets[i].size() > upper) {
                std::cout << "  partition_id = " << partition_id
                          << ": num_kmers belonging to buckets of size > " << lower
                          << " and <= " << upper << ": " << num_kmers_in_partition[partition_id]
                          << std::endl;
                num_kmers_in_skew_index += num_kmers_in_partition[partition_id];

                if (i == buckets.size()) break;

                lower = upper;
                upper = 2 * lower;
                partition_id += 1;
                if (partition_id == num_partitions - 1) upper = max_bucket_size;
            }

            if (i == buckets.size()) break;

            assert(buckets[i].size() > lower and buckets[i].size() <= upper);
            for (auto mt : buckets[i]) {
                num_kmers_in_partition[partition_id] += mt.num_kmers_in_super_kmer;
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
        mphf_build_config.lambda =
            build_config.lambda + 2.0; /* Use higher lambda here since we have less keys. */
        mphf_build_config.alpha = 0.94;
        mphf_build_config.seed = util::get_seed_for_hash_function(build_config);
        mphf_build_config.verbose = false;
        mphf_build_config.num_threads = build_config.num_threads;
        mphf_build_config.avg_partition_size = constants::avg_partition_size;

        uint64_t partition_id = 0;
        uint64_t lower = min_size;
        uint64_t upper = 2 * lower;
        uint64_t num_bits_per_pos = constants::min_l + 1;

        /* Temporary storage for kmers and positions within a partition. */
        std::vector<Kmer> kmers;
        std::vector<uint32_t> positions_in_bucket;
        bits::compact_vector::builder cvb_positions;
        kmers.reserve(num_kmers_in_partition[partition_id]);
        positions_in_bucket.reserve(num_kmers_in_partition[partition_id]);
        cvb_positions.resize(num_kmers_in_partition[partition_id], num_bits_per_pos);

        const uint64_t k = build_config.k;

        for (uint64_t i = buckets.size() - num_buckets_in_skew_index; i <= buckets.size(); ++i)  //
        {
            while (i == buckets.size() or buckets[i].size() > upper) {
                std::cout << "  lower = " << lower << "; upper = " << upper
                          << "; num_bits_per_pos = " << num_bits_per_pos
                          << "; num_kmers_in_partition = " << kmers.size() << std::endl;
                assert(num_kmers_in_partition[partition_id] == kmers.size());
                assert(num_kmers_in_partition[partition_id] == positions_in_bucket.size());

                if (num_kmers_in_partition[partition_id] > 0)  //
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

                    auto& mphf = d.m_ssi.m_skew_index.mphfs[partition_id];
                    mphf.build_in_internal_memory(kmers.begin(), kmers.size(), mphf_build_config);

                    std::cout << "    built mphs[" << partition_id << "] for " << kmers.size()
                              << " kmers; bits/key = "
                              << static_cast<double>(mphf.num_bits()) / mphf.num_keys()
                              << std::endl;

                    for (uint64_t i = 0; i != kmers.size(); ++i) {
                        Kmer kmer = kmers[i];
                        uint64_t pos = mphf(kmer);
                        uint32_t pos_in_bucket = positions_in_bucket[i];
                        cvb_positions.set(pos, pos_in_bucket);
                    }
                    auto& positions = d.m_ssi.m_skew_index.positions[partition_id];
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
                    upper = max_bucket_size;
                    num_bits_per_pos = log2_max_bucket_size;
                }

                kmers.clear();
                positions_in_bucket.clear();
                kmers.reserve(num_kmers_in_partition[partition_id]);
                positions_in_bucket.reserve(num_kmers_in_partition[partition_id]);
                cvb_positions.resize(num_kmers_in_partition[partition_id], num_bits_per_pos);
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

                mt.pos_in_seq = d.m_spss.strings_offsets.decode(mt.pos_in_seq).absolute_offset;

                const uint64_t starting_pos_of_super_kmer = mt.pos_in_seq - mt.pos_in_kmer;
                kmer_iterator<Kmer, bits::bit_vector> it(
                    d.m_spss.strings, k, Kmer::bits_per_char * starting_pos_of_super_kmer);
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

    timer.stop();
    std::cout << "computing skew index took: " << timer.elapsed() / 1000000 << " [sec]"
              << std::endl;

    return buckets_stats;
}

}  // namespace sshash