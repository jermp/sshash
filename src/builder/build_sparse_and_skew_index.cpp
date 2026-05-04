#include "include/buckets_statistics.hpp"

namespace sshash {

/*
    A request to extract `num_kmers_in_super_kmer` consecutive k-mers from
    `strings` starting at base position `starting_pos`. After requests are
    externally sorted by `starting_pos`, k-mer extraction reduces to a single
    forward scan over `strings`.
*/
#pragma pack(push, 4)
struct kmer_extraction_request {
    kmer_extraction_request() {}
    kmer_extraction_request(uint64_t starting_pos, uint32_t partition_id,
                            uint32_t pos_in_bucket, uint32_t num_kmers_in_super_kmer)
        : starting_pos(starting_pos)
        , partition_id(partition_id)
        , pos_in_bucket(pos_in_bucket)
        , num_kmers_in_super_kmer(num_kmers_in_super_kmer) {}

    bool operator<(kmer_extraction_request const& o) const {
        return starting_pos < o.starting_pos;
    }
    bool operator>(kmer_extraction_request const& o) const {
        return starting_pos > o.starting_pos;
    }

    static kmer_extraction_request max() {
        return kmer_extraction_request(uint64_t(-1), uint32_t(-1), uint32_t(-1), uint32_t(-1));
    }

    uint64_t starting_pos;
    uint32_t partition_id;
    uint32_t pos_in_bucket;
    uint32_t num_kmers_in_super_kmer;
};
#pragma pack(pop)

template <typename Kmer, typename Offsets>
void dictionary_builder<Kmer, Offsets>::build_sparse_and_skew_index(
    dictionary<Kmer, Offsets>& d)  //
{
    essentials::timer_type timer;
    timer.start();

    const uint64_t num_minimizer_positions = minimizers.num_minimizer_positions();
    const uint64_t num_minimizers = minimizers.num_minimizers();
    const uint64_t min_size = 1ULL << constants::min_l;
    const uint64_t num_bits_per_offset = strings_offsets_builder.num_bits_per_offset();

    mm::file_source<minimizer_tuple> input(minimizers.get_minimizers_filename(),
                                           mm::advice::sequential);

    buckets_statistics buckets_stats(num_minimizers, num_kmers, num_minimizer_positions);

    uint64_t num_buckets_larger_than_1_not_in_skew_index = 0;
    uint64_t num_buckets_in_skew_index = 0;
    uint64_t num_super_kmers_in_buckets_larger_than_1 = 0;
    uint64_t num_minimizer_positions_of_buckets_larger_than_1 = 0;
    uint64_t num_minimizer_positions_of_buckets_in_skew_index = 0;

    /*
        First pass: collect bucket statistics to compute tighter bound.
    */
    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size());  //
         it.has_next(); it.next())                                                  //
    {
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

        for (auto mt : bucket) {
            buckets_stats.add_num_kmers_in_super_kmer(bucket_size, mt.num_kmers_in_super_kmer);
        }
    }

    assert(buckets_stats.num_buckets() == num_minimizers);

    /*
        Calculate bits needed for control codewords encoding.
        Encoding format:
            ((list_id << min_l) | (bucket_size - 2)) << 2 | status_code
        We need: 2 bits (status) + min_l bits (bucket_size) + bits for list_id.
        list_id is bounded by the maximum number of buckets sharing the same size.
    */
    const uint64_t bits_for_list_id =
        std::ceil(std::log2(buckets_stats.max_sparse_buckets_per_size() + 1));
    const uint64_t num_bits_for_control =
        std::max(num_bits_per_offset + 1, 2 + constants::min_l + bits_for_list_id);

    if (build_config.verbose) {
        std::cout << "num_bits_per_offset = " << num_bits_per_offset << std::endl;
        std::cout << "max_list_id = " << buckets_stats.max_sparse_buckets_per_size() << std::endl;
        std::cout << "bits_for_list_id = " << bits_for_list_id << std::endl;
        std::cout << "num_bits_for_control = " << num_bits_for_control << std::endl;
    }

    bits::compact_vector::builder control_codewords_builder;
    control_codewords_builder.resize(num_minimizers, num_bits_for_control);

    strings_offsets_builder.build(d.m_spss.strings_offsets);
    strings_builder.build(d.m_spss.strings);

    /* step 1. build sparse index */
    assert(buckets_stats.num_buckets() == num_minimizers);

    const uint64_t max_bucket_size = buckets_stats.max_bucket_size();
    const uint64_t log2_max_bucket_size = std::ceil(std::log2(max_bucket_size));

    if (build_config.verbose) {
        std::cout << "num_buckets_larger_than_1_not_in_skew_index "
                  << num_buckets_larger_than_1_not_in_skew_index << "/"
                  << buckets_stats.num_buckets() << " ("
                  << (num_buckets_larger_than_1_not_in_skew_index * 100.0) /
                         buckets_stats.num_buckets()
                  << "%)" << std::endl;
        std::cout << "num_buckets_in_skew_index " << num_buckets_in_skew_index << "/"
                  << buckets_stats.num_buckets() << " ("
                  << (num_buckets_in_skew_index * 100.0) / buckets_stats.num_buckets() << "%)"
                  << std::endl;
        std::cout << "max_bucket_size " << max_bucket_size << std::endl;
        std::cout << "log2_max_bucket_size " << log2_max_bucket_size << std::endl;
    }

    std::vector<bucket_type> buckets;
    buckets.reserve(num_buckets_larger_than_1_not_in_skew_index + num_buckets_in_skew_index);
    std::vector<minimizer_tuple> tuples;  // backed memory
    tuples.reserve(num_super_kmers_in_buckets_larger_than_1);

    // Second pass: collect buckets > 1 for sorting AND handle size-1 buckets
    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size());  //
         it.has_next(); it.next())                                                  //
    {
        const uint64_t bucket_id = it.minimizer();
        auto bucket = it.bucket();
        const uint64_t bucket_size = bucket.size();
        if (bucket_size == 1) {
            // Handle size-1 buckets: encode directly into control codewords
            uint64_t prev_pos_in_seq = constants::invalid_uint64;
            for (auto mt : bucket) {
                if (mt.pos_in_seq != prev_pos_in_seq) {
                    /*
                        For minimizers occurring once, store a (log(N)+1)-bit
                        code, as follows: |offset|0|, i.e., the LSB is 0.
                    */
                    uint64_t code = mt.pos_in_seq << 1;  // first LS bit encodes status code: 0
                    assert(code < (uint64_t(1) << num_bits_for_control));
                    control_codewords_builder.set(bucket_id, code);
                    prev_pos_in_seq = mt.pos_in_seq;
                }
            }
        } else {
            // Collect buckets > 1 for later processing
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
    assert(num_partitions <= 8);  // so that we need 3 bits to encode a partition_id

    if (build_config.verbose) {
        std::cout << "num_partitions in skew index " << num_partitions << std::endl;
        std::cout << "num_minimizer_positions_of_buckets_larger_than_1 "
                  << num_minimizer_positions_of_buckets_larger_than_1 << "/"
                  << num_minimizer_positions << " ("
                  << (num_minimizer_positions_of_buckets_larger_than_1 * 100.0) /
                         num_minimizer_positions
                  << "%)" << std::endl;
        std::cout << "num_minimizer_positions_of_buckets_in_skew_index "
                  << num_minimizer_positions_of_buckets_in_skew_index << "/"
                  << num_minimizer_positions << " ("
                  << (num_minimizer_positions_of_buckets_in_skew_index * 100.0) /
                         num_minimizer_positions
                  << "%)" << std::endl;
    }

    {
        bits::compact_vector::builder mid_load_buckets_builder;
        bits::compact_vector::builder heavy_load_buckets_builder;
        mid_load_buckets_builder.resize(num_minimizer_positions_of_buckets_larger_than_1,
                                        num_bits_per_offset);
        heavy_load_buckets_builder.resize(num_minimizer_positions_of_buckets_in_skew_index,
                                          num_bits_per_offset);

        std::vector<uint32_t> begin_buckets_of_size;
        begin_buckets_of_size.resize(min_size + 1, 0);

        uint64_t curr_bucket_size = 2;
        uint64_t list_id = 0;
        uint64_t mid_load_buckets_size = 0;
        uint64_t heavy_load_buckets_size = 0;

        uint64_t partition_id = 0;
        uint64_t lower = min_size;
        uint64_t upper = 2 * lower;

        for (auto bucket : buckets) {
            const uint64_t bucket_size = bucket.size();
            assert(bucket_size >= 2);

            if (bucket_size > curr_bucket_size) {
                while (bucket_size > curr_bucket_size) ++curr_bucket_size;
                if (curr_bucket_size <= min_size) {
                    begin_buckets_of_size[curr_bucket_size] = mid_load_buckets_size;
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
                        assert(code < (uint64_t(1) << num_bits_for_control));
                        control_codewords_builder.set(mt.minimizer, code);
                    }
                    if (mt.pos_in_seq != prev_pos_in_seq) {
                        mid_load_buckets_builder.push_back(mt.pos_in_seq);
                        prev_pos_in_seq = mt.pos_in_seq;
                        mid_load_buckets_size += 1;
                    }
                }
                ++list_id;
            } else {
                uint64_t prev_pos_in_seq = constants::invalid_uint64;
                for (auto mt : bucket) {
                    if (prev_pos_in_seq == constants::invalid_uint64) {  // only once
                        assert(partition_id < 8);
                        uint64_t p = (heavy_load_buckets_size << 3) | partition_id;
                        uint64_t code = (p << 2) | 3;  // first two LS bits encode status code: 11
                        assert(code < (uint64_t(1) << num_bits_for_control));
                        control_codewords_builder.set(mt.minimizer, code);
                    }
                    if (mt.pos_in_seq != prev_pos_in_seq) {
                        heavy_load_buckets_builder.push_back(mt.pos_in_seq);
                        prev_pos_in_seq = mt.pos_in_seq;
                        heavy_load_buckets_size += 1;
                    }
                }
            }
        }

        d.m_ssi.begin_buckets_of_size = std::move(begin_buckets_of_size);

        control_codewords_builder.build(d.m_ssi.codewords.control_codewords);
        mid_load_buckets_builder.build(d.m_ssi.mid_load_buckets);
        heavy_load_buckets_builder.build(d.m_ssi.ski.heavy_load_buckets);
    }

    timer.stop();

    build_stats.add("step 7.1 (build sparse index)", uint64_t(timer.elapsed()));

    if (build_config.verbose) {
        print_time(uint64_t(timer.elapsed()), "step 7.1 (build sparse index)");
    }

    timer.reset();

    if (num_buckets_in_skew_index == 0) {
        if (build_config.verbose) buckets_stats.print_less();
        return;
    }

    /*
        step 2. build skew index

        We do this in three sub-steps:
        (A) walk the heavy-load buckets in size-sorted order, decode each
            super-kmer's absolute starting position in `strings` and emit a
            `kmer_extraction_request`. Requests are sort+flushed to disk in
            chunks (external sort by `starting_pos`).
        (B) merge the sorted runs and walk `strings` in a single forward
            sequential pass, extracting the requested k-mers. For each k-mer
            we append `(kmer.bits, pos_in_bucket)` to a per-partition tmp file.
        (C) for each partition, read its tmp file, build the MPHF, then build
            the positions compact vector. The skew index is assembled
            partition by partition.

        Avoiding the random access pattern over `strings` in (B) is the
        precondition for moving `strings` itself out of RAM in a later step.
    */
    timer.start();

    std::vector<uint64_t> num_kmers_in_partition(num_partitions, 0);

    /* unique run identifier for the tmp files produced by this step */
    const uint64_t skew_run_id = pthash::clock_type::now().time_since_epoch().count();
    auto request_run_filename = [&](uint64_t id) {
        std::stringstream ss;
        ss << build_config.tmp_dirname << "/sshash.tmp.run_" << skew_run_id
           << ".kmer_requests." << id << ".bin";
        return ss.str();
    };
    auto skew_partition_filename = [&](uint64_t pid) {
        std::stringstream ss;
        ss << build_config.tmp_dirname << "/sshash.tmp.run_" << skew_run_id
           << ".skew_kmers." << pid << ".bin";
        return ss.str();
    };

    /* (A) emit kmer-extraction requests, externally sorted by `starting_pos` */
    std::atomic<uint64_t> num_request_runs{0};
    {
        const uint64_t request_buffer_capacity = std::max<uint64_t>(
            uint64_t(1) << 16,
            (build_config.ram_limit_in_GiB * essentials::GiB) /
                (4 * sizeof(kmer_extraction_request)));

        std::vector<kmer_extraction_request> request_buffer;
        request_buffer.reserve(request_buffer_capacity);

        auto flush_request_buffer = [&]() {
            if (request_buffer.empty()) return;
            parallel_sort(request_buffer, build_config.num_threads,
                          [](kmer_extraction_request const& a,
                             kmer_extraction_request const& b) {
                              return a.starting_pos < b.starting_pos;
                          });
            const uint64_t id = num_request_runs.fetch_add(1);
            const std::string fn = request_run_filename(id);
            if (build_config.verbose) {
                std::cout << "saving to file '" << fn << "'..." << std::endl;
            }
            std::ofstream out(fn, std::ofstream::binary);
            if (!out.is_open()) throw std::runtime_error("cannot open file");
            out.write(reinterpret_cast<char const*>(request_buffer.data()),
                      request_buffer.size() * sizeof(kmer_extraction_request));
            out.close();
            request_buffer.clear();
        };

        uint64_t partition_id = 0;
        uint64_t lower = min_size;
        uint64_t upper = 2 * lower;

        for (uint64_t i = buckets.size() - num_buckets_in_skew_index; i < buckets.size(); ++i)  //
        {
            auto const& bucket = buckets[i];
            const uint64_t bucket_size = bucket.size();
            while (bucket_size > upper)  //
            {
                lower = upper;
                upper = 2 * lower;
                partition_id += 1;
                if (partition_id == num_partitions - 1) upper = max_bucket_size;
            }
            assert(bucket_size > lower and bucket_size <= upper);
            assert(partition_id < num_partitions);

            uint32_t pos_in_bucket = uint32_t(-1);
            uint64_t prev_pos_in_seq = constants::invalid_uint64;
            for (auto mt : bucket)  //
            {
                num_kmers_in_partition[partition_id] += mt.num_kmers_in_super_kmer;
                if (mt.pos_in_seq != prev_pos_in_seq) {
                    prev_pos_in_seq = mt.pos_in_seq;
                    ++pos_in_bucket;
                }
                assert(mt.pos_in_seq >= mt.pos_in_kmer);
                const uint64_t abs_offset =
                    d.m_spss.strings_offsets.decode(mt.pos_in_seq).absolute_offset;
                const uint64_t starting_pos = abs_offset - mt.pos_in_kmer;
                if (request_buffer.size() == request_buffer_capacity) flush_request_buffer();
                request_buffer.emplace_back(starting_pos,                          //
                                            uint32_t(partition_id),                //
                                            pos_in_bucket,                         //
                                            uint32_t(mt.num_kmers_in_super_kmer)); //
            }
        }
        flush_request_buffer();
        assert(partition_id == num_partitions - 1);
    }

    if (build_config.verbose) {
        uint64_t total_kmers_in_skew = 0;
        for (uint64_t p = 0; p != num_partitions; ++p) {
            total_kmers_in_skew += num_kmers_in_partition[p];
            std::cout << "  partition = " << p
                      << ": num kmers in partition = " << num_kmers_in_partition[p] << std::endl;
        }
        std::cout << "num kmers in skew index = " << total_kmers_in_skew << " ("
                  << (total_kmers_in_skew * 100.0) / buckets_stats.num_kmers() << "%)" << std::endl;
    }

    /* (B) sequential extraction over `strings` -> per-partition kmer tmp files */
    {
        struct request_run_names_iterator {
            request_run_names_iterator(std::string const& tmp_dirname, uint64_t skew_run_id)
                : i(0), skew_run_id(skew_run_id), tmp_dirname(tmp_dirname) {}

            std::string operator*() const {
                std::stringstream ss;
                ss << tmp_dirname << "/sshash.tmp.run_" << skew_run_id
                   << ".kmer_requests." << i << ".bin";
                return ss.str();
            }
            void operator++() { ++i; }

            uint64_t i;
            uint64_t skew_run_id;
            std::string tmp_dirname;
        };

        request_run_names_iterator names_it(build_config.tmp_dirname, skew_run_id);
        file_merging_iterator<kmer_extraction_request> merger(names_it, num_request_runs.load());

        std::vector<std::ofstream> partition_writers(num_partitions);
        for (uint64_t p = 0; p != num_partitions; ++p) {
            if (num_kmers_in_partition[p] == 0) continue;
            partition_writers[p].open(skew_partition_filename(p),
                                      std::ofstream::binary | std::ofstream::trunc);
            if (!partition_writers[p].is_open()) {
                throw std::runtime_error("cannot open skew-partition tmp file");
            }
        }

        const uint64_t k = build_config.k;
        const bool canonical = build_config.canonical;
        kmer_iterator<Kmer, bits::bit_vector> kmer_it(d.m_spss.strings, k);

        while (merger.has_next())  //
        {
            const kmer_extraction_request req = *merger;
            kmer_it.at(Kmer::bits_per_char * req.starting_pos);
            for (uint32_t i = 0; i != req.num_kmers_in_super_kmer; ++i) {
                Kmer kmer = kmer_it.get();
                if (canonical) {
                    Kmer kmer_rc = kmer;
                    kmer_rc.reverse_complement_inplace(k);
                    kmer = std::min(kmer, kmer_rc);
                }
                auto& w = partition_writers[req.partition_id];
                /* write only `kmer.bits` (avoids serializing the vptr that
                   `uint_kmer_t` carries due to its virtual destructor) */
                w.write(reinterpret_cast<char const*>(&kmer.bits), sizeof(kmer.bits));
                w.write(reinterpret_cast<char const*>(&req.pos_in_bucket),
                        sizeof(req.pos_in_bucket));
                kmer_it.next();
            }
            merger.next();
        }
        merger.close();

        for (auto& w : partition_writers) {
            if (w.is_open()) w.close();
        }

        for (uint64_t i = 0; i != num_request_runs.load(); ++i) {
            std::remove(request_run_filename(i).c_str());
        }
    }

    /* (C) per-partition MPHF + positions build */
    {
        std::vector<kmers_pthash_type<Kmer>> mphfs;
        std::vector<bits::compact_vector> positions;
        mphfs.resize(num_partitions);
        positions.resize(num_partitions);

        pthash::build_configuration mphf_build_config;
        mphf_build_config.lambda =
            build_config.lambda + 2.0; /* Use higher lambda here since we have less keys. */
        mphf_build_config.alpha = 0.94;
        mphf_build_config.seed = util::get_seed_for_hash_function(build_config);
        mphf_build_config.verbose = false;
        mphf_build_config.num_threads = build_config.num_threads;
        mphf_build_config.avg_partition_size = constants::avg_partition_size;

        uint64_t lower = min_size;
        uint64_t upper = 2 * lower;
        uint64_t num_bits_per_pos = constants::min_l + 1;
        if (num_partitions == 1) {
            upper = max_bucket_size;
            num_bits_per_pos = log2_max_bucket_size;
        }

        for (uint64_t partition_id = 0; partition_id != num_partitions; ++partition_id)  //
        {
            const uint64_t n = num_kmers_in_partition[partition_id];

            if (build_config.verbose) {
                std::cout << "  lower = " << lower << "; upper = " << upper
                          << "; num_bits_per_pos = " << num_bits_per_pos
                          << "; num_kmers_in_partition = " << n << std::endl;
            }

            if (n > 0)  //
            {
                std::vector<Kmer> kmers;
                std::vector<uint32_t> positions_in_bucket;
                kmers.reserve(n);
                positions_in_bucket.reserve(n);

                {
                    const std::string fn = skew_partition_filename(partition_id);
                    std::ifstream in(fn, std::ifstream::binary);
                    if (!in.is_open()) {
                        throw std::runtime_error("cannot open skew-partition tmp file");
                    }
                    for (uint64_t i = 0; i != n; ++i) {
                        Kmer kmer;
                        in.read(reinterpret_cast<char*>(&kmer.bits), sizeof(kmer.bits));
                        uint32_t pib;
                        in.read(reinterpret_cast<char*>(&pib), sizeof(pib));
                        kmers.push_back(kmer);
                        positions_in_bucket.push_back(pib);
                    }
                    in.close();
                    std::remove(fn.c_str());
                }

                bits::compact_vector::builder cvb_positions;
                cvb_positions.resize(n, num_bits_per_pos);

                if (build_config.verbose) {
                    const uint64_t avg_partition_size =
                        pthash::compute_avg_partition_size(kmers.size(), mphf_build_config);
                    const uint64_t pthash_num_partitions =
                        pthash::compute_num_partitions(kmers.size(), avg_partition_size);
                    assert(pthash_num_partitions > 0);
                    std::cout << "    building MPHF with " << mphf_build_config.num_threads
                              << " threads and " << pthash_num_partitions
                              << " partitions (avg. partition size = " << avg_partition_size
                              << ")..." << std::endl;
                }

                auto& F = mphfs[partition_id];
                F.build_in_internal_memory(kmers.begin(), kmers.size(), mphf_build_config);

                if (build_config.verbose) {
                    std::cout << "    built mphs[" << partition_id << "] for " << kmers.size()
                              << " kmers; bits/key = "
                              << static_cast<double>(F.num_bits()) / F.num_keys() << std::endl;
                }

                for (uint64_t i = 0; i != kmers.size(); ++i) {
                    Kmer kmer = kmers[i];
                    uint64_t pos = F(kmer);
                    uint32_t pos_in_bucket = positions_in_bucket[i];
                    cvb_positions.set(pos, pos_in_bucket);
                }
                auto& P = positions[partition_id];
                cvb_positions.build(P);

                if (build_config.verbose) {
                    std::cout << "    built positions[" << partition_id << "] for " << P.size()
                              << " kmers; bits/key = " << (P.num_bytes() * 8.0) / P.size()
                              << std::endl;
                }
            }

            /* advance partition state for the next iteration */
            lower = upper;
            upper = 2 * lower;
            num_bits_per_pos += 1;
            if (partition_id + 1 == num_partitions - 1) {
                upper = max_bucket_size;
                num_bits_per_pos = log2_max_bucket_size;
            }
        }

        d.m_ssi.ski.mphfs = std::move(mphfs);
        d.m_ssi.ski.positions = std::move(positions);
    }

    timer.stop();

    build_stats.add("step 7.2 (build skew index)", uint64_t(timer.elapsed()));

    if (build_config.verbose) {
        print_time(uint64_t(timer.elapsed()), "step 7.2 (build skew index)");
        buckets_stats.print_less();
    }
}

}  // namespace sshash
