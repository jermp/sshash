#pragma once

#include "include/dictionary.hpp"
#include "essentials.hpp"
#include "include/builder/util.hpp"
#include "include/buckets_statistics.hpp"

#include <numeric>  // for std::accumulate

namespace sshash {

template <class kmer_t, class Endpoints>
struct dictionary_builder  //
{
    dictionary_builder(build_configuration const& build_config)
        : build_config(build_config), num_kmers(0), minimizers(build_config) {}

    void build(dictionary<kmer_t, Endpoints>& d, std::string const& filename)  //
    {
        d.m_k = build_config.k;
        d.m_m = build_config.m;
        d.m_canonical = build_config.canonical;
        d.m_hasher.seed(build_config.seed);

        std::vector<double> timings;
        timings.reserve(6);
        essentials::timer_type timer;

        /*
            step 1: parse the input file, encode sequences (1.1), and compute minimizer tuples (1.2)
        */
        timer.start();
        parse_file(filename);
        d.m_num_kmers = num_kmers;
        assert(strings_endpoints_builder.size() >= 2);
        d.m_num_strings = strings_endpoints_builder.size() - 1;
        if (build_config.weighted) {
            essentials::timer_type timer;
            timer.start();
            weights_builder.build(d.m_weights);
            timer.stop();
            print_time(timer.elapsed(), num_kmers, "step 1.3: 'build weights'");
            if (build_config.verbose) {
                double entropy_weights = weights_builder.print_info(num_kmers);
                double avg_bits_per_weight =
                    static_cast<double>(d.m_weights.num_bits()) / num_kmers;
                std::cout << "weights: " << avg_bits_per_weight << " [bits/kmer]" << std::endl;
                std::cout << "  (" << entropy_weights / avg_bits_per_weight
                          << "x smaller than the empirical entropy)" << std::endl;
            }
        }
        timer.stop();
        timings.push_back(timer.elapsed());
        print_time(timings.back(), num_kmers, "step 1: 'parse file'");
        timer.reset();
        /******/

        /* step 2: merge minimizer tuples and build MPHF ***/
        {
            timer.start();
            minimizers.merge();
            timer.stop();
            timings.push_back(timer.elapsed());
            print_time(timings.back(), num_kmers, "step 2.1: 'merging minimizers tuples'");

            std::cout << "num_minimizers = " << minimizers.num_minimizers() << std::endl;
            std::cout << "num_minimizer_positions = " << minimizers.num_minimizer_positions()
                      << std::endl;
            std::cout << "num_super_kmers = " << minimizers.num_super_kmers() << std::endl;

            timer.reset();

            timer.start();
            const uint64_t num_minimizers = minimizers.num_minimizers();
            mm::file_source<minimizer_tuple> input(minimizers.get_minimizers_filename(),
                                                   mm::advice::sequential);
            minimizers_tuples_iterator iterator(input.data(), input.data() + input.size());
            d.m_minimizers.build(iterator, num_minimizers, build_config);
            input.close();
            assert(d.m_minimizers.size() == num_minimizers);
            timer.stop();
            timings.push_back(timer.elapsed());
            print_time(timings.back(), num_kmers, "step 2.2: 'build minimizers mphf'");

            timer.reset();
        }

        {
            if (build_config.verbose) std::cout << "re-sorting minimizer tuples..." << std::endl;

            timer.start();

            std::string filename = minimizers.get_minimizers_filename();
            std::ifstream input(filename, std::ifstream::binary);

            auto const& f = d.m_minimizers;
            const uint64_t num_threads = build_config.num_threads;
            const uint64_t num_files_to_merge = minimizers.num_files_to_merge();

            minimizers.init();

            const uint64_t num_super_kmers = minimizers.num_super_kmers();
            const uint64_t buffer_size = num_files_to_merge == 1
                                             ? num_super_kmers
                                             : ((build_config.ram_limit_in_GiB * essentials::GiB) /
                                                (2 * sizeof(minimizer_tuple)));
            const uint64_t num_blocks = (num_super_kmers + buffer_size - 1) / buffer_size;
            assert(num_super_kmers > (num_blocks - 1) * buffer_size);

            std::vector<std::thread> threads;
            threads.reserve(num_threads);

            std::vector<minimizer_tuple> buffer;
            for (uint64_t i = 0; i != num_blocks; ++i) {
                const uint64_t n = (i == num_blocks - 1)
                                       ? num_super_kmers - (num_blocks - 1) * buffer_size
                                       : buffer_size;
                buffer.resize(n);
                input.read(reinterpret_cast<char*>(buffer.data()),
                           buffer.size() * sizeof(minimizer_tuple));
                const uint64_t chunk_size = (n + num_threads - 1) / num_threads;
                for (uint64_t t = 0; t * chunk_size < n; ++t) {
                    uint64_t begin = t * chunk_size;
                    uint64_t end = std::min(n, begin + chunk_size);
                    threads.emplace_back([begin, end, &buffer, &f]() {
                        for (uint64_t i = begin; i < end; ++i) {
                            buffer[i].minimizer = f.lookup(buffer[i].minimizer);
                        }
                    });
                }
                for (auto& t : threads) {
                    if (t.joinable()) t.join();
                }
                threads.clear();
                minimizers.sort_and_flush(buffer);
            }
            assert(buffer.empty());

            timer.stop();
            timings.push_back(timer.elapsed());
            print_time(timings.back(), num_kmers,
                       "step 2.3: 'replacing minimizer values with MPHF hashes'");
            timer.reset();

            timer.start();
            minimizers.merge();
            input.close();
            timer.stop();
            timings.push_back(timer.elapsed());
            print_time(timings.back(), num_kmers, "step 2.4: 'merging minimizers tuples '");
            timer.reset();
        }
        /******/

        /* step 3: build sparse and skew index ***/
        timer.start();
        auto buckets_stats = build_sparse_and_skew_index(d.m_buckets, d.m_skew_index);
        assert(strings_endpoints_builder.size() == 0);
        timer.stop();
        timings.push_back(timer.elapsed());
        print_time(timings.back(), num_kmers, "step 3: 'build sparse and skew index'");
        timer.reset();
        /******/

        assert(timings.size() == 6);
        double total_time = std::accumulate(timings.begin(), timings.end(), 0.0);
        print_time(total_time, num_kmers, "total_time");

        d.print_space_breakdown();

        if (build_config.verbose) buckets_stats.print_less();

        minimizers.remove_tmp_file();
    }

    build_configuration build_config;
    uint64_t num_kmers;
    minimizers_tuples minimizers;
    typename Endpoints::builder strings_endpoints_builder;
    bits::bit_vector::builder strings_builder;
    weights::builder weights_builder;

private:
    void parse_file(std::string const& filename);

    void parse_file(std::istream& is, const input_file_type fmt);

    buckets_statistics build_sparse_and_skew_index(buckets<kmer_t, Endpoints>& buckets,
                                                   skew_index<kmer_t>& skew_index);
};

}  // namespace sshash