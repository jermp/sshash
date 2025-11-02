#pragma once

#include "essentials.hpp"
#include "include/dictionary.hpp"
#include "include/offsets.hpp"
#include "include/builder/util.hpp"
#include "include/buckets_statistics.hpp"

#include <numeric>  // for std::accumulate

namespace sshash {

template <typename Kmer, typename Offsets>
struct dictionary_builder  //
{
    dictionary_builder(build_configuration const& build_config)
        : build_config(build_config), num_kmers(0), minimizers(build_config) {}

    void build(dictionary<Kmer, Offsets>& d, std::string const& filename)  //
    {
        d.m_k = build_config.k;
        d.m_m = build_config.m;
        d.m_spss.k = build_config.k;
        d.m_spss.m = build_config.m;
        d.m_canonical = build_config.canonical;
        d.m_hasher.seed(build_config.seed);

        std::vector<double> timings;
        timings.reserve(7);

        essentials::timer_type timer;
        // TODO: json_lines

        /*
            step 1: encode strings
        */
        timer.start();
        encode_strings(filename);
        d.m_num_kmers = num_kmers;
        assert(strings_offsets_builder.size() >= 2);
        d.m_num_strings = strings_offsets_builder.size() - 1;
        timer.stop();
        timings.push_back(timer.elapsed());
        if (build_config.verbose) {
            print_time(timings.back(), num_kmers, "step 1: 'encode strings'");
        }
        timer.reset();

        /*
            step 1.1: build weights (if dictionary is weighted)
        */
        if (build_config.weighted) {
            timer.start();
            weights_builder.build(d.m_weights);
            timer.stop();
            if (build_config.verbose) {
                print_time(timings.back(), num_kmers, "step 1.1: 'build weights'");
            }
            timer.reset();
        }

        /*
            step 2: compute minimizer tuples
        */
        timer.start();
        compute_minimizer_tuples();
        timer.stop();
        timings.push_back(timer.elapsed());
        if (build_config.verbose) {
            print_time(timings.back(), num_kmers, "step 2: 'compute minimizer tuples'");
        }
        timer.reset();

        /*
            step 3: merge minimizer tuples
        */
        timer.start();
        minimizers.merge();
        timer.stop();
        timings.push_back(timer.elapsed());
        if (build_config.verbose) {
            print_time(timings.back(), num_kmers, "step 3: 'merging minimizers tuples'");
            std::cout << "num_minimizers = " << minimizers.num_minimizers() << std::endl;
            std::cout << "num_minimizer_positions = " << minimizers.num_minimizer_positions()
                      << std::endl;
            std::cout << "num_super_kmers = " << minimizers.num_super_kmers() << std::endl;
        }
        timer.reset();

        /*
            step 4: build mphf
        */
        timer.start();
        build_mphf(d);
        timer.stop();
        timings.push_back(timer.elapsed());
        if (build_config.verbose) { print_time(timings.back(), num_kmers, "step 4: 'build mphf'"); }
        timer.reset();

        /*
            step 5: hash minimizers
        */
        timer.start();
        hash_minimizers(d);
        timer.stop();
        timings.push_back(timer.elapsed());
        if (build_config.verbose) {
            print_time(timings.back(), num_kmers,
                       "step 5: 'replacing minimizer values with MPHF hashes'");
        }
        timer.reset();

        /*
            step 6: merge minimizer tuples
        */
        timer.start();
        minimizers.merge();
        timer.stop();
        timings.push_back(timer.elapsed());
        if (build_config.verbose) {
            print_time(timings.back(), num_kmers, "step 6: 'merging minimizers tuples '");
        }
        timer.reset();

        /*
            step 7: build sparse and skew index
        */
        timer.start();
        build_sparse_and_skew_index(d);
        assert(strings_offsets_builder.size() == 0);
        timer.stop();
        timings.push_back(timer.elapsed());
        if (build_config.verbose) {
            print_time(timings.back(), num_kmers, "step 7: 'build sparse and skew index'");
        }
        timer.reset();

        double total_time = std::accumulate(timings.begin(), timings.end(), 0.0);
        if (build_config.verbose) {
            print_time(total_time, num_kmers, "total_time");
            d.print_space_breakdown();
        }

        minimizers.remove_tmp_file();
    }

    build_configuration build_config;
    uint64_t num_kmers;
    minimizers_tuples minimizers;
    typename Offsets::builder strings_offsets_builder;
    bits::bit_vector::builder strings_builder;
    weights::builder weights_builder;

    essentials::json_lines build_stats;

private:
    void encode_strings(std::string const& filename);
    void encode_strings(std::istream& is, const input_file_t fmt);
    void compute_minimizer_tuples();
    void build_sparse_and_skew_index(dictionary<Kmer, Offsets>& d);

    void build_mphf(dictionary<Kmer, Offsets>& d) {
        const uint64_t num_minimizers = minimizers.num_minimizers();
        mm::file_source<minimizer_tuple> input(minimizers.get_minimizers_filename(),
                                               mm::advice::sequential);
        minimizers_tuples_iterator iterator(input.data(), input.data() + input.size());
        d.m_ssi.codewords.build(iterator, num_minimizers, build_config);
        input.close();
        assert(d.m_ssi.codewords.size() == num_minimizers);
    }

    void hash_minimizers(dictionary<Kmer, Offsets>& d) {
        std::string filename = minimizers.get_minimizers_filename();
        std::ifstream input(filename, std::ifstream::binary);

        auto const& f = d.m_ssi.codewords.mphf;
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
                        buffer[i].minimizer = f(buffer[i].minimizer);
                    }
                });
            }
            for (auto& t : threads) {
                if (t.joinable()) t.join();
            }
            threads.clear();
            minimizers.sort_and_flush(buffer);
        }

        input.close();
    }
};

}  // namespace sshash