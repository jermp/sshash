#include "include/dictionary.hpp"
#include "essentials.hpp"
#include "include/builder/util.hpp"

#include "include/builder/parse_file.hpp"
#include "include/builder/build_sparse_index.hpp"
#include "include/builder/build_skew_index.hpp"

#include <numeric>  // for std::accumulate

namespace sshash {

template <class kmer_t>
void dictionary<kmer_t>::build(std::string const& filename,
                               build_configuration const& build_config) {
    /* Validate the build configuration. */
    if (build_config.k == 0) throw std::runtime_error("k must be > 0");
    if (build_config.k > kmer_t::max_k) {
        throw std::runtime_error("k must be less <= " + std::to_string(kmer_t::max_k) +
                                 " but got k = " + std::to_string(build_config.k));
    }
    if (build_config.m == 0) throw std::runtime_error("m must be > 0");
    if (build_config.m > kmer_t::max_m) {
        throw std::runtime_error("m must be less <= " + std::to_string(kmer_t::max_m) +
                                 " but got m = " + std::to_string(build_config.m));
    }
    if (build_config.m > build_config.k) throw std::runtime_error("m must be <= k");
    if (build_config.l >= constants::max_l) {
        throw std::runtime_error("l must be < " + std::to_string(constants::max_l));
    }

    m_k = build_config.k;
    m_m = build_config.m;
    m_canonical = build_config.canonical;
    m_skew_index.min_log2 = build_config.l;
    m_hasher.seed(build_config.seed);

    std::vector<double> timings;
    timings.reserve(7);
    essentials::timer_type timer;

    /* step 1: parse the input file, encode sequences (1.1), and compute minimizer tuples (1.2) ***/
    timer.start();
    parse_data<kmer_t> data(build_config);
    parse_file<kmer_t>(filename, data, build_config);
    m_size = data.num_kmers;
    if (build_config.weighted) {
        essentials::timer_type timer;
        timer.start();
        data.weights_builder.build(m_weights);
        timer.stop();
        print_time(timer.elapsed(), data.num_kmers, "step 1.3: 'build_weights'");
        if (build_config.verbose) {
            double entropy_weights = data.weights_builder.print_info(data.num_kmers);
            double avg_bits_per_weight = static_cast<double>(m_weights.num_bits()) / data.num_kmers;
            std::cout << "weights: " << avg_bits_per_weight << " [bits/kmer]" << std::endl;
            std::cout << "  (" << entropy_weights / avg_bits_per_weight
                      << "x smaller than the empirical entropy)" << std::endl;
        }
    }
    timer.stop();
    timings.push_back(timer.elapsed());
    print_time(timings.back(), data.num_kmers, "step 1: 'parse_file'");
    timer.reset();
    /******/

    /* step 2: merge minimizer tuples and build MPHF ***/
    {
        timer.start();
        data.minimizers.merge();
        timer.stop();
        timings.push_back(timer.elapsed());
        print_time(timings.back(), data.num_kmers, "step 2.1: 'merging_minimizers_tuples'");

        timer.reset();

        timer.start();
        const uint64_t num_minimizers = data.minimizers.num_minimizers();
        mm::file_source<minimizer_tuple> input(data.minimizers.get_minimizers_filename(),
                                               mm::advice::sequential);
        minimizers_tuples_iterator iterator(input.data(), input.data() + input.size());
        m_minimizers.build(iterator, num_minimizers, build_config);
        input.close();
        assert(m_minimizers.size() == num_minimizers);
        timer.stop();
        timings.push_back(timer.elapsed());
        print_time(timings.back(), data.num_kmers, "step 2.2: 'build_minimizers_mphf'");

        timer.reset();
    }

    {
        if (build_config.verbose) std::cout << "re-sorting minimizer tuples..." << std::endl;

        timer.start();

        std::string filename = data.minimizers.get_minimizers_filename();
        std::ifstream input(filename, std::ifstream::binary);

        auto const& f = m_minimizers;
        const uint64_t num_threads = build_config.num_threads;
        const uint64_t num_files_to_merge = data.minimizers.num_files_to_merge();

        data.minimizers.init();

        const uint64_t num_super_kmers = data.minimizers.num_super_kmers();
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
            for (uint64_t t = 0; t != num_threads; ++t) {
                uint64_t begin = t * chunk_size;
                uint64_t end = (t == num_threads - 1) ? n : begin + chunk_size;
                threads.emplace_back([begin, end, &buffer, &f]() {
                    for (uint64_t i = begin; i != end; ++i) {
                        buffer[i].minimizer = f.lookup(buffer[i].minimizer);
                    }
                });
            }
            for (auto& t : threads) {
                if (t.joinable()) t.join();
            }
            threads.clear();
            data.minimizers.sort_and_flush(buffer);
        }
        assert(buffer.empty());

        timer.stop();
        timings.push_back(timer.elapsed());
        print_time(timings.back(), data.num_kmers,
                   "step 2.3: 'replacing minimizer values with MPHF hashes'");
        timer.reset();

        timer.start();
        data.minimizers.merge();
        input.close();
        timer.stop();
        timings.push_back(timer.elapsed());
        print_time(timings.back(), data.num_kmers, "step 2.4: 'merging_minimizers_tuples '");
        timer.reset();
    }
    /******/

    /* step 3: build sparse index ***/
    timer.start();
    auto buckets_stats = build_sparse_index(data, m_buckets, build_config);
    timer.stop();
    timings.push_back(timer.elapsed());
    print_time(timings.back(), data.num_kmers, "step 3: 'build_sparse_index'");
    timer.reset();
    /******/

    /* step 4: build skew index ***/
    timer.start();
    build_skew_index(m_skew_index, data, m_buckets, build_config, buckets_stats);
    timer.stop();
    timings.push_back(timer.elapsed());
    print_time(timings.back(), data.num_kmers, "step 4: 'build_skew_index'");
    timer.reset();
    /******/

    assert(timings.size() == 7);
    double total_time = std::accumulate(timings.begin(), timings.end(), 0.0);
    print_time(total_time, data.num_kmers, "total_time");

    print_space_breakdown();

    if (build_config.verbose) buckets_stats.print_less();

    data.minimizers.remove_tmp_file();
}

}  // namespace sshash
