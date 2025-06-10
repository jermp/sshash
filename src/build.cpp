#include "include/dictionary.hpp"
#include "essentials.hpp"
#include "include/builder/util.hpp"

#include "include/builder/parse_file.hpp"
#include "include/builder/build_index.hpp"
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
    m_seed = build_config.seed;
    m_canonical = build_config.canonical;
    m_skew_index.min_log2 = build_config.l;

    std::vector<double> timings;
    timings.reserve(6);
    essentials::timer_type timer;

    /* step 1: parse the input file and build compact string pool ***/
    timer.start();
    parse_data<kmer_t> data = parse_file<kmer_t>(filename, build_config);
    m_size = data.num_kmers;
    timer.stop();
    timings.push_back(timer.elapsed());
    print_time(timings.back(), data.num_kmers, "step 1: 'parse_file'");
    timer.reset();
    /******/

    if (build_config.weighted) {
        /* step 1.1: compress weights ***/
        timer.start();
        data.weights_builder.build(m_weights);
        timer.stop();
        timings.push_back(timer.elapsed());
        print_time(timings.back(), data.num_kmers, "step 1.1.: 'build_weights'");
        timer.reset();
        /******/
        if (build_config.verbose) {
            double entropy_weights = data.weights_builder.print_info(data.num_kmers);
            double avg_bits_per_weight = static_cast<double>(m_weights.num_bits()) / data.num_kmers;
            std::cout << "weights: " << avg_bits_per_weight << " [bits/kmer]" << std::endl;
            std::cout << "  (" << entropy_weights / avg_bits_per_weight
                      << "x smaller than the empirical entropy)" << std::endl;
        }
    }

    /* step 2: merge minimizers and build MPHF ***/
    timer.start();
    data.minimizers.merge();
    const uint64_t num_minimizers = data.minimizers.num_minimizers();
    const uint64_t num_super_kmers = data.strings.num_super_kmers();
    {
        mm::file_source<minimizer_tuple> input(data.minimizers.get_minimizers_filename(),
                                               mm::advice::sequential);
        minimizers_tuples_iterator iterator(input.data(), input.data() + input.size());
        assert(input.size() == num_super_kmers);
        if (build_config.verbose) std::cout << "num_minimizers " << num_minimizers << std::endl;
        m_minimizers.build(iterator, num_minimizers, build_config);
        input.close();
        assert(m_minimizers.size() == num_minimizers);
    }
    timer.stop();
    timings.push_back(timer.elapsed());
    print_time(timings.back(), data.num_kmers, "step 2: 'build_minimizers'");
    timer.reset();
    /******/

    /* step 2.1: resort minimizers by MPHF order ***/
    timer.start();
    {
        if (build_config.verbose) std::cout << "re-sorting minimizer tuples..." << std::endl;
        std::ifstream input(data.minimizers.get_minimizers_filename(), std::ifstream::binary);

        auto const& f = m_minimizers;
        const uint64_t num_threads = build_config.num_threads;
        const uint64_t num_files_to_merge = data.minimizers.num_files_to_merge();

        data.minimizers.init();

        const uint64_t buffer_size =
            num_files_to_merge == 1 ? num_super_kmers : data.minimizers.buffer_size();
        const uint64_t num_blocks = (num_super_kmers + buffer_size - 1) / buffer_size;
        assert(num_super_kmers > (num_blocks - 1) * buffer_size);

        std::vector<std::thread> threads;
        threads.reserve(num_threads);

        auto& buff = data.minimizers.buffer();

        for (uint64_t i = 0; i != num_blocks; ++i) {
            const uint64_t n = (i == num_blocks - 1)
                                   ? num_super_kmers - (num_blocks - 1) * buffer_size
                                   : buffer_size;
            buff.resize(n);
            input.read(reinterpret_cast<char*>(buff.data()), buff.size() * sizeof(minimizer_tuple));
            const uint64_t chunk_size = (n + num_threads - 1) / num_threads;
            for (uint64_t t = 0; t != num_threads; ++t) {
                uint64_t begin = t * chunk_size;
                uint64_t end = (t == num_threads - 1) ? n : begin + chunk_size;
                threads.emplace_back([begin, end, &buff, &f]() {
                    for (uint64_t i = begin; i != end; ++i) {
                        buff[i].minimizer = f.lookup(buff[i].minimizer);
                    }
                });
            }
            for (auto& t : threads) {
                if (t.joinable()) t.join();
            }
            threads.clear();
            data.minimizers.sort_and_flush();
        }

        data.minimizers.finalize();
        data.minimizers.merge();
        input.close();
    }
    timer.stop();
    timings.push_back(timer.elapsed());
    print_time(timings.back(), data.num_kmers, "step 2.1: 're-sorting minimizers tuples'");
    timer.reset();
    /******/

    /* step 3: build index ***/
    timer.start();
    auto buckets_stats = build_index(data, m_buckets, num_minimizers, build_config);
    timer.stop();
    timings.push_back(timer.elapsed());
    print_time(timings.back(), data.num_kmers, "step 3: 'build_index'");
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

    double total_time = std::accumulate(timings.begin(), timings.end(), 0.0);
    print_time(total_time, data.num_kmers, "total_time");

    print_space_breakdown();

    if (build_config.verbose) buckets_stats.print_less();

    data.minimizers.remove_tmp_file();
}

}  // namespace sshash
