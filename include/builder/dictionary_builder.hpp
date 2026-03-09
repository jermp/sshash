#pragma once

#include "essentials.hpp"
#include "include/dictionary.hpp"
#include "include/offsets.hpp"
#include "include/builder/util.hpp"
#include "include/buckets_statistics.hpp"

namespace sshash {

template <typename Kmer, typename Offsets>
struct dictionary_builder  //
{
    dictionary_builder(build_configuration const& build_config)
        : build_config(build_config), num_kmers(0), minimizers(build_config), total_time_musec(0) {}

    void build(dictionary<Kmer, Offsets>& d, std::string const& filename)  //
    {
        d.m_k = build_config.k;
        d.m_m = build_config.m;
        d.m_spss.k = build_config.k;
        d.m_spss.m = build_config.m;
        d.m_canonical = build_config.canonical;
        d.m_hasher.seed(build_config.seed);

        build_stats.add("input_filename", filename.c_str());
        build_stats.add("k", d.m_k);
        build_stats.add("m", d.m_m);
        build_stats.add("canonical", d.m_canonical ? "true" : "false");
        build_stats.add("seed", build_config.seed);
        build_stats.add("num_threads", build_config.num_threads);

        total_time_musec = 0;

        do_step("step 1 (encode strings)", [&]() {
            encode_strings(filename);
            d.m_num_kmers = num_kmers;
            assert(strings_offsets_builder.size() >= 2);
            d.m_num_strings = strings_offsets_builder.size() - 1;
        });

        if (build_config.weighted) {
            do_step("step 1.1 (build weights)", [&]() { weights_builder.build(d.m_weights); });
        }

        do_step("step 2 (compute minimizer tuples)", [&]() { compute_minimizer_tuples(); });

        do_step("step 3 (merging minimizer tuples)", [&]() { minimizers.merge(); });
        if (build_config.verbose) {
            std::cout << "num_minimizers = " << minimizers.num_minimizers() << std::endl;
            std::cout << "num_minimizer_positions = " << minimizers.num_minimizer_positions()
                      << std::endl;
            std::cout << "num_super_kmers = " << minimizers.num_super_kmers() << std::endl;
        }

        do_step("step 4 (build mphf)", [&]() { build_mphf(d); });

        do_step("step 5 (replacing minimizer values with MPHF hashes)",
                [&]() { hash_minimizers(d); });

        do_step("step 6 (merging minimizers tuples)", [&]() { minimizers.merge(); });

        do_step("step 7 (build sparse and skew index)", [&]() {
            build_sparse_and_skew_index(d);
            minimizers.remove_tmp_file();
            assert(strings_offsets_builder.size() == 0);
        });

        if (build_config.verbose) {
            print_time(total_time_musec, "total time");
            d.print_space_breakdown();
        }

        build_stats.add("total_build_time_in_microsec", total_time_musec);
        build_stats.add("index_size_in_bytes", (d.num_bits() + 7) / 8);
        build_stats.add("num_kmers", d.num_kmers());

        if (build_config.verbose) build_stats.print();
    }

    build_configuration build_config;
    uint64_t num_kmers;
    minimizers_tuples minimizers;
    typename Offsets::builder strings_offsets_builder;
    bits::bit_vector::builder strings_builder;
    weights::builder weights_builder;

    essentials::timer_type timer;
    essentials::json_lines build_stats;
    uint64_t total_time_musec;

private:
    void print_time(double time_in_musec, std::string const& message) {
        std::cout << "=== " << message << ": " << time_in_musec / 1'000'000 << " [sec] ("
                  << (time_in_musec * 1000) / num_kmers << " [ns/kmer])" << std::endl;
    }

    template <typename Callback>
    void do_step(std::string const& step, Callback const& f) {
        timer.start();
        f();
        timer.stop();
        uint64_t step_elapsed_time_musec = timer.elapsed();
        total_time_musec += step_elapsed_time_musec;
        if (build_config.verbose) print_time(step_elapsed_time_musec, step);
        build_stats.add(step, step_elapsed_time_musec);
        timer.reset();
    }

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