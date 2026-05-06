#pragma once

#include <unordered_map>

#include "essentials.hpp"
#include "include/dictionary.hpp"
#include "include/offsets.hpp"
#include "include/builder/util.hpp"
#include "include/builder/disk_backed_strings.hpp"
#include "include/builder/disk_backed_offsets_builder.hpp"
#include "include/builder/streaming_compact_vector_writer.hpp"
#include "include/builder/streaming_save.hpp"
#include "include/buckets_statistics.hpp"

namespace sshash {

/*
    Helper: load a serialized bits::compact_vector back from a tmp file
    into the given in-RAM compact_vector. Used by the materializing build
    flow (after step 7) so that --check / queries can run.
*/
inline void materialize_compact_vector_from_file(bits::compact_vector& cv,
                                                 std::string const& filename) {
    essentials::loader loader(filename.c_str());
    loader.visit(cv);
}

/*
    Tmp file paths for the compact_vectors that step 7 spills to disk.
    Populated by build_sparse_and_skew_index; consumed by step 8 (either
    materialized back into RAM for `build()`, or injected into the output
    by `build_streaming_save()`).
*/
struct spilled_components {
    std::string control_codewords_path;
    std::string mid_load_buckets_path;
    std::string heavy_load_buckets_path;
    std::vector<std::string> skew_positions_paths;  // one entry per skew partition
    std::string codewords_mphf_path;                // step-4 minimizers MPHF
    std::vector<std::string> skew_mphfs_paths;      // one entry per skew partition

    void clear_files() {
        if (!control_codewords_path.empty()) std::remove(control_codewords_path.c_str());
        if (!mid_load_buckets_path.empty()) std::remove(mid_load_buckets_path.c_str());
        if (!heavy_load_buckets_path.empty()) std::remove(heavy_load_buckets_path.c_str());
        if (!codewords_mphf_path.empty()) std::remove(codewords_mphf_path.c_str());
        for (auto const& p : skew_positions_paths) {
            if (!p.empty()) std::remove(p.c_str());
        }
        for (auto const& p : skew_mphfs_paths) {
            if (!p.empty()) std::remove(p.c_str());
        }
    }
};

template <typename Kmer, typename Offsets>
struct dictionary_builder  //
{
    dictionary_builder(build_configuration const& build_config)
        : build_config(build_config)
        , num_kmers(0)
        , minimizers(build_config)
        , strings_run_id(pthash::clock_type::now().time_since_epoch().count())
        , total_time_musec(0) {}

    ~dictionary_builder() {
        strings_builder.remove_file();
        strings_offsets_builder.remove_file();
        spilled.clear_files();
    }

    /*
        Build a query-ready dictionary in `d`. After this returns, all
        spilled components and `d.m_spss.strings` are materialized in RAM
        (peak briefly equals the index size). Use this when the caller
        needs to query `d` post-build (e.g., `--check`).
    */
    void build(dictionary<Kmer, Offsets>& d, std::string const& filename) {
        run_steps_1_through_7(d, filename);
        do_step("step 8 (materialize spilled components to RAM)", [&]() {
            materialize_spilled_into(d);
            strings_builder.load_into(d.m_spss.strings);
            strings_builder.remove_file();
            spilled.clear_files();
        });
        finalize_stats(d);
    }

    /*
        Build the dictionary and stream-save it to `output_filename` without
        ever materializing the spilled components or `strings` in RAM.
        After this returns, `d` is *not* query-ready. Use this when the
        caller only needs the on-disk index file and wants to keep peak RAM
        bounded by the build phase.
    */
    void build_streaming_save(dictionary<Kmer, Offsets>& d,                  //
                              std::string const& filename,                   //
                              std::string const& output_filename)            //
    {
        run_steps_1_through_7(d, filename);
        do_step("step 8 (stream-save dictionary to disk)", [&]() {
            /* Address+type-keyed substitution map. The saver replaces the
               visit byte content of any registered (address, type) pair
               with the bytes of the corresponding tmp file. Type matching
               disambiguates the case where a struct's address coincides
               with the address of its first member. */
            std::unordered_map<void const*, typed_address_sub> subs;
            if (!spilled.control_codewords_path.empty()) {
                register_sub(subs, &d.m_ssi.codewords.control_codewords,
                             spilled.control_codewords_path);
            }
            if (!spilled.mid_load_buckets_path.empty()) {
                register_sub(subs, &d.m_ssi.mid_load_buckets, spilled.mid_load_buckets_path);
            }
            if (!spilled.heavy_load_buckets_path.empty()) {
                register_sub(subs, &d.m_ssi.ski.heavy_load_buckets,
                             spilled.heavy_load_buckets_path);
            }
            if (!spilled.codewords_mphf_path.empty()) {
                register_sub(subs, &d.m_ssi.codewords.mphf, spilled.codewords_mphf_path);
            }
            /* Skew positions / mphfs: populate the owning_spans with
               placeholders so the visit walks the right number of entries
               and we can take their addresses for substitution. */
            const std::size_t num_part = std::max(spilled.skew_positions_paths.size(),
                                                  spilled.skew_mphfs_paths.size());
            if (num_part > 0) {
                std::vector<bits::compact_vector> position_placeholders(num_part);
                std::vector<kmers_pthash_type<Kmer>> mphf_placeholders(num_part);
                d.m_ssi.ski.positions = std::move(position_placeholders);
                d.m_ssi.ski.mphfs = std::move(mphf_placeholders);
                for (std::size_t i = 0; i != spilled.skew_positions_paths.size(); ++i) {
                    if (!spilled.skew_positions_paths[i].empty()) {
                        register_sub(subs, &d.m_ssi.ski.positions[i],
                                     spilled.skew_positions_paths[i]);
                    }
                }
                for (std::size_t i = 0; i != spilled.skew_mphfs_paths.size(); ++i) {
                    if (!spilled.skew_mphfs_paths[i].empty()) {
                        register_sub(subs, &d.m_ssi.ski.mphfs[i],
                                     spilled.skew_mphfs_paths[i]);
                    }
                }
            }
            save_streaming(d, output_filename.c_str(), &d.m_spss.strings, strings_builder,
                           std::move(subs));
            strings_builder.remove_file();
            spilled.clear_files();
        });
        finalize_stats(d);
    }

    build_configuration build_config;
    uint64_t num_kmers;
    minimizers_tuples minimizers;
    disk_backed_offsets_builder<Offsets> strings_offsets_builder;
    disk_backed_strings strings_builder;
    weights::builder weights_builder;
    spilled_components spilled;

    uint64_t strings_run_id;

    essentials::timer_type timer;
    essentials::json_lines build_stats;
    uint64_t total_time_musec;

private:
    /* Load each spilled compact_vector tmp file back into the corresponding
       in-RAM compact_vector inside `d`. Used by the materializing build
       flow so queries can run against `d` (e.g., during --check). */
    void materialize_spilled_into(dictionary<Kmer, Offsets>& d) {
        if (!spilled.control_codewords_path.empty()) {
            materialize_compact_vector_from_file(d.m_ssi.codewords.control_codewords,
                                                 spilled.control_codewords_path);
        }
        if (!spilled.mid_load_buckets_path.empty()) {
            materialize_compact_vector_from_file(d.m_ssi.mid_load_buckets,
                                                 spilled.mid_load_buckets_path);
        }
        if (!spilled.heavy_load_buckets_path.empty()) {
            materialize_compact_vector_from_file(d.m_ssi.ski.heavy_load_buckets,
                                                 spilled.heavy_load_buckets_path);
        }
        /* Reload the spilled MPHFs back into RAM so queries work. */
        if (!spilled.codewords_mphf_path.empty()) {
            essentials::loader loader(spilled.codewords_mphf_path.c_str());
            loader.visit(d.m_ssi.codewords.mphf);
        }
        const std::size_t num_part = std::max(spilled.skew_positions_paths.size(),
                                              spilled.skew_mphfs_paths.size());
        if (num_part > 0) {
            std::vector<bits::compact_vector> positions_vec(num_part);
            std::vector<kmers_pthash_type<Kmer>> mphfs_vec(num_part);
            for (std::size_t i = 0; i != spilled.skew_positions_paths.size(); ++i) {
                if (!spilled.skew_positions_paths[i].empty()) {
                    materialize_compact_vector_from_file(positions_vec[i],
                                                         spilled.skew_positions_paths[i]);
                }
            }
            for (std::size_t i = 0; i != spilled.skew_mphfs_paths.size(); ++i) {
                if (!spilled.skew_mphfs_paths[i].empty()) {
                    essentials::loader loader(spilled.skew_mphfs_paths[i].c_str());
                    loader.visit(mphfs_vec[i]);
                }
            }
            d.m_ssi.ski.positions = std::move(positions_vec);
            d.m_ssi.ski.mphfs = std::move(mphfs_vec);
        }
    }

    void run_steps_1_through_7(dictionary<Kmer, Offsets>& d, std::string const& filename) {
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

        {
            std::stringstream ss_strings;
            ss_strings << build_config.tmp_dirname << "/sshash.tmp.run_" << strings_run_id
                       << ".strings.bin";
            strings_builder.open_for_writing(ss_strings.str());
            std::stringstream ss_offsets;
            ss_offsets << build_config.tmp_dirname << "/sshash.tmp.run_" << strings_run_id
                       << ".strings_offsets.bin";
            strings_offsets_builder.open_for_writing(ss_offsets.str());
        }

        do_step("step 1 (encode strings)", [&]() {
            encode_strings(filename);
            strings_builder.freeze();
            strings_offsets_builder.freeze();
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
    }

    void finalize_stats(dictionary<Kmer, Offsets>& d) {
        if (build_config.verbose) {
            print_time(total_time_musec, "total time");
            /* `print_space_breakdown` reads d.m_spss.strings; only safe in
               the materialize-to-RAM flow. */
            if (d.m_spss.strings.num_bits() > 0) d.print_space_breakdown();
        }

        build_stats.add("total_build_time_in_microsec", total_time_musec);
        build_stats.add("index_size_in_bytes", (d.num_bits() + 7) / 8);
        build_stats.add("num_kmers", d.num_kmers());

        if (build_config.verbose) build_stats.print();
    }

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
        /* Stream minimizers from disk via std::ifstream (no mmap); the
           iterator yields each distinct minimizer once, matching what
           `minimizers_tuples_iterator` did over the mmap'd file. */
        streaming_minimizers_iterator iterator;
        iterator.open(minimizers.get_minimizers_filename());
        d.m_ssi.codewords.build(iterator, num_minimizers, build_config);
        iterator.close();
        assert(d.m_ssi.codewords.size() == num_minimizers);
    }

    void hash_minimizers(dictionary<Kmer, Offsets>& d) {
        std::string filename = minimizers.get_minimizers_filename();
        std::ifstream input(filename, std::ifstream::binary);

        auto& f_mut = d.m_ssi.codewords.mphf;
        auto const& f = f_mut;
        const uint64_t num_threads = build_config.num_threads;
        const uint64_t num_files_to_merge = minimizers.num_files_to_merge();

        minimizers.init();

        uint64_t RAM_available_in_bytes = essentials::GiB / 2;  // at least 0.5 GB
        {
            /* `strings_builder` is now disk-backed; its in-RAM footprint is
               bounded by its window size, not by the strings size. */
            const uint64_t RAM_taken_in_bytes =
                f.num_bits() / 8 + strings_offsets_builder.num_bytes();
            const uint64_t RAM_limit_in_bytes = build_config.ram_limit_in_GiB * essentials::GiB;
            if (RAM_limit_in_bytes > RAM_taken_in_bytes) {
                RAM_available_in_bytes = std::max<uint64_t>(RAM_limit_in_bytes - RAM_taken_in_bytes,
                                                            RAM_available_in_bytes);
            }
        }

        const uint64_t num_super_kmers = minimizers.num_super_kmers();
        /* Cap the in-RAM buffer at ram_limit/8 worth of tuples so that
           even when subsequent steps fragment the heap, step 5's lingering
           pages don't blow past the budget when stacked with later step's
           allocations. */
        const uint64_t buffer_cap_bytes =
            (build_config.ram_limit_in_GiB * essentials::GiB) / 8;
        const uint64_t buffer_cap_records =
            std::max<uint64_t>(uint64_t(1) << 16, buffer_cap_bytes / sizeof(minimizer_tuple));
        const uint64_t buffer_size_unbounded =
            num_files_to_merge == 1
                ? num_super_kmers
                : (RAM_available_in_bytes / (3 * sizeof(minimizer_tuple)));
        const uint64_t buffer_size = std::min(buffer_size_unbounded, buffer_cap_records);
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

        /* The codewords MPHF is no longer needed during build (step 6 onward
           reads minimizer values that step 5 has already replaced with
           mphf hashes; step 7 references mphf hashes only as bucket ids).
           Spill it to disk and free its in-RAM footprint. */
        {
            std::stringstream ss;
            ss << build_config.tmp_dirname << "/sshash.tmp.run_" << strings_run_id
               << ".codewords_mphf.bin";
            spilled.codewords_mphf_path = ss.str();
            essentials::save(f_mut, spilled.codewords_mphf_path.c_str());
            f_mut = minimizers_pthash_type{};
        }
    }
};

}  // namespace sshash
