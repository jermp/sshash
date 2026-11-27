#pragma once

#include "include/batched_lookup.hpp"
#include "include/sorted_lookup.hpp"
#include "include/hugepages.hpp"

namespace sshash {

/*
    The empirical crossover between the two batched random-lookup engines:
    lookup_query_batched hides the latency of random accesses with a small
    in-flight batch, lookup_query_sorted eliminates the randomness by
    sorting a large batch along each data structure's own layout. The first
    is insensitive to batch size, the second amortizes its sorts better the
    larger the batch is relative to the index; this driver measures both on
    the same query sets, at batch sizes from 10^4 up to a given maximum.

    Queries are packed kmers (no string parsing), so both engines run their
    uint-kmer API head-to-head and 10^8 queries fit in memory.
*/
template <typename Dict>
void crossover_bench(Dict const& dict, const uint64_t max_n,
                     essentials::json_lines& perf_stats)  //
{
    using kmer_t = typename Dict::kmer_type;
    const uint64_t k = dict.k();

    std::cout << "generating " << max_n << " positive and " << max_n << " negative queries..."
              << std::endl;
    essentials::uniform_int_rng<uint64_t> distr(0, dict.num_kmers() - 1,
                                                essentials::get_random_seed());
    std::string kmer(k, 0);
    std::vector<kmer_t> positives(max_n), negatives(max_n);
    for (uint64_t i = 0; i != max_n; ++i) {
        dict.access(distr.gen(), kmer.data());
        kmer_t x = lookup_kernels::pack_string<kmer_t>(kmer.data(), k);
        if ((i & 1) == 0) x.reverse_complement_inplace(k);
        positives[i] = x;
    }
    for (uint64_t i = 0; i != max_n; ++i) {
        random_kmer(kmer.data(), k);
        negatives[i] = lookup_kernels::pack_string<kmer_t>(kmer.data(), k);
    }
    std::vector<lookup_result> results(max_n);

    lookup_query_batched<Dict> batched(&dict);
    lookup_query_sorted<Dict> sorted(&dict);

    auto measure = [&](auto& engine, kmer_t const* queries, const uint64_t n) {
        const uint64_t runs = std::max<uint64_t>(1, std::min<uint64_t>(5, 10'000'000 / n));
        perf::timer_type t;
        t.start();
        for (uint64_t r = 0; r != runs; ++r) {
            engine.lookup_batch(queries, n, results.data());
            essentials::do_not_optimize_away(results[r % n].kmer_id);
        }
        t.stop();
        return t.elapsed() / (runs * n);
    };

    for (uint64_t n = 10'000; n <= max_n; n *= 10) {
        for (int w = 0; w != 2; ++w) {
            kmer_t const* queries = (w == 0 ? positives : negatives).data();
            char const* wl = w == 0 ? "positive" : "negative";
            const double b = measure(batched, queries, n);
            const double s = measure(sorted, queries, n);
            std::cout << "n = " << n << ", " << wl  //
                      << ": batched = " << b << " ns/kmer, sorted = " << s << " ns/kmer"
                      << std::endl;
            perf_stats.add(("batched " + std::string(wl) + " n=" + std::to_string(n)).c_str(), b);
            perf_stats.add(("sorted " + std::string(wl) + " n=" + std::to_string(n)).c_str(), s);
        }
    }
}

}  // namespace sshash

int crossover(int argc, char** argv) {
    using namespace sshash;
    cmd_line_parser::parser parser(argc, argv);
    parser.add("index_filename", "Must be a file generated with the tool 'build'.", "-i", true);
    parser.add("max_n",
               "Maximum batch size (default 10^7). Memory grows with it: two query arrays and "
               "one result array of this many entries.",
               "-n", false);
    parser.add("hugepages",
               "Advise the kernel to back the index with huge pages (Linux only), removing most "
               "TLB misses of large indexes.",
               "--hugepages", false, true);
    parser.add("verbose", "Verbose output.", "--verbose", false, true);
    if (!parser.parse()) return 0;
    auto index_filename = parser.get<std::string>("index_filename");
    bool verbose = parser.get<bool>("verbose");
    uint64_t max_n = 10'000'000;
    if (parser.parsed("max_n")) max_n = parser.get<uint64_t>("max_n");
    if (max_n < 10'000) max_n = 10'000;

    dictionary_type dict;
    open_dictionary(dict, index_filename, false, verbose);

    if (parser.get<bool>("hugepages")) {
        const auto advisor = advise_hugepages(dict);
        std::cout << "huge pages: " << advisor.collapsed_bytes()
                  << " bytes collapsed synchronously, " << advisor.lazy_bytes()
                  << " bytes advised for background collapse" << std::endl;
    }

    essentials::json_lines perf_stats;
    perf_stats.add("index_filename", index_filename.c_str());
    perf_stats.add("k", dict.k());
    perf_stats.add("m", dict.m());
    crossover_bench(dict, max_n, perf_stats);
    perf_stats.print();
    return 0;
}
