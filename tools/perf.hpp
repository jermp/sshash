#pragma once

namespace sshash {

namespace perf {
using timer_type = essentials::timer<std::chrono::high_resolution_clock, std::chrono::nanoseconds>;
}

template <typename Dict>
void perf_test_iterator(Dict const& dict, essentials::json_lines& perf_stats) {
    perf::timer_type t;
    t.start();
    auto it = dict.begin();
    uint64_t n = std::min<uint64_t>(dict.num_kmers(), 100'000'000);
    for (uint64_t i = 0; i != n; ++i) {
        auto [kmer_id, kmer] = it.next();
        essentials::do_not_optimize_away(kmer_id);
        essentials::do_not_optimize_away(kmer.at(0));
    }
    t.stop();
    double avg_nanosec = t.elapsed() / n;
    std::cout << "iterator (avg_nanosec_per_kmer) = " << avg_nanosec << std::endl;
    perf_stats.add("iterator (avg_nanosec_per_kmer)", avg_nanosec);
}

template <typename Dict>
void perf_test_lookup_access(Dict const& dict, essentials::json_lines& perf_stats)  //
{
    using kmer_t = typename Dict::kmer_type;
    constexpr uint64_t num_queries = 1'000'000;
    constexpr uint64_t runs = 5;
    essentials::uniform_int_rng<uint64_t> distr(0, dict.num_kmers() - 1,
                                                essentials::get_random_seed());
    const uint64_t k = dict.k();
    std::string kmer(k, 0);
    std::string kmer_rc(k, 0);

    {
        std::vector<std::string> lookup_queries;
        lookup_queries.reserve(num_queries);
        for (uint64_t i = 0; i != num_queries; ++i) {
            uint64_t id = distr.gen();
            dict.access(id, kmer.data());
            if ((i & 1) == 0) {
                /* transform 50% of the kmers into their reverse complements */
                kmer_t::compute_reverse_complement(kmer.data(), kmer_rc.data(), k);
                lookup_queries.push_back(kmer_rc);
            } else {
                lookup_queries.push_back(kmer);
            }
        }

        perf::timer_type t;
        t.start();
        for (uint64_t r = 0; r != runs; ++r) {
            for (auto const& string : lookup_queries) {
                auto res = dict.lookup(string.c_str());
                essentials::do_not_optimize_away(res.kmer_id);
            }
        }
        t.stop();
        double nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
        std::cout << "positive lookup (avg_nanosec_per_kmer) = " << nanosec_per_lookup << std::endl;
        perf_stats.add("positive lookup (avg_nanosec_per_kmer)", nanosec_per_lookup);
    }

    {
        // perf test negative lookup
        std::vector<std::string> lookup_queries;
        lookup_queries.reserve(num_queries);
        for (uint64_t i = 0; i != num_queries; ++i) {
            random_kmer(kmer.data(), k);
            lookup_queries.push_back(kmer);
        }
        perf::timer_type t;
        t.start();
        for (uint64_t r = 0; r != runs; ++r) {
            for (auto const& string : lookup_queries) {
                auto res = dict.lookup(string.c_str());
                essentials::do_not_optimize_away(res.kmer_id);
            }
        }
        t.stop();
        double nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
        std::cout << "negative lookup (avg_nanosec_per_kmer) " << nanosec_per_lookup << std::endl;
        perf_stats.add("negative lookup (avg_nanosec_per_kmer)", nanosec_per_lookup);
    }

    {
        // perf test access
        std::vector<uint64_t> access_queries;
        access_queries.reserve(num_queries);
        for (uint64_t i = 0; i != num_queries; ++i) access_queries.push_back(distr.gen());
        perf::timer_type t;
        t.start();
        for (uint64_t r = 0; r != runs; ++r) {
            for (auto id : access_queries) {
                dict.access(id, kmer.data());
                essentials::do_not_optimize_away(kmer[0]);
            }
        }
        t.stop();
        double nanosec_per_access = t.elapsed() / static_cast<double>(runs * access_queries.size());
        std::cout << "access (avg_nanosec_per_kmer) = " << nanosec_per_access << std::endl;
        perf_stats.add("access (avg_nanosec_per_kmer)", nanosec_per_access);
    }
}

template <typename Dict>
void perf_test_lookup_weight(Dict const& dict, essentials::json_lines& perf_stats)  //
{
    using kmer_t = typename Dict::kmer_type;

    if (!dict.weighted()) {
        std::cerr << "ERROR: the dictionary does not store weights" << std::endl;
        return;
    }

    constexpr uint64_t num_queries = 1'000'000;
    constexpr uint64_t runs = 5;
    essentials::uniform_int_rng<uint64_t> distr(0, dict.num_kmers() - 1,
                                                essentials::get_random_seed());
    uint64_t k = dict.k();
    std::string kmer(k, 0);
    std::string kmer_rc(k, 0);

    std::vector<std::string> lookup_queries;
    lookup_queries.reserve(num_queries);
    for (uint64_t i = 0; i != num_queries; ++i) {
        uint64_t id = distr.gen();
        dict.access(id, kmer.data());
        if ((i & 1) == 0) {
            /* transform 50% of the kmers into their reverse complements */
            kmer_t::compute_reverse_complement(kmer.data(), kmer_rc.data(), k);
            lookup_queries.push_back(kmer_rc);
        } else {
            lookup_queries.push_back(kmer);
        }
    }

    perf::timer_type t;
    t.start();
    for (uint64_t r = 0; r != runs; ++r) {
        for (auto const& string : lookup_queries) {
            auto res = dict.lookup(string.c_str());
            auto w = dict.weight(res.kmer_id);
            essentials::do_not_optimize_away(w);
        }
    }
    t.stop();
    double nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
    std::cout << "positive lookup + weight (avg_nanosec_per_kmer) = " << nanosec_per_lookup
              << std::endl;
    perf_stats.add("positive lookup + weight (avg_nanosec_per_kmer)", nanosec_per_lookup);
}

}  // namespace sshash