#pragma once

namespace sshash {

void perf_test_iterator(dictionary const& dict) {
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::nanoseconds> t;
    t.start();
    auto it = dict.begin();
    while (it.has_next()) {
        auto [kmer_id, kmer] = it.next();
        essentials::do_not_optimize_away(kmer_id);
        essentials::do_not_optimize_away(kmer[0]);
    }
    t.stop();
    double avg_nanosec = t.elapsed() / dict.size();
    std::cout << "iterator: avg_nanosec_per_kmer " << avg_nanosec << std::endl;
}

void perf_test_lookup_access(dictionary const& dict) {
    constexpr uint64_t num_queries = 1000000;
    constexpr uint64_t runs = 5;
    essentials::uniform_int_rng<uint64_t> distr(0, dict.size() - 1, essentials::get_random_seed());
    uint64_t k = dict.k();
    std::string kmer(k, 0);
    std::string kmer_rc(k, 0);

    {
        // perf test positive lookup
        std::vector<std::string> lookup_queries;
        lookup_queries.reserve(num_queries);
        for (uint64_t i = 0; i != num_queries; ++i) {
            uint64_t id = distr.gen();
            dict.access(id, kmer.data());
            if ((i & 1) == 0) {
                /* transform 50% of the kmers into their reverse complements */
                util::compute_reverse_complement(kmer.data(), kmer_rc.data(), k);
                lookup_queries.push_back(kmer_rc);
            } else {
                lookup_queries.push_back(kmer);
            }
        }
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::nanoseconds> t;
        t.start();
        for (uint64_t r = 0; r != runs; ++r) {
            for (auto const& string : lookup_queries) {
                auto id = dict.lookup(string.c_str());
                essentials::do_not_optimize_away(id);
            }
        }
        t.stop();
        double nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
        std::cout << "avg_nanosec_per_positive_lookup " << nanosec_per_lookup << std::endl;
    }
    {
        // perf test negative lookup
        std::vector<std::string> lookup_queries;
        lookup_queries.reserve(num_queries);
        for (uint64_t i = 0; i != num_queries; ++i) {
            random_kmer(kmer.data(), k);
            lookup_queries.push_back(kmer);
        }
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::nanoseconds> t;
        t.start();
        for (uint64_t r = 0; r != runs; ++r) {
            for (auto const& string : lookup_queries) {
                auto id = dict.lookup(string.c_str());
                essentials::do_not_optimize_away(id);
            }
        }
        t.stop();
        double nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
        std::cout << "avg_nanosec_per_negative_lookup " << nanosec_per_lookup << std::endl;
    }
    {
        // perf test positive lookup_advanced
        std::vector<std::string> lookup_queries;
        lookup_queries.reserve(num_queries);
        for (uint64_t i = 0; i != num_queries; ++i) {
            uint64_t id = distr.gen();
            dict.access(id, kmer.data());
            if ((i & 1) == 0) {
                /* transform 50% of the kmers into their reverse complements */
                util::compute_reverse_complement(kmer.data(), kmer_rc.data(), k);
                lookup_queries.push_back(kmer_rc);
            } else {
                lookup_queries.push_back(kmer);
            }
        }
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::nanoseconds> t;
        t.start();
        for (uint64_t r = 0; r != runs; ++r) {
            for (auto const& string : lookup_queries) {
                auto res = dict.lookup_advanced(string.c_str());
                essentials::do_not_optimize_away(res.kmer_id);
            }
        }
        t.stop();
        double nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
        std::cout << "avg_nanosec_per_positive_lookup_advanced " << nanosec_per_lookup << std::endl;
    }
    {
        // perf test negative lookup_advanced
        std::vector<std::string> lookup_queries;
        lookup_queries.reserve(num_queries);
        for (uint64_t i = 0; i != num_queries; ++i) {
            random_kmer(kmer.data(), k);
            lookup_queries.push_back(kmer);
        }
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::nanoseconds> t;
        t.start();
        for (uint64_t r = 0; r != runs; ++r) {
            for (auto const& string : lookup_queries) {
                auto res = dict.lookup_advanced(string.c_str());
                essentials::do_not_optimize_away(res.kmer_id);
            }
        }
        t.stop();
        double nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
        std::cout << "avg_nanosec_per_negative_lookup_advanced " << nanosec_per_lookup << std::endl;
    }
    {
        // perf test access
        std::vector<uint64_t> access_queries;
        access_queries.reserve(num_queries);

        for (uint64_t i = 0; i != num_queries; ++i) access_queries.push_back(distr.gen());
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::nanoseconds> t;
        t.start();
        for (uint64_t r = 0; r != runs; ++r) {
            for (auto id : access_queries) {
                dict.access(id, kmer.data());
                essentials::do_not_optimize_away(kmer[0]);
            }
        }
        t.stop();
        double nanosec_per_access = t.elapsed() / static_cast<double>(runs * access_queries.size());
        std::cout << "avg_nanosec_per_access " << nanosec_per_access << std::endl;
    }
}

void perf_test_lookup_weight(dictionary const& dict) {
    if (!dict.weighted()) {
        std::cerr << "ERROR: the dictionary does not store weights" << std::endl;
        return;
    }

    constexpr uint64_t num_queries = 1000000;
    constexpr uint64_t runs = 5;
    essentials::uniform_int_rng<uint64_t> distr(0, dict.size() - 1, essentials::get_random_seed());
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
            util::compute_reverse_complement(kmer.data(), kmer_rc.data(), k);
            lookup_queries.push_back(kmer_rc);
        } else {
            lookup_queries.push_back(kmer);
        }
    }

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::nanoseconds> t;
    t.start();
    for (uint64_t r = 0; r != runs; ++r) {
        for (auto const& string : lookup_queries) {
            auto id = dict.lookup(string.c_str());
            auto w = dict.weight(id);
            essentials::do_not_optimize_away(w);
        }
    }
    t.stop();
    double nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
    std::cout << "avg_nanosec_per_positive_lookup_with_weight " << nanosec_per_lookup << std::endl;
}

}  // namespace sshash