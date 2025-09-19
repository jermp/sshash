#pragma once

namespace sshash {

namespace perf {
using timer_type = essentials::timer<std::chrono::high_resolution_clock, std::chrono::nanoseconds>;
}

template <class kmer_t>
void perf_test_iterator(dictionary<kmer_t> const& dict) {
    perf::timer_type t;
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

template <class kmer_t>
void perf_test_lookup_access(dictionary<kmer_t> const& dict) {
    constexpr uint64_t num_queries = 1000000;
    constexpr uint64_t runs = 5;
    essentials::uniform_int_rng<uint64_t> distr(0, dict.size() - 1, essentials::get_random_seed());
    const uint64_t k = dict.k();
    // const uint64_t m = dict.m();
    std::string kmer(k, 0);
    std::string kmer_rc(k, 0);

    // std::string kmers;
    {
        // perf test positive lookup, using a std::vector<std::string>
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
                auto id = dict.lookup(string.c_str());
                essentials::do_not_optimize_away(id);
            }
        }
        t.stop();
        double nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
        std::cout << "lookup: avg_nanosec_per_positive_lookup " << nanosec_per_lookup << std::endl;

        std::vector<kmer_t> lookup_queries_uint;
        lookup_queries_uint.reserve(num_queries);
        for (auto const& kmer : lookup_queries) {
            kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(kmer.c_str(), k);
            lookup_queries_uint.push_back(uint_kmer);
        }
        t.reset();
        t.start();
        for (uint64_t r = 0; r != runs; ++r) {
            for (auto uint_kmer : lookup_queries_uint) {
                auto id = dict.lookup_uint(uint_kmer);
                essentials::do_not_optimize_away(id);
            }
        }
        t.stop();
        nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
        std::cout << "lookup_uint: avg_nanosec_per_positive_lookup " << nanosec_per_lookup
                  << std::endl;

        // std::vector<std::pair<kmer_t, minimizer_info>> lookup_queries_uint_minimizer;
        // lookup_queries_uint_minimizer.reserve(num_queries);
        // for (auto uint_kmer : lookup_queries_uint) {
        //     lookup_queries_uint_minimizer.push_back(
        //         {uint_kmer, util::compute_minimizer(uint_kmer, k, m, dict.hasher())});
        // }
        // t.reset();
        // t.start();
        // for (uint64_t r = 0; r != runs; ++r) {
        //     for (auto [uint_kmer, mini_info] : lookup_queries_uint_minimizer) {
        //         auto res = dict.lookup_uint_regular(uint_kmer, mini_info);
        //         essentials::do_not_optimize_away(res.kmer_id);
        //     }
        // }
        // t.stop();
        // nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
        // std::cout << "lookup_uint no-minimizer: avg_nanosec_per_positive_lookup "
        //           << nanosec_per_lookup << std::endl;

        // kmers.resize(num_queries * k);
        // uint64_t pos = 0;
        // for (auto const& string : lookup_queries) {
        //     kmers.replace(pos, k, string);
        //     pos += k;
        // }
    }
    // {
    // perf test positive lookup, using a single std::string with all kmers contatenated
    // perf::timer_type t;
    // t.start();
    // uint64_t pos = 0;
    // for (uint64_t r = 0; r != runs; ++r) {
    //     for (uint64_t i = 0; i != num_queries; ++i, pos += k) {
    //         auto id = dict.lookup(kmers.data() + pos);
    //         essentials::do_not_optimize_away(id);
    //     }
    //     pos = 0;

    //         /*
    //             loop-unrolling
    //         */
    //         for (uint64_t i = 0; i < num_queries; i += 8) {
    //             auto id0 = dict.lookup(kmers.data() + pos + 8 * 0);
    //             essentials::do_not_optimize_away(id0);

    //             auto id1 = dict.lookup(kmers.data() + pos + 8 * 1);
    //             essentials::do_not_optimize_away(id1);

    //             auto id2 = dict.lookup(kmers.data() + pos + 8 * 2);
    //             essentials::do_not_optimize_away(id2);

    //             auto id3 = dict.lookup(kmers.data() + pos + 8 * 3);
    //             essentials::do_not_optimize_away(id3);

    //             auto id4 = dict.lookup(kmers.data() + pos + 8 * 4);
    //             essentials::do_not_optimize_away(id4);

    //             auto id5 = dict.lookup(kmers.data() + pos + 8 * 5);
    //             essentials::do_not_optimize_away(id5);

    //             auto id6 = dict.lookup(kmers.data() + pos + 8 * 6);
    //             essentials::do_not_optimize_away(id6);

    //             auto id7 = dict.lookup(kmers.data() + pos + 8 * 7);
    //             essentials::do_not_optimize_away(id7);

    //             pos += 8 * 8;
    //         }
    //         pos = 0;
    //     }
    //     t.stop();
    //     double nanosec_per_lookup = t.elapsed() / (runs * num_queries);
    //     std::cout << "2. avg_nanosec_per_positive_lookup " << nanosec_per_lookup << std::endl;
    // }

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
                auto id = dict.lookup(string.c_str());
                essentials::do_not_optimize_away(id);
            }
        }
        t.stop();
        double nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
        std::cout << "avg_nanosec_per_negative_lookup " << nanosec_per_lookup << std::endl;
    }
    // {
    //     // perf test positive lookup_advanced
    //     std::vector<std::string> lookup_queries;
    //     lookup_queries.reserve(num_queries);
    //     for (uint64_t i = 0; i != num_queries; ++i) {
    //         uint64_t id = distr.gen();
    //         dict.access(id, kmer.data());
    //         if ((i & 1) == 0) {
    //             /* transform 50% of the kmers into their reverse complements */
    //             kmer_t::compute_reverse_complement(kmer.data(), kmer_rc.data(), k);
    //             lookup_queries.push_back(kmer_rc);
    //         } else {
    //             lookup_queries.push_back(kmer);
    //         }
    //     }
    //     perf::timer_type t;
    //     t.start();
    //     for (uint64_t r = 0; r != runs; ++r) {
    //         for (auto const& string : lookup_queries) {
    //             auto res = dict.lookup_advanced(string.c_str());
    //             essentials::do_not_optimize_away(res.kmer_id);
    //         }
    //     }
    //     t.stop();
    //     double nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
    //     std::cout << "avg_nanosec_per_positive_lookup_advanced " << nanosec_per_lookup <<
    //     std::endl;
    // }
    // {
    //     // perf test negative lookup_advanced
    //     std::vector<std::string> lookup_queries;
    //     lookup_queries.reserve(num_queries);
    //     for (uint64_t i = 0; i != num_queries; ++i) {
    //         random_kmer(kmer.data(), k);
    //         lookup_queries.push_back(kmer);
    //     }
    //     perf::timer_type t;
    //     t.start();
    //     for (uint64_t r = 0; r != runs; ++r) {
    //         for (auto const& string : lookup_queries) {
    //             auto res = dict.lookup_advanced(string.c_str());
    //             essentials::do_not_optimize_away(res.kmer_id);
    //         }
    //     }
    //     t.stop();
    //     double nanosec_per_lookup = t.elapsed() / (runs * lookup_queries.size());
    //     std::cout << "avg_nanosec_per_negative_lookup_advanced " << nanosec_per_lookup <<
    //     std::endl;
    // }
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
        std::cout << "avg_nanosec_per_access " << nanosec_per_access << std::endl;
    }
}  // namespace sshash

template <class kmer_t>
void perf_test_lookup_weight(dictionary<kmer_t> const& dict) {
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