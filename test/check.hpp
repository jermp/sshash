#pragma once

#include <mutex>

namespace sshash {

template <typename Dict>
bool check_dictionary(Dict const& dict) {
    const uint64_t k = dict.k();
    const uint64_t n = dict.num_kmers();

    const uint64_t num_threads = std::thread::hardware_concurrency();
    std::cout << "checking correctness of access and positive lookup using " << num_threads
              << " threads..." << std::endl;

    std::mutex print_mutex;

    auto worker = [&](uint64_t start, uint64_t end, size_t thread_id) {
        std::string kmer(k, 0);
        for (uint64_t id = start; id != end; ++id)  //
        {
            uint64_t count = id - start;
            if (count != 0 and count % 15'000'000 == 0) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[Thread " << thread_id << "] Checked " << count
                          << " kmers (local progress)" << std::endl;
            }

            dict.access(id, kmer.data());
            uint64_t got_id = dict.lookup(kmer.c_str()).kmer_id;

            if (got_id == constants::invalid_uint64) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cerr << "[Thread " << thread_id << "] kmer '" << kmer << "' not found!\n";
                return;
            }
            if (got_id >= n) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cerr << "[Thread " << thread_id << "] ERROR: id out of range " << got_id << "/"
                          << n << "\n";
                return;
            }
            if (got_id != id) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cerr << "[Thread " << thread_id << "] expected id " << id << " but got id "
                          << got_id << "\n";
                return;
            }
        }
        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[Thread " << thread_id << "] Finished range [" << start << ", " << end
                      << ")\n";
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (size_t t = 0, chunk_size = (n + num_threads - 1) / num_threads; t != num_threads; ++t) {
        uint64_t start = t * chunk_size;
        uint64_t end = std::min(n, start + chunk_size);
        threads.emplace_back(worker, start, end, t);
    }

    for (auto& th : threads) th.join();

    std::cout << "EVERYTHING OK!" << std::endl;

    return check_correctness_negative_lookup(dict);
}

template <typename Dict>
bool check_correctness_negative_lookup(Dict const& dict) {
    std::cout << "checking correctness of negative lookup with random kmers..." << std::endl;
    const uint64_t num_lookups = std::min<uint64_t>(1000000, dict.num_kmers());
    std::string kmer(dict.k(), 0);
    for (uint64_t i = 0; i != num_lookups; ++i) {
        random_kmer(kmer.data(), dict.k());
        /*
            We could use a std::unordered_set to check if kmer is really absent,
            but that would take much more memory...
        */
        auto res = dict.lookup(kmer.c_str());
        if (res.kmer_id != constants::invalid_uint64) {
            std::cout << "kmer '" << kmer << "' found!" << std::endl;
        }
    }
    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

template <typename Dict>
bool check_correctness_navigational_contig_query(Dict const& dict)  //
{
    std::cout << "checking correctness of navigational queries for contigs..." << std::endl;
    using kmer_t = typename Dict::kmer_type;
    const uint64_t num_strings = dict.num_strings();
    const uint64_t k = dict.k();
    uint64_t kmer_id = 0;
    std::string kmer(k, 0);
    uint64_t contig_id = 0;
    for (; contig_id != num_strings; ++contig_id) {
        if (contig_id != 0 and contig_id % 1000000 == 0) {
            std::cout << "checked " << contig_id << "/" << num_strings << " contigs" << std::endl;
        }

        auto res = dict.contig_neighbours(contig_id);
        uint64_t contig_size = dict.contig_size(contig_id);

        uint64_t begin_kmer_id = kmer_id;
        dict.access(begin_kmer_id, kmer.data());
        auto backward = dict.kmer_backward_neighbours(kmer.data());
        for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
            equal_lookup_result(backward.backward[i], res.backward[i]);
        }

        uint64_t end_kmer_id = kmer_id + contig_size - 1;
        dict.access(end_kmer_id, kmer.data());
        auto forward = dict.kmer_forward_neighbours(kmer.data());
        for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
            equal_lookup_result(forward.forward[i], res.forward[i]);
        }
        kmer_id += contig_size;
    }
    std::cout << "checked " << contig_id << " contigs" << std::endl;
    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

template <typename Dict>
bool check_correctness_kmer_iterator(Dict const& dict) {
    std::cout << "checking correctness of kmer iterator..." << std::endl;
    const uint64_t k = dict.k();
    std::string read_kmer(k, 0);
    std::string expected_kmer(k, 0);
    constexpr uint64_t runs = 4;
    essentials::uniform_int_rng<uint64_t> distr(0, dict.num_kmers() - 1,
                                                essentials::get_random_seed());
    for (uint64_t run = 0; run != runs; ++run) {
        uint64_t from_kmer_id = run == 0 ? 0 : distr.gen();
        auto it = dict.at_kmer_id(from_kmer_id);
        while (it.has_next()) {
            auto [kmer_id, kmer] = it.next();
            util::uint_kmer_to_string<typename Dict::kmer_type>(kmer, read_kmer.data(), k);
            dict.access(kmer_id, expected_kmer.data());
            if (read_kmer != expected_kmer or kmer_id != from_kmer_id) {
                std::cout << "got (" << kmer_id << ",'" << read_kmer << "')";
                std::cout << " but ";
                std::cout << "expected (" << from_kmer_id << ",'" << expected_kmer << "')"
                          << std::endl;
                return false;
            }
            ++from_kmer_id;
        }
        assert(from_kmer_id == dict.num_kmers());
    }
    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

template <typename Dict>
bool check_correctness_contig_iterator(Dict const& dict) {
    std::cout << "checking correctness of contig iterator..." << std::endl;
    const uint64_t k = dict.k();
    std::string read_kmer(k, 0);
    std::string expected_kmer(k, 0);
    for (uint64_t contig_id = 0; contig_id != dict.num_strings(); ++contig_id) {
        auto [begin, _] = dict.contig_offsets(contig_id);
        uint64_t from_kmer_id = begin - contig_id * (dict.k() - 1);
        auto it = dict.at_contig_id(contig_id);
        while (it.has_next()) {
            auto [kmer_id, kmer] = it.next();
            util::uint_kmer_to_string<typename Dict::kmer_type>(kmer, read_kmer.data(), k);
            dict.access(kmer_id, expected_kmer.data());
            if (read_kmer != expected_kmer or kmer_id != from_kmer_id) {
                std::cout << "got (" << kmer_id << ",'" << read_kmer << "')";
                std::cout << " but ";
                std::cout << "expected (" << from_kmer_id << ",'" << expected_kmer << "')"
                          << std::endl;
                return false;
            }
            ++from_kmer_id;
        }
    }
    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

}  // namespace sshash