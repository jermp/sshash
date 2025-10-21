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

    auto worker = [&](uint64_t start, uint64_t end, uint64_t thread_id) {
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
    for (uint64_t t = 0, chunk_size = (n + num_threads - 1) / num_threads; t != num_threads; ++t) {
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
bool check_correctness_navigational_string_query(Dict const& dict)  //
{
    using kmer_t = typename Dict::kmer_type;
    const uint64_t num_strings = dict.num_strings();
    const uint64_t k = dict.k();

    const uint64_t num_threads = std::thread::hardware_concurrency();
    std::cout << "checking correctness of navigational queries for strings using " << num_threads
              << " threads ..." << std::endl;

    std::mutex print_mutex;

    auto worker = [&](uint64_t start, uint64_t end, uint64_t start_kmer_id, uint64_t thread_id) {
        std::string kmer(k, 0);
        uint64_t kmer_id = start_kmer_id;

        for (uint64_t string_id = start; string_id < end; ++string_id) {
            if (string_id != start && (string_id - start) % 1'000'000 == 0) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[Thread " << thread_id << "] checked " << (string_id - start)
                          << " strings (local progress)\n";
            }

            auto res = dict.string_neighbours(string_id);
            uint64_t string_size = dict.string_size(string_id);

            // Check backward neighbours at beginning k-mer
            uint64_t begin_kmer_id = kmer_id;
            dict.access(begin_kmer_id, kmer.data());
            auto backward = dict.kmer_backward_neighbours(kmer.data());
            for (uint64_t i = 0; i < kmer_t::alphabet_size; i++) {
                equal_lookup_result(backward.backward[i], res.backward[i]);
            }

            // Check forward neighbours at end k-mer
            uint64_t end_kmer_id = kmer_id + string_size - 1;
            dict.access(end_kmer_id, kmer.data());
            auto forward = dict.kmer_forward_neighbours(kmer.data());
            for (uint64_t i = 0; i < kmer_t::alphabet_size; i++) {
                equal_lookup_result(forward.forward[i], res.forward[i]);
            }

            kmer_id += string_size;
        }
        {
            std::lock_guard<std::mutex> lock(print_mutex);
            std::cout << "[Thread " << thread_id << "] Finished range [" << start << ", " << end
                      << ")\n";
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    for (uint64_t t = 0, current_start = 0, current_kmer_id = 0,
                  chunk_size = (num_strings + num_threads - 1) / num_threads;
         t < num_threads && current_start < num_strings; ++t)  //
    {
        uint64_t start = current_start;
        uint64_t end = std::min(num_strings, start + chunk_size);

        // compute starting kmer_id for this thread
        uint64_t start_kmer_id = current_kmer_id;
        for (uint64_t i = start; i < end; ++i) current_kmer_id += dict.string_size(i);

        threads.emplace_back(worker, start, end, start_kmer_id, t);
        current_start = end;
    }

    for (auto& th : threads) th.join();

    std::cout << "checked " << num_strings << " strings" << std::endl;
    std::cout << "EVERYTHING OK!" << std::endl;

    return true;
}

template <typename Dict>
bool check_correctness_kmer_iterator(Dict const& dict)  //
{
    const uint64_t num_kmers = dict.num_kmers();
    const uint64_t k = dict.k();
    const uint64_t num_threads = std::thread::hardware_concurrency();
    std::cout << "checking correctness of kmer iterator using " << num_threads << " threads ..."
              << std::endl;

    std::mutex print_mutex;

    auto worker = [&](uint64_t start, uint64_t end, uint64_t thread_id) {
        assert(end > start);
        std::string read_kmer(k, 0);
        std::string expected_kmer(k, 0);
        for (auto it = dict.at_kmer_id(start); start != end; ++start) {
            uint64_t count = end - start;
            if (count != 0 and count % 100'000'000 == 0) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[Thread " << thread_id << "] Checked " << count
                          << " kmers (local progress)" << std::endl;
            }
            auto [kmer_id, kmer] = it.next();
            util::uint_kmer_to_string<typename Dict::kmer_type>(kmer, read_kmer.data(), k);
            dict.access(kmer_id, expected_kmer.data());
            if (read_kmer != expected_kmer or kmer_id != start) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cerr << "[Thread " << thread_id << "] ";
                std::cerr << "got (" << kmer_id << ",'" << read_kmer << "')";
                std::cerr << " but ";
                std::cerr << "expected (" << start << ",'" << expected_kmer << "')" << std::endl;
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
    for (uint64_t t = 0, chunk_size = (num_kmers + num_threads - 1) / num_threads; t != num_threads;
         ++t) {
        uint64_t start = t * chunk_size;
        uint64_t end = std::min(num_kmers, start + chunk_size);
        threads.emplace_back(worker, start, end, t);
    }

    for (auto& th : threads) th.join();

    std::cout << "EVERYTHING OK!" << std::endl;

    return true;
}

template <typename Dict>
bool check_correctness_string_iterator(Dict const& dict) {
    const uint64_t k = dict.k();
    const uint64_t num_strings = dict.num_strings();

    const uint64_t num_threads = std::thread::hardware_concurrency();
    std::cout << "checking correctness of string iterator using " << num_threads << " threads..."
              << std::endl;

    std::mutex print_mutex;

    auto worker = [&](uint64_t start, uint64_t end, uint64_t thread_id) {
        std::string read_kmer(k, 0);
        std::string expected_kmer(k, 0);

        for (uint64_t string_id = start; string_id < end; ++string_id) {
            auto [begin, _] = dict.string_offsets(string_id);
            uint64_t from_kmer_id = begin - string_id * (dict.k() - 1);
            auto it = dict.at_string_id(string_id);

            while (it.has_next()) {
                auto [kmer_id, kmer] = it.next();
                util::uint_kmer_to_string<typename Dict::kmer_type>(kmer, read_kmer.data(), k);
                dict.access(kmer_id, expected_kmer.data());

                if (read_kmer != expected_kmer || kmer_id != from_kmer_id) {
                    std::lock_guard<std::mutex> lock(print_mutex);
                    std::cerr << "[Thread " << thread_id << "] ERROR at string_id " << string_id
                              << ": got (" << kmer_id << ", '" << read_kmer << "') but expected ("
                              << from_kmer_id << ", '" << expected_kmer << "')\n";
                    return;
                }
                ++from_kmer_id;
            }

            if ((string_id - start) % 1'000'000 == 0 && string_id != start) {
                std::lock_guard<std::mutex> lock(print_mutex);
                std::cout << "[Thread " << thread_id << "] checked " << (string_id - start)
                          << " strings (local progress)\n";
            }
        }

        std::lock_guard<std::mutex> lock(print_mutex);
        std::cout << "[Thread " << thread_id << "] Finished range [" << start << ", " << end
                  << ")\n";
    };

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (uint64_t t = 0, chunk_size = (num_strings + num_threads - 1) / num_threads;
         t < num_threads; ++t)  //
    {
        uint64_t start = t * chunk_size;
        uint64_t end = std::min(num_strings, start + chunk_size);
        if (start >= end) break;
        threads.emplace_back(worker, start, end, t);
    }

    for (auto& th : threads) th.join();

    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

}  // namespace sshash