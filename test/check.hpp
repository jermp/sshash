#pragma once

namespace sshash {

template <class kmer_t>
bool check_dictionary(dictionary<kmer_t> const& dict) {
    const uint64_t k = dict.k();
    const uint64_t n = dict.size();
    std::cout << "checking correctness of access and positive lookup..." << std::endl;
    uint64_t id = 0;
    std::string kmer(k, 0);
    for (; id != n; ++id) {
        if (id != 0 and id % 5000000 == 0) std::cout << "checked " << id << " kmers" << std::endl;
        dict.access(id, kmer.data());
        uint64_t got_id = dict.lookup(kmer.c_str());
        if (got_id == constants::invalid_uint64) {
            std::cout << "kmer '" << kmer << "' not found!" << std::endl;
            return false;
        }
        if (got_id >= n) {
            std::cout << "ERROR: id out of range " << got_id << "/" << n << std::endl;
            return false;
        }
        if (got_id != id) {
            std::cout << "expected id " << id << " but got id " << got_id << std::endl;
            return false;
        }
    }
    std::cout << "checked " << id << " kmers" << std::endl;
    std::cout << "EVERYTHING OK!" << std::endl;
    return check_correctness_negative_lookup(dict);
}

template <class kmer_t>
bool check_correctness_negative_lookup(dictionary<kmer_t> const& dict) {
    std::cout << "checking correctness of negative lookup with random kmers..." << std::endl;
    const uint64_t num_lookups = std::min<uint64_t>(1000000, dict.size());
    std::string kmer(dict.k(), 0);
    for (uint64_t i = 0; i != num_lookups; ++i) {
        random_kmer(kmer.data(), dict.k());
        /*
            We could use a std::unordered_set to check if kmer is really absent,
            but that would take much more memory...
        */
        uint64_t id = dict.lookup(kmer.c_str());
        if (id != constants::invalid_uint64) {
            std::cout << "kmer '" << kmer << "' found!" << std::endl;
        }
    }
    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

template <class kmer_t>
bool check_correctness_navigational_contig_query(dictionary<kmer_t> const& dict) {
    std::cout << "checking correctness of navigational queries for contigs..." << std::endl;
    const uint64_t num_contigs = dict.num_contigs();
    const uint64_t k = dict.k();
    uint64_t kmer_id = 0;
    std::string kmer(k, 0);
    uint64_t contig_id = 0;
    for (; contig_id != num_contigs; ++contig_id) {
        if (contig_id != 0 and contig_id % 1000000 == 0) {
            std::cout << "checked " << contig_id << "/" << num_contigs << " contigs" << std::endl;
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

template <class kmer_t>
bool check_correctness_kmer_iterator(dictionary<kmer_t> const& dict) {
    std::cout << "checking correctness of kmer iterator..." << std::endl;
    std::string expected_kmer(dict.k(), 0);
    constexpr uint64_t runs = 3;
    essentials::uniform_int_rng<uint64_t> distr(0, dict.size() - 1, essentials::get_random_seed());
    for (uint64_t run = 0; run != runs; ++run) {
        uint64_t from_kmer_id = distr.gen();
        auto it = dict.at_kmer_id(from_kmer_id);
        while (it.has_next()) {
            auto [kmer_id, kmer] = it.next();
            dict.access(kmer_id, expected_kmer.data());
            if (kmer != expected_kmer or kmer_id != from_kmer_id) {
                std::cout << "got (" << kmer_id << ",'" << kmer << "')";
                std::cout << " but ";
                std::cout << "expected (" << from_kmer_id << ",'" << expected_kmer << "')"
                          << std::endl;
                return false;
            }
            ++from_kmer_id;
        }
        assert(from_kmer_id == dict.size());
    }
    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

template <class kmer_t>
bool check_correctness_contig_iterator(dictionary<kmer_t> const& dict) {
    std::cout << "checking correctness of contig iterator..." << std::endl;
    std::string expected_kmer(dict.k(), 0);
    for (uint64_t contig_id = 0; contig_id != dict.num_contigs(); ++contig_id) {
        auto [begin, _] = dict.contig_offsets(contig_id);
        uint64_t from_kmer_id = begin - contig_id * (dict.k() - 1);
        auto it = dict.at_contig_id(contig_id);
        while (it.has_next()) {
            auto [kmer_id, kmer] = it.next();
            dict.access(kmer_id, expected_kmer.data());
            if (kmer != expected_kmer or kmer_id != from_kmer_id) {
                std::cout << "got (" << kmer_id << ",'" << kmer << "')";
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