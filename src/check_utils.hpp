#pragma once

#include "../include/gz/zip_stream.hpp"

namespace sshash {

bool check_correctness_lookup_access(std::istream& is, dictionary const& dict) {
    uint64_t k = dict.k();
    uint64_t n = dict.size();

    std::string line;
    uint64_t pos = 0;
    uint64_t num_kmers = 0;
    lookup_result prev;
    prev.contig_id = 0;

    std::string got_kmer_str(k, 0);
    std::string expected_kmer_str(k, 0);

    std::cout << "checking correctness of access and positive lookup..." << std::endl;

    while (appendline(is, line)) {
        if (line.size() == pos || line[pos] == '>' || line[pos] == ';') {
            // comment or empty line restart the term buffer
            line.clear();
            continue;
        }
        for (uint64_t i = 0; i + k <= line.size(); ++i) {
            assert(util::is_valid(line.data() + i, k));
            kmer_t uint_kmer = util::string_to_uint_kmer_no_reverse(line.data() + i, k);
            bool orientation = constants::forward_orientation;

            if (num_kmers != 0 and num_kmers % 5000000 == 0) {
                std::cout << "checked " << num_kmers << " kmers" << std::endl;
            }

            if ((num_kmers & 1) == 0) {
                /* transform 50% of the kmers into their reverse complements */
                uint_kmer = util::compute_reverse_complement(uint_kmer, k);
                orientation = constants::backward_orientation;
            }

            util::uint_kmer_to_string_no_reverse(uint_kmer, expected_kmer_str.data(), k);
            uint64_t id = dict.lookup(expected_kmer_str.c_str());

            /*
                Since we assume that we stream through the file from which the index was built,
                ids are assigned sequentially to kmers, so it must be id == num_kmers.
            */
            if (id != num_kmers) std::cout << "wrong id assigned" << std::endl;

            if (id == constants::invalid_uint64) {
                std::cout << "kmer '" << expected_kmer_str << "' not found!" << std::endl;
            }
            assert(id != constants::invalid_uint64);

            auto curr = dict.lookup_advanced(expected_kmer_str.c_str());
            assert(curr.kmer_id == id);

            if (curr.kmer_orientation != orientation) {
                std::cout << "ERROR: got orientation " << int(curr.kmer_orientation)
                          << " but expected " << int(orientation) << std::endl;
            }
            assert(curr.kmer_orientation == orientation);

            if (num_kmers == 0) {
                if (curr.contig_id != 0) {
                    std::cout << "contig_id " << curr.contig_id << " but expected 0" << std::endl;
                }
                assert(curr.contig_id == 0);  // at the beginning, contig_id must be 0
            } else {
                if (curr.kmer_id != prev.kmer_id + 1) {
                    std::cout << "ERROR: got curr.kmer_id " << curr.kmer_id << " but expected "
                              << prev.kmer_id + 1 << std::endl;
                }
                assert(curr.kmer_id == prev.kmer_id + 1);  // kmer_id must be sequential

                if (curr.kmer_id_in_contig >= curr.contig_size) {
                    std::cout << "ERROR: got curr.kmer_id_in_contig " << curr.kmer_id_in_contig
                              << " but expected something < " << curr.contig_size << std::endl;
                }
                assert(curr.kmer_id_in_contig <
                       curr.contig_size);  // kmer_id_in_contig must always be < contig_size

                if (curr.contig_id == prev.contig_id) {
                    /* same contig */
                    if (curr.contig_size != prev.contig_size) {
                        std::cout << "ERROR: got curr.contig_size " << curr.contig_size
                                  << " but expected " << prev.contig_size << std::endl;
                    }
                    assert(curr.contig_size == prev.contig_size);  // contig_size must be same
                    if (curr.kmer_id_in_contig != prev.kmer_id_in_contig + 1) {
                        std::cout << "ERROR: got curr.kmer_id_in_contig " << curr.kmer_id_in_contig
                                  << " but expected " << prev.kmer_id_in_contig + 1 << std::endl;
                    }
                    assert(curr.kmer_id_in_contig ==
                           prev.kmer_id_in_contig + 1);  // kmer_id_in_contig must be sequential
                } else {
                    /* we have changed contig */
                    if (curr.contig_id != prev.contig_id + 1) {
                        std::cout << "ERROR: got curr.contig_id " << curr.contig_id
                                  << " but expected " << prev.contig_id + 1 << std::endl;
                    }
                    assert(curr.contig_id ==
                           prev.contig_id + 1);  // contig_id must be sequential since we stream
                    if (curr.kmer_id_in_contig != 0) {
                        std::cout << "ERROR: got curr.kmer_id_in_contig " << curr.kmer_id_in_contig
                                  << " but expected 0" << std::endl;
                    }
                    assert(curr.kmer_id_in_contig ==
                           0);  // kmer_id_in_contig must be 0 when we change contig
                }
            }

            /* check also contig_size() */
            uint64_t contig_size = dict.contig_size(curr.contig_id);
            if (contig_size != curr.contig_size) {
                std::cout << "ERROR: got contig_size " << contig_size << " but expected "
                          << curr.contig_size << std::endl;
            }
            assert(contig_size == curr.contig_size);

            prev = curr;

            // check access
            dict.access(id, got_kmer_str.data());
            kmer_t got_uint_kmer = util::string_to_uint_kmer_no_reverse(got_kmer_str.data(), k);
            kmer_t got_uint_kmer_rc = util::compute_reverse_complement(got_uint_kmer, k);
            if (got_uint_kmer != uint_kmer and got_uint_kmer_rc != uint_kmer) {
                std::cout << "ERROR: got '" << got_kmer_str << "' but expected '"
                          << expected_kmer_str << "'" << std::endl;
            }
            ++num_kmers;
        }
        if (line.size() > k - 1) {
            std::copy(line.data() + line.size() - (k - 1), line.data() + line.size(), line.data());
            line.resize(k - 1);
            pos = line.size();
        } else {
            pos = 0;
        }
    }
    std::cout << "checked " << num_kmers << " kmers" << std::endl;

    std::cout << "EVERYTHING OK!" << std::endl;

    std::cout << "checking correctness of negative lookup with random kmers..." << std::endl;
    uint64_t num_lookups = std::min<uint64_t>(1000000, n);
    for (uint64_t i = 0; i != num_lookups; ++i) {
        random_kmer(got_kmer_str.data(), k);
        /*
            We could use a std::unordered_set to check if kmer is really absent,
            but that would take much more memory...
        */
        uint64_t id = dict.lookup(got_kmer_str.c_str());
        if (id != constants::invalid_uint64) {
            std::cout << "kmer '" << got_kmer_str << "' found!" << std::endl;
        }
    }

    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

bool check_correctness_navigational_kmer_query(std::istream& is, dictionary const& dict) {
    uint64_t k = dict.k();
    std::string line;
    uint64_t pos = 0;
    uint64_t num_kmers = 0;

    std::cout << "checking correctness of navigational queries for kmers..." << std::endl;
    while (appendline(is, line)) {
        if (line.size() == pos || line[pos] == '>' || line[pos] == ';') {
            // comment or empty line restart the term buffer
            line.clear();
            continue;
        }
        for (uint64_t i = 0; i + k <= line.size(); ++i) {
            assert(util::is_valid(line.data() + i, k));
            if (num_kmers != 0 and num_kmers % 5000000 == 0) {
                std::cout << "checked " << num_kmers << " kmers" << std::endl;
            }

            neighbourhood curr = dict.kmer_neighbours(line.data() + i);

            char next_nuc = line[i + k];
            switch (next_nuc) {
                case 'A':
                    if (curr.forward_A.kmer_id == constants::invalid_uint64) {
                        std::cout << "expected forward_A" << std::endl;
                    }
                    assert(curr.forward_A.kmer_id != constants::invalid_uint64);
                    break;
                case 'C':
                    if (curr.forward_C.kmer_id == constants::invalid_uint64) {
                        std::cout << "expected forward_C" << std::endl;
                    }
                    assert(curr.forward_C.kmer_id != constants::invalid_uint64);
                    break;
                case 'G':
                    if (curr.forward_G.kmer_id == constants::invalid_uint64) {
                        std::cout << "expected forward_G" << std::endl;
                    }
                    assert(curr.forward_G.kmer_id != constants::invalid_uint64);
                    break;
                case 'T':
                    if (curr.forward_T.kmer_id == constants::invalid_uint64) {
                        std::cout << "expected forward_T" << std::endl;
                    }
                    assert(curr.forward_T.kmer_id != constants::invalid_uint64);
                    break;
            }

            if (i != 0) {
                char prev_nuc = line[i - 1];
                switch (prev_nuc) {
                    case 'A':
                        if (curr.backward_A.kmer_id == constants::invalid_uint64) {
                            std::cout << "expected backward_A" << std::endl;
                        }
                        assert(curr.backward_A.kmer_id != constants::invalid_uint64);
                        break;
                    case 'C':
                        if (curr.backward_C.kmer_id == constants::invalid_uint64) {
                            std::cout << "expected backward_C" << std::endl;
                        }
                        assert(curr.backward_C.kmer_id != constants::invalid_uint64);
                        break;
                    case 'G':
                        if (curr.backward_G.kmer_id == constants::invalid_uint64) {
                            std::cout << "expected backward_G" << std::endl;
                        }
                        assert(curr.backward_G.kmer_id != constants::invalid_uint64);
                        break;
                    case 'T':
                        if (curr.backward_T.kmer_id == constants::invalid_uint64) {
                            std::cout << "expected backward_T" << std::endl;
                        }
                        assert(curr.backward_T.kmer_id != constants::invalid_uint64);
                        break;
                }
            }

            ++num_kmers;
        }
        if (line.size() > k - 1) {
            std::copy(line.data() + line.size() - (k - 1), line.data() + line.size(), line.data());
            line.resize(k - 1);
            pos = line.size();
        } else {
            pos = 0;
        }
    }
    std::cout << "checked " << num_kmers << " kmers" << std::endl;

    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

bool check_correctness_navigational_contig_query(dictionary const& dict) {
    std::cout << "checking correctness of navigational queries for contigs..." << std::endl;
    uint64_t num_contigs = dict.num_contigs();
    uint64_t k = dict.k();
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
        equal_lookup_result(backward.backward_A, res.backward_A);
        equal_lookup_result(backward.backward_C, res.backward_C);
        equal_lookup_result(backward.backward_G, res.backward_G);
        equal_lookup_result(backward.backward_T, res.backward_T);

        uint64_t end_kmer_id = kmer_id + contig_size - 1;
        dict.access(end_kmer_id, kmer.data());
        auto forward = dict.kmer_forward_neighbours(kmer.data());
        equal_lookup_result(forward.forward_A, res.forward_A);
        equal_lookup_result(forward.forward_C, res.forward_C);
        equal_lookup_result(forward.forward_G, res.forward_G);
        equal_lookup_result(forward.forward_T, res.forward_T);

        kmer_id += contig_size;
    }
    std::cout << "checked " << contig_id << " contigs" << std::endl;
    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

bool check_correctness_weights(std::istream& is, dictionary const& dict) {
    uint64_t k = dict.k();
    std::string line;
    uint64_t kmer_id = 0;

    if (!dict.weighted()) {
        std::cerr << "ERROR: the dictionary does not store weights" << std::endl;
        return false;
    }

    std::cout << "checking correctness of weights..." << std::endl;

    while (!is.eof()) {
        std::getline(is, line);  // header line
        if (line.empty()) break;

        uint64_t i = 0;
        i = line.find_first_of(' ', i);
        assert(i != std::string::npos);

        i += 1;
        i += 5;  // skip "LN:i:"
        uint64_t j = line.find_first_of(' ', i);
        assert(j != std::string::npos);

        char* end;
        uint64_t seq_len = std::strtoull(line.data() + i, &end, 10);
        i = j + 1;
        i += 5;  // skip "ab:Z:"

        for (uint64_t j = 0; j != seq_len - k + 1; ++j, ++kmer_id) {
            uint64_t expected = std::strtoull(line.data() + i, &end, 10);
            i = line.find_first_of(' ', i) + 1;
            uint64_t got = dict.weight(kmer_id);
            if (expected != got) {
                std::cout << "ERROR for kmer_id " << kmer_id << ": expected " << expected
                          << " but got " << got << std::endl;
            }
            if (kmer_id != 0 and kmer_id % 5000000 == 0) {
                std::cout << "checked " << kmer_id << " weights" << std::endl;
            }
        }

        std::getline(is, line);  // skip DNA sequence
    }

    std::cout << "checked " << kmer_id << " weights" << std::endl;
    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

/*
   The input file must be the one the index was built from.
   Throughout the code, we assume the input does not contain any duplicate.
*/
bool check_correctness_lookup_access(dictionary const& dict, std::string const& filename) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    bool good = true;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        good = check_correctness_lookup_access(zis, dict);
    } else {
        good = check_correctness_lookup_access(is, dict);
    }
    is.close();
    return good;
}

/*
   The input file must be the one the index was built from.
   Throughout the code, we assume the input does not contain any duplicate.
*/
bool check_correctness_navigational_kmer_query(dictionary const& dict,
                                               std::string const& filename) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    bool good = true;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        good = check_correctness_navigational_kmer_query(zis, dict);
    } else {
        good = check_correctness_navigational_kmer_query(is, dict);
    }
    is.close();
    return good;
}

/*
   The input file must be the one the index was built from.
*/
bool check_correctness_weights(dictionary const& dict, std::string const& filename) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    bool good = true;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        good = check_correctness_weights(zis, dict);
    } else {
        good = check_correctness_weights(is, dict);
    }
    is.close();
    return good;
}

bool check_dictionary(dictionary const& dict) {
    uint64_t k = dict.k();
    uint64_t n = dict.size();
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
    return true;
}

bool check_correctness_iterator(dictionary const& dict) {
    std::cout << "checking correctness of iterator..." << std::endl;
    std::string expected_kmer(dict.k(), 0);
    constexpr uint64_t runs = 3;
    essentials::uniform_int_rng<uint64_t> distr(0, dict.size() - 1, essentials::get_random_seed());
    for (uint64_t run = 0; run != runs; ++run) {
        uint64_t from_kmer_id = distr.gen();
        auto it = dict.at(from_kmer_id);
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

}  // namespace sshash