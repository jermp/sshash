#pragma once

#include <algorithm>  // for std::transform

#include "external/gz/zip_stream.hpp"

namespace sshash {

template <typename Dict, input_file_t fmt>
bool check_correctness_lookup_access(std::istream& is, Dict const& dict)  //
{
    using kmer_t = typename Dict::kmer_type;
    const uint64_t k = dict.k();
    std::string sequence;
    uint64_t num_kmers = 0;
    uint64_t num_sequences = 0;
    lookup_result prev;
    prev.string_id = 0;

    std::string got_kmer_str(k, 0);
    std::string expected_kmer_str(k, 0);

    std::cout << "checking correctness of access and positive lookup..." << std::endl;

    while (!is.eof())  //
    {
        if constexpr (fmt == input_file_t::cf_seg) {
            std::getline(is, sequence, '\t');  // skip '\t'
            std::getline(is, sequence);        // DNA sequence
        } else {
            static_assert(fmt == input_file_t::fasta);
            std::getline(is, sequence);  // header sequence
            std::getline(is, sequence);  // DNA sequence
        }

        if (sequence.length() < k) continue;

        /* transform 50% of the read nucleotides into lower-case letters
           (assuming the input is upper-case):
           lower-case kmers must be found anyway in the index */
        if ((num_sequences & 1) == 0) {
            std::transform(sequence.begin(), sequence.end(), sequence.begin(),
                           [](char c) { return std::tolower(c); });
        }
        ++num_sequences;

        for (uint64_t i = 0; i + k <= sequence.length(); ++i) {
            assert(util::is_valid<kmer_t>(sequence.data() + i, k));

            // std::cout << "kmer = '" << std::string(sequence.data() + i, k) << "'" << std::endl;

            kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(sequence.data() + i, k);
            auto orientation = constants::forward_orientation;

            if (num_kmers != 0 and num_kmers % 5000000 == 0) {
                std::cout << "checked " << num_kmers << " kmers" << std::endl;
            }

            /* transform 50% of the kmers into their reverse complements */
            if ((num_kmers & 1) == 0) {
                uint_kmer.reverse_complement_inplace(k);
                orientation = constants::backward_orientation;
            }

            util::uint_kmer_to_string(uint_kmer, expected_kmer_str.data(), k);
            auto curr = dict.lookup(expected_kmer_str.c_str());

            /*
                Since we assume that we stream through the file from which the index was built,
                ids are assigned sequentially to kmers, so it must be curr.kmer_id == num_kmers.
            */
            if (curr.kmer_id != num_kmers) std::cout << "wrong id assigned" << std::endl;

            if (curr.kmer_id == constants::invalid_uint64) {
                std::cout << "kmer '" << expected_kmer_str << "' not found!" << std::endl;
            }
            assert(curr.kmer_id != constants::invalid_uint64);

            if (curr.kmer_orientation != orientation) {
                std::cout << "ERROR: got orientation " << int(curr.kmer_orientation)
                          << " but expected " << int(orientation) << std::endl;
            }
            assert(curr.kmer_orientation == orientation);

            uint64_t curr_string_size = curr.string_end - curr.string_begin - k + 1;

            if (num_kmers == 0) {
                if (curr.string_id != 0) {
                    std::cout << "string_id " << curr.string_id << " but expected 0" << std::endl;
                }
                assert(curr.string_id == 0);  // at the beginning, string_id must be 0
            } else {
                if (curr.kmer_id != prev.kmer_id + 1) {
                    std::cout << "ERROR: got curr.kmer_id " << curr.kmer_id << " but expected "
                              << prev.kmer_id + 1 << std::endl;
                }
                assert(curr.kmer_id == prev.kmer_id + 1);  // kmer_id must be sequential

                if (curr.kmer_id_in_string >= curr_string_size) {
                    std::cout << "ERROR: got curr.kmer_id_in_string " << curr.kmer_id_in_string
                              << " but expected something < " << curr_string_size << std::endl;
                }
                assert(curr.kmer_id_in_string <
                       curr_string_size);  // kmer_id_in_string must always be < string_size

                if (curr.string_id == prev.string_id) {
                    /* same string */
                    uint64_t prev_string_size = prev.string_end - prev.string_begin - k + 1;
                    if (curr_string_size != prev_string_size) {
                        std::cout << "ERROR: got curr_string_size " << curr_string_size
                                  << " but expected " << prev_string_size << std::endl;
                    }
                    if (curr.kmer_id_in_string != prev.kmer_id_in_string + 1) {
                        std::cout << "ERROR: got curr.kmer_id_in_string " << curr.kmer_id_in_string
                                  << " but expected " << prev.kmer_id_in_string + 1 << std::endl;
                    }
                    assert(curr.kmer_id_in_string ==
                           prev.kmer_id_in_string + 1);  // kmer_id_in_string must be sequential
                } else {
                    /* we have changed string */
                    if (curr.string_id != prev.string_id + 1) {
                        std::cout << "ERROR: got curr.string_id " << curr.string_id
                                  << " but expected " << prev.string_id + 1 << std::endl;
                    }
                    assert(curr.string_id ==
                           prev.string_id + 1);  // string_id must be sequential since we stream
                    if (curr.kmer_id_in_string != 0) {
                        std::cout << "ERROR: got curr.kmer_id_in_string " << curr.kmer_id_in_string
                                  << " but expected 0" << std::endl;
                    }
                    assert(curr.kmer_id_in_string ==
                           0);  // kmer_id_in_string must be 0 when we change string
                }
            }

            /* check also string_size() */
            uint64_t string_size = dict.string_size(curr.string_id);
            if (string_size != curr_string_size) {
                std::cout << "ERROR: got string_size " << string_size << " but expected "
                          << curr_string_size << std::endl;
            }
            assert(string_size == curr_string_size);

            prev = curr;

            // check access
            dict.access(curr.kmer_id, got_kmer_str.data());
            kmer_t got_uint_kmer = util::string_to_uint_kmer<kmer_t>(got_kmer_str.data(), k);
            kmer_t got_uint_kmer_rc = got_uint_kmer;
            got_uint_kmer_rc.reverse_complement_inplace(k);
            if (got_uint_kmer != uint_kmer and got_uint_kmer_rc != uint_kmer) {
                std::cout << "ERROR: got '" << got_kmer_str << "' but expected '"
                          << expected_kmer_str << "'" << std::endl;
                return false;
            }

            ++num_kmers;
        }
    }
    std::cout << "checked " << num_kmers << " kmers" << std::endl;
    std::cout << "EVERYTHING OK!" << std::endl;

    return check_correctness_negative_lookup(dict);
}

template <typename Dict, input_file_t fmt>
bool check_correctness_navigational_kmer_query(std::istream& is, Dict const& dict)  //
{
    using kmer_t = typename Dict::kmer_type;
    const uint64_t k = dict.k();
    std::string sequence;
    uint64_t num_kmers = 0;

    std::cout << "checking correctness of navigational queries for kmers..." << std::endl;

    while (!is.eof())  //
    {
        if constexpr (fmt == input_file_t::cf_seg) {
            std::getline(is, sequence, '\t');  // skip '\t'
            std::getline(is, sequence);        // DNA sequence
        } else {
            static_assert(fmt == input_file_t::fasta);
            std::getline(is, sequence);  // header sequence
            std::getline(is, sequence);  // DNA sequence
        }
        for (uint64_t i = 0; i + k <= sequence.length(); ++i) {
            assert(util::is_valid<kmer_t>(sequence.data() + i, k));
            if (num_kmers != 0 and num_kmers % 5000000 == 0) {
                std::cout << "checked " << num_kmers << " kmers" << std::endl;
            }
            neighbourhood<kmer_t> curr = dict.kmer_neighbours(sequence.data() + i);
            if (i + k < sequence.length()) {
                char next_nuc = sequence[i + k];
                bool next_nuc_not_found = curr.forward[kmer_t::char_to_uint(next_nuc)].kmer_id ==
                                          constants::invalid_uint64;
                if (next_nuc_not_found) {
                    std::cout << "expected forward[" << next_nuc << "]" << std::endl;
                }
                assert(!next_nuc_not_found);
            }

            if (i != 0) {
                char prev_nuc = sequence[i - 1];
                bool prev_nuc_not_found = curr.backward[kmer_t::char_to_uint(prev_nuc)].kmer_id ==
                                          constants::invalid_uint64;
                if (prev_nuc_not_found) {
                    std::cout << "expected backward[" << prev_nuc << "]" << std::endl;
                }
                assert(!prev_nuc_not_found);
            }

            ++num_kmers;
        }
    }
    std::cout << "checked " << num_kmers << " kmers" << std::endl;

    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

template <typename Dict>
bool check_correctness_weights(std::istream& is, Dict const& dict) {
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
template <typename Dict>
bool check_correctness_lookup_access(Dict const& dict, std::string const& filename)  //
{
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    bool good = true;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        if (util::ends_with(filename, ".cf_seg.gz")) {
            good = check_correctness_lookup_access<Dict, input_file_t::cf_seg>(zis, dict);
        } else {
            good = check_correctness_lookup_access<Dict, input_file_t::fasta>(zis, dict);
        }
    } else {
        if (util::ends_with(filename, ".cf_seg")) {
            good = check_correctness_lookup_access<Dict, input_file_t::cf_seg>(is, dict);
        } else {
            good = check_correctness_lookup_access<Dict, input_file_t::fasta>(is, dict);
        }
    }
    is.close();
    return good;
}

/*
   The input file must be the one the index was built from.
   Throughout the code, we assume the input does not contain any duplicate.
*/
template <typename Dict>
bool check_correctness_navigational_kmer_query(Dict const& dict, std::string const& filename) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    bool good = true;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        if (util::ends_with(filename, ".cf_seg.gz")) {
            good = check_correctness_navigational_kmer_query<Dict, input_file_t::cf_seg>(zis, dict);
        } else {
            good = check_correctness_navigational_kmer_query<Dict, input_file_t::fasta>(zis, dict);
        }
    } else {
        if (util::ends_with(filename, ".cf_seg")) {
            good = check_correctness_navigational_kmer_query<Dict, input_file_t::cf_seg>(is, dict);
        } else {
            good = check_correctness_navigational_kmer_query<Dict, input_file_t::fasta>(is, dict);
        }
    }
    is.close();
    return good;
}

/*
   The input file must be the one the index was built from.
   Only for FASTA files since CUTTLEFISH does not store kmer weights.
*/
template <typename Dict>
bool check_correctness_weights(Dict const& dict, std::string const& filename) {
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

}  // namespace sshash