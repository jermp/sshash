#pragma once

#include "external/gz/zip_stream.hpp"
#include "include/minimizer_iterator.hpp"

namespace sshash {

template <class kmer_t>
struct parse_data {
    parse_data(build_configuration const& build_config) : num_kmers(0), minimizers(build_config) {}

    uint64_t num_kmers;
    minimizers_tuples minimizers;
    compact_string_pool<kmer_t> strings;
    weights::builder weights_builder;
};

template <class kmer_t, input_file_type fmt>
void parse_file(std::istream& is, parse_data<kmer_t>& data,
                build_configuration const& build_config)  //
{
    essentials::timer_type timer;
    timer.start();

    const uint64_t k = build_config.k;
    const uint64_t m = build_config.m;
    const uint64_t max_num_kmers_in_super_kmer = k - m + 1;
    assert(k > 0 and m <= k);

    if (max_num_kmers_in_super_kmer >= (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8))) {
        throw std::runtime_error(
            "max_num_kmers_in_super_kmer " + std::to_string(max_num_kmers_in_super_kmer) +
            " does not fit into " + std::to_string(sizeof(num_kmers_in_super_kmer_uint_type) * 8) +
            " bits");
    }

    /* fit into the wanted number of bits */
    assert(max_num_kmers_in_super_kmer < (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8)));

    typename compact_string_pool<kmer_t>::builder builder;

    std::string sequence;
    uint64_t num_sequences = 0;
    uint64_t num_bases = 0;

    hasher_type hasher(build_config.seed);
    minimizer_iterator<kmer_t> minimizer_it(k, m, hasher);
    minimizer_iterator_rc<kmer_t> minimizer_it_rc(k, m, hasher);
    uint64_t seq_len = 0;
    uint64_t sum_of_weights = 0;
    data.weights_builder.init();

    /* intervals of weights */
    uint64_t weight_value = constants::invalid_uint64;
    uint64_t weight_length = 0;

    while (true)  //
    {
        if constexpr (fmt == input_file_type::cf_seg) {
            std::getline(is, sequence, '\t');  // skip '\t'
            std::getline(is, sequence);        // DNA sequence
        } else {
            static_assert(fmt == input_file_type::fasta);
            std::getline(is, sequence);   // header sequence
            if (build_config.weighted) {  // parse header
                if (sequence.empty()) return;

                /*
                    Heder format:
                    >[id] LN:i:[seq_len] ab:Z:[weight_seq]
                    where [weight_seq] is a space-separated sequence of integer counters (the
                   weights), whose length is equal to [seq_len]-k+1
                */

                // example header: '>12 LN:i:41 ab:Z:2 2 2 2 2 2 2 2 2 2 2'

                expect(sequence[0], '>');
                uint64_t i = 0;
                i = sequence.find_first_of(' ', i);
                if (i == std::string::npos) throw parse_runtime_error();

                i += 1;
                expect(sequence[i + 0], 'L');
                expect(sequence[i + 1], 'N');
                expect(sequence[i + 2], ':');
                expect(sequence[i + 3], 'i');
                expect(sequence[i + 4], ':');
                i += 5;
                uint64_t j = sequence.find_first_of(' ', i);
                if (j == std::string::npos) throw parse_runtime_error();

                seq_len = std::strtoull(sequence.data() + i, nullptr, 10);
                i = j + 1;
                expect(sequence[i + 0], 'a');
                expect(sequence[i + 1], 'b');
                expect(sequence[i + 2], ':');
                expect(sequence[i + 3], 'Z');
                expect(sequence[i + 4], ':');
                i += 5;

                for (uint64_t j = 0; j != seq_len - k + 1; ++j) {
                    uint64_t weight = std::strtoull(sequence.data() + i, nullptr, 10);
                    i = sequence.find_first_of(' ', i) + 1;

                    data.weights_builder.eat(weight);
                    sum_of_weights += weight;

                    if (weight == weight_value) {
                        weight_length += 1;
                    } else {
                        if (weight_value != constants::invalid_uint64) {
                            data.weights_builder.push_weight_interval(weight_value, weight_length);
                        }
                        weight_value = weight;
                        weight_length = 1;
                    }
                }
            }
            std::getline(is, sequence);  // DNA sequence
        }

        if (is.eof()) {
            assert(sequence.empty());
            break;
        }

        assert(sequence.length() >= k);

        if (++num_sequences % 100000 == 0) {
            std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
                      << data.num_kmers << " kmers" << std::endl;
        }

        builder.new_piece();
        num_bases += sequence.length();

        if (build_config.weighted and seq_len != sequence.length()) {
            std::cout << "ERROR: expected a sequence of length " << seq_len
                      << " but got one of length " << sequence.length() << std::endl;
            throw std::runtime_error("file is malformed");
        }

        data.num_kmers += sequence.length() - k + 1;
        for (uint64_t i = 0; i != sequence.length(); ++i) {
            assert(kmer_t::is_valid(sequence[i]));
            builder.append(sequence[i]);
        }
    }

    builder.finalize();
    builder.build(data.strings);

    assert(data.strings.pieces.front() == 0);
    assert(data.strings.pieces.size() == num_sequences + 1);

    timer.stop();
    print_time(timer.elapsed(), data.num_kmers, "step 1.1: 'encoding_input'");

    std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
              << data.num_kmers << " kmers" << std::endl;
    std::cout << "num_kmers " << data.num_kmers << std::endl;
    std::cout << "cost: 2.0 + "
              << static_cast<double>(kmer_t::bits_per_char * num_sequences * (k - 1)) /
                     data.num_kmers
              << " [bits/kmer]" << std::endl;

    timer.reset();
    timer.start();

    /*
        TODO:
        We could distribute the sequences to multiple threads
        and do this computation in parallel.
    */
    for (uint64_t i = 0; i != num_sequences; ++i) {
        const uint64_t begin = data.strings.pieces[i];
        const uint64_t end = data.strings.pieces[i + 1];
        const uint64_t sequence_len = end - begin;
        assert(end > begin);

        kmer_iterator<kmer_t> it(data.strings.strings, k, kmer_t::bits_per_char * begin);
        minimizer_info prev_mini_info;
        assert(prev_mini_info.minimizer == constants::invalid_uint64);
        bool prev_mini_rc_is_min = false;
        uint64_t num_kmers_in_super_kmer = 0;
        minimizer_it.set_position(begin);
        minimizer_it_rc.set_position(begin);

        for (uint64_t j = 0; j != sequence_len - k + 1; ++j) {
            auto uint_kmer = it.get();
            auto mini_info = minimizer_it.next(uint_kmer);
            assert(mini_info.pos_in_seq < end - m + 1);
            assert(mini_info.pos_in_kmer < k - m + 1);

            bool mini_rc_is_min = false;
            if (build_config.canonical) {
                auto uint_kmer_rc = uint_kmer;
                uint_kmer_rc.reverse_complement_inplace(k);
                auto mini_info_rc = minimizer_it_rc.next(uint_kmer_rc);
                assert(mini_info_rc.pos_in_seq < end - m + 1);
                assert(mini_info_rc.pos_in_kmer < k - m + 1);
                if (mini_info_rc.minimizer < mini_info.minimizer) {
                    mini_info = mini_info_rc;
                    mini_info.pos_in_kmer = k - m - mini_info.pos_in_kmer;
                    mini_rc_is_min = true;
                }
            }

            if (prev_mini_info.minimizer == constants::invalid_uint64) {
                prev_mini_info = mini_info;
                prev_mini_rc_is_min = mini_rc_is_min;
            }

            if (mini_info.minimizer != prev_mini_info.minimizer or
                mini_info.pos_in_seq != prev_mini_info.pos_in_seq) {
                assert(num_kmers_in_super_kmer <= max_num_kmers_in_super_kmer);
                data.minimizers.emplace_back(prev_mini_info, num_kmers_in_super_kmer,
                                             prev_mini_rc_is_min);
                prev_mini_info = mini_info;
                prev_mini_rc_is_min = mini_rc_is_min;
                num_kmers_in_super_kmer = 0;
            }

            num_kmers_in_super_kmer += 1;
            it.next();
        }

        assert(num_kmers_in_super_kmer <= max_num_kmers_in_super_kmer);
        data.minimizers.emplace_back(prev_mini_info, num_kmers_in_super_kmer, prev_mini_rc_is_min);
    }

    data.minimizers.finalize();

    timer.stop();
    print_time(timer.elapsed(), data.num_kmers, "step 1.2: 'computing_minimizers_tuples'");

    if (build_config.weighted) {
        std::cout << "sum_of_weights " << sum_of_weights << std::endl;
        data.weights_builder.push_weight_interval(weight_value, weight_length);
        data.weights_builder.finalize(data.num_kmers);
    }
}

template <class kmer_t>
parse_data<kmer_t> parse_file(std::string const& filename,
                              build_configuration const& build_config) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    std::cout << "reading file '" << filename << "'..." << std::endl;
    parse_data<kmer_t> data(build_config);
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        if (util::ends_with(filename, ".cf_seg.gz")) {
            parse_file<kmer_t, input_file_type::cf_seg>(zis, data, build_config);
        } else {
            parse_file<kmer_t, input_file_type::fasta>(zis, data, build_config);
        }
    } else {
        if (util::ends_with(filename, ".cf_seg")) {
            parse_file<kmer_t, input_file_type::cf_seg>(is, data, build_config);
        } else {
            parse_file<kmer_t, input_file_type::fasta>(is, data, build_config);
        }
    }
    is.close();
    return data;
}

}  // namespace sshash