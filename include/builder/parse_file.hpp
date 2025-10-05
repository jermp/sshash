#pragma once

#include "external/gz/zip_stream.hpp"
#include "include/minimizer_iterator.hpp"

namespace sshash {

template <class kmer_t>
struct parse_data {
    parse_data(build_configuration const& build_config) : num_kmers(0), minimizers(build_config) {}

    uint64_t num_kmers;
    minimizers_tuples minimizers;
    // compact_string_pool<kmer_t> strings;
    std::vector<uint64_t> pieces;
    bits::bit_vector strings;
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
    assert(k > 0 and k >= m);
    const uint64_t max_num_kmers_in_super_kmer = k - m + 1;

    if (max_num_kmers_in_super_kmer >= (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8))) {
        throw std::runtime_error(
            "max_num_kmers_in_super_kmer " + std::to_string(max_num_kmers_in_super_kmer) +
            " does not fit into " + std::to_string(sizeof(num_kmers_in_super_kmer_uint_type) * 8) +
            " bits");
    }

    /* fit into the wanted number of bits */
    assert(max_num_kmers_in_super_kmer < (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8)));

    const uint64_t num_bits = 8 * 8 * essentials::GB;  // reserve 8 GB of memory
    bits::bit_vector::builder bvb_strings;
    bvb_strings.reserve(num_bits);

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
            std::getline(is, sequence, '\t');  // skip until and consume '\t'
        } else {
            static_assert(fmt == input_file_type::fasta);
            if (build_config.weighted) {     // parse header
                std::getline(is, sequence);  // header sequence
                if (sequence.empty()) break;

                /*
                    Heder format:
                    >[id] LN:i:[seq_len] ab:Z:[weight_seq]
                    where [weight_seq] is a space-separated sequence of integer counters
                    (the weights), whose length is equal to [seq_len]-k+1.
                    Example: '>12 LN:i:41 ab:Z:2 2 2 2 2 2 2 2 2 2 2'
                */

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
            } else {
                // skip header sequence
                is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }

        std::getline(is, sequence);  // DNA sequence

        if (is.eof()) break;

        const uint64_t n = sequence.length();
        assert(n >= k);

        ++num_sequences;
        if (num_sequences % 100000 == 0) {
            std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
                      << data.num_kmers << " kmers" << std::endl;
        }

        assert(bvb_strings.num_bits() % kmer_t::bits_per_char == 0);
        data.pieces.push_back(bvb_strings.num_bits() / kmer_t::bits_per_char);

        num_bases += n;

        if (build_config.weighted and seq_len != n) {
            std::cout << "ERROR: expected a sequence of length " << seq_len
                      << " but got one of length " << n << std::endl;
            throw std::runtime_error("file is malformed");
        }

        data.num_kmers += n - k + 1;

        uint64_t i = 0;
        if constexpr (kmer_t::bits_per_char == 2) {
#if !defined(SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING) and defined(__AVX2__)

            /* process 32 bytes at a time */
            for (; i + 32 <= n; i += 32) {
                __m256i v = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(&sequence[i]));
                uint64_t word = pack2bits_shift1(v);
                bvb_strings.append_bits(word, 64);
            }
#endif
        }
        for (; i < n; ++i) {
            bvb_strings.append_bits(kmer_t::char_to_uint(sequence[i]), kmer_t::bits_per_char);
        }
    }

    /*
        So pieces will be of size p+1, where p is the number of DNA sequences
        in the input file.
    */
    data.pieces.push_back(bvb_strings.num_bits() / kmer_t::bits_per_char);
    assert(data.pieces.front() == 0);

    /* Push a final sentinel (dummy) value to avoid bounds' checking in
       kmer_iterator::fill_buff(). */
    bvb_strings.append_bits(0, kmer_t::uint_kmer_bits);

    bvb_strings.build(data.strings);

    assert(data.pieces.front() == 0);
    assert(data.pieces.size() == num_sequences + 1);

    timer.stop();
    print_time(timer.elapsed(), data.num_kmers, "step 1.1: 'encoding_input'");

    std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
              << data.num_kmers << " kmers" << std::endl;
    std::cout << "num_kmers " << data.num_kmers << std::endl;
    std::cout << "cost: 2.0 + "
              << static_cast<double>(kmer_t::bits_per_char * num_sequences * (k - 1)) /
                     data.num_kmers
              << " [bits/kmer]" << std::endl;

    /*
        The parameter m (minimizer length) should be at least
            ceil(log_s(N))+1
        where N is the number of nucleotides in the input and s is the alphabet size.
        We warn the user if the used m is less than this lower bound.
    */
    const uint64_t s = uint64_t(1) << kmer_t::bits_per_char;
    const uint64_t lower_bound_on_m = std::ceil(std::log(num_bases) / std::log(s)) + 1;
    if (m < lower_bound_on_m) {
        std::cout << "\n--> WARNING: using minimizer length " << m
                  << " but it should be at least ceil(log_" << s << "(" << num_bases
                  << "))+1 = " << lower_bound_on_m << '\n'
                  << std::endl;
    }

    timer.reset();
    timer.start();

    const uint64_t num_threads = build_config.num_threads;
    const uint64_t num_sequences_per_thread = (num_sequences + num_threads - 1) / num_threads;
    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    for (uint64_t t = 0; t * num_sequences_per_thread < num_sequences; ++t)  //
    {
        threads.emplace_back([&, t] {
            std::vector<minimizer_tuple> buffer;
            const uint64_t buffer_size = (build_config.ram_limit_in_GiB * essentials::GiB) /
                                         (2 * sizeof(minimizer_tuple) * num_threads);
            buffer.reserve(buffer_size);

            auto save = [&](minimizer_info mini_info,
                            uint64_t num_kmers_in_super_kmer)  //
            {
                assert(num_kmers_in_super_kmer <= max_num_kmers_in_super_kmer);
                if (!buffer.empty() and                                   //
                    buffer.back().minimizer == mini_info.minimizer and    //
                    buffer.back().pos_in_seq == mini_info.pos_in_seq and  //
                    buffer.back().pos_in_kmer == mini_info.pos_in_kmer)   //
                {
                    buffer.back().num_kmers_in_super_kmer += num_kmers_in_super_kmer;
                    return;
                }
                if (buffer.size() == buffer_size) {
                    data.minimizers.sort_and_flush(buffer);
                    buffer.clear();
                }
                buffer.emplace_back(mini_info, num_kmers_in_super_kmer);
            };

            const uint64_t index_begin = t * num_sequences_per_thread;
            const uint64_t index_end =
                std::min<uint64_t>(index_begin + num_sequences_per_thread, num_sequences);

            kmer_iterator<kmer_t> it(data.strings, k);
            minimizer_iterator<kmer_t> minimizer_it(k, m, hasher);
            minimizer_iterator_rc<kmer_t> minimizer_it_rc(k, m, hasher);

            for (uint64_t i = index_begin; i < index_end; ++i)  //
            {
                const uint64_t begin = data.pieces[i];
                const uint64_t end = data.pieces[i + 1];
                const uint64_t sequence_len = end - begin;
                assert(sequence_len >= k);

                minimizer_info prev_mini_info;
                assert(prev_mini_info.minimizer == constants::invalid_uint64);
                uint64_t num_kmers_in_super_kmer = 0;

                it.at(kmer_t::bits_per_char * begin);
                minimizer_it.set_position(begin);
                minimizer_it_rc.set_position(begin);

                for (uint64_t j = 0; j != sequence_len - k + 1; ++j) {
                    auto uint_kmer = it.get();
                    auto mini_info = minimizer_it.next(uint_kmer);
                    assert(mini_info.pos_in_seq < end - m + 1);
                    assert(mini_info.pos_in_kmer < k - m + 1);

                    if (build_config.canonical) {
                        auto uint_kmer_rc = uint_kmer;
                        uint_kmer_rc.reverse_complement_inplace(k);
                        auto mini_info_rc = minimizer_it_rc.next(uint_kmer_rc);
                        assert(mini_info_rc.pos_in_seq < end - m + 1);
                        assert(mini_info_rc.pos_in_kmer < k - m + 1);
                        if (mini_info_rc.minimizer < mini_info.minimizer) {
                            mini_info = mini_info_rc;
                            mini_info.pos_in_kmer = k - m - mini_info.pos_in_kmer;
                        }
                    }

                    if (prev_mini_info.minimizer == constants::invalid_uint64) {
                        prev_mini_info = mini_info;
                    }

                    if (mini_info.minimizer != prev_mini_info.minimizer or
                        mini_info.pos_in_seq != prev_mini_info.pos_in_seq)  //
                    {
                        save(prev_mini_info, num_kmers_in_super_kmer);
                        prev_mini_info = mini_info;
                        num_kmers_in_super_kmer = 0;
                    }

                    num_kmers_in_super_kmer += 1;
                    it.next();
                }

                save(prev_mini_info, num_kmers_in_super_kmer);
            }

            /* flush leftover */
            if (!buffer.empty()) data.minimizers.sort_and_flush(buffer);
        });
    }

    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }

    timer.stop();
    print_time(timer.elapsed(), data.num_kmers, "step 1.2: 'computing_minimizers_tuples'");

    if (build_config.weighted) {
        std::cout << "sum_of_weights " << sum_of_weights << std::endl;
        data.weights_builder.push_weight_interval(weight_value, weight_length);
        data.weights_builder.finalize(data.num_kmers);
    }
}

template <class kmer_t>
void parse_file(std::string const& filename, parse_data<kmer_t>& data,
                build_configuration const& build_config)  //
{
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    std::cout << "reading file '" << filename << "'..." << std::endl;
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
}

}  // namespace sshash