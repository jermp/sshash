#pragma once

#include "external/gz/zip_stream.hpp"
#include "include/minimizer_enumerator.hpp"

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
    const uint64_t k = build_config.k;
    const uint64_t m = build_config.m;
    const uint64_t max_num_kmers_in_super_kmer = k - m + 1;
    const uint64_t block_size = 2 * k - m;  // max_num_kmers_in_super_kmer + k - 1
    hasher_type hasher(build_config.seed);

    if (max_num_kmers_in_super_kmer >= (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8))) {
        throw std::runtime_error(
            "max_num_kmers_in_super_kmer " + std::to_string(max_num_kmers_in_super_kmer) +
            " does not fit into " + std::to_string(sizeof(num_kmers_in_super_kmer_uint_type) * 8) +
            " bits");
    }

    /* fit into the wanted number of bits */
    assert(max_num_kmers_in_super_kmer < (1ULL << (sizeof(num_kmers_in_super_kmer_uint_type) * 8)));

    typename compact_string_pool<kmer_t>::builder builder(k);

    std::string sequence;
    uint64_t prev_minimizer = constants::invalid_uint64;

    uint64_t begin = 0;  // begin of parsed super_kmer in sequence
    uint64_t end = 0;    // end of parsed super_kmer in sequence
    uint64_t num_sequences = 0;
    uint64_t num_bases = 0;
    bool glue = false;

    auto append_super_kmer = [&]() {
        if (sequence.empty() or prev_minimizer == constants::invalid_uint64 or begin == end) {
            return;
        }

        assert(end > begin);
        char const* super_kmer = sequence.data() + begin;
        uint64_t size = (end - begin) + k - 1;
        assert(util::is_valid<kmer_t>(super_kmer, size));

        /* if num_kmers_in_super_kmer > k - m + 1, then split the super_kmer into blocks */
        uint64_t num_kmers_in_super_kmer = end - begin;
        uint64_t num_blocks = num_kmers_in_super_kmer / max_num_kmers_in_super_kmer +
                              (num_kmers_in_super_kmer % max_num_kmers_in_super_kmer != 0);
        assert(num_blocks > 0);
        for (uint64_t i = 0; i != num_blocks; ++i) {
            uint64_t n = block_size;
            if (i == num_blocks - 1) n = size;
            uint64_t num_kmers_in_block = n - k + 1;
            assert(num_kmers_in_block <= max_num_kmers_in_super_kmer);
            data.minimizers.emplace_back(prev_minimizer, builder.offset, num_kmers_in_block);
            builder.append(super_kmer + i * max_num_kmers_in_super_kmer, n, glue);
            if (glue) {
                assert(data.minimizers.back().offset > k - 1);
                data.minimizers.back().offset -= k - 1;
            }
            size -= max_num_kmers_in_super_kmer;
            glue = true;
        }
    };

    minimizer_enumerator<kmer_t> minimizer_enum(k, m, hasher);
    minimizer_enumerator<kmer_t> minimizer_enum_rc(k, m, hasher);
    uint64_t seq_len = 0;
    uint64_t sum_of_weights = 0;
    data.weights_builder.init();

    /* intervals of weights */
    uint64_t weight_value = constants::invalid_uint64;
    uint64_t weight_length = 0;

    while (!is.eof())  //
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

        if (sequence.length() < k) continue;

        if (++num_sequences % 100000 == 0) {
            std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
                      << data.num_kmers << " kmers" << std::endl;
        }

        begin = 0;
        end = 0;
        glue = false;  // start a new piece
        prev_minimizer = constants::invalid_uint64;
        num_bases += sequence.length();

        if (build_config.weighted and seq_len != sequence.length()) {
            std::cout << "ERROR: expected a sequence of length " << seq_len
                      << " but got one of length " << sequence.length() << std::endl;
            throw std::runtime_error("file is malformed");
        }

        bool start = true;
        kmer_t uint_kmer = 0;
        while (end != sequence.length() - k + 1) {
            char const* kmer = sequence.data() + end;
            assert(util::is_valid<kmer_t>(kmer, k));

            if (!start) {
                uint_kmer.drop_char();
                uint_kmer.set(k - 1, kmer_t::char_to_uint(kmer[k - 1]));
                assert(uint_kmer == util::string_to_uint_kmer<kmer_t>(kmer, k));
            } else {
                uint_kmer = util::string_to_uint_kmer<kmer_t>(kmer, k);
            }

            uint64_t minimizer = minimizer_enum.template next<false>(uint_kmer, start);
            if (build_config.canonical) {
                kmer_t uint_kmer_rc = uint_kmer;
                uint_kmer_rc.reverse_complement_inplace(k);
                uint64_t minimizer_rc = minimizer_enum_rc.template next<true>(uint_kmer_rc, start);
                minimizer = std::min(minimizer, minimizer_rc);
            }

            if (prev_minimizer == constants::invalid_uint64) prev_minimizer = minimizer;
            if (minimizer != prev_minimizer) {
                append_super_kmer();
                begin = end;
                prev_minimizer = minimizer;
                glue = true;
            }

            ++data.num_kmers;
            ++end;
            start = false;
        }

        append_super_kmer();
    }

    data.minimizers.finalize();
    builder.finalize();
    builder.build(data.strings);

    std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
              << data.num_kmers << " kmers" << std::endl;
    std::cout << "num_kmers " << data.num_kmers << std::endl;
    std::cout << "num_super_kmers " << data.strings.num_super_kmers() << std::endl;
    std::cout << "num_pieces " << data.strings.pieces.size() << " (+"
              << (2.0 * data.strings.pieces.size() * (k - 1)) / data.num_kmers << " [bits/kmer])"
              << std::endl;
    assert(data.strings.pieces.size() == num_sequences + 1);

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