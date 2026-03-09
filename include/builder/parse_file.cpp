#include "dictionary_builder.hpp"
#include "util.hpp"
#include "external/gz/zip_stream.hpp"

#if defined(__AVX2__)
#include <immintrin.h>
#include <x86intrin.h>
#endif

namespace sshash {

namespace util {

#if defined(__AVX2__)
/*
    This function takes 32 bytes and packs the two bits
    in positions 1 and 2 (from right) of each byte into
    a single 64-bit word.

    This works with the map:
    A -> 00; C -> 01; G -> 11; T -> 10.
*/
inline uint64_t pack2bits_shift1(__m256i v) {
    // shift >> 1, then mask by 3 to isolate the relevant bits
    __m256i shifted = _mm256_srli_epi16(v, 1);
    __m256i values = _mm256_and_si256(shifted, _mm256_set1_epi8(3));

    // collect bit-0 plane
    __m256i bit0 = _mm256_slli_epi16(values, 7);
    uint32_t mask0 = _mm256_movemask_epi8(bit0);

    // collect bit-1 plane
    __m256i bit1 = _mm256_slli_epi16(values, 6);
    uint32_t mask1 = _mm256_movemask_epi8(bit1);

    // interleave into the 64-bit result
    uint64_t even = _pdep_u64(mask0, 0x5555555555555555ULL);  // 010101...
    uint64_t odd = _pdep_u64(mask1, 0xAAAAAAAAAAAAAAAAULL);   // 101010...
    return even | odd;
}
#endif

}  // namespace util

template <typename Kmer, typename Offsets>
void dictionary_builder<Kmer, Offsets>::encode_strings(std::istream& is,
                                                       const input_file_t fmt)  //
{
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

    {
        const uint64_t num_bits_for_strings = 8 * 8 * essentials::GB;  // reserve 8 GB of memory
        const uint64_t num_sequences = 100'000'000;
        strings_builder.reserve(num_bits_for_strings);
        strings_offsets_builder.reserve(num_sequences);
    }

    std::string sequence;
    uint64_t num_bases = 0;
    uint64_t max_len = 0;
    uint64_t seq_len = 0;
    weights_builder.init();

    /* intervals of weights */
    uint64_t weight_value = constants::invalid_uint64;
    uint64_t weight_length = 0;

    while (true)  //
    {
        if (fmt == input_file_t::cf_seg) {
            std::getline(is, sequence, '\t');  // skip until '\t' and consume it
        } else {
            assert(fmt == input_file_t::fasta);
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
                    weights_builder.eat(weight);
                    if (weight == weight_value) {
                        weight_length += 1;
                    } else {
                        if (weight_value != constants::invalid_uint64) {
                            weights_builder.push_weight_interval(weight_value, weight_length);
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
        max_len = n > max_len ? n : max_len;
        assert(strings_builder.num_bits() % Kmer::bits_per_char == 0);
        strings_offsets_builder.push_back(strings_builder.num_bits() / Kmer::bits_per_char);
        num_kmers += n - k + 1;
        num_bases += n;

        if (build_config.weighted and seq_len != n) {
            std::cout << "ERROR: expected a sequence of length " << seq_len
                      << " but got one of length " << n << std::endl;
            throw std::runtime_error("file is malformed");
        }

        if (build_config.verbose and strings_offsets_builder.size() % 1'000'000 == 0) {
            std::cout << "read " << strings_offsets_builder.size() << " sequences, " << num_bases
                      << " bases, " << num_kmers << " kmers" << std::endl;
        }

        uint64_t i = 0;
        if constexpr (Kmer::bits_per_char == 2) {
#if !defined(SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING) and defined(__AVX2__)

            /* process 32 bytes at a time */
            for (; i + 32 <= n; i += 32) {
                __m256i v = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(&sequence[i]));
                uint64_t word = util::pack2bits_shift1(v);
                strings_builder.append_bits(word, 64);
            }
#endif
        }
        for (; i < n; ++i) {
            strings_builder.append_bits(Kmer::char_to_uint(sequence[i]), Kmer::bits_per_char);
        }
    }

    strings_offsets_builder.push_back(strings_builder.num_bits() / Kmer::bits_per_char);
    assert(strings_offsets_builder.front() == 0);
    assert(strings_offsets_builder.size() >= 2);

    /* Push a final sentinel (dummy) value to avoid bounds' checking in
       kmer_iterator::fill_buff(). */
    static_assert(Kmer::uint_kmer_bits % 64 == 0);
    for (int dummy_bits = Kmer::uint_kmer_bits; dummy_bits; dummy_bits -= 64) {
        strings_builder.append_bits(0, 64);
    }

    const uint64_t num_sequences = strings_offsets_builder.size() - 1;

    if (build_config.verbose) {
        std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases, "
                  << num_kmers << " kmers" << std::endl;
        std::cout << "num_kmers " << num_kmers << std::endl;
        std::cout << "cost: 2.0 + "
                  << static_cast<double>(Kmer::bits_per_char * num_sequences * (k - 1)) / num_kmers
                  << " [bits/kmer]" << std::endl;
    }

    /*
        The parameter m (minimizer length) should be at least
            ceil(log_s(N))+1
        where N is the number of nucleotides in the input and s is the alphabet size.
        We warn the user if the used m is less than this lower bound.
    */
    const uint64_t s = uint64_t(1) << Kmer::bits_per_char;
    const uint64_t lower_bound_on_m = std::ceil(std::log(num_bases) / std::log(s)) + 1;
    if (build_config.verbose and m < lower_bound_on_m) {
        std::cout << "\n--> WARNING: using minimizer length " << m
                  << " but it should be at least ceil(log_" << s << "(" << num_bases
                  << "))+1 = " << lower_bound_on_m << '\n'
                  << std::endl;
    }

    if (build_config.weighted) {
        weights_builder.push_weight_interval(weight_value, weight_length);
        weights_builder.finalize(num_kmers);
    }

    num_bits nb;
    nb.per_absolute_offset = std::ceil(std::log2(strings_offsets_builder.back()));
    nb.per_relative_offset = std::ceil(std::log2(max_len - m + 1));
    nb.per_string_id = std::ceil(std::log2(num_sequences));

    if (build_config.verbose) {
        std::cout << "max string length = " << max_len << std::endl;
        std::cout << "num bits per_absolute_offset = " << nb.per_absolute_offset << std::endl;
        std::cout << "num bits per_relative_offset = " << nb.per_relative_offset << std::endl;
        std::cout << "num bits per_string_id = " << nb.per_string_id << std::endl;
    }

    if (nb.per_string_id + nb.per_relative_offset > 64) {
        throw std::runtime_error("minimizer offset does not fit within 64 bits");
    }

    strings_offsets_builder.set_num_bits(nb);
}

template <typename Kmer, typename Offsets>
void dictionary_builder<Kmer, Offsets>::encode_strings(std::string const& filename)  //
{
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    if (build_config.verbose) std::cout << "reading file '" << filename << "'..." << std::endl;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        if (util::ends_with(filename, ".cf_seg.gz")) {
            encode_strings(zis, input_file_t::cf_seg);
        } else {
            encode_strings(zis, input_file_t::fasta);
        }
    } else {
        if (util::ends_with(filename, ".cf_seg")) {
            encode_strings(is, input_file_t::cf_seg);
        } else {
            encode_strings(is, input_file_t::fasta);
        }
    }
    is.close();
}

}  // namespace sshash