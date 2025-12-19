#include "dictionary_builder.hpp"
#include "util.hpp"
#include "include/kmer_iterator.hpp"
#include "include/minimizer_iterator.hpp"

namespace sshash {

template <typename Kmer, typename Offsets>
void dictionary_builder<Kmer, Offsets>::compute_minimizer_tuples()  //
{
    const uint64_t num_threads = build_config.num_threads;
    const uint64_t num_sequences = strings_offsets_builder.size() - 1;
    const uint64_t num_sequences_per_thread = (num_sequences + num_threads - 1) / num_threads;
    const uint64_t k = build_config.k;
    const uint64_t m = build_config.m;

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
                assert(num_kmers_in_super_kmer <= k - m + 1 /* max num kmers in super-kmer */);
                if (!buffer.empty() and                                   //
                    buffer.back().minimizer == mini_info.minimizer and    //
                    buffer.back().pos_in_seq == mini_info.pos_in_seq and  //
                    buffer.back().pos_in_kmer == mini_info.pos_in_kmer)   //
                {
                    buffer.back().num_kmers_in_super_kmer += num_kmers_in_super_kmer;
                    return;
                }
                if (buffer.size() == buffer_size) {
                    minimizers.sort_and_flush(buffer);
                    buffer.clear();
                }
                buffer.emplace_back(mini_info, num_kmers_in_super_kmer);
            };

            const uint64_t index_begin = t * num_sequences_per_thread;
            const uint64_t index_end =
                std::min<uint64_t>(index_begin + num_sequences_per_thread, num_sequences);

            kmer_iterator<Kmer, bits::bit_vector::builder> kmer_it(strings_builder, k);
            hasher_type hasher(build_config.seed);
            minimizer_iterator<Kmer> minimizer_it(k, m, hasher);
            minimizer_iterator_rc<Kmer> minimizer_it_rc(k, m, hasher);

            for (uint64_t i = index_begin; i < index_end; ++i)  //
            {
                const uint64_t begin = strings_offsets_builder[i];
                const uint64_t end = strings_offsets_builder[i + 1];
                const uint64_t sequence_len = end - begin;
                assert(sequence_len >= k);

                minimizer_info prev_mini_info;
                assert(prev_mini_info.minimizer == constants::invalid_uint64);
                uint64_t num_kmers_in_super_kmer = 0;

                kmer_it.at(Kmer::bits_per_char * begin);
                minimizer_it.set_position(begin);
                minimizer_it_rc.set_position(begin);

                for (uint64_t j = 0; j != sequence_len - k + 1; ++j) {
                    auto uint_kmer = kmer_it.get();
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

                    mini_info.pos_in_seq =
                        strings_offsets_builder.encode(mini_info.pos_in_seq, begin, i);

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
                    kmer_it.next();
                }

                save(prev_mini_info, num_kmers_in_super_kmer);
            }

            /* flush leftover */
            if (!buffer.empty()) minimizers.sort_and_flush(buffer);
        });
    }

    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }
}

}  // namespace sshash
