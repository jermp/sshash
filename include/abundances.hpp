#pragma once

#include <vector>
#include <unordered_map>  // count the distinct abundances
#include "ef_sequence.hpp"

namespace sshash {

struct abundances {
    struct builder {
        builder() : m_most_frequent_abundance(0) {}

        void init(uint64_t most_frequent_abundance) {
            m_most_frequent_abundance = most_frequent_abundance;
            m_kmer_id_interval_lengths.push_back(0);
            m_abundance_interval_lengths.push_back(0);
        }

        void eat(uint64_t abundance) {
            assert(abundance > 0);
            auto it = m_abundances_map.find(abundance);
            if (it != m_abundances_map.cend()) {  // found
                (*it).second += 1;
            } else {
                m_abundances_map[abundance] = 1;
            }
        }

        void push_kmer_id_interval(uint64_t value, uint64_t length) {
            m_kmer_id_interval_values.push_back(value);
            m_kmer_id_interval_lengths.push_back(m_kmer_id_interval_lengths.back() + length);
        }

        void push_abundance_interval(uint64_t value, uint64_t length) {
            m_abundance_interval_values.push_back(value);
            m_abundance_interval_lengths.push_back(m_abundance_interval_lengths.back() + length);
        }

        uint64_t num_kmer_id_intervals() const { return m_kmer_id_interval_values.size(); }
        uint64_t num_abundance_intervals() const { return m_abundance_interval_values.size(); }

        void finalize(uint64_t num_kmers) {
            assert(
                std::is_sorted(m_kmer_id_interval_values.begin(), m_kmer_id_interval_values.end()));
            assert(std::is_sorted(m_kmer_id_interval_lengths.begin(),
                                  m_kmer_id_interval_lengths.end()));
            assert(std::is_sorted(m_abundance_interval_lengths.begin(),
                                  m_abundance_interval_lengths.end()));

            std::cout << "num_kmer_id_intervals " << num_kmer_id_intervals() << std::endl;
            std::cout << "num_abundance_intervals " << num_abundance_intervals() << std::endl;

            uint64_t num_distinct_abundances = m_abundances_map.size();

            std::cout << "found " << num_distinct_abundances << " distint abundances (ceil(log2("
                      << num_distinct_abundances
                      << ")) = " << std::ceil(std::log2(num_distinct_abundances)) << ")"
                      << std::endl;

            m_abundances.reserve(num_distinct_abundances);
            uint64_t n = 0;
            uint64_t largest_ab = 0;
            for (auto p : m_abundances_map) {
                if (p.first > largest_ab) largest_ab = p.first;
                n += p.second;
                m_abundances.push_back(p);
            }
            assert(largest_ab > 0);

            std::cout << "largest_ab+1 = " << largest_ab + 1 << " (ceil(log2(" << largest_ab + 1
                      << ")) = " << std::ceil(std::log2(largest_ab + 1)) << ")" << std::endl;

            if (n != num_kmers) {
                std::cout << "ERROR: expected " << num_kmers << " kmers but got " << n << std::endl;
                throw std::runtime_error("file is malformed");
            }

            std::sort(m_abundances.begin(), m_abundances.end(), [](auto const& x, auto const& y) {
                if (x.second != y.second) return x.second > y.second;
                return x.first < y.first;
            });

            uint64_t rest = num_kmers - m_abundances.front().second;
            std::cout << "kmers that do not have the most frequent ab: " << rest << " ("
                      << (rest * 100.0) / num_kmers << "%)" << std::endl;
            std::cout << "cumulative_kmer_id_interval_lengths " << m_kmer_id_interval_lengths.back()
                      << '/' << rest << std::endl;
            std::cout << "cumulative_abundance_interval_lengths "
                      << m_abundance_interval_lengths.back() << '/' << rest << std::endl;

            m_abundance_dictionary_builder.resize(num_distinct_abundances,
                                                  std::ceil(std::log2(largest_ab + 1)));
            for (uint64_t id = 0; id != num_distinct_abundances; ++id) {
                uint64_t ab_value = m_abundances[id].first;
                m_abundance_dictionary_builder.set(id, ab_value);
                m_abundances_map[ab_value] = id;
            }
        }

        void build(abundances& index) {
            std::swap(index.m_most_frequent_abundance, m_most_frequent_abundance);

            index.m_kmer_id_interval_values.encode(m_kmer_id_interval_values.begin(),
                                                   m_kmer_id_interval_values.size());
            index.m_kmer_id_interval_lengths.encode(m_kmer_id_interval_lengths.begin(),
                                                    m_kmer_id_interval_lengths.size());

            uint64_t num_distinct_abundances = m_abundance_dictionary_builder.size();
            pthash::compact_vector::builder abundance_interval_values;
            abundance_interval_values.resize(
                m_abundance_interval_values.size(),
                num_distinct_abundances == 1 ? 1 : std::ceil(std::log2(num_distinct_abundances)));
            for (uint64_t i = 0; i != m_abundance_interval_values.size(); ++i) {
                uint64_t abundance = m_abundance_interval_values[i];
                uint64_t id = m_abundances_map[abundance];
                assert(id < num_distinct_abundances);
                abundance_interval_values.set(i, id);
            }
            abundance_interval_values.build(index.m_abundance_interval_values);
            index.m_abundance_interval_lengths.encode(m_abundance_interval_lengths.begin(),
                                                      m_abundance_interval_lengths.size());

            m_abundance_dictionary_builder.build(index.m_abundance_dictionary);
        }

        // return the average empirical entropy per abundance
        double print_info(uint64_t num_kmers) {
            assert(!m_abundances.empty());
            double expected_ab_value = 0.0;
            double entropy_ab = 0.0;
            uint64_t print = 0;
            for (auto p : m_abundances) {
                double prob = static_cast<double>(p.second) / num_kmers;
                expected_ab_value += p.first * prob;
                entropy_ab += prob * std::log2(1.0 / prob);
                print += 1;
                if (print <= 10) {
                    std::cout << "ab:" << p.first << " freq:" << p.second << " ("
                              << (p.second * 100.0) / num_kmers << "%)" << std::endl;
                }
            }
            std::cout << "expected_ab_value " << expected_ab_value << std::endl;
            std::cout << "entropy_ab " << entropy_ab << " [bits/kmer]" << std::endl;
            return entropy_ab;
        }

    private:
        uint64_t m_most_frequent_abundance;

        /* (abundance,frequency) pairs during construction, then (abundance,id) after sorting */
        std::unordered_map<uint64_t, uint64_t> m_abundances_map;
        std::vector<std::pair<uint64_t, uint64_t>> m_abundances;  // (abundance,frequency)

        std::vector<uint64_t> m_kmer_id_interval_values;
        std::vector<uint64_t> m_kmer_id_interval_lengths;

        std::vector<uint64_t> m_abundance_interval_values;
        std::vector<uint64_t> m_abundance_interval_lengths;

        pthash::compact_vector::builder m_abundance_dictionary_builder;
    };

    bool empty() const { return m_abundance_dictionary.size() == 0; }

    uint64_t abundance(uint64_t kmer_id) const {
        bool is_present = false;
        uint64_t rank = 0;

        auto [pos, val] = m_kmer_id_interval_values.next_geq(kmer_id);
        if (val == kmer_id) {
            is_present = true;
            rank = m_kmer_id_interval_lengths.access(pos);
        } else {
            if (pos > 0) {
                --pos;
                val = m_kmer_id_interval_values.access(pos);
                rank = m_kmer_id_interval_lengths.access(pos);
                uint64_t length = m_kmer_id_interval_lengths.access(pos + 1) - rank;
                if (kmer_id < val + length) {
                    is_present = true;
                    assert(kmer_id >= val);
                    rank += kmer_id - val;
                }
            }
        }

        if (!is_present) return m_most_frequent_abundance;

        uint64_t i = m_abundance_interval_lengths.prev_leq(rank);
        uint64_t id = m_abundance_interval_values.access(i);
        uint64_t abundance = m_abundance_dictionary.access(id);

        return abundance;
    }

    uint64_t num_bits() const {
        return sizeof(m_most_frequent_abundance) * 8 + m_kmer_id_interval_values.num_bits() +
               m_kmer_id_interval_lengths.num_bits() + m_abundance_interval_values.bytes() * 8 +
               m_abundance_interval_lengths.num_bits() + m_abundance_dictionary.bytes() * 8;
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_most_frequent_abundance);
        visitor.visit(m_kmer_id_interval_values);
        visitor.visit(m_kmer_id_interval_lengths);
        visitor.visit(m_abundance_interval_values);
        visitor.visit(m_abundance_interval_lengths);
        visitor.visit(m_abundance_dictionary);
    }

private:
    uint64_t m_most_frequent_abundance;

    /*****
       We model abundances as two lists of intervals.
       Each interval is a pair (value,length).
       - First list: kmer_ids list. In this case, a pair (value,length)
       represents all kmer_ids = value, value+1, value+2, ..., value+length-1.
       - Second list: abundance list. In this case, a pair (value,length)
       represents that the abundance [value] repeats for [length] times.
       */

    ef_sequence<true> m_kmer_id_interval_values;
    ef_sequence<false> m_kmer_id_interval_lengths;

    pthash::compact_vector m_abundance_interval_values;
    ef_sequence<true> m_abundance_interval_lengths;
    /***/

    /* Abundance dictionary, listing all distinct abundances sorted by decreasing frequency. */
    pthash::compact_vector m_abundance_dictionary;
};

}  // namespace sshash