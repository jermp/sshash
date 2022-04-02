#pragma once

#include <vector>
#include <unordered_map>  // count the distinct abundances with freq information

#include "ef_sequence.hpp"

namespace sshash {

struct abundances {
    struct builder {
        builder() {}

        void init() { m_abundance_interval_lengths.push_back(0); }

        void eat(uint64_t abundance) {
            assert(abundance > 0);
            auto it = m_abundances_map.find(abundance);
            if (it != m_abundances_map.cend()) {  // found
                (*it).second += 1;
            } else {
                m_abundances_map[abundance] = 1;
            }
        }

        void push_abundance_interval(uint64_t value, uint64_t length) {
            m_abundance_interval_values.push_back(value);
            m_abundance_interval_lengths.push_back(m_abundance_interval_lengths.back() + length);
        }

        uint64_t num_abundance_intervals() const { return m_abundance_interval_values.size(); }

        void finalize(uint64_t num_kmers) {
            assert(std::is_sorted(m_abundance_interval_lengths.begin(),
                                  m_abundance_interval_lengths.end()));

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

            m_abundance_dictionary_builder.resize(num_distinct_abundances,
                                                  std::ceil(std::log2(largest_ab + 1)));
            for (uint64_t id = 0; id != num_distinct_abundances; ++id) {
                uint64_t ab_value = m_abundances[id].first;
                m_abundance_dictionary_builder.set(id, ab_value);
                m_abundances_map[ab_value] = id;
            }
        }

        void build(abundances& index) {
            uint64_t num_distinct_abundances = m_abundance_dictionary_builder.size();
            pthash::compact_vector::builder abundance_interval_values;
            abundance_interval_values.resize(
                m_abundance_interval_values.size(),
                num_distinct_abundances == 1 ? 1 : std::ceil(std::log2(num_distinct_abundances)));
            uint64_t prev_abundance = constants::invalid;
            for (uint64_t i = 0; i != m_abundance_interval_values.size(); ++i) {
                uint64_t abundance = m_abundance_interval_values[i];
                if (i != 0) {
                    if (abundance == prev_abundance) {
                        std::cerr << "Error at " << i << "/" << m_abundance_interval_values.size()
                                  << ": cannot have two consecutive intervals with the same "
                                     "abundance value"
                                  << std::endl;
                        throw std::runtime_error("abundance intervals are malformed");
                    }
                }
                prev_abundance = abundance;
                uint64_t id = m_abundances_map[abundance];
                assert(id < num_distinct_abundances);
                abundance_interval_values.set(i, id);
            }
            abundance_interval_values.build(index.m_abundance_interval_values);
            index.m_abundance_interval_lengths.encode(m_abundance_interval_lengths.begin(),
                                                      m_abundance_interval_lengths.size());

            m_abundance_dictionary_builder.build(index.m_abundance_dictionary);
        }

        /* return the average empirical entropy per abundance */
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
        /* (abundance,frequency) pairs during construction, then (abundance,id) after sorting */
        std::unordered_map<uint64_t, uint64_t> m_abundances_map;
        std::vector<std::pair<uint64_t, uint64_t>> m_abundances;  // (abundance,frequency)

        std::vector<uint64_t> m_abundance_interval_values;
        std::vector<uint64_t> m_abundance_interval_lengths;

        pthash::compact_vector::builder m_abundance_dictionary_builder;
    };

    bool empty() const { return m_abundance_dictionary.size() == 0; }

    uint64_t abundance(uint64_t kmer_id) const {
        uint64_t i = m_abundance_interval_lengths.prev_leq(kmer_id);
        uint64_t id = m_abundance_interval_values.access(i);
        uint64_t abundance = m_abundance_dictionary.access(id);
        return abundance;
    }

    uint64_t num_bits() const {
        return m_abundance_interval_values.bytes() * 8 + m_abundance_interval_lengths.num_bits() +
               m_abundance_dictionary.bytes() * 8;
    }

    void print_space_breakdown(uint64_t num_kmers) const {
        std::cout << "    abundance_interval_values: "
                  << static_cast<double>(m_abundance_interval_values.bytes() * 8) / num_kmers
                  << " [bits/kmer]\n";
        std::cout << "    abundance_interval_lengths: "
                  << static_cast<double>(m_abundance_interval_lengths.num_bits()) / num_kmers
                  << " [bits/kmer]\n";
        std::cout << "    abundance_dictionary: "
                  << static_cast<double>(m_abundance_dictionary.bytes() * 8) / num_kmers
                  << " [bits/kmer]\n";
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_abundance_interval_values);
        visitor.visit(m_abundance_interval_lengths);
        visitor.visit(m_abundance_dictionary);
    }

private:
    pthash::compact_vector m_abundance_interval_values;
    ef_sequence<true> m_abundance_interval_lengths;
    pthash::compact_vector m_abundance_dictionary;
};

}  // namespace sshash