#pragma once

namespace sshash {

struct buckets_statistics {
    static const uint64_t max_bucket_size = 4 * 1024;
    static const uint64_t max_string_size = 256;

    buckets_statistics(uint64_t num_buckets, uint64_t num_kmers, uint64_t num_super_kmers)
        : m_num_buckets(num_buckets)
        , m_num_kmers(num_kmers)
        , m_num_super_kmers(num_super_kmers)
        , m_max_num_kmers_in_super_kmer(0)
        , m_max_num_super_kmers_in_bucket(0) {
        m_bucket_sizes.resize(max_bucket_size + 1, 0);
        m_total_kmers.resize(max_bucket_size + 1, 0);
        m_string_sizes.resize(max_string_size + 1, 0);
    }

    void add_num_super_kmers_in_bucket(uint64_t num_super_kmers_in_bucket) {
        if (num_super_kmers_in_bucket < max_bucket_size + 1) {
            m_bucket_sizes[num_super_kmers_in_bucket] += 1;
        }
    }

    void add_num_kmers_in_super_kmer(uint64_t num_super_kmers_in_bucket,
                                     uint64_t num_kmers_in_super_kmer) {
        if (num_super_kmers_in_bucket < max_bucket_size + 1) {
            m_total_kmers[num_super_kmers_in_bucket] += num_kmers_in_super_kmer;
        }
        if (num_kmers_in_super_kmer > m_max_num_kmers_in_super_kmer) {
            m_max_num_kmers_in_super_kmer = num_kmers_in_super_kmer;
        }
        if (num_super_kmers_in_bucket > m_max_num_super_kmers_in_bucket) {
            m_max_num_super_kmers_in_bucket = num_super_kmers_in_bucket;
        }
        if (num_kmers_in_super_kmer < max_string_size + 1)
            m_string_sizes[num_kmers_in_super_kmer] += 1;
    }

    uint64_t num_kmers() const { return m_num_kmers; }
    uint64_t num_buckets() const { return m_num_buckets; }
    uint64_t max_num_super_kmers_in_bucket() const { return m_max_num_super_kmers_in_bucket; }

    void print_full() const {
        std::cout << " === bucket statistics (full) === \n";
        for (uint64_t bucket_size = 1, prev_bucket_size = 0, prev_kmers_in_buckets = 0,
                      kmers_in_buckets = 0;
             bucket_size != max_bucket_size + 1; ++bucket_size) {
            if (m_bucket_sizes[bucket_size] > 0) {
                std::cout << "buckets with " << bucket_size
                          << " super_kmers = " << m_bucket_sizes[bucket_size] << " ("
                          << (m_bucket_sizes[bucket_size] * 100.0) / m_num_buckets
                          << "%) | num_kmers = " << m_total_kmers[bucket_size] << " ("
                          << (m_total_kmers[bucket_size] * 100.0) / m_num_kmers
                          << "%)"
                          // << "|avg_num_kmers_per_bucket="
                          // << static_cast<double>(m_total_kmers[bucket_size]) /
                          //        m_bucket_sizes[bucket_size]
                          // << "|avg_num_kmers_per_string="
                          // << static_cast<double>(m_total_kmers[bucket_size]) /
                          //        (m_bucket_sizes[bucket_size] * bucket_size)
                          << std::endl;
                kmers_in_buckets += m_total_kmers[bucket_size];
            }
            if (bucket_size == 4 or bucket_size == 8 or bucket_size == 16 or bucket_size == 32 or
                bucket_size == 64 or bucket_size == 128 or bucket_size == 256 or
                bucket_size == 512 or bucket_size == 1024 or bucket_size == max_bucket_size) {
                assert(kmers_in_buckets >= prev_kmers_in_buckets);

                std::cout << " *** num_kmers in buckets of size > " << prev_bucket_size
                          << " and <= " << bucket_size << ": "
                          << kmers_in_buckets - prev_kmers_in_buckets << " ("
                          << (100.0 * (kmers_in_buckets - prev_kmers_in_buckets)) / m_num_kmers
                          << "%)" << std::endl;
                std::cout << " *** num_kmers in buckets of size <= " << bucket_size << ": "
                          << kmers_in_buckets << " (" << (100.0 * kmers_in_buckets) / m_num_kmers
                          << "%)" << std::endl;

                prev_bucket_size = bucket_size;
                prev_kmers_in_buckets = kmers_in_buckets;
            }
        }

        std::cout << " === super_kmer statistics === \n";
        uint64_t total_super_kmers = 0;
        uint64_t total_kmers = 0;
        for (uint64_t string_size = 1; string_size != max_string_size + 1; ++string_size) {
            if (m_string_sizes[string_size] > 0) {
                std::cout << "super_kmers with " << string_size
                          << " kmer = " << m_string_sizes[string_size] << " ("
                          << (m_string_sizes[string_size] * 100.0) / m_num_super_kmers
                          << "%) | num_kmers = " << (string_size * m_string_sizes[string_size])
                          << " ("
                          << (string_size * m_string_sizes[string_size] * 100.0) / m_num_kmers
                          << "%)" << std::endl;
                total_super_kmers += m_string_sizes[string_size];
                total_kmers += string_size * m_string_sizes[string_size];
            }
        }
        std::cout << "total_super_kmers " << total_super_kmers << "/" << m_num_super_kmers << "("
                  << (total_super_kmers * 100.0) / m_num_super_kmers << "%)" << std::endl;
        std::cout << "total_kmers " << total_kmers << "/" << m_num_kmers << " ("
                  << (total_kmers * 100.0) / m_num_kmers << "%)" << std::endl;
        std::cout << "max_num_kmers_in_super_kmer " << m_max_num_kmers_in_super_kmer << std::endl;
    }

    void print_less() const {
        std::cout << " === bucket statistics (less) === \n";
        for (uint64_t bucket_size = 1; bucket_size != 16 + 1; ++bucket_size) {
            if (m_bucket_sizes[bucket_size] > 0) {
                std::cout << "buckets with " << bucket_size << " super_kmers = "
                          << (m_bucket_sizes[bucket_size] * 100.0) / m_num_buckets << "%"
                          << std::endl;
            }
        }
        std::cout << "max_num_super_kmers_in_bucket " << m_max_num_super_kmers_in_bucket
                  << std::endl;
    }

private:
    uint64_t m_num_buckets;
    uint64_t m_num_kmers;
    uint64_t m_num_super_kmers;
    uint64_t m_max_num_kmers_in_super_kmer;
    uint64_t m_max_num_super_kmers_in_bucket;
    std::vector<uint64_t> m_bucket_sizes;
    std::vector<uint64_t> m_total_kmers;
    std::vector<uint64_t> m_string_sizes;
};

}  // namespace sshash