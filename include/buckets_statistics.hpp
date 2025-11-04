#pragma once

namespace sshash {

struct buckets_statistics {
    static const uint64_t MAX_BUCKET_SIZE = 4 * 1024;
    static const uint64_t MAX_STRING_SIZE = 256;

    buckets_statistics()
        : m_num_buckets(0)
        , m_num_kmers(0)
        , m_num_minimizer_positions(0)
        , m_max_num_kmers_in_super_kmer(0)
        , m_max_bucket_size(0) {}

    buckets_statistics(uint64_t num_buckets, uint64_t num_kmers, uint64_t num_minimizer_positions)
        : m_num_buckets(num_buckets)
        , m_num_kmers(num_kmers)
        , m_num_minimizer_positions(num_minimizer_positions)
        , m_max_num_kmers_in_super_kmer(0)
        , m_max_bucket_size(0)  //
    {
        m_bucket_sizes.resize(MAX_BUCKET_SIZE + 1, 0);
        m_total_kmers.resize(MAX_BUCKET_SIZE + 1, 0);
        m_super_kmer_sizes.resize(MAX_STRING_SIZE + 1, 0);
    }

    void add_bucket_size(uint64_t bucket_size) {
        if (bucket_size < MAX_BUCKET_SIZE + 1) { m_bucket_sizes[bucket_size] += 1; }
        if (bucket_size > m_max_bucket_size) { m_max_bucket_size = bucket_size; }
    }

    void add_num_kmers_in_super_kmer(uint64_t bucket_size,
                                     uint64_t num_kmers_in_super_kmer)  //
    {
        if (bucket_size < MAX_BUCKET_SIZE + 1) {
            m_total_kmers[bucket_size] += num_kmers_in_super_kmer;
        }
        if (num_kmers_in_super_kmer > m_max_num_kmers_in_super_kmer) {
            m_max_num_kmers_in_super_kmer = num_kmers_in_super_kmer;
        }
        if (num_kmers_in_super_kmer < MAX_STRING_SIZE + 1) {
            m_super_kmer_sizes[num_kmers_in_super_kmer] += 1;
        }
    }

    uint64_t num_buckets() const { return m_num_buckets; }
    uint64_t num_kmers() const { return m_num_kmers; }
    uint64_t num_minimizer_positions() const { return m_num_minimizer_positions; }
    uint64_t max_num_kmers_in_super_kmer() const { return m_max_num_kmers_in_super_kmer; }
    uint64_t max_bucket_size() const { return m_max_bucket_size; }

    void print_full() const {
        std::cout << "=== bucket statistics (full) === \n";
        for (uint64_t bucket_size = 1, prev_bucket_size = 0, prev_kmers_in_buckets = 0,
                      kmers_in_buckets = 0;
             bucket_size != MAX_BUCKET_SIZE + 1; ++bucket_size) {
            if (m_bucket_sizes[bucket_size] > 0) {
                std::cout << "buckets with " << bucket_size
                          << " minimizer positions = " << m_bucket_sizes[bucket_size] << " ("
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
                bucket_size == 512 or bucket_size == 1024 or bucket_size == MAX_BUCKET_SIZE) {
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

        std::cout << "=== super_kmer statistics === \n";
        uint64_t total_super_kmers = 0;
        uint64_t total_kmers = 0;
        for (uint64_t string_size = 1; string_size != MAX_STRING_SIZE + 1; ++string_size) {
            if (m_super_kmer_sizes[string_size] > 0) {
                std::cout << "super_kmers with " << string_size
                          << " kmer = " << m_super_kmer_sizes[string_size] << " ("
                          << (m_super_kmer_sizes[string_size] * 100.0) / m_num_minimizer_positions
                          << "%) | num_kmers = " << (string_size * m_super_kmer_sizes[string_size])
                          << " ("
                          << (string_size * m_super_kmer_sizes[string_size] * 100.0) / m_num_kmers
                          << "%)" << std::endl;
                total_super_kmers += m_super_kmer_sizes[string_size];
                total_kmers += string_size * m_super_kmer_sizes[string_size];
            }
        }
        assert(total_kmers == m_num_kmers);
        std::cout << "total_super_kmers " << total_super_kmers << "/" << m_num_minimizer_positions
                  << "(" << (total_super_kmers * 100.0) / m_num_minimizer_positions << "%)"
                  << std::endl;
        std::cout << "total_kmers " << total_kmers << "/" << m_num_kmers << " ("
                  << (total_kmers * 100.0) / m_num_kmers << "%)" << std::endl;
        std::cout << "max_num_kmers_in_super_kmer " << m_max_num_kmers_in_super_kmer << std::endl;
    }

    void print_less() const {
        std::cout << "=== bucket statistics (less) === \n";
        for (uint64_t bucket_size = 1; bucket_size != 16 + 1; ++bucket_size) {
            if (m_bucket_sizes[bucket_size] > 0) {
                std::cout << "buckets with " << bucket_size << " minimizer positions = "
                          << (m_bucket_sizes[bucket_size] * 100.0) / m_num_buckets << "%"
                          << std::endl;
            }
        }
        std::cout << "max_bucket_size = " << m_max_bucket_size << std::endl;
    }

    void operator+=(buckets_statistics const& rhs) {
        assert(m_num_buckets == rhs.num_buckets());
        assert(m_num_kmers == rhs.num_kmers());
        assert(m_num_minimizer_positions == rhs.num_minimizer_positions());

        if (rhs.max_num_kmers_in_super_kmer() > m_max_num_kmers_in_super_kmer) {
            m_max_num_kmers_in_super_kmer = rhs.max_num_kmers_in_super_kmer();
        }
        if (rhs.max_bucket_size() > m_max_bucket_size) {
            m_max_bucket_size = rhs.max_bucket_size();
        }

        assert(m_bucket_sizes.size() == rhs.m_bucket_sizes.size());
        for (uint64_t i = 0; i != m_bucket_sizes.size(); ++i) {
            m_bucket_sizes[i] += rhs.m_bucket_sizes[i];
        }
        assert(m_total_kmers.size() == rhs.m_total_kmers.size());
        for (uint64_t i = 0; i != m_total_kmers.size(); ++i) {
            m_total_kmers[i] += rhs.m_total_kmers[i];
        }
        assert(m_super_kmer_sizes.size() == rhs.m_super_kmer_sizes.size());
        for (uint64_t i = 0; i != m_super_kmer_sizes.size(); ++i) {
            m_super_kmer_sizes[i] += rhs.m_super_kmer_sizes[i];
        }
    }

private:
    uint64_t m_num_buckets;
    uint64_t m_num_kmers;
    uint64_t m_num_minimizer_positions;

    uint64_t m_max_num_kmers_in_super_kmer;
    uint64_t m_max_bucket_size;

    std::vector<uint64_t> m_bucket_sizes;
    std::vector<uint64_t> m_total_kmers;
    std::vector<uint64_t> m_super_kmer_sizes;
};

}  // namespace sshash