#include <numeric>  // for std::partial_sum

namespace sshash {

buckets_statistics build_index(parse_data& data, minimizers const& m_minimizers,
                               buckets& m_buckets) {
    uint64_t num_buckets = m_minimizers.size();
    uint64_t num_kmers = data.num_kmers;
    uint64_t num_super_kmers = data.strings.num_super_kmers();
    std::vector<uint64_t> num_super_kmers_before_bucket(num_buckets + 1, 0);
    pthash::compact_vector::builder offsets;
    offsets.resize(num_super_kmers, std::ceil(std::log2(data.strings.num_bits() / 2)));

    std::cout << "bits_per_offset = ceil(log2(" << data.strings.num_bits() / 2
              << ")) = " << std::ceil(std::log2(data.strings.num_bits() / 2)) << std::endl;

    for (auto it = data.minimizers.begin(); it.has_next(); it.next()) {
        assert(it.list().size() > 0);
        if (it.list().size() != 1) {
            uint64_t bucket_id = m_minimizers.lookup(it.minimizer());
            num_super_kmers_before_bucket[bucket_id + 1] = it.list().size() - 1;
        }
        // else: it.list().size() = 1 and num_super_kmers_before_bucket[bucket_id + 1] is already 0
    }
    std::partial_sum(num_super_kmers_before_bucket.begin(), num_super_kmers_before_bucket.end(),
                     num_super_kmers_before_bucket.begin());

    buckets_statistics buckets_stats(num_buckets, num_kmers, num_super_kmers);

    uint64_t num_singletons = 0;
    for (auto it = data.minimizers.begin(); it.has_next(); it.next()) {
        uint64_t bucket_id = m_minimizers.lookup(it.minimizer());
        uint64_t base = num_super_kmers_before_bucket[bucket_id] + bucket_id;
        uint64_t num_super_kmers_in_bucket =
            (num_super_kmers_before_bucket[bucket_id + 1] + bucket_id + 1) - base;
        assert(num_super_kmers_in_bucket > 0);
        if (num_super_kmers_in_bucket == 1) ++num_singletons;
        buckets_stats.add_num_super_kmers_in_bucket(num_super_kmers_in_bucket);
        uint64_t offset_pos = 0;
        auto list = it.list();
        for (auto [offset, num_kmers_in_super_kmer] : list) {
            offsets.set(base + offset_pos++, offset);
            buckets_stats.add_num_kmers_in_super_kmer(num_super_kmers_in_bucket,
                                                      num_kmers_in_super_kmer);
        }
        assert(offset_pos == num_super_kmers_in_bucket);
    }

    std::cout << "num_singletons " << num_singletons << "/" << num_buckets << " ("
              << (num_singletons * 100.0) / num_buckets << "%)" << std::endl;

    m_buckets.pieces.encode(data.strings.pieces.begin(), data.strings.pieces.size());
    m_buckets.num_super_kmers_before_bucket.encode(num_super_kmers_before_bucket.begin(),
                                                   num_super_kmers_before_bucket.size());
    offsets.build(m_buckets.offsets);
    m_buckets.strings.swap(data.strings.strings);

    return buckets_stats;
}

}  // namespace sshash