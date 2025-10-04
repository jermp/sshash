#pragma once

#include "util.hpp"

namespace sshash {

struct minimizers {
    template <typename ForwardIterator>
    void build(ForwardIterator begin, uint64_t size, build_configuration const& build_config)  //
    {
        pthash::build_configuration mphf_build_config;
        mphf_build_config.lambda = 5.0;
        mphf_build_config.alpha = 0.94;
        mphf_build_config.seed = util::get_seed_for_hash_function(build_config);
        mphf_build_config.verbose = false;
        mphf_build_config.num_threads = build_config.num_threads;
        mphf_build_config.avg_partition_size = constants::avg_partition_size;
        mphf_build_config.ram = (build_config.ram_limit_in_GiB * essentials::GiB) / 2;
        mphf_build_config.tmp_dir = build_config.tmp_dirname;

        if (build_config.verbose) {
            const uint64_t avg_partition_size =
                pthash::compute_avg_partition_size(size, mphf_build_config);
            const uint64_t num_partitions =
                pthash::compute_num_partitions(size, avg_partition_size);
            assert(num_partitions > 0);
            std::cout << "building minimizers MPHF with " << mphf_build_config.num_threads
                      << " threads and " << num_partitions
                      << " partitions (avg. partition size = " << avg_partition_size << ")..."
                      << std::endl;
        }

        m_mphf.build_in_external_memory(begin, size, mphf_build_config);
    }

    uint64_t lookup(uint64_t uint64_minimizer) const {
        uint64_t bucket_id = m_mphf(uint64_minimizer);
        return bucket_id;
    }

    uint64_t size() const { return m_mphf.num_keys(); }
    uint64_t num_bits() const { return m_mphf.num_bits(); }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visitor.visit(m_mphf);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_mphf);
    }

private:
    minimizers_pthash_type m_mphf;
};

}  // namespace sshash
