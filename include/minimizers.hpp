#pragma once

#include "util.hpp"

namespace sshash {

struct minimizers {
    template <typename ForwardIterator>
    void build(ForwardIterator begin, uint64_t size, build_configuration const& build_config) {
        pthash::build_configuration mphf_config;
        mphf_config.c = 6.0;
        mphf_config.alpha = 0.94;
        mphf_config.seed = util::get_seed_for_hash_function(build_config);
        mphf_config.minimal_output = true;
        mphf_config.verbose_output = false;
        mphf_config.num_threads = std::thread::hardware_concurrency();
        mphf_config.num_partitions = 4 * mphf_config.num_threads;

        if (size / mphf_config.num_partitions < pthash::constants::min_partition_size) {
            mphf_config.num_partitions =
                std::max<uint64_t>(1, size / (2 * pthash::constants::min_partition_size));
        }

        if (build_config.verbose) {
            std::cout << "building minimizers MPHF with " << mphf_config.num_threads
                      << " threads and " << mphf_config.num_partitions << " partitions..."
                      << std::endl;
        }

        mphf_config.ram = 4 * essentials::GB;
        mphf_config.tmp_dir = build_config.tmp_dirname;
        m_mphf.build_in_external_memory(begin, size, mphf_config);
    }

    uint64_t lookup(uint64_t uint64_minimizer) const {
        uint64_t bucket_id = m_mphf(uint64_minimizer);
        return bucket_id;
    }

    uint64_t size() const { return m_mphf.num_keys(); }
    uint64_t num_bits() const { return m_mphf.num_bits(); }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_mphf);
    }

private:
    minimizers_pthash_type m_mphf;
};

}  // namespace sshash