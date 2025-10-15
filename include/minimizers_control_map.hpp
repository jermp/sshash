#pragma once

namespace sshash {

struct minimizers_control_map  //
{
    template <typename ForwardIterator>
    void build(ForwardIterator begin, const uint64_t size,  //
               build_configuration const& build_config)     //
    {
        pthash::build_configuration mphf_build_config;
        mphf_build_config.lambda = build_config.lambda;
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

        mphf.build_in_external_memory(begin, size, mphf_build_config);
    }

    uint64_t lookup(uint64_t minimizer) const {
        uint64_t minimizer_id = mphf(minimizer);
        return control_codewords.access(minimizer_id);
    }

    uint64_t size() const { return mphf.num_keys(); }

    uint64_t num_bits() const { return mphf.num_bits() + 8 * control_codewords.num_bytes(); }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

    minimizers_pthash_type mphf;
    bits::compact_vector control_codewords;

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.mphf);
        visitor.visit(t.control_codewords);
    }
};

}  // namespace sshash