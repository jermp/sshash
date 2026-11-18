#pragma once

namespace sshash {

struct minimizers_control_map  //
{
    template <typename ForwardIterator>
    void build(ForwardIterator begin, const uint64_t size,  //
               build_configuration const& build_config)     //
    {
        phastr::build_configuration mphf_build_config;
        mphf_build_config.seed = util::get_seed_for_hash_function(build_config);
        mphf_build_config.num_threads = build_config.num_threads;
        mphf_build_config.tmp_dir = build_config.tmp_dirname;

        if (build_config.verbose) {
            std::cout << "building minimizers MPHF (PHastR) with " << mphf_build_config.num_threads
                      << " threads in external memory..." << std::endl;
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

    minimizers_mphf_type mphf;
    bits::compact_vector control_codewords;

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.mphf);
        visitor.visit(t.control_codewords);
    }
};

}  // namespace sshash
