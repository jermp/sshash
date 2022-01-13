#pragma once

#include "util.hpp"

namespace sshash {

struct minimizers {
    template <typename ForwardIterator>
    void build(ForwardIterator begin, uint64_t size) {
        util::check_hash_collision_probability(size);
        pthash::build_configuration mphf_config;
        mphf_config.c = 6.0;
        mphf_config.alpha = 0.94;
        mphf_config.seed = 1234567890;  // my favourite seed
        mphf_config.minimal_output = true;
        mphf_config.verbose_output = false;
        /*
            We use one thread here because the keys iterator is not
            a random-access iterator, just forward
        */
        mphf_config.num_threads = 1;
        m_mphf.build_in_internal_memory(begin, size, mphf_config);
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
    pthash_mphf_type m_mphf;
};

}  // namespace sshash