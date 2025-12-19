#include "include/dictionary.hpp"

#include "include/builder/util.hpp"

namespace sshash {

template <typename Kmer, typename Offsets>
lookup_result dictionary<Kmer, Offsets>::lookup_regular(const Kmer uint_kmer) const {
    auto mini_info = util::compute_minimizer(uint_kmer, m_k, m_m, m_hasher);
    return lookup_regular(uint_kmer, mini_info);
}

template <typename Kmer, typename Offsets>
lookup_result dictionary<Kmer, Offsets>::lookup_regular(const Kmer uint_kmer,                  //
                                                        const minimizer_info mini_info) const  //
{
    assert(minimizer_info(mini_info.minimizer, mini_info.pos_in_kmer) ==
           util::compute_minimizer(uint_kmer, m_k, m_m, m_hasher));

    auto it = m_ssi.lookup(uint_kmer, mini_info);
    return m_spss.lookup_regular(it, uint_kmer, mini_info);
}

template <typename Kmer, typename Offsets>
lookup_result dictionary<Kmer, Offsets>::lookup_canonical(Kmer uint_kmer) const  //
{
    Kmer uint_kmer_rc = uint_kmer;
    uint_kmer_rc.reverse_complement_inplace(m_k);
    auto mini_info = util::compute_minimizer(uint_kmer, m_k, m_m, m_hasher);
    auto mini_info_rc = util::compute_minimizer(uint_kmer_rc, m_k, m_m, m_hasher);
    if (mini_info.minimizer < mini_info_rc.minimizer) {
        return lookup_canonical(uint_kmer, uint_kmer_rc, mini_info);
    } else if (mini_info_rc.minimizer < mini_info.minimizer) {
        return lookup_canonical(uint_kmer, uint_kmer_rc, mini_info_rc);
    } else {
        auto res = lookup_canonical(uint_kmer, uint_kmer_rc, mini_info);
        if (res.kmer_id == constants::invalid_uint64) {
            res = lookup_canonical(uint_kmer, uint_kmer_rc, mini_info_rc);
        }
        return res;
    }
}

template <typename Kmer, typename Offsets>
lookup_result dictionary<Kmer, Offsets>::lookup_canonical(const Kmer uint_kmer,                  //
                                                          const Kmer uint_kmer_rc,               //
                                                          const minimizer_info mini_info) const  //
{
    assert(mini_info.minimizer ==
           std::min(util::compute_minimizer(uint_kmer, m_k, m_m, m_hasher).minimizer,
                    util::compute_minimizer(uint_kmer_rc, m_k, m_m, m_hasher).minimizer));

    const Kmer uint_kmer_canon = std::min(uint_kmer, uint_kmer_rc);
    auto it = m_ssi.lookup(uint_kmer_canon, mini_info);
    return m_spss.lookup_canonical(it, uint_kmer, uint_kmer_rc, mini_info);
}

template <typename Kmer, typename Offsets>
lookup_result dictionary<Kmer, Offsets>::lookup(char const* string_kmer,
                                                bool check_reverse_complement) const {
    Kmer uint_kmer = util::string_to_uint_kmer<Kmer>(string_kmer, m_k);
    return lookup(uint_kmer, check_reverse_complement);
}
template <typename Kmer, typename Offsets>
lookup_result dictionary<Kmer, Offsets>::lookup(Kmer uint_kmer,
                                                bool check_reverse_complement) const  //
{
    if (m_canonical) return lookup_canonical(uint_kmer);
    auto res = lookup_regular(uint_kmer);
    assert(res.kmer_orientation == constants::forward_orientation);
    if (check_reverse_complement and res.kmer_id == constants::invalid_uint64) {
        Kmer uint_kmer_rc = uint_kmer;
        uint_kmer_rc.reverse_complement_inplace(m_k);
        res = lookup_regular(uint_kmer_rc);
        res.kmer_orientation = constants::backward_orientation;
    }
    return res;
}

template <typename Kmer, typename Offsets>
bool dictionary<Kmer, Offsets>::is_member(char const* string_kmer,
                                          bool check_reverse_complement) const {
    return lookup(string_kmer, check_reverse_complement) != constants::invalid_uint64;
}
template <typename Kmer, typename Offsets>
bool dictionary<Kmer, Offsets>::is_member(Kmer uint_kmer, bool check_reverse_complement) const {
    return lookup(uint_kmer, check_reverse_complement) != constants::invalid_uint64;
}

template <typename Kmer, typename Offsets>
void dictionary<Kmer, Offsets>::access(uint64_t kmer_id, char* string_kmer) const {
    assert(kmer_id < num_kmers());
    m_spss.access(kmer_id, string_kmer);
}

template <typename Kmer, typename Offsets>
uint64_t dictionary<Kmer, Offsets>::weight(uint64_t kmer_id) const {
    assert(kmer_id < num_kmers());
    return m_weights.weight(kmer_id);
}

template <typename Kmer, typename Offsets>
uint64_t dictionary<Kmer, Offsets>::string_size(uint64_t string_id) const {
    assert(string_id < num_strings());
    auto [begin, end] = m_spss.string_offsets(string_id);
    uint64_t string_length = end - begin;
    assert(string_length >= m_k);
    return string_length - m_k + 1;
}

template <typename Kmer, typename Offsets>
void dictionary<Kmer, Offsets>::forward_neighbours(Kmer suffix, neighbourhood<Kmer>& res,
                                                   bool check_reverse_complement) const {
    for (size_t i = 0; i < Kmer::alphabet_size; i++) {
        Kmer new_kmer = suffix;
        new_kmer.set(m_k - 1, Kmer::char_to_uint(Kmer::alphabet[i]));
        res.forward[i] = lookup(new_kmer, check_reverse_complement);
    }
}
template <typename Kmer, typename Offsets>
void dictionary<Kmer, Offsets>::backward_neighbours(Kmer prefix, neighbourhood<Kmer>& res,
                                                    bool check_reverse_complement) const {
    for (size_t i = 0; i < Kmer::alphabet_size; i++) {
        Kmer new_kmer = prefix;
        new_kmer.set(0, Kmer::char_to_uint(Kmer::alphabet[i]));
        res.backward[i] = lookup(new_kmer, check_reverse_complement);
    }
}

template <typename Kmer, typename Offsets>
neighbourhood<Kmer> dictionary<Kmer, Offsets>::kmer_forward_neighbours(
    char const* string_kmer, bool check_reverse_complement) const {
    Kmer uint_kmer = util::string_to_uint_kmer<Kmer>(string_kmer, m_k);
    return kmer_forward_neighbours(uint_kmer, check_reverse_complement);
}

template <typename Kmer, typename Offsets>
Kmer dictionary<Kmer, Offsets>::get_suffix(Kmer kmer) const {
    Kmer suffix = kmer;
    suffix.drop_char();
    return suffix;
}
template <typename Kmer, typename Offsets>
neighbourhood<Kmer> dictionary<Kmer, Offsets>::kmer_forward_neighbours(
    Kmer uint_kmer, bool check_reverse_complement) const {
    neighbourhood<Kmer> res;
    forward_neighbours(get_suffix(uint_kmer), res, check_reverse_complement);
    return res;
}

template <typename Kmer, typename Offsets>
neighbourhood<Kmer> dictionary<Kmer, Offsets>::kmer_backward_neighbours(
    char const* string_kmer, bool check_reverse_complement) const {
    Kmer uint_kmer = util::string_to_uint_kmer<Kmer>(string_kmer, m_k);
    return kmer_backward_neighbours(uint_kmer, check_reverse_complement);
}

template <typename Kmer, typename Offsets>
Kmer dictionary<Kmer, Offsets>::get_prefix(Kmer kmer) const {
    Kmer prefix = kmer;
    prefix.pad_char();
    prefix.take_chars(m_k);
    return prefix;
}

template <typename Kmer, typename Offsets>
neighbourhood<Kmer> dictionary<Kmer, Offsets>::kmer_backward_neighbours(
    Kmer uint_kmer, bool check_reverse_complement) const {
    neighbourhood<Kmer> res;
    backward_neighbours(get_prefix(uint_kmer), res, check_reverse_complement);
    return res;
}

template <typename Kmer, typename Offsets>
neighbourhood<Kmer> dictionary<Kmer, Offsets>::kmer_neighbours(
    char const* string_kmer, bool check_reverse_complement) const {
    Kmer uint_kmer = util::string_to_uint_kmer<Kmer>(string_kmer, m_k);
    return kmer_neighbours(uint_kmer, check_reverse_complement);
}

template <typename Kmer, typename Offsets>
neighbourhood<Kmer> dictionary<Kmer, Offsets>::kmer_neighbours(
    Kmer uint_kmer, bool check_reverse_complement) const {
    neighbourhood<Kmer> res;
    forward_neighbours(get_suffix(uint_kmer), res, check_reverse_complement);
    backward_neighbours(get_prefix(uint_kmer), res, check_reverse_complement);
    return res;
}

template <typename Kmer, typename Offsets>
neighbourhood<Kmer> dictionary<Kmer, Offsets>::string_neighbours(
    uint64_t string_id, bool check_reverse_complement) const {
    assert(string_id < num_strings());
    neighbourhood<Kmer> res;
    Kmer suffix = m_spss.string_suffix(string_id);
    forward_neighbours(suffix, res, check_reverse_complement);
    Kmer prefix = m_spss.string_prefix(string_id);
    prefix.pad_char();
    backward_neighbours(prefix, res, check_reverse_complement);
    return res;
}

template <typename Kmer, typename Offsets>
uint64_t dictionary<Kmer, Offsets>::num_bits() const {
    return 8 * (sizeof(m_vnum) + sizeof(m_num_kmers) + sizeof(m_num_strings) + sizeof(m_k) +
                sizeof(m_m) + sizeof(m_canonical) + sizeof(m_hasher)) +
           m_spss.num_bits() + m_ssi.num_bits() + m_weights.num_bits();
}

}  // namespace sshash
