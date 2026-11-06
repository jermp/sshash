#include "include/dictionary.hpp"

#include "include/builder/util.hpp"

namespace sshash {

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
    return m_ssi.strand_specific
        ? lookup_mode_impl<true>(uint_kmer.bits, check_reverse_complement)
        : lookup_mode_impl<false>(uint_kmer.bits, check_reverse_complement);
}

template <typename Kmer, typename Offsets>
template <bool Oriented>
[[gnu::flatten]] lookup_result dictionary<Kmer, Offsets>::lookup_mode_impl(
    decltype(Kmer::bits) bits, bool check_reverse_complement) const {
    Kmer uint_kmer;
    uint_kmer.bits = bits;
    Kmer uint_kmer_rc = uint_kmer;
    uint_kmer_rc.reverse_complement_inplace(m_k);
    /* the reverse complement is needed anyway to tell the two orientations
       apart, so computing the minimizer from it costs nothing extra */
    auto mini_info = util::compute_minimizer(uint_kmer, uint_kmer_rc, m_k, m_m, m_hasher);
    auto res = lookup_bucket_impl<Oriented>(uint_kmer, uint_kmer_rc, mini_info, nullptr, !check_reverse_complement);
    if (!check_reverse_complement and res.kmer_orientation == constants::backward_orientation) {  //
        return lookup_result();
    }
    return res;
}
template <typename Kmer, typename Offsets>
lookup_result dictionary<Kmer, Offsets>::lookup(const Kmer uint_kmer,                           //
                                                const Kmer uint_kmer_rc,                        //
                                                const minimizer_info mini_info,                 //
                                                typename spss_type::bucket_cache* cache, bool forward_only) const  //
{
    return m_ssi.strand_specific
        ? lookup_bucket_impl<true>(uint_kmer, uint_kmer_rc, mini_info, cache, forward_only)
        : lookup_bucket_impl<false>(uint_kmer, uint_kmer_rc, mini_info, cache, forward_only);
}

template <typename Kmer, typename Offsets>
template <bool Oriented>
lookup_result dictionary<Kmer, Offsets>::lookup_bucket_impl(
    Kmer uint_kmer, Kmer uint_kmer_rc, minimizer_info mini_info,
    typename spss_type::bucket_cache* cache, bool forward_only) const {
    assert(minimizer_info(mini_info.minimizer, mini_info.pos_in_kmer) ==
           util::compute_minimizer(uint_kmer, uint_kmer_rc, m_k, m_m, m_hasher));
    auto it = m_ssi.template lookup_impl<Oriented>(uint_kmer, uint_kmer_rc, mini_info);
    if constexpr (Oriented) {
        auto res = m_spss.lookup(it, uint_kmer, uint_kmer, mini_info, cache);
        if (res.kmer_id != constants::invalid_uint64) {
            res.kmer_orientation = constants::forward_orientation;
            return res;
        }
        if (forward_only) return res;
        auto reverse_it = m_ssi.template lookup_impl<Oriented>(uint_kmer_rc, uint_kmer, mini_info);
        return m_spss.lookup(reverse_it, uint_kmer, uint_kmer_rc, mini_info, cache);
    }
    return m_spss.lookup(it, uint_kmer, uint_kmer_rc, mini_info, cache);
}

template <typename Kmer, typename Offsets>
bool dictionary<Kmer, Offsets>::is_member(char const* string_kmer,
                                          bool check_reverse_complement) const {
    return lookup(string_kmer, check_reverse_complement).kmer_id != constants::invalid_uint64;
}
template <typename Kmer, typename Offsets>
bool dictionary<Kmer, Offsets>::is_member(Kmer uint_kmer, bool check_reverse_complement) const {
    return lookup(uint_kmer, check_reverse_complement).kmer_id != constants::invalid_uint64;
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
                sizeof(m_m) + sizeof(m_hasher) + sizeof(m_ssi.strand_specific)) +
           m_spss.num_bits() + m_ssi.num_bits() + m_weights.num_bits();
}

}  // namespace sshash
