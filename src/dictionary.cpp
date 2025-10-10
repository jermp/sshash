#include "include/dictionary.hpp"

#include "include/builder/util.hpp"

namespace sshash {

template <class kmer_t>
lookup_result dictionary<kmer_t>::lookup_uint_regular(kmer_t uint_kmer) const {
    auto mini_info = util::compute_minimizer(uint_kmer, m_k, m_m, m_hasher);
    return lookup_uint_regular(uint_kmer, mini_info);
}

template <class kmer_t>
lookup_result dictionary<kmer_t>::lookup_uint_regular(kmer_t uint_kmer,                //
                                                      minimizer_info mini_info) const  //
{
    assert(minimizer_info(mini_info.minimizer, mini_info.pos_in_kmer) ==
           util::compute_minimizer(uint_kmer, m_k, m_m, m_hasher));

    const uint64_t minimizer_id = m_minimizers.lookup(mini_info.minimizer);

    uint64_t code = m_buckets.offsets.access(minimizer_id);
    uint64_t status = code & 1;
    if (status == 0) {  // minimizer occurs once
        uint64_t offset = code >> 1;
        auto res = m_buckets.lookup_at_offset(offset, uint_kmer, mini_info, m_k, m_m);
        // res.list_size = 1;
        return res;
    }

    status = code & 3;
    code >>= 2;
    if (status == 1) {  // minimizer occurs more than once, but is not part of the skew index
        constexpr uint64_t mask = (uint64_t(1) << constants::min_l) - 1;
        uint64_t list_size = (code & mask) + 2;
        uint64_t list_id = code >> constants::min_l;
        assert(list_size < m_buckets.start_lists_of_size.size());
        uint64_t begin = m_buckets.start_lists_of_size[list_size] + list_id * list_size;
        uint64_t end = begin + list_size;
        auto res = m_buckets.lookup(begin, end, uint_kmer, mini_info, m_k, m_m);
        // res.list_size = list_size;
        return res;
    }

    // minimizer is part of the skew index
    assert(status == 3);
    uint64_t partition_id = code & 7;
    uint64_t begin = code >> 3;
    uint64_t pos_in_bucket = m_skew_index.lookup(uint_kmer, partition_id);
    uint64_t offset = m_buckets.offsets3.access(begin + pos_in_bucket);
    auto res = m_buckets.lookup_at_offset(offset, uint_kmer, mini_info, m_k, m_m);
    /*
        The function `lookup_at_offset` determines if the minimizer is found at the given
        `offset`, not whether the minimizer does not appear at all.
        It can happen that the `mini_info.minimizer` appears somewhere but not at the
        computed `offset` because `pos_in_bucket` might be larger than the size of
        the bucket (which we do not know for minimizers in the skew index).
        Since for streaming queries we keep track of the presence
        of minimizers, only in this special case we set `res.minimizer_found` to true
        to indicate that we do not know whether the minimizer appears in the index or not.
    */
    res.minimizer_found = true;
    return res;
}

template <class kmer_t>
lookup_result dictionary<kmer_t>::lookup_uint_canonical(kmer_t uint_kmer) const  //
{
    kmer_t uint_kmer_rc = uint_kmer;
    uint_kmer_rc.reverse_complement_inplace(m_k);
    auto mini_info = util::compute_minimizer(uint_kmer, m_k, m_m, m_hasher);
    auto mini_info_rc = util::compute_minimizer(uint_kmer_rc, m_k, m_m, m_hasher);
    if (mini_info.minimizer < mini_info_rc.minimizer) {
        return lookup_uint_canonical(uint_kmer, uint_kmer_rc, mini_info);
    } else if (mini_info_rc.minimizer < mini_info.minimizer) {
        return lookup_uint_canonical(uint_kmer, uint_kmer_rc, mini_info_rc);
    } else {
        auto res = lookup_uint_canonical(uint_kmer, uint_kmer_rc, mini_info);
        if (res.kmer_id == constants::invalid_uint64) {
            res = lookup_uint_canonical(uint_kmer, uint_kmer_rc, mini_info_rc);
        }
        return res;
    }
}

template <class kmer_t>
lookup_result dictionary<kmer_t>::lookup_uint_canonical(kmer_t uint_kmer, kmer_t uint_kmer_rc,
                                                        minimizer_info mini_info) const  //
{
    assert(mini_info.minimizer ==
           std::min(util::compute_minimizer(uint_kmer, m_k, m_m, m_hasher).minimizer,
                    util::compute_minimizer(uint_kmer_rc, m_k, m_m, m_hasher).minimizer));

    const uint64_t minimizer_id = m_minimizers.lookup(mini_info.minimizer);

    uint64_t code = m_buckets.offsets.access(minimizer_id);
    uint64_t status = code & 1;
    if (status == 0) {  // minimizer occurs once
        uint64_t offset = code >> 1;
        auto res = m_buckets.lookup_canonical_at_offset(          //
            offset, uint_kmer, uint_kmer_rc, mini_info, m_k, m_m  //
        );
        // res.list_size = 1;
        return res;
    }

    status = code & 3;
    code >>= 2;
    if (status == 1) {  // minimizer occurs more than once, but is not part of the skew index
        constexpr uint64_t mask = (uint64_t(1) << constants::min_l) - 1;
        uint64_t list_size = (code & mask) + 2;
        uint64_t list_id = code >> constants::min_l;
        assert(list_size < m_buckets.start_lists_of_size.size());
        uint64_t begin = m_buckets.start_lists_of_size[list_size] + list_id * list_size;
        uint64_t end = begin + list_size;
        auto res = m_buckets.lookup_canonical(                        //
            begin, end, uint_kmer, uint_kmer_rc, mini_info, m_k, m_m  //
        );
        // res.list_size = list_size;
        return res;
    }

    // minimizer is part of the skew index
    assert(status == 3);
    uint64_t partition_id = code & 7;
    uint64_t begin = code >> 3;
    const auto uint_kmer_canon = std::min(uint_kmer, uint_kmer_rc);
    uint64_t pos_in_bucket = m_skew_index.lookup(uint_kmer_canon, partition_id);
    uint64_t offset = m_buckets.offsets3.access(begin + pos_in_bucket);
    auto res = m_buckets.lookup_canonical_at_offset(          //
        offset, uint_kmer, uint_kmer_rc, mini_info, m_k, m_m  //
    );
    res.minimizer_found = true;  // see note above in the function `lookup_uint_regular`
    return res;
}

template <class kmer_t>
lookup_result dictionary<kmer_t>::lookup(char const* string_kmer,
                                         bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(string_kmer, m_k);
    /*
        SIMD here does not help, as expected, because it is only used at the
        beginning of each query. To be useful, we would need to process a
        batch of random lookup queries and execute this preliminary step
        for all queries in a first pass, then invoke
        `lookup_uint(uint_kmer, check_reverse_complement)` directly.
    */
    // __m256i v = _mm256_loadu_si256(reinterpret_cast<__m256i const*>(string_kmer));
    // uint64_t word = pack2bits_shift1(v);
    // kmer_t uint_kmer(word);
    return lookup_uint(uint_kmer, check_reverse_complement);
}
template <class kmer_t>
lookup_result dictionary<kmer_t>::lookup_uint(kmer_t uint_kmer,
                                              bool check_reverse_complement) const {
    if (m_canonical) return lookup_uint_canonical(uint_kmer);
    auto res = lookup_uint_regular(uint_kmer);
    assert(res.kmer_orientation == constants::forward_orientation);
    if (check_reverse_complement and res.kmer_id == constants::invalid_uint64) {
        kmer_t uint_kmer_rc = uint_kmer;
        uint_kmer_rc.reverse_complement_inplace(m_k);
        res = lookup_uint_regular(uint_kmer_rc);
        res.kmer_orientation = constants::backward_orientation;
    }
    return res;
}

template <class kmer_t>
bool dictionary<kmer_t>::is_member(char const* string_kmer, bool check_reverse_complement) const {
    return lookup(string_kmer, check_reverse_complement) != constants::invalid_uint64;
}
template <class kmer_t>
bool dictionary<kmer_t>::is_member_uint(kmer_t uint_kmer, bool check_reverse_complement) const {
    return lookup_uint(uint_kmer, check_reverse_complement) != constants::invalid_uint64;
}

template <class kmer_t>
void dictionary<kmer_t>::access(uint64_t kmer_id, char* string_kmer) const {
    assert(kmer_id < num_kmers());
    m_buckets.access(kmer_id, string_kmer, m_k);
}

template <class kmer_t>
uint64_t dictionary<kmer_t>::weight(uint64_t kmer_id) const {
    assert(kmer_id < num_kmers());
    return m_weights.weight(kmer_id);
}

template <class kmer_t>
uint64_t dictionary<kmer_t>::contig_size(uint64_t contig_id) const {
    assert(contig_id < num_strings());
    auto [begin, end] = m_buckets.contig_offsets(contig_id);
    uint64_t contig_length = end - begin;
    assert(contig_length >= m_k);
    return contig_length - m_k + 1;
}

template <class kmer_t>
void dictionary<kmer_t>::forward_neighbours(kmer_t suffix, neighbourhood<kmer_t>& res,
                                            bool check_reverse_complement) const {
    for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
        kmer_t new_kmer = suffix;
        new_kmer.set(m_k - 1, kmer_t::char_to_uint(kmer_t::alphabet[i]));
        res.forward[i] = lookup_uint(new_kmer, check_reverse_complement);
    }
}
template <class kmer_t>
void dictionary<kmer_t>::backward_neighbours(kmer_t prefix, neighbourhood<kmer_t>& res,
                                             bool check_reverse_complement) const {
    for (size_t i = 0; i < kmer_t::alphabet_size; i++) {
        kmer_t new_kmer = prefix;
        new_kmer.set(0, kmer_t::char_to_uint(kmer_t::alphabet[i]));
        res.backward[i] = lookup_uint(new_kmer, check_reverse_complement);
    }
}

template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::kmer_forward_neighbours(
    char const* string_kmer, bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(string_kmer, m_k);
    return kmer_forward_neighbours(uint_kmer, check_reverse_complement);
}

template <class kmer_t>
kmer_t dictionary<kmer_t>::get_suffix(kmer_t kmer) const {
    kmer_t suffix = kmer;
    suffix.drop_char();
    return suffix;
}
template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::kmer_forward_neighbours(
    kmer_t uint_kmer, bool check_reverse_complement) const {
    neighbourhood<kmer_t> res;
    forward_neighbours(get_suffix(uint_kmer), res, check_reverse_complement);
    return res;
}

template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::kmer_backward_neighbours(
    char const* string_kmer, bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(string_kmer, m_k);
    return kmer_backward_neighbours(uint_kmer, check_reverse_complement);
}

template <class kmer_t>
kmer_t dictionary<kmer_t>::get_prefix(kmer_t kmer) const {
    kmer_t prefix = kmer;
    prefix.pad_char();
    prefix.take_chars(m_k);
    return prefix;
}

template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::kmer_backward_neighbours(
    kmer_t uint_kmer, bool check_reverse_complement) const {
    neighbourhood<kmer_t> res;
    backward_neighbours(get_prefix(uint_kmer), res, check_reverse_complement);
    return res;
}

template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::kmer_neighbours(char const* string_kmer,
                                                          bool check_reverse_complement) const {
    kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(string_kmer, m_k);
    return kmer_neighbours(uint_kmer, check_reverse_complement);
}

template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::kmer_neighbours(kmer_t uint_kmer,
                                                          bool check_reverse_complement) const {
    neighbourhood<kmer_t> res;
    forward_neighbours(get_suffix(uint_kmer), res, check_reverse_complement);
    backward_neighbours(get_prefix(uint_kmer), res, check_reverse_complement);
    return res;
}

template <class kmer_t>
neighbourhood<kmer_t> dictionary<kmer_t>::contig_neighbours(uint64_t contig_id,
                                                            bool check_reverse_complement) const {
    assert(contig_id < num_strings());
    neighbourhood<kmer_t> res;
    kmer_t suffix = m_buckets.contig_suffix(contig_id, m_k);
    forward_neighbours(suffix, res, check_reverse_complement);
    kmer_t prefix = m_buckets.contig_prefix(contig_id, m_k);
    prefix.pad_char();
    backward_neighbours(prefix, res, check_reverse_complement);
    return res;
}

template <class kmer_t>
uint64_t dictionary<kmer_t>::num_bits() const {
    return 8 * (sizeof(m_vnum) + sizeof(m_num_kmers) + sizeof(m_num_strings) + sizeof(m_hasher) +
                sizeof(m_k) + sizeof(m_m) + sizeof(m_canonical)) +
           m_minimizers.num_bits() + m_buckets.num_bits() + m_skew_index.num_bits() +
           m_weights.num_bits();
}

}  // namespace sshash
