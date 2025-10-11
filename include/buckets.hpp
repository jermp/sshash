#pragma once

#include "util.hpp"
#include "kmer_iterator.hpp"
#include "endpoints.hpp"

namespace sshash {

template <class kmer_t, class Endpoints>
struct buckets  //
{
    /* Return where the contig begins and ends in strings. */
    std::pair<uint64_t, uint64_t>  // [begin, end)
    contig_offsets(const uint64_t contig_id) const {
        uint64_t begin = strings_endpoints.access(contig_id);
        uint64_t end = strings_endpoints.access(contig_id + 1);
        assert(end > begin);
        return {begin, end};
    }

    kmer_t contig_prefix(const uint64_t contig_id, const uint64_t k) const {
        uint64_t contig_begin = strings_endpoints.access(contig_id);
        return util::read_kmer_at<kmer_t>(strings, k - 1, kmer_t::bits_per_char * contig_begin);
    }

    kmer_t contig_suffix(const uint64_t contig_id, const uint64_t k) const {
        uint64_t contig_end = strings_endpoints.access(contig_id + 1);
        return util::read_kmer_at<kmer_t>(strings, k - 1,
                                          kmer_t::bits_per_char * (contig_end - k + 1));
    }

    lookup_result lookup(const uint64_t begin, const uint64_t end,           //
                         const kmer_t kmer, const minimizer_info mini_info,  //
                         const uint64_t k, const uint64_t m) const           //
    {
        /* check minimizer first */
        uint64_t minimizer_offset = offsets2.access(begin);
        auto p = strings_endpoints.decode(minimizer_offset);
        uint64_t read_mmer = uint64_t(
            util::read_kmer_at<kmer_t>(strings, m, kmer_t::bits_per_char * p.absolute_offset));
        if (read_mmer != mini_info.minimizer) return lookup_result(false);

        auto res = lookup_at_offset_no_check_minimizer(p, kmer, mini_info, k);
        if (res.kmer_id != constants::invalid_uint64) {
            assert(is_valid(res));
            return res;
        }

        for (uint64_t i = begin + 1; i < end; ++i) {
            minimizer_offset = offsets2.access(i);
            p = strings_endpoints.decode(minimizer_offset);
            res = lookup_at_offset_no_check_minimizer(p, kmer, mini_info, k);
            if (res.kmer_id != constants::invalid_uint64) {
                assert(is_valid(res));
                return res;
            }
        }

        return lookup_result();
    }

    lookup_result lookup_at_offset_no_check_minimizer(endpoints::decoded_offset p,     //
                                                      const kmer_t kmer,               //
                                                      const minimizer_info mini_info,  //
                                                      const uint64_t k) const          //
    {
        auto res = strings_endpoints.offset_to_id(p, mini_info.pos_in_kmer, k);
        if (res.kmer_id != constants::invalid_uint64) {
            uint64_t kmer_offset = res.kmer_offset(k);
            if (kmer_offset + k - 1 < res.contig_end(k)) {
                auto read_kmer =
                    util::read_kmer_at<kmer_t>(strings, k, kmer_t::bits_per_char * kmer_offset);
                if (read_kmer == kmer) return res;
            }
        }
        return lookup_result();
    }

    lookup_result lookup_at_offset(const uint64_t minimizer_offset,  //
                                   const kmer_t kmer,                //
                                   const minimizer_info mini_info,   //
                                   const uint64_t k,                 //
                                   const uint64_t m) const           //
    {
        /* check minimizer first */
        auto p = strings_endpoints.decode(minimizer_offset);
        uint64_t read_mmer = uint64_t(
            util::read_kmer_at<kmer_t>(strings, m, kmer_t::bits_per_char * p.absolute_offset));
        if (read_mmer != mini_info.minimizer) return lookup_result(false);
        return lookup_at_offset_no_check_minimizer(p, kmer, mini_info, k);
    }

    lookup_result lookup_canonical(const uint64_t begin, const uint64_t end,  //
                                   const kmer_t kmer, const kmer_t kmer_rc,   //
                                   const minimizer_info mini_info,            //
                                   const uint64_t k,                          //
                                   const uint64_t m) const                    //
    {
        /* check minimizer first */
        uint64_t minimizer_offset = offsets2.access(begin);
        auto p = strings_endpoints.decode(minimizer_offset);
        uint64_t read_mmer = uint64_t(
            util::read_kmer_at<kmer_t>(strings, m, kmer_t::bits_per_char * p.absolute_offset));
        uint64_t minimizer_rc = 0;
        {
            auto tmp = kmer_t(mini_info.minimizer);
            tmp.reverse_complement_inplace(m);
            minimizer_rc = uint64_t(tmp);
        }
        if (read_mmer != mini_info.minimizer and read_mmer != minimizer_rc) {
            return lookup_result(false);
        }

        auto res = lookup_canonical_at_offset_no_check_minimizer(p, kmer, kmer_rc, mini_info, k, m);
        if (res.kmer_id != constants::invalid_uint64) {
            assert(is_valid(res));
            return res;
        }

        for (uint64_t i = begin + 1; i < end; ++i) {
            minimizer_offset = offsets2.access(i);
            p = strings_endpoints.decode(minimizer_offset);
            res = lookup_canonical_at_offset_no_check_minimizer(p, kmer, kmer_rc, mini_info, k, m);
            if (res.kmer_id != constants::invalid_uint64) {
                assert(is_valid(res));
                return res;
            }
        }
        return lookup_result();
    }

    lookup_result lookup_canonical_at_offset_no_check_minimizer(endpoints::decoded_offset p,  //
                                                                const kmer_t kmer,
                                                                const kmer_t kmer_rc,            //
                                                                const minimizer_info mini_info,  //
                                                                const uint64_t k,                //
                                                                const uint64_t m) const          //
    {
        uint64_t pos_in_kmer = mini_info.pos_in_kmer;
        auto res = check_offset(p, pos_in_kmer, kmer, kmer_rc, k);
        if (res.kmer_id != constants::invalid_uint64) {
            assert(is_valid(res));
            return res;
        }
        pos_in_kmer = k - m - mini_info.pos_in_kmer;
        return check_offset(p, pos_in_kmer, kmer, kmer_rc, k);
    }

    lookup_result lookup_canonical_at_offset(const uint64_t minimizer_offset,  //
                                             const kmer_t kmer,                //
                                             const kmer_t kmer_rc,             //
                                             const minimizer_info mini_info,   //
                                             const uint64_t k,                 //
                                             const uint64_t m) const           //
    {
        /* check minimizer first */
        auto p = strings_endpoints.decode(minimizer_offset);
        uint64_t read_mmer = uint64_t(
            util::read_kmer_at<kmer_t>(strings, m, kmer_t::bits_per_char * p.absolute_offset));
        uint64_t minimizer_rc = 0;
        {
            auto tmp = kmer_t(mini_info.minimizer);
            tmp.reverse_complement_inplace(m);
            minimizer_rc = uint64_t(tmp);
        }
        if (read_mmer != mini_info.minimizer and read_mmer != minimizer_rc) {
            return lookup_result(false);
        }
        return lookup_canonical_at_offset_no_check_minimizer(p, kmer, kmer_rc, mini_info, k, m);
    }

    lookup_result check_offset(endpoints::decoded_offset p,              //
                               const uint64_t pos_in_kmer,               //
                               const kmer_t kmer, const kmer_t kmer_rc,  //
                               const uint64_t k) const                   //
    {
        auto res = strings_endpoints.offset_to_id(p, pos_in_kmer, k);
        if (res.kmer_id != constants::invalid_uint64) {
            uint64_t kmer_offset = res.kmer_offset(k);
            if (kmer_offset + k - 1 < res.contig_end(k)) {
                auto read_kmer =
                    util::read_kmer_at<kmer_t>(strings, k, kmer_t::bits_per_char * kmer_offset);
                if (read_kmer == kmer) return res;
                if (read_kmer == kmer_rc) {
                    res.kmer_orientation = constants::backward_orientation;
                    return res;
                }
            }
        }
        return lookup_result();
    }

    void access(const uint64_t kmer_id, char* string_kmer, const uint64_t k) const {
        auto [_, offset] = strings_endpoints.id_to_offset(kmer_id, k);
        auto read_kmer = util::read_kmer_at<kmer_t>(strings, k, kmer_t::bits_per_char * offset);
        util::uint_kmer_to_string(read_kmer, string_kmer, k);
    }

    struct iterator {
        iterator() {}

        iterator(buckets const* ptr,                                        //
                 const uint64_t begin_kmer_id, const uint64_t end_kmer_id,  // [begin,end)
                 const uint64_t k)
            : m_buckets(ptr)
            , m_begin_kmer_id(begin_kmer_id)
            , m_end_kmer_id(end_kmer_id)
            , m_k(k)
            , m_it(ptr->strings, m_k)  //
        {
            auto [pos, val] = m_buckets->strings_endpoints.id_to_offset(m_begin_kmer_id, k);
            m_offset = val;
            m_strings_endpoints_it = m_buckets->strings_endpoints.get_iterator_at(pos);
            assert(m_strings_endpoints_it.value() > m_offset);
            next_piece();
            m_ret.second.resize(m_k, 0);
        }

        bool has_next() const { return m_begin_kmer_id != m_end_kmer_id; }

        std::pair<uint64_t, std::string> next() {
            if (m_offset == m_next_offset - m_k + 1) {
                m_offset = m_next_offset;
                next_piece();
            }
            m_ret.first = m_begin_kmer_id;
            if (m_clear) {
                util::uint_kmer_to_string(m_it.get(), m_ret.second.data(), m_k);
                assert(kmer_t::bits_per_char * m_offset == m_it.position());
                m_it.at(kmer_t::bits_per_char * (m_offset + m_k));
            } else {
                memmove(m_ret.second.data(), m_ret.second.data() + 1, m_k - 1);
                m_ret.second[m_k - 1] = kmer_t::uint64_to_char(m_it.get_next_char());
            }
            m_clear = false;
            ++m_begin_kmer_id;
            ++m_offset;
            return m_ret;
        }

    private:
        std::pair<uint64_t, std::string> m_ret;
        buckets const* m_buckets;
        uint64_t m_begin_kmer_id, m_end_kmer_id;
        uint64_t m_k;
        uint64_t m_offset;
        uint64_t m_next_offset;
        kmer_iterator<kmer_t, bits::bit_vector> m_it;
        typename Endpoints::iterator m_strings_endpoints_it;
        bool m_clear;

        void next_piece() {
            m_it.at(kmer_t::bits_per_char * m_offset);
            m_next_offset = m_strings_endpoints_it.value();
            assert(m_next_offset > m_offset);
            m_clear = true;
            m_strings_endpoints_it.next();
        }
    };

    iterator at(const uint64_t begin_kmer_id, const uint64_t end_kmer_id, const uint64_t k) const {
        return iterator(this, begin_kmer_id, end_kmer_id, k);
    }

    uint64_t num_bits() const {
        return 8 * (strings_endpoints.num_bytes() + essentials::vec_bytes(start_lists_of_size) +  //
                    offsets.num_bytes() + offsets2.num_bytes() + offsets3.num_bytes() +           //
                    strings.num_bytes());                                                         //
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

    Endpoints strings_endpoints;
    // bits::endpoints_sequence<> strings_endpoints;
    // bits::compact_vector strings_endpoints;

    std::vector<uint32_t> start_lists_of_size;
    bits::compact_vector offsets;
    bits::compact_vector offsets2;
    bits::compact_vector offsets3;

    bits::bit_vector strings;

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.strings_endpoints);
        visitor.visit(t.start_lists_of_size);
        visitor.visit(t.offsets);
        visitor.visit(t.offsets2);
        visitor.visit(t.offsets3);
        visitor.visit(t.strings);
    }

    bool is_valid(lookup_result res) const {
        return res.contig_size != constants::invalid_uint64 and  //
               res.kmer_id_in_contig < res.contig_size and       //
               res.contig_id < strings_endpoints.size();         //
    }
};

}  // namespace sshash