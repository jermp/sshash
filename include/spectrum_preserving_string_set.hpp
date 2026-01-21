#pragma once

#include "kmer_iterator.hpp"

namespace sshash {

template <typename Kmer, typename Offsets>
struct spectrum_preserving_string_set  //
{
    /* Return where the string begins and ends in `strings`. */
    std::pair<uint64_t, uint64_t>  // [begin, end)
    string_offsets(const uint64_t string_id) const {
        uint64_t begin = strings_offsets.access(string_id);
        uint64_t end = strings_offsets.access(string_id + 1);
        assert(end > begin);
        return {begin, end};
    }

    Kmer string_prefix(const uint64_t string_id) const {
        uint64_t string_begin = strings_offsets.access(string_id);
        return util::read_kmer_at<Kmer>(strings, k - 1, Kmer::bits_per_char * string_begin);
    }

    Kmer string_suffix(const uint64_t string_id) const {
        uint64_t string_end = strings_offsets.access(string_id + 1);
        return util::read_kmer_at<Kmer>(strings, k - 1, Kmer::bits_per_char * (string_end - k + 1));
    }

    template <typename Iterator>
    lookup_result lookup_regular(Iterator it,                           //
                                 const Kmer kmer,                       //
                                 const minimizer_info mini_info) const  //
    {
        const uint64_t size = it.size();
        assert(size > 0);

        static thread_local  //
            std::array<typename Offsets::decoded_offset, 1ULL << constants::min_l>
                v;

        for (uint64_t i = 0; i != size; ++i, ++it) {
            uint64_t minimizer_offset = *it;
            v[i] = strings_offsets.decode(minimizer_offset);
        }

        /* check minimizer first */
        if (uint64_t read_mmer = uint64_t(
                util::read_kmer_at<Kmer>(strings, m, Kmer::bits_per_char * v[0].absolute_offset));
            read_mmer != mini_info.minimizer)  //
        {
            /*
               The function `lookup_regular` determines if the minimizer is found at the
               offset `Kmer::bits_per_char * p.absolute_offset`, not whether the minimizer
               does not appear at all. In fact, it can happen that the minimizer appear but
               not at the specified offset, so it would be wrong to set `res.minimizer_found`
               to `false`. This can happen for HEAVYLOAD buckets only because their lookup is
               resolved via the skew index and `pos_in_bucket` might be larger than the size
               of the bucket (which we do not know for a HEAVYLOAD bucket). Since for streaming
               queries we keep track of the presence of minimizers (i.e., whether they appear
               in the index or not), only in this special case we set
               `res.minimizer_found` to `true` to indicate that we do not know whether the
               minimizer appears in the index or not.
            */
            return lookup_result(it.bucket_type() != bucket_t::HEAVYLOAD ? false : true);
        }

        lookup_result res;
        for (uint64_t i = 0; i != size; ++i) {
            if (_lookup_regular(res, v[i], kmer, mini_info)) return res;
        }

        return lookup_result();
    }

    template <typename Iterator>
    lookup_result lookup_canonical(Iterator it,                           //
                                   const Kmer kmer, const Kmer kmer_rc,   //
                                   const minimizer_info mini_info) const  //
    {
        const uint64_t size = it.size();
        assert(size > 0);

        static thread_local  //
            std::array<typename Offsets::decoded_offset, 1ULL << constants::min_l>
                v;

        for (uint64_t i = 0; i != size; ++i, ++it) {
            uint64_t minimizer_offset = *it;
            v[i] = strings_offsets.decode(minimizer_offset);
        }

        /* check minimizer first */
        if (uint64_t read_mmer = uint64_t(
                util::read_kmer_at<Kmer>(strings, m, Kmer::bits_per_char * v[0].absolute_offset));
            read_mmer != mini_info.minimizer)  //
        {
            Kmer tmp = mini_info.minimizer;
            tmp.reverse_complement_inplace(m);
            uint64_t minimizer_rc = uint64_t(tmp);
            if (read_mmer != minimizer_rc) {
                /* Same note as for the function `lookup_regular`. */
                return lookup_result(it.bucket_type() != bucket_t::HEAVYLOAD ? false : true);
            }
        }

        lookup_result res;
        for (uint64_t i = 0; i != size; ++i) {
            if (_lookup_canonical(res, v[i], kmer, kmer_rc, mini_info)) return res;
        }

        return lookup_result();
    }

    void access(const uint64_t kmer_id, char* string_kmer) const {
        auto [_, offset] = strings_offsets.id_to_offset(kmer_id, k);
        auto read_kmer = util::read_kmer_at<Kmer>(strings, k, Kmer::bits_per_char * offset);
        util::uint_kmer_to_string(read_kmer, string_kmer, k);
    }

    struct iterator {
        iterator() {}

        iterator(spectrum_preserving_string_set const* ptr,                 //
                 const uint64_t begin_kmer_id, const uint64_t end_kmer_id,  // [begin,end)
                 const uint64_t k)
            : m_ptr(ptr)
            , m_begin_kmer_id(begin_kmer_id)
            , m_end_kmer_id(end_kmer_id)
            , k(k)
            , m_it(ptr->strings, k)  //
        {
            auto [pos, val] = m_ptr->strings_offsets.id_to_offset(m_begin_kmer_id, k);
            m_offset = val;
            m_strings_offsets_it = m_ptr->strings_offsets.get_iterator_at(pos);
            assert(m_strings_offsets_it.value() > m_offset);
            next_piece();
        }

        bool has_next() const { return m_begin_kmer_id != m_end_kmer_id; }

        std::pair<uint64_t, Kmer>  // (kmer-id, encoded kmer)
        next() {
            if (m_offset == m_next_offset - k + 1) {
                m_offset = m_next_offset;
                next_piece();
            }
            m_ret.first = m_begin_kmer_id;
            if (m_clear) {
                m_ret.second = m_it.get();
                assert(Kmer::bits_per_char * m_offset == m_it.position());
                m_it.at(Kmer::bits_per_char * (m_offset + k));
            } else {
                m_ret.second.drop_char();
                m_ret.second.set(k - 1, m_it.get_next_char());
            }
            m_clear = false;
            ++m_begin_kmer_id;
            ++m_offset;
            return m_ret;
        }

    private:
        std::pair<uint64_t, Kmer> m_ret;
        spectrum_preserving_string_set const* m_ptr;
        uint64_t m_begin_kmer_id, m_end_kmer_id;
        uint64_t k;
        uint64_t m_offset, m_next_offset;
        kmer_iterator<Kmer, bits::bit_vector> m_it;
        typename Offsets::iterator m_strings_offsets_it;
        bool m_clear;

        void next_piece() {
            m_it.at(Kmer::bits_per_char * m_offset);
            m_next_offset = m_strings_offsets_it.value();
            assert(m_next_offset > m_offset);
            m_clear = true;
            m_strings_offsets_it.next();
        }
    };

    iterator at(const uint64_t begin_kmer_id, const uint64_t end_kmer_id) const {
        return iterator(this, begin_kmer_id, end_kmer_id, k);
    }

    uint64_t num_bits() const {
        return 8 * (sizeof(k) + sizeof(m) + strings_offsets.num_bytes() + strings.num_bytes());
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

    uint16_t k;
    uint16_t m;
    Offsets strings_offsets;
    bits::bit_vector strings;

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.k);
        visitor.visit(t.m);
        visitor.visit(t.strings_offsets);
        visitor.visit(t.strings);
    }

    bool _lookup_regular(lookup_result& res,                    //
                         typename Offsets::decoded_offset p,    //
                         const Kmer kmer,                       //
                         const minimizer_info mini_info) const  //
    {
        if (p.absolute_offset < mini_info.pos_in_kmer) return false;

        res.kmer_offset = p.absolute_offset - mini_info.pos_in_kmer;

        if (kmer != util::read_kmer_at<Kmer>(strings, k, Kmer::bits_per_char * res.kmer_offset)) {
            return false;
        }

        if (res.kmer_offset >= res.string_begin and res.kmer_offset < res.string_end - k + 1) {
            res.kmer_id = res.kmer_offset - res.string_id * (k - 1);     // absolute kmer id
            res.kmer_id_in_string = res.kmer_offset - res.string_begin;  // relative kmer id
        } else {
            strings_offsets.offset_to_id(res, p, k);
        }

        if (res.kmer_offset < res.string_end - k + 1) return true;
        return false;
    }

    bool _lookup_canonical(lookup_result& res,                    //
                           typename Offsets::decoded_offset p,    //
                           const Kmer kmer,                       //
                           const Kmer kmer_rc,                    //
                           const minimizer_info mini_info) const  //
    {
        uint64_t pos_in_kmer = mini_info.pos_in_kmer;
        if (__lookup_canonical(res, p, kmer, kmer_rc, pos_in_kmer)) return true;
        pos_in_kmer = k - m - mini_info.pos_in_kmer;
        return __lookup_canonical(res, p, kmer, kmer_rc, pos_in_kmer);
    }

    bool __lookup_canonical(lookup_result& res,                  //
                            typename Offsets::decoded_offset p,  //
                            const Kmer kmer,                     //
                            const Kmer kmer_rc,                  //
                            const uint64_t pos_in_kmer) const    //
    {
        if (p.absolute_offset < pos_in_kmer) return false;

        res.kmer_offset = p.absolute_offset - pos_in_kmer;

        auto read_kmer =
            util::read_kmer_at<Kmer>(strings, k, Kmer::bits_per_char * res.kmer_offset);
        if (read_kmer != kmer and read_kmer != kmer_rc) return false;

        res.kmer_orientation =
            read_kmer == kmer_rc ? constants::backward_orientation : constants::forward_orientation;

        if (res.kmer_offset >= res.string_begin and res.kmer_offset < res.string_end - k + 1) {
            res.kmer_id = res.kmer_offset - res.string_id * (k - 1);     // absolute kmer id
            res.kmer_id_in_string = res.kmer_offset - res.string_begin;  // relative kmer id
        } else {
            strings_offsets.offset_to_id(res, p, k);
        }

        if (res.kmer_offset < res.string_end - k + 1) return true;
        return false;
    }
};

}  // namespace sshash