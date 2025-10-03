#pragma once

namespace sshash {

struct endpoints  //
{
    struct string_group_info {
        uint32_t strings_length;  // all strings in the group have the same length
        uint32_t first_id;        // absolute id of the first string in the group
        uint64_t group_offset;    // absolute offset where this group of strings begins
    };

    struct builder  //
    {
        uint64_t build_from(std::vector<uint64_t> const& se) {
            reserve(128);  // reserve space for super groups
            assert(se.front() == 0);
            uint64_t prev_len = se[1] - se[0];
            uint64_t log2_prev_len = bits::util::ceil_log2_uint64(prev_len);

            // std::cout << "log2_prev_len = " << log2_prev_len << std::endl;

            uint64_t first_id = 0;
            uint64_t group_offset = 0;
            new_super_group(1ULL << 15);

            uint64_t prev_group_offset = 0;
            uint64_t num_bits_per_offset = 0;

            for (uint64_t i = 1; i != se.size(); ++i) {
                uint64_t curr_len = se[i] - se[i - 1];
                uint64_t log2_curr_len = bits::util::ceil_log2_uint64(curr_len);

                // std::cout << "log2_curr_len = " << log2_curr_len << std::endl;

                if (curr_len >= (1ULL << 32)) {
                    std::cerr
                        << "A string is longer than 2^32-1 chars...a u32 is not enough to hold a "
                           "string length"
                        << std::endl;
                }
                assert(curr_len >= prev_len);
                if (curr_len > prev_len)  //
                {
                    if (first_id >= (1ULL << 32)) {
                        std::cerr << "More than 2^32-1 strings...a u32 is not enough to hold a "
                                     "string id"
                                  << std::endl;
                    }
                    // std::cout << "len = " << prev_len << ", first_id = " << first_id
                    //           << ", offset = " << group_offset << std::endl;

                    push_back({static_cast<uint32_t>(prev_len),  //
                               static_cast<uint32_t>(first_id),  //
                               group_offset});                   //

                    if (group_offset > 0) {
                        // std::cout << "cumulative_len = " << (group_offset - prev_group_offset)
                        //           << std::endl;
                        uint64_t log2_group_offset =
                            bits::util::ceil_log2_uint64(group_offset - prev_group_offset) +
                            (m_data.back().size() == 1
                                 ? 1
                                 : bits::util::ceil_log2_uint64(m_data.back().size()));
                        if (log2_group_offset > num_bits_per_offset) {
                            num_bits_per_offset = log2_group_offset;
                        }
                        // std::cout << "log2_group_offset = " << log2_group_offset << std::endl;
                        prev_group_offset = group_offset;
                    }

                    if (log2_curr_len > log2_prev_len) {
                        new_super_group(1ULL << 15);
                        log2_prev_len = log2_curr_len;
                    }

                    prev_len = curr_len;
                    first_id = i - 1;
                    group_offset = se[i - 1];
                }
            }
            // std::cout << "len = " << prev_len << ", first_id = " << first_id
            //           << ", offset = " << group_offset << std::endl;
            push_back({static_cast<uint32_t>(prev_len),  //
                       static_cast<uint32_t>(first_id),  //
                       group_offset});                   //

            if (group_offset > 0) {
                // std::cout << "cumulative_len = " << (se.back() - prev_group_offset) << std::endl;
                uint64_t log2_group_offset =
                    bits::util::ceil_log2_uint64(se.back() - prev_group_offset) +
                    (m_data.back().size() == 1
                         ? 1
                         : bits::util::ceil_log2_uint64(m_data.back().size()));
                if (log2_group_offset > num_bits_per_offset) {
                    num_bits_per_offset = log2_group_offset;
                }
                // std::cout << "log2_group_offset = " << log2_group_offset << std::endl;
            }

            num_bits_per_offset +=
                m_data.size() == 1 ? 1 : bits::util::ceil_log2_uint64(m_data.size());
            return num_bits_per_offset;
        }

        /* Reserve space for n super groups. */
        void reserve(uint64_t n = 0) { m_data.reserve(n); }

        /* Reserve space for n group in the new group. */
        void new_super_group(uint64_t n = 0) {
            std::vector<string_group_info> v;
            v.reserve(n);
            m_data.emplace_back(v);
        }

        void push_back(string_group_info sgi) { m_data.back().push_back(sgi); }

        void build(endpoints& e) { e.m_data.swap(m_data); }

    private:
        std::vector<std::vector<string_group_info>> m_data;
    };

    /* (absolute offset, string_group_info) */
    std::pair<uint64_t, string_group_info> decode(const uint64_t offset) const  //
    {
        const uint64_t p = num_bits_per_super_group();
        // std::cout << "num_bits_per_super_group = " << p << std::endl;

        uint64_t super_group = offset & ((1ULL << p) - 1);
        // std::cout << "super_group = " << super_group << std::endl;

        assert(super_group < m_data.size());
        auto const& v = m_data[super_group];
        const uint64_t q = num_bits_per_group(v.size());
        // std::cout << "num_bits_per_group = " << q << std::endl;

        uint64_t group = (offset >> p) & ((1ULL << q) - 1);
        // std::cout << "group = " << group << std::endl;

        assert(group < v.size());
        return {(offset >> (p + q)) + v[group].group_offset, v[group]};
    }

    uint64_t num_bits_per_super_group() const {
        uint64_t b = m_data.size() == 1 ? 1 : bits::util::ceil_log2_uint64(m_data.size());
        assert(b > 0);
        return b;
    }
    uint64_t num_bits_per_group(uint64_t group_size) const {
        uint64_t b = group_size == 1 ? 1 : bits::util::ceil_log2_uint64(group_size);
        assert(b > 0);
        return b;
    }

    uint64_t id_to_offset(const uint64_t kmer_id, const uint64_t k) const {
        return 0;
        // TODO

        // constexpr uint64_t linear_scan_threshold = 32;
        // uint64_t lo = 0;
        // uint64_t hi = strings_endpoints.size() - 1;
        // assert(strings_endpoints.access(0) == 0);
        // while (hi - lo > linear_scan_threshold) {
        //     uint64_t mid = lo + (hi - lo) / 2;
        //     uint64_t val = strings_endpoints.access(mid);
        //     assert(val >= strings_groups_info[mid].string_id * (k - 1));
        //     uint64_t id = val - strings_groups_info[mid].string_id * (k - 1);
        //     if (kmer_id <= id) {
        //         hi = mid;
        //     } else {
        //         lo = mid + 1;
        //     }
        // }
        // assert(lo < hi);
        // assert(hi < strings_endpoints.size());
        // for (auto it = strings_endpoints.get_iterator_at(lo); lo < hi; ++lo, it.next()) {
        //     uint64_t val = it.value() - strings_groups_info[lo].string_id * (k - 1);
        //     if (val > kmer_id) break;
        // }
        // assert(lo > 0);
        // auto const& p = strings_groups_info[lo - 1];
        // uint64_t id = strings_endpoints.access(lo - 1) - p.string_id * (k - 1);
        // assert(id <= kmer_id);
        // uint64_t offset =
        //     kmer_id + ((kmer_id - id) / (p.string_length - k + 1) + p.string_id) * (k - 1);
        // return offset;
    }

    std::vector<std::vector<string_group_info>> const& data() const { return m_data; }

    uint64_t num_bytes() const {
        uint64_t n = sizeof(uint64_t);  // for std::vector::size
        for (auto const& v : m_data) n += sizeof(uint64_t) + essentials::vec_bytes(v);
        return n;
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visitor.visit(m_data);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_data);
    }

private:
    std::vector<std::vector<string_group_info>> m_data;
};

/*

    If we sort the strings by increasing lengths,
    we can avoid storing an offset *per string* and the
    resulting locate query.
    Rather, we can keep an offset *per group of lists*
    having the same length. Call A such sequence.
    We can then use a locate query on A to determine the
    group of lists to which an offset belongs.

    The common API is, given an offset:
    return (contig_id, contig_begin, contig_end).

*/

struct locate_result {
    uint64_t string_id, string_begin, string_end;
};

/*
    This solution stores an array of triples, T, one triple
    per group of strings having the same length:
    (string length,
     absolute string id of the first string in the group,
     offset of group).

    We use 32-bit ints for the first two components, and
    a 64-bit int for the third component.

    Additionally, we store a sequence of locate query "hints", H,
    of size g, where g := next_power_of_two(U/n).
    Then for a query, we do i = x >> log(g) and use this as an index
    into H: H[i] indicates the position in T of the first element of T
    that has its high bits equal to x >> log(g), and start scanning T
    from there.
*/
struct endpoints1 {};

/*
    This solution stores an array of pairs, P, one pair
    per group of strings having the same length:
    (string length,
     absolute string id of the first string in the group).

    We use 32-bit ints for both components.

    The offsets of the groups are encoded using `bits::endpoints_sequence`.
*/
struct endpoints2 {
    endpoints2() {}

    void build(std::vector<uint64_t> const& se) {
        std::vector<uint64_t> v;
        v.reserve(se.size());
        strings_groups_info.reserve(se.size());

        assert(se.front() == 0);
        uint64_t prev_len = se[1] - se[0];
        uint64_t string_id = 0;
        uint64_t string_group_begin = 0;
        for (uint64_t i = 1; i != se.size(); ++i) {
            uint64_t curr_len = se[i] - se[i - 1];
            if (curr_len >= 1ULL << 32) {
                std::cerr << "A string is longer than 2^32-1 chars...a u32 is not enough to hold a "
                             "string length"
                          << std::endl;
            }
            assert(curr_len >= prev_len);
            if (curr_len > prev_len) {
                if (string_id >= 1ULL << 32) {
                    std::cerr
                        << "More than 2^32-1 strings...a u32 is not enough to hold a string_id"
                        << std::endl;
                }
                // std::cout << "len = " << prev_len << ", string_id = " << string_id
                //           << ", offset = " << string_group_begin << std::endl;
                strings_groups_info.push_back({static_cast<uint32_t>(prev_len),  //
                                               static_cast<uint32_t>(string_id)});
                v.push_back(string_group_begin);
                prev_len = curr_len;
                string_id = i - 1;
                string_group_begin = se[i - 1];
            }
        }
        // std::cout << "len = " << prev_len << ", string_id = " << string_id
        //           << ", offset = " << string_group_begin << std::endl;
        strings_groups_info.push_back({static_cast<uint32_t>(prev_len),  //
                                       static_cast<uint32_t>(string_id)});
        v.push_back(string_group_begin);
        v.push_back(string_group_begin + prev_len);
        strings_endpoints.encode(v.begin(), v.size(), v.back());
    }

    locate_result locate(const uint64_t kmer_offset) const {
        auto p = strings_endpoints.locate(kmer_offset);
        uint64_t i = p.first.pos;
        uint64_t string_group_begin = p.first.val;

        assert(i < strings_groups_info.size());
        auto q = strings_groups_info[i];
        uint64_t string_id_relative_to_group = (kmer_offset - string_group_begin) / q.string_length;
        uint64_t string_id = q.string_id + string_id_relative_to_group;
        uint64_t string_begin = string_group_begin + string_id_relative_to_group * q.string_length;
        uint64_t string_end = string_begin + q.string_length;

        return {string_id, string_begin, string_end};
    }

    uint64_t id_to_offset(const uint64_t kmer_id, const uint64_t k) const {
        constexpr uint64_t linear_scan_threshold = 32;
        uint64_t lo = 0;
        uint64_t hi = strings_endpoints.size() - 1;
        assert(strings_endpoints.access(0) == 0);
        while (hi - lo > linear_scan_threshold) {
            uint64_t mid = lo + (hi - lo) / 2;
            uint64_t val = strings_endpoints.access(mid);
            assert(val >= strings_groups_info[mid].string_id * (k - 1));
            uint64_t id = val - strings_groups_info[mid].string_id * (k - 1);
            if (kmer_id <= id) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        assert(lo < hi);
        assert(hi < strings_endpoints.size());
        for (auto it = strings_endpoints.get_iterator_at(lo); lo < hi; ++lo, it.next()) {
            uint64_t val = it.value() - strings_groups_info[lo].string_id * (k - 1);
            if (val > kmer_id) break;
        }
        assert(lo > 0);
        auto const& p = strings_groups_info[lo - 1];
        uint64_t id = strings_endpoints.access(lo - 1) - p.string_id * (k - 1);
        assert(id <= kmer_id);
        uint64_t offset =
            kmer_id + ((kmer_id - id) / (p.string_length - k + 1) + p.string_id) * (k - 1);
        return offset;
    }

    uint64_t num_bytes() const {
        return strings_endpoints.num_bytes() + essentials::vec_bytes(strings_groups_info);
    }

    struct string_group_info {
        uint32_t string_length;  // all strings in the group have the same length
        uint32_t string_id;      // first id in the group
    };

    bits::endpoints_sequence<> strings_endpoints;
    std::vector<string_group_info> strings_groups_info;

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.strings_endpoints);
        visitor.visit(t.strings_groups_info);
    }
};

}  // namespace sshash