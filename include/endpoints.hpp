#pragma once

namespace sshash {

struct endpoints  //
{
    struct string_group_info {
        uint32_t strings_length;  // all strings in the group have the same length
        uint32_t first_id;        // absolute id of the first string in the group
        uint64_t group_offset;    // absolute offset where this group of strings begins
    };

    struct string_super_group_info {
        uint8_t num_bits_per_group;
        std::vector<string_group_info> string_groups_info;

        template <typename Visitor>
        void visit(Visitor& visitor) const {
            visitor.visit(num_bits_per_group);
            visitor.visit(string_groups_info);
        }

        template <typename Visitor>
        void visit(Visitor& visitor) {
            visitor.visit(num_bits_per_group);
            visitor.visit(string_groups_info);
        }
    };

    struct builder  //
    {
        builder() : m_num_bits_per_super_group(0), m_num_bits_per_offset(0) {}

        void build_from(std::vector<uint64_t> const& se)  //
        {
            m_string_super_groups_info.reserve(128);  // reserve space for super groups

            assert(se.front() == 0);
            uint64_t prev_len = se[1] - se[0];
            uint64_t log2_prev_len = bits::util::ceil_log2_uint64(prev_len);

            uint64_t first_id = 0;
            uint64_t group_offset = 0;

            /* Reserve space for n groups in the new super group. */
            auto new_super_group = [&](uint64_t n = 0) {
                string_super_group_info ssgi;
                ssgi.num_bits_per_group = 0;
                ssgi.string_groups_info.reserve(n);
                m_string_super_groups_info.emplace_back(ssgi);
            };

            new_super_group(1ULL << 15);
            uint64_t prev_group_offset = 0;
            m_num_bits_per_offset = 0;

            for (uint64_t i = 1; i != se.size(); ++i) {
                uint64_t curr_len = se[i] - se[i - 1];
                uint64_t log2_curr_len = bits::util::ceil_log2_uint64(curr_len);
                assert(curr_len >= prev_len);

                if (curr_len > prev_len)  //
                {
                    // std::cout << "len = " << prev_len << ", first_id = " << first_id
                    //           << ", offset = " << group_offset << std::endl;

                    m_string_super_groups_info.back().string_groups_info.push_back(
                        {static_cast<uint32_t>(prev_len),  //
                         static_cast<uint32_t>(first_id),  //
                         group_offset});                   //

                    if (group_offset > 0) {
                        // std::cout << "cumulative_len = " << (group_offset - prev_group_offset)
                        //           << std::endl;
                        uint64_t log2_group_offset =
                            bits::util::ceil_log2_uint64(group_offset - prev_group_offset) +
                            (m_string_super_groups_info.back().string_groups_info.size() == 1
                                 ? 1
                                 : bits::util::ceil_log2_uint64(m_string_super_groups_info.back()
                                                                    .string_groups_info.size()));
                        if (log2_group_offset > m_num_bits_per_offset) {
                            m_num_bits_per_offset = log2_group_offset;
                        }
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

            m_string_super_groups_info.back().string_groups_info.push_back(
                {static_cast<uint32_t>(prev_len),  //
                 static_cast<uint32_t>(first_id),  //
                 group_offset});                   //

            if (group_offset > 0) {
                // std::cout << "cumulative_len = " << (se.back() - prev_group_offset) << std::endl;
                uint64_t log2_group_offset =
                    bits::util::ceil_log2_uint64(se.back() - prev_group_offset) +
                    (m_string_super_groups_info.back().string_groups_info.size() == 1
                         ? 1
                         : bits::util::ceil_log2_uint64(
                               m_string_super_groups_info.back().string_groups_info.size()));
                if (log2_group_offset > m_num_bits_per_offset) {
                    m_num_bits_per_offset = log2_group_offset;
                }
                // std::cout << "log2_group_offset = " << log2_group_offset << std::endl;
            }

            m_num_bits_per_super_group =
                m_string_super_groups_info.size() == 1
                    ? 1
                    : bits::util::ceil_log2_uint64(m_string_super_groups_info.size());
            assert(m_num_bits_per_super_group > 0);
            std::cout << "num_bits_per_super_group = " << m_num_bits_per_super_group << std::endl;

            m_num_bits_per_offset += m_num_bits_per_super_group;

            /* set num. bits per group for all super groups */
            for (auto& g : m_string_super_groups_info) {
                g.num_bits_per_group =
                    g.string_groups_info.size() == 1
                        ? 1
                        : bits::util::ceil_log2_uint64(g.string_groups_info.size());
            }

            {
                // print
                for (uint64_t i = 0; i != m_string_super_groups_info.size(); ++i) {
                    std::cout << "super_group " << i << ":\n\t";
                    auto const& g = m_string_super_groups_info[i];
                    for (auto sg_info : g.string_groups_info) {
                        std::cout << "(" << sg_info.strings_length << "," << sg_info.first_id << ","
                                  << sg_info.group_offset << ")" << ' ';
                    }
                    std::cout << std::endl;
                    std::cout << "num_bits_per_group = " << int(g.num_bits_per_group) << std::endl;
                }
            }
        }

        uint64_t num_bits_per_offset() const { return m_num_bits_per_offset; }

        uint64_t encode(uint64_t offset, uint64_t super_group_id, uint64_t group_id,
                        uint64_t group_offset)  //
        {
            assert(offset >= group_offset);
            offset -= group_offset;  // relative to beginning of group
            auto const& g = m_string_super_groups_info[super_group_id];
            assert(group_id < g.string_groups_info.size());
            offset <<= g.num_bits_per_group;
            offset += group_id;
            assert(super_group_id < (1ULL << m_num_bits_per_super_group));
            offset <<= m_num_bits_per_super_group;
            offset += super_group_id;
            return offset;
        }

        void build(endpoints& e) {
            std::swap(e.m_num_bits_per_super_group, m_num_bits_per_super_group);
            e.m_string_super_groups_info.swap(m_string_super_groups_info);
        }

    private:
        uint8_t m_num_bits_per_super_group;
        std::vector<string_super_group_info> m_string_super_groups_info;

        uint64_t m_num_bits_per_offset;
    };

    struct decoded_offset {
        uint64_t offset;  // absolute offset
        string_group_info group_info;
    };
    decoded_offset decode(const uint64_t offset) const  //
    {
        const uint64_t p = m_num_bits_per_super_group;
        uint64_t super_group_id = offset & ((1ULL << p) - 1);
        assert(super_group_id < m_string_super_groups_info.size());
        auto const& g = m_string_super_groups_info[super_group_id];
        const uint64_t q = g.num_bits_per_group;
        uint64_t group_id = (offset >> p) & ((1ULL << q) - 1);
        assert(group_id < g.string_groups_info.size());
        auto sgi = g.string_groups_info[group_id];
        return {(offset >> (p + q)) + sgi.group_offset, sgi};
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

    uint64_t num_bytes() const {
        uint64_t n =
            sizeof(m_num_bits_per_super_group) + sizeof(uint64_t);  // for std::vector::size
        for (auto const& g : m_string_super_groups_info) {
            n += sizeof(uint64_t) + sizeof(g.num_bits_per_group) +
                 essentials::vec_bytes(g.string_groups_info);
        }
        return n;
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

private:
    uint8_t m_num_bits_per_super_group;
    std::vector<string_super_group_info> m_string_super_groups_info;

    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_num_bits_per_super_group);
        visitor.visit(t.m_string_super_groups_info);
    }
};

}  // namespace sshash