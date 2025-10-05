#pragma once

namespace sshash {

struct endpoints  //
{
    struct string_group_info {
        uint32_t strings_length;  // all strings in the group have the same length
        uint32_t first_id;        // absolute id of the first string in the group
        uint64_t group_offset;    // where the group of strings begins
    };

    struct string_super_group_info {
        uint8_t num_bits_per_group;
        uint32_t super_group_offset;  // where the super group info begins
    };

    struct builder  //
    {
        builder() : m_num_bits_per_super_group(0), m_num_bits_per_offset(0) {}

        void build_from(std::vector<uint64_t> const& se)  //
        {
            m_string_super_groups_info.reserve(1ULL << 7);
            m_string_groups_info.reserve(1ULL << 16);

            assert(se.front() == 0);
            uint64_t prev_len = se[1] - se[0];
            uint64_t log2_prev_len = bits::util::ceil_log2_uint64(prev_len);

            uint64_t first_id = 0;
            uint64_t group_offset = 0;

            auto new_super_group = [&]() {
                m_string_super_groups_info.push_back(
                    {0, static_cast<uint32_t>(m_string_groups_info.size())});
            };

            new_super_group();

            for (uint64_t i = 1; i != se.size(); ++i) {
                uint64_t curr_len = se[i] - se[i - 1];
                uint64_t log2_curr_len = bits::util::ceil_log2_uint64(curr_len);
                assert(curr_len >= prev_len);
                if (curr_len > prev_len)  //
                {
                    m_string_groups_info.push_back({static_cast<uint32_t>(prev_len),  //
                                                    static_cast<uint32_t>(first_id),  //
                                                    group_offset});                   //

                    if (log2_curr_len > log2_prev_len) {
                        new_super_group();
                        log2_prev_len = log2_curr_len;
                    }

                    prev_len = curr_len;
                    first_id = i - 1;
                    group_offset = se[i - 1];
                }
            }

            m_string_groups_info.push_back({static_cast<uint32_t>(prev_len),  //
                                            static_cast<uint32_t>(first_id),  //
                                            group_offset});                   //

            m_num_bits_per_super_group =
                m_string_super_groups_info.size() == 1
                    ? 1
                    : bits::util::ceil_log2_uint64(m_string_super_groups_info.size());
            assert(m_num_bits_per_super_group > 0);
            std::cout << "num_bits_per_super_group = " << int(m_num_bits_per_super_group)
                      << std::endl;

            /*
                Push a last dummy super-group so that
                super_group_size = m_string_super_groups_info[i+1].super_group_offset -
                m_string_super_groups_info[i].super_group_offset;
                for any super_group_id i
            */
            new_super_group();

            /* set num. bits per group for all super groups */
            for (uint64_t i = 0, prev_offset = 0; i != m_string_super_groups_info.size() - 1; ++i) {
                std::cout << "super_group " << i << ": ";
                uint64_t begin = m_string_super_groups_info[i].super_group_offset;
                uint64_t end = m_string_super_groups_info[i + 1].super_group_offset;
                uint64_t super_group_size = end - begin;
                m_string_super_groups_info[i].num_bits_per_group =
                    super_group_size == 1 ? 1 : bits::util::ceil_log2_uint64(super_group_size);
                std::cout << "super_group_size = " << super_group_size << " "
                          << "num_bits_per_group = "
                          << int(m_string_super_groups_info[i].num_bits_per_group) << " ";

                uint64_t longest_cumulative_length = 0;
                for (uint64_t j = begin; j != end; ++j)  //
                {
                    uint64_t next_offset = j + 1 < m_string_groups_info.size()
                                               ? m_string_groups_info[j + 1].group_offset
                                               : se.back();
                    assert(next_offset > prev_offset);
                    uint64_t cumulative_length = next_offset - prev_offset;
                    if (cumulative_length > longest_cumulative_length) {
                        longest_cumulative_length = cumulative_length;
                    }
                    prev_offset = next_offset;
                }

                assert(longest_cumulative_length > 1);
                uint64_t num_bits_per_relative_offset =
                    bits::util::ceil_log2_uint64(longest_cumulative_length);
                uint64_t num_bits_per_offset =
                    m_string_super_groups_info[i].num_bits_per_group + num_bits_per_relative_offset;
                std::cout << "num_bits_per_relative_offset = " << num_bits_per_relative_offset
                          << " num_bits_per_offset = " << num_bits_per_offset << std::endl;

                if (num_bits_per_offset > m_num_bits_per_offset) {
                    m_num_bits_per_offset = num_bits_per_offset;
                }
            }

            m_num_bits_per_offset += m_num_bits_per_super_group;
            std::cout << "num_bits_per_offset = " << m_num_bits_per_offset << std::endl;
            assert(m_num_bits_per_offset > 0);
        }

        uint64_t num_bits_per_offset() const { return m_num_bits_per_offset; }

        uint64_t encode(uint64_t offset, uint64_t super_group_id, uint64_t group_id,
                        uint64_t group_offset)  //
        {
            assert(offset >= group_offset);
            offset -= group_offset;  // relative to beginning of group
            auto g = m_string_super_groups_info[super_group_id];
            assert(group_id < m_string_super_groups_info[super_group_id + 1].super_group_offset -
                                  g.super_group_offset);
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
            e.m_string_groups_info.swap(m_string_groups_info);
        }

    private:
        uint8_t m_num_bits_per_super_group;
        std::vector<string_super_group_info> m_string_super_groups_info;
        std::vector<string_group_info> m_string_groups_info;

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
        auto g = m_string_super_groups_info[super_group_id];
        const uint64_t q = g.num_bits_per_group;
        uint64_t group_id = (offset >> p) & ((1ULL << q) - 1);
        assert(group_id < m_string_super_groups_info[super_group_id + 1].super_group_offset -
                              g.super_group_offset);
        auto sgi = m_string_groups_info[g.super_group_offset + group_id];
        return {(offset >> (p + q)) + sgi.group_offset, sgi};
    }

    uint64_t id_to_offset(const uint64_t kmer_id, const uint64_t k) const  //
    {
        auto it =
            std::upper_bound(m_string_groups_info.begin(), m_string_groups_info.end(), kmer_id,
                             [&](const uint64_t x, string_group_info const& sgi) {
                                 assert(sgi.group_offset >= sgi.first_id * (k - 1));
                                 return x < sgi.group_offset - sgi.first_id * (k - 1);
                             }) -
            1;
        assert(it != m_string_groups_info.end());
        auto sgi = *it;
        uint64_t id = sgi.group_offset - sgi.first_id * (k - 1);
        assert(id <= kmer_id);
        uint64_t offset =
            kmer_id + ((kmer_id - id) / (sgi.strings_length - k + 1) + sgi.first_id) * (k - 1);
        return offset;
    }

    uint64_t num_bytes() const {
        return sizeof(m_num_bits_per_super_group) + 2 * sizeof(uint64_t)  // for std::vector::size
               + essentials::vec_bytes(m_string_super_groups_info) +
               essentials::vec_bytes(m_string_groups_info);
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
    std::vector<string_group_info> m_string_groups_info;

    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_num_bits_per_super_group);
        visitor.visit(t.m_string_super_groups_info);
        visitor.visit(t.m_string_groups_info);
    }
};

}  // namespace sshash