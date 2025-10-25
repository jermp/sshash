#pragma once

#include <functional>

namespace sshash {

template <typename T>
struct file_merging_iterator  //
{
    const uint64_t scan_threshold = 32;

    template <typename FileNamesIterator>
    file_merging_iterator(FileNamesIterator file_names_iterator, uint64_t num_files_to_merge)
        : m_mm_files(num_files_to_merge)  //
    {
        if (num_files_to_merge == 0) return;

        /* create the input iterators and make the heap */
        m_iterators.reserve(num_files_to_merge);
        for (uint64_t i = 0; i != num_files_to_merge; ++i, ++file_names_iterator) {
            m_mm_files[i].open(*file_names_iterator, mm::advice::sequential);
            m_iterators.push_back(
                {m_mm_files[i].data(), m_mm_files[i].data() + m_mm_files[i].size()});
        }

        m_min_idx = 0;
        if (m_iterators.size() <= scan_threshold) {
            compute_min();
        } else {
            std::make_heap(m_iterators.begin(), m_iterators.end(), heap_comparator);
        }
    }

    bool has_next() { return m_iterators.size() != 0; }
    void next() { advance_heap_head(); }
    T operator*() const { return *(m_iterators[m_min_idx].begin); }

    void close() {
        for (uint64_t i = 0; i != m_mm_files.size(); ++i) m_mm_files[i].close();
        m_iterators.clear();
        m_mm_files.clear();
    }

private:
    struct pointer_pair {
        T const* begin;
        T const* end;
    };
    std::vector<pointer_pair> m_iterators;
    std::vector<mm::file_source<T const>> m_mm_files;
    uint64_t m_min_idx;

    std::function<bool(pointer_pair, pointer_pair)> heap_comparator =
        [](pointer_pair x, pointer_pair y) { return *(x.begin) > *(y.begin); };

    void advance_heap_head() {
        auto& it = m_iterators[m_min_idx];
        it.begin += 1;
        if (m_iterators.size() <= scan_threshold) {  // compute min with a linear scan
            if (it.begin == it.end) {
                m_iterators.erase(m_iterators.begin() + m_min_idx);
                m_min_idx = 0;
                if (m_iterators.size() == 0) return;
            }
            compute_min();
        } else {  // update the min-heap
            if (it.begin != it.end) {
                uint64_t pos = 0;
                uint64_t size = m_iterators.size();
                while (2 * pos + 1 < size) {
                    uint64_t i = 2 * pos + 1;
                    if (i + 1 < size and heap_comparator(m_iterators[i], m_iterators[i + 1])) ++i;
                    if (heap_comparator(m_iterators[i], m_iterators[pos])) break;
                    std::swap(m_iterators[pos], m_iterators[i]);
                    pos = i;
                }
            } else {
                std::pop_heap(m_iterators.begin(), m_iterators.end(), heap_comparator);
                m_iterators.pop_back();
            }
            assert(m_min_idx == 0);
        }
    };

    void compute_min() {
        m_min_idx = 0;
        auto min_val = *m_iterators.front().begin;
        for (uint64_t i = 1; i != m_iterators.size(); ++i) {
            assert(m_iterators[i].begin != m_iterators[i].end);
            auto val = *m_iterators[i].begin;
            if (val < min_val) {
                min_val = val;
                m_min_idx = i;
            }
        }
    }
};

}  // namespace sshash