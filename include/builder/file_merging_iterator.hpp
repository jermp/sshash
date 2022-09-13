#pragma once

#include <functional>

namespace sshash {

template <typename FileIterator>
struct file_merging_iterator {
    template <typename FileNamesIterator>
    file_merging_iterator(FileNamesIterator file_names_iterator, uint64_t num_files_to_merge)
        : m_mm_files(num_files_to_merge) {
        if (num_files_to_merge == 0) return;
        assert(num_files_to_merge > 0);

        m_iterators.reserve(num_files_to_merge);
        m_idx_heap.reserve(num_files_to_merge);

        /* create the input iterators and make the heap */
        for (uint64_t i = 0; i != num_files_to_merge; ++i, ++file_names_iterator) {
            m_mm_files[i].open(*file_names_iterator, mm::advice::sequential);
            m_iterators.emplace_back(m_mm_files[i].data(),
                                     m_mm_files[i].data() + m_mm_files[i].size());
            m_idx_heap.push_back(i);
        }
        std::make_heap(m_idx_heap.begin(), m_idx_heap.end(), heap_idx_comparator);
    }

    bool has_next() { return !m_idx_heap.empty(); }
    void next() { advance_heap_head(); }
    FileIterator const& operator*() const { return m_iterators[m_idx_heap.front()]; }

    void close() {
        for (uint64_t i = 0; i != m_mm_files.size(); ++i) m_mm_files[i].close();
        m_iterators.clear();
        m_idx_heap.clear();
        m_mm_files.clear();
    }

private:
    std::vector<FileIterator> m_iterators;
    std::vector<uint32_t> m_idx_heap;
    std::vector<mm::file_source<uint8_t>> m_mm_files;

    std::function<bool(uint32_t, uint32_t)> heap_idx_comparator = [&](uint32_t i, uint32_t j) {
        return (*m_iterators[i]) > (*m_iterators[j]);
    };

    void advance_heap_head() {
        uint32_t idx = m_idx_heap.front();
        m_iterators[idx].next();
        if (m_iterators[idx].has_next()) {  // percolate down the head
            uint64_t pos = 0;
            uint64_t size = m_idx_heap.size();
            while (2 * pos + 1 < size) {
                uint64_t i = 2 * pos + 1;
                if (i + 1 < size and heap_idx_comparator(m_idx_heap[i], m_idx_heap[i + 1])) ++i;
                if (heap_idx_comparator(m_idx_heap[i], m_idx_heap[pos])) break;
                std::swap(m_idx_heap[pos], m_idx_heap[i]);
                pos = i;
            }
        } else {
            std::pop_heap(m_idx_heap.begin(), m_idx_heap.end(), heap_idx_comparator);
            m_idx_heap.pop_back();
        }
    };
};

}  // namespace sshash