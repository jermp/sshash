#pragma once

#include <algorithm>
#include <functional>

namespace sshash {

template <typename T, typename Compare>
void parallel_sort(std::vector<T>& data, const uint64_t num_threads, Compare comp) {
    if (num_threads <= 1) {
        std::sort(data.begin(), data.end(), comp);
        return;
    }

    const uint64_t data_size = data.size();
    const uint64_t chunk_size = (data_size + num_threads - 1) / num_threads;
    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    std::vector<std::vector<T>> sorted_chunks(num_threads);

    for (uint64_t i = 0; i != num_threads; ++i) {
        uint64_t begin = i * chunk_size;
        uint64_t end = (i == num_threads - 1) ? data_size : begin + chunk_size;
        sorted_chunks[i].resize(end - begin);
        std::copy(data.begin() + begin, data.begin() + end, sorted_chunks[i].begin());
        threads.emplace_back([&sorted_chunks, i, comp]() {
            std::sort(sorted_chunks[i].begin(), sorted_chunks[i].end(), comp);
        });
    }
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }

    using heap_element = std::pair<T, uint32_t>;  // (value, index in heap)
    std::vector<heap_element> heap;
    heap.reserve(num_threads);
    std::vector<typename std::vector<T>::iterator> iterators(num_threads);

    for (uint64_t i = 0; i != num_threads; ++i) {
        if (!sorted_chunks[i].empty()) {
            heap.emplace_back(sorted_chunks[i][0], i);
            iterators[i] = sorted_chunks[i].begin();
        }
    }

    auto neg_comp = [&](heap_element const& x, heap_element const& y) {
        return !comp(x.first, y.first);
    };

    std::make_heap(heap.begin(), heap.end(), neg_comp);
    data.clear();
    data.reserve(data_size);

    while (!heap.empty()) {
        auto [min, i] = heap.front();
        data.push_back(min);
        ++iterators[i];
        if (iterators[i] != sorted_chunks[i].end()) {  // percolate down the head
            heap[0].first = *iterators[i];
            uint64_t pos = 0;
            uint64_t size = heap.size();
            while (2 * pos + 1 < size) {
                uint64_t i = 2 * pos + 1;
                if (i + 1 < size and neg_comp(heap[i], heap[i + 1])) ++i;
                if (neg_comp(heap[i], heap[pos])) break;
                std::swap(heap[i], heap[pos]);
                pos = i;
            }
        } else {
            std::pop_heap(heap.begin(), heap.end(), neg_comp);
            heap.pop_back();
        }
    }
}

}  // namespace sshash