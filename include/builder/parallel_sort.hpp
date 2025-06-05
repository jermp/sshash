#pragma once

#include <algorithm>
// #include <functional>

namespace sshash {

template <typename Iterator, typename Compare>
void parallel_merge(Iterator A_begin, Iterator A_end,                                  //
                    Iterator B_begin, Iterator B_end,                                  //
                    std::vector<typename Iterator::value_type>& output, Compare comp)  //
{
    using T = typename Iterator::value_type;
    assert(std::is_sorted(A_begin, A_end, comp));
    assert(std::is_sorted(B_begin, B_end, comp));

    const uint64_t size_A = A_end - A_begin;
    const uint64_t size_B = B_end - B_begin;

    assert(output.size() == size_A + size_B);

    constexpr uint64_t sequential_merge_threshold = 1024 * 1024;
    if (size_A + size_B <= sequential_merge_threshold) {  // sequential merge for small inputs
        std::merge(A_begin, A_end, B_begin, B_end, output.begin(), comp);
        return;
    }

    Iterator larger_begin = A_begin;
    Iterator larger_end = A_end;
    Iterator smaller_begin = B_begin;
    Iterator smaller_end = B_end;
    uint64_t larger_size = size_A;
    if (size_A < size_B) {
        larger_begin = B_begin;
        larger_end = B_end;
        smaller_begin = A_begin;
        smaller_end = A_end;
        larger_size = size_B;
    }

    uint64_t pos_A = larger_size / 2;
    T mid_val = *(larger_begin + pos_A);

    auto it = std::lower_bound(smaller_begin, smaller_end, mid_val, comp);
    size_t pos_B = std::distance(smaller_begin, it);

    std::vector<T> output1(pos_A + pos_B);
    std::vector<T> output2(size_A + size_B - (pos_A + pos_B));
    std::thread thread1([&]() {
        parallel_merge(larger_begin, larger_begin + pos_A,    //
                       smaller_begin, smaller_begin + pos_B,  //
                       output1, comp);
        std::copy(output1.begin(), output1.end(), output.begin());
    });
    std::thread thread2([&]() {
        parallel_merge(larger_begin + pos_A, larger_end,    //
                       smaller_begin + pos_B, smaller_end,  //
                       output2, comp);
        std::copy(output2.begin(), output2.end(), output.begin() + output1.size());
    });
    thread1.join();
    thread2.join();

    assert(std::is_sorted(output.begin(), output.end(), comp));
}

template <typename T, typename Compare>
void parallel_merge(std::vector<std::vector<T>>& sorted_ranges, std::vector<T>& output,
                    Compare comp)  //
{
    if (sorted_ranges.empty()) return;

    if (sorted_ranges.size() == 1) {
        output.swap(sorted_ranges[0]);
        return;
    }

    std::vector<std::vector<T>> current_ranges;
    current_ranges.swap(sorted_ranges);
    std::vector<std::vector<T>> next_ranges;

    while (current_ranges.size() > 1) {
        next_ranges.clear();
        for (size_t i = 0; i < current_ranges.size(); i += 2) {
            if (i + 1 < current_ranges.size()) {
                std::vector<T> merged_range(current_ranges[i].size() +
                                            current_ranges[i + 1].size());
                parallel_merge(current_ranges[i].begin(), current_ranges[i].end(),
                               current_ranges[i + 1].begin(), current_ranges[i + 1].end(),
                               merged_range, comp);
                next_ranges.push_back(std::move(merged_range));
            } else {
                next_ranges.push_back(std::move(current_ranges[i]));
            }
        }
        current_ranges.swap(next_ranges);
    }

    output.swap(current_ranges[0]);
}

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

    data.clear();
    parallel_merge(sorted_chunks, data, comp);

    // using heap_element = std::pair<T, uint32_t>;  // (value, index in heap)
    // std::vector<heap_element> heap;
    // heap.reserve(num_threads);
    // std::vector<typename std::vector<T>::iterator> iterators(num_threads);

    // for (uint64_t i = 0; i != num_threads; ++i) {
    //     if (!sorted_chunks[i].empty()) {
    //         heap.emplace_back(sorted_chunks[i][0], i);
    //         iterators[i] = sorted_chunks[i].begin();
    //     }
    // }

    // auto neg_comp = [&](heap_element const& x, heap_element const& y) {
    //     return !comp(x.first, y.first);
    // };

    // std::make_heap(heap.begin(), heap.end(), neg_comp);
    // data.clear();
    // data.reserve(data_size);

    // while (!heap.empty()) {
    //     auto [min, i] = heap.front();
    //     data.push_back(min);
    //     ++iterators[i];
    //     if (iterators[i] != sorted_chunks[i].end()) {  // percolate down the head
    //         heap[0].first = *iterators[i];
    //         uint64_t pos = 0;
    //         uint64_t size = heap.size();
    //         while (2 * pos + 1 < size) {
    //             uint64_t i = 2 * pos + 1;
    //             if (i + 1 < size and neg_comp(heap[i], heap[i + 1])) ++i;
    //             if (neg_comp(heap[i], heap[pos])) break;
    //             std::swap(heap[i], heap[pos]);
    //             pos = i;
    //         }
    //     } else {
    //         std::pop_heap(heap.begin(), heap.end(), neg_comp);
    //         heap.pop_back();
    //     }
    // }
}

}  // namespace sshash