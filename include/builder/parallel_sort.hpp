#pragma once

#include <algorithm>

namespace sshash {

template <typename Iterator, typename Compare>
void parallel_merge(Iterator A_begin, Iterator A_end,                   //
                    Iterator B_begin, Iterator B_end,                   //
                    Iterator output_begin,                              //
                    Compare comp, uint64_t sequential_merge_threshold)  //
{
    assert(std::is_sorted(A_begin, A_end, comp));
    assert(std::is_sorted(B_begin, B_end, comp));

    const uint64_t size_A = std::distance(A_begin, A_end);
    const uint64_t size_B = std::distance(B_begin, B_end);
    if (size_A + size_B <= sequential_merge_threshold) {
        std::merge(A_begin, A_end, B_begin, B_end, output_begin, comp);
        return;
    }

    Iterator larger_begin = A_begin;
    Iterator larger_end = A_end;
    Iterator smaller_begin = B_begin;
    Iterator smaller_end = B_end;
    if (size_A < size_B) {
        std::swap(larger_begin, smaller_begin);
        std::swap(larger_end, smaller_end);
    }

    Iterator larger_mid = larger_begin + (std::distance(larger_begin, larger_end) / 2);
    auto larger_mid_val = *larger_mid;
    Iterator smaller_mid = std::lower_bound(smaller_begin, smaller_end, larger_mid_val, comp);

    std::thread thread1([&]() {
        parallel_merge(larger_begin, larger_mid,    //
                       smaller_begin, smaller_mid,  //
                       output_begin,                //
                       comp, sequential_merge_threshold);
    });
    std::thread thread2([&]() {
        parallel_merge(larger_mid, larger_end,    //
                       smaller_mid, smaller_end,  //
                       output_begin + std::distance(larger_begin, larger_mid) +
                           std::distance(smaller_begin, smaller_mid),  //
                       comp, sequential_merge_threshold);
    });
    thread1.join();
    thread2.join();
}

/*
    Data to sort is in *data*. The vector *temp_data* is the temporary
    working memory, which must be of same size as data.
*/
template <typename T, typename Compare>
void parallel_sort(std::vector<T>& data, const uint64_t num_threads, Compare comp)  //
{
    std::vector<T> temp_data;
    temp_data.resize(data.size());

    if (num_threads <= 1) {
        std::sort(data.begin(), data.end(), comp);
        return;
    }

    const uint64_t data_size = data.size();
    const uint64_t chunk_size = (data_size + num_threads - 1) / num_threads;
    const uint64_t sequential_merge_threshold = data_size / uint64_t(std::log2(num_threads));
    assert(sequential_merge_threshold > 0);

    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    using iterator_t = typename std::vector<T>::iterator;
    std::vector<std::pair<iterator_t, iterator_t>> ranges(num_threads);

    for (uint64_t i = 0; i != num_threads; ++i) {
        uint64_t begin = i * chunk_size;
        uint64_t end = (i == num_threads - 1) ? data_size : begin + chunk_size;
        ranges[i] = {data.begin() + begin, data.begin() + end};
        threads.emplace_back(
            [&, begin, end]() { std::sort(data.begin() + begin, data.begin() + end, comp); });
    }
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }

    bool swap = false;
    std::vector<std::pair<iterator_t, iterator_t>> next_ranges;

    while (ranges.size() != 1)  //
    {
        next_ranges.clear();
        for (uint64_t i = 0; i != ranges.size(); i += 2) {
            if (i + 1 < ranges.size()) {
                auto [begin1, end1] = ranges[i];
                auto [begin2, end2] = ranges[i + 1];
                uint64_t output_size = (end1 - begin1) + (end2 - begin2);

                auto input = data.begin();
                auto output = temp_data.begin();
                if (swap) std::swap(input, output);
                uint64_t offset = std::distance(input, begin1);
                auto output_iterator = output + offset;
                assert(offset <= data_size);

                parallel_merge(begin1, end1, begin2, end2, output_iterator, comp,
                               sequential_merge_threshold);
                assert(std::is_sorted(output_iterator, output_iterator + output_size, comp));

                auto merged_begin = output_iterator;
                auto merged_end = merged_begin + output_size;
                next_ranges.push_back({merged_begin, merged_end});
            } else {
                // Odd range out: carry it forward to the next iteration
                next_ranges.push_back(ranges[i]);
            }
        }
        ranges.swap(next_ranges);
        swap = !swap;
    }

    if (swap) data.swap(temp_data);
}

}  // namespace sshash