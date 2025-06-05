#pragma once

#include <algorithm>

namespace sshash {

template <typename Iterator, typename Compare>
void parallel_merge(Iterator A_begin, Iterator A_end,                   //
                    Iterator B_begin, Iterator B_end,                   //
                    Iterator output_begin,                              //
                    Compare comp, uint64_t sequential_merge_threshold)  //
{
    using T = typename std::iterator_traits<Iterator>::value_type;
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
    T larger_mid_val = *larger_mid;
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

template <typename T, typename Compare>
void parallel_sort(std::vector<T>& data, const uint64_t num_threads, Compare comp)  //
{
    if (num_threads <= 1) {
        std::sort(data.begin(), data.end(), comp);
        return;
    }

    assert((num_threads & (num_threads - 1)) == 0);
    const uint64_t data_size = data.size();
    const uint64_t chunk_size = (data_size + num_threads - 1) / num_threads;
    const uint64_t sequential_merge_threshold = data_size / uint64_t(std::log2(num_threads));

    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    std::vector<std::pair<typename std::vector<T>::iterator,  //
                          typename std::vector<T>::iterator>>
        ranges(num_threads);

    for (uint64_t i = 0; i != num_threads; ++i) {
        size_t begin = i * chunk_size;
        size_t end = (i == num_threads - 1) ? data_size : begin + chunk_size;
        // std::cout << "[" << begin << "," << end << ")" << std::endl;
        ranges[i] = {data.begin() + begin, data.begin() + end};
        threads.emplace_back(
            [&, begin, end]() { std::sort(data.begin() + begin, data.begin() + end, comp); });
    }
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }

    std::vector<T> temp_data;
    temp_data.resize(data.size());

    bool use_data = true;
    std::vector<std::pair<typename std::vector<T>::iterator,  //
                          typename std::vector<T>::iterator>>
        next_ranges;

    while (ranges.size() != 1)  //
    {
        next_ranges.clear();
        for (uint64_t i = 0; i != ranges.size(); i += 2) {
            if (i + 1 < ranges.size()) {
                auto begin1 = ranges[i].first;
                auto end1 = ranges[i].second;
                auto begin2 = ranges[i + 1].first;
                auto end2 = ranges[i + 1].second;
                uint64_t output_size = (end1 - begin1) + (end2 - begin2);

                uint64_t offset = 0;
                typename std::vector<T>::iterator output_iterator;

                if (use_data) {
                    offset = std::distance(data.begin(), begin1);
                    output_iterator = temp_data.begin() + offset;
                } else {
                    offset = std::distance(temp_data.begin(), begin1);
                    output_iterator = data.begin() + offset;
                }
                assert(offset <= data_size);

                parallel_merge(begin1, end1, begin2, end2, output_iterator, comp,
                               sequential_merge_threshold);
                assert(std::is_sorted(output_iterator, output_iterator + output_size, comp));

                auto merged_begin = output_iterator;
                auto merged_end = merged_begin + output_size;
                next_ranges.push_back({merged_begin, merged_end});
            }
        }
        ranges.swap(next_ranges);
        use_data = !use_data;
    }

    if (!use_data) data.swap(temp_data);
}

}  // namespace sshash