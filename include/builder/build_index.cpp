namespace sshash {

#pragma pack(push, 4)
struct bucket_pair {
    bucket_pair(uint64_t id, uint32_t size) : id(id), size(size) {}
    uint64_t id;
    uint32_t size;
};
#pragma pack(pop)

struct bucket_pairs_iterator : std::forward_iterator_tag {
    bucket_pairs_iterator(bucket_pair const* begin, bucket_pair const* end)
        : m_begin(begin)
        , m_end(end)
        , m_i(0)
        , m_val(0)  // first returned val is always 0
    {}

    uint64_t operator*() const { return m_val; }
    void operator++() {
        ++m_i;
        if (m_i == (*m_begin).id) {
            m_val += (*m_begin).size;
            ++m_begin;
            assert(m_begin <= m_end);
        }
    }

private:
    bucket_pair const* m_begin;
    bucket_pair const* m_end;
    uint64_t m_i;
    uint64_t m_val;
};

struct bucket_pairs {
    static constexpr uint64_t ram_limit = 0.25 * essentials::GB;

    bucket_pairs(std::string tmp_dirname = constants::default_tmp_dirname)
        : m_buffer_size(0)
        , m_num_files_to_merge(0)
        , m_run_identifier(pthash::clock_type::now().time_since_epoch().count())
        , m_tmp_dirname(tmp_dirname) {
        m_buffer_size = ram_limit / sizeof(bucket_pair);
        std::cout << "m_buffer_size " << m_buffer_size << std::endl;
    }

    void emplace_back(uint64_t id, uint32_t size) {
        if (m_buffer.size() == m_buffer_size) sort_and_flush();
        m_buffer.emplace_back(id, size);
    }

    void sort_and_flush() {
        std::cout << "sorting buffer..." << std::endl;
        std::sort(m_buffer.begin(), m_buffer.end(),
                  [](bucket_pair const& x, bucket_pair const& y) { return x.id < y.id; });

        auto tmp_output_filename = get_tmp_output_filename(m_num_files_to_merge);
        std::cout << "saving to file '" << tmp_output_filename << "'..." << std::endl;
        std::ofstream out(tmp_output_filename.c_str(), std::ofstream::binary);
        if (!out.is_open()) throw std::runtime_error("cannot open file");
        out.write(reinterpret_cast<char const*>(m_buffer.data()),
                  m_buffer.size() * sizeof(bucket_pair));
        out.close();

        m_buffer.clear();
        ++m_num_files_to_merge;
    }

    void finalize() {
        if (!m_buffer.empty()) sort_and_flush();
    }

    std::string get_bucket_pairs_filename() const {
        std::stringstream filename;
        filename << m_tmp_dirname << "/sshash.tmp.run_" << m_run_identifier << ".bucket_pairs.bin";
        return filename.str();
    }

    void merge() {
        if (m_num_files_to_merge == 0) return;

        assert(m_num_files_to_merge > 0);
        std::cout << "files to merge = " << m_num_files_to_merge << std::endl;

        struct iterator_type {
            iterator_type(bucket_pair const* b, bucket_pair const* e) : begin(b), end(e) {}
            bucket_pair const* begin;
            bucket_pair const* end;
        };
        std::vector<iterator_type> iterators;
        std::vector<uint32_t> idx_heap;
        iterators.reserve(m_num_files_to_merge);
        idx_heap.reserve(m_num_files_to_merge);
        std::vector<mm::file_source<bucket_pair>> mm_files(m_num_files_to_merge);

        auto heap_idx_comparator = [&](uint32_t i, uint32_t j) {
            bucket_pair const* begin_i = iterators[i].begin;
            bucket_pair const* begin_j = iterators[j].begin;
            return (*begin_i).id > (*begin_j).id;
        };

        auto advance_heap_head = [&]() {
            uint32_t idx = idx_heap.front();
            iterators[idx].begin += 1;
            if (iterators[idx].begin != iterators[idx].end) {  // percolate down the head
                uint64_t pos = 0;
                uint64_t size = idx_heap.size();
                while (2 * pos + 1 < size) {
                    uint64_t i = 2 * pos + 1;
                    if (i + 1 < size and heap_idx_comparator(idx_heap[i], idx_heap[i + 1])) ++i;
                    if (heap_idx_comparator(idx_heap[i], idx_heap[pos])) break;
                    std::swap(idx_heap[pos], idx_heap[i]);
                    pos = i;
                }
            } else {
                std::pop_heap(idx_heap.begin(), idx_heap.end(), heap_idx_comparator);
                idx_heap.pop_back();
            }
        };

        /* create the input iterators and make the heap */
        for (uint64_t i = 0; i != m_num_files_to_merge; ++i) {
            auto tmp_output_filename = get_tmp_output_filename(i);
            mm_files[i].open(tmp_output_filename, mm::advice::sequential);
            iterators.emplace_back(mm_files[i].data(), mm_files[i].data() + mm_files[i].size());
            idx_heap.push_back(i);
        }
        std::make_heap(idx_heap.begin(), idx_heap.end(), heap_idx_comparator);

        std::ofstream out(get_bucket_pairs_filename().c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open file");

        uint64_t num_written_pairs = 0;
        while (!idx_heap.empty()) {
            bucket_pair const* begin = iterators[idx_heap.front()].begin;
            out.write(reinterpret_cast<char const*>(begin), sizeof(bucket_pair));
            num_written_pairs += 1;
            if (num_written_pairs % 50000000 == 0) {
                std::cout << "num_written_pairs = " << num_written_pairs << std::endl;
            }
            advance_heap_head();
        }
        std::cout << "num_written_pairs = " << num_written_pairs << std::endl;
        out.close();

        /* remove tmp files */
        for (uint64_t i = 0; i != m_num_files_to_merge; ++i) {
            mm_files[i].close();
            auto tmp_output_filename = get_tmp_output_filename(i);
            std::remove(tmp_output_filename.c_str());
        }

        std::vector<bucket_pair>().swap(m_buffer);
        m_num_files_to_merge = 0;  // any other call to merge() will do nothing
    }

    void remove_tmp_file() { std::remove(get_bucket_pairs_filename().c_str()); }

private:
    uint64_t m_buffer_size;
    uint64_t m_num_files_to_merge;
    uint64_t m_run_identifier;
    std::string m_tmp_dirname;
    std::vector<bucket_pair> m_buffer;

    std::string get_tmp_output_filename(uint64_t id) {
        std::stringstream filename;
        filename << m_tmp_dirname << "/sshash.tmp.run_" << m_run_identifier << ".bucket_pairs."
                 << id << ".bin";
        return filename.str();
    }
};

buckets_statistics build_index(parse_data& data, minimizers const& m_minimizers,
                               buckets& m_buckets) {
    uint64_t num_buckets = m_minimizers.size();
    uint64_t num_kmers = data.num_kmers;
    uint64_t num_super_kmers = data.strings.num_super_kmers();

    pthash::compact_vector::builder offsets;
    offsets.resize(num_super_kmers, std::ceil(std::log2(data.strings.num_bits() / 2)));

    std::cout << "bits_per_offset = ceil(log2(" << data.strings.num_bits() / 2
              << ")) = " << std::ceil(std::log2(data.strings.num_bits() / 2)) << std::endl;

    mm::file_source<minimizer_tuple> input(data.minimizers.get_minimizers_filename(),
                                           mm::advice::sequential);

    bucket_pairs bucket_pairs_manager;
    uint64_t num_singletons = 0;
    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size()); it.has_next();
         it.next()) {
        uint32_t list_size = it.list().size();
        assert(list_size > 0);
        if (list_size != 1) {
            uint64_t bucket_id = m_minimizers.lookup(it.minimizer());
            bucket_pairs_manager.emplace_back(bucket_id + 1, list_size - 1);
        } else {
            ++num_singletons;
        }
    }
    bucket_pairs_manager.finalize();

    std::cout << "num_singletons " << num_singletons << "/" << num_buckets << " ("
              << (num_singletons * 100.0) / num_buckets << "%)" << std::endl;

    bucket_pairs_manager.merge();

    {
        mm::file_source<bucket_pair> bucket_pairs_file(
            bucket_pairs_manager.get_bucket_pairs_filename(), mm::advice::sequential);
        bucket_pairs_iterator iterator(bucket_pairs_file.data(),
                                       bucket_pairs_file.data() + bucket_pairs_file.size());
        m_buckets.num_super_kmers_before_bucket.encode(iterator, num_buckets + 1,
                                                       num_super_kmers - num_buckets);
        bucket_pairs_file.close();
    }

    bucket_pairs_manager.remove_tmp_file();

    buckets_statistics buckets_stats(num_buckets, num_kmers, num_super_kmers);

    for (minimizers_tuples_iterator it(input.data(), input.data() + input.size()); it.has_next();
         it.next()) {
        uint64_t bucket_id = m_minimizers.lookup(it.minimizer());
        uint64_t base = m_buckets.num_super_kmers_before_bucket.access(bucket_id) + bucket_id;
        uint64_t num_super_kmers_in_bucket =
            (m_buckets.num_super_kmers_before_bucket.access(bucket_id + 1) + bucket_id + 1) - base;
        assert(num_super_kmers_in_bucket > 0);
        buckets_stats.add_num_super_kmers_in_bucket(num_super_kmers_in_bucket);
        uint64_t offset_pos = 0;
        auto list = it.list();
        for (auto [offset, num_kmers_in_super_kmer] : list) {
            offsets.set(base + offset_pos++, offset);
            buckets_stats.add_num_kmers_in_super_kmer(num_super_kmers_in_bucket,
                                                      num_kmers_in_super_kmer);
        }
        assert(offset_pos == num_super_kmers_in_bucket);
    }

    m_buckets.pieces.encode(data.strings.pieces.begin(), data.strings.pieces.size(),
                            data.strings.pieces.back());
    offsets.build(m_buckets.offsets);
    m_buckets.strings.swap(data.strings.strings);

    input.close();

    return buckets_stats;
}

}  // namespace sshash