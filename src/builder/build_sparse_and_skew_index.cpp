#include <functional>

#include "include/buckets_statistics.hpp"
#include "include/builder/streaming_compact_vector_writer.hpp"

namespace sshash {

/*
    A request to extract `num_kmers_in_super_kmer` consecutive k-mers from
    `strings` starting at base position `starting_pos`. After requests are
    externally sorted by `starting_pos`, k-mer extraction reduces to a single
    forward scan over `strings`.
*/
#pragma pack(push, 4)
struct kmer_extraction_request {
    kmer_extraction_request() {}
    kmer_extraction_request(uint64_t starting_pos, uint32_t partition_id, uint32_t pos_in_bucket,
                            uint32_t num_kmers_in_super_kmer)
        : starting_pos(starting_pos)
        , partition_id(partition_id)
        , pos_in_bucket(pos_in_bucket)
        , num_kmers_in_super_kmer(num_kmers_in_super_kmer) {}

    bool operator<(kmer_extraction_request const& o) const { return starting_pos < o.starting_pos; }
    bool operator>(kmer_extraction_request const& o) const { return starting_pos > o.starting_pos; }

    static kmer_extraction_request max() {
        return kmer_extraction_request(uint64_t(-1), uint32_t(-1), uint32_t(-1), uint32_t(-1));
    }

    uint64_t starting_pos;
    uint32_t partition_id;
    uint32_t pos_in_bucket;
    uint32_t num_kmers_in_super_kmer;
};
#pragma pack(pop)

/*
    A (mphf_pos, pos_in_bucket) record used to spill the per-skew-partition
    `cvb_positions` to disk. We external-sort these by mphf_pos so the
    streaming compact_vector writer can pack the final cvb_positions file in
    a single forward pass.
*/
#pragma pack(push, 4)
struct position_tuple {
    position_tuple() {}
    position_tuple(uint64_t mphf_pos, uint32_t pib) : mphf_pos(mphf_pos), pib(pib) {}

    bool operator<(position_tuple const& o) const { return mphf_pos < o.mphf_pos; }
    bool operator>(position_tuple const& o) const { return mphf_pos > o.mphf_pos; }
    static position_tuple max() { return {uint64_t(-1), uint32_t(-1)}; }

    uint64_t mphf_pos;
    uint32_t pib;
};
#pragma pack(pop)

/*
    Per-skew-partition tmp file record (written by step 7.2 phase (B),
    consumed by phase (C)): a kmer's bit pattern + the pos_in_bucket
    we'll later pack into the partition's positions compact_vector.
*/
#pragma pack(push, 4)
template <typename Kmer>
struct skew_kmer_record_t {
    using kmer_bits_t = decltype(Kmer{}.bits);
    kmer_bits_t kmer_bits;
    uint32_t pib;
};
#pragma pack(pop)

/*
    Forward iterator over a per-skew-partition tmp file produced by phase
    (B). Yields successive Kmer values via the minimal interface (`*it`,
    `++it`) that pthash's external-memory partitioned PHF builder
    consumes.

    pthash takes the iterator by value, so it must be copyable. The
    underlying buffered_record_stream is held via shared_ptr and shared
    between copies; pthash's copy advances the shared stream state, which
    is fine because the original at the call site is unused after the
    build returns.
*/
template <typename Kmer>
struct skew_partition_kmer_iterator {
    using iterator_category = std::forward_iterator_tag;
    using value_type = Kmer;
    using difference_type = std::ptrdiff_t;
    using reference = Kmer const&;
    using pointer = Kmer const*;

    skew_partition_kmer_iterator() = default;

    void open(std::string const& filename) {
        m_stream = std::make_shared<buffered_record_stream<skew_kmer_record_t<Kmer>>>();
        m_stream->open(filename);
        if (!m_stream->empty()) m_current.bits = m_stream->current().kmer_bits;
    }

    void close() {
        if (m_stream) m_stream->close();
        m_stream.reset();
    }

    Kmer const& operator*() const { return m_current; }
    skew_partition_kmer_iterator& operator++() {
        if (!m_stream->empty()) {
            m_stream->advance();
            if (!m_stream->empty()) m_current.bits = m_stream->current().kmer_bits;
        }
        return *this;
    }

private:
    std::shared_ptr<buffered_record_stream<skew_kmer_record_t<Kmer>>> m_stream;
    Kmer m_current;
};

template <typename Kmer, typename Offsets>
void dictionary_builder<Kmer, Offsets>::build_sparse_and_skew_index(
    dictionary<Kmer, Offsets>& d)  //
{
    essentials::timer_type timer;
    timer.start();

    const uint64_t num_minimizer_positions = minimizers.num_minimizer_positions();
    const uint64_t num_minimizers = minimizers.num_minimizers();
    const uint64_t min_size = 1ULL << constants::min_l;
    const uint64_t num_bits_per_offset = strings_offsets_builder.num_bits_per_offset();

    const std::string minimizers_filename = minimizers.get_minimizers_filename();

    buckets_statistics buckets_stats(num_minimizers, num_kmers, num_minimizer_positions);

    uint64_t num_buckets_larger_than_1_not_in_skew_index = 0;
    uint64_t num_buckets_in_skew_index = 0;
    uint64_t num_minimizer_positions_of_buckets_larger_than_1 = 0;
    uint64_t num_minimizer_positions_of_buckets_in_skew_index = 0;

    /*
        Pass 1: streaming statistics over the merged minimizers file. Buckets
        are accumulated one at a time via std::ifstream-backed reads (no
        mmap), so RAM usage is bounded by max_bucket_size * sizeof(tuple).
    */
    {
        streaming_minimizer_bucket_reader reader;
        reader.open(minimizers_filename);
        std::vector<minimizer_tuple> bucket_buf;
        while (reader.has_next_bucket()) {
            reader.next_bucket(bucket_buf);
            uint64_t bucket_size = 0;
            {
                uint64_t prev = constants::invalid_uint64;
                for (auto const& mt : bucket_buf) {
                    if (mt.pos_in_seq != prev) {
                        ++bucket_size;
                        prev = mt.pos_in_seq;
                    }
                }
            }
            buckets_stats.add_bucket_size(bucket_size);
            if (bucket_size > 1) {
                if (bucket_size <= min_size) {
                    ++num_buckets_larger_than_1_not_in_skew_index;
                    num_minimizer_positions_of_buckets_larger_than_1 += bucket_size;
                } else {
                    ++num_buckets_in_skew_index;
                    num_minimizer_positions_of_buckets_in_skew_index += bucket_size;
                }
            }
            for (auto const& mt : bucket_buf) {
                buckets_stats.add_num_kmers_in_super_kmer(bucket_size, mt.num_kmers_in_super_kmer);
            }
        }
        reader.close();
    }

    assert(buckets_stats.num_buckets() == num_minimizers);

    /*
        Calculate bits needed for control codewords encoding.
        Encoding format:
            ((list_id << min_l) | (bucket_size - 2)) << 2 | status_code
    */
    const uint64_t bits_for_list_id =
        std::ceil(std::log2(buckets_stats.max_sparse_buckets_per_size() + 1));
    const uint64_t num_bits_for_control =
        std::max(num_bits_per_offset + 1, 2 + constants::min_l + bits_for_list_id);

    if (build_config.verbose) {
        std::cout << "num_bits_per_offset = " << num_bits_per_offset << std::endl;
        std::cout << "max_list_id = " << buckets_stats.max_sparse_buckets_per_size() << std::endl;
        std::cout << "bits_for_list_id = " << bits_for_list_id << std::endl;
        std::cout << "num_bits_for_control = " << num_bits_for_control << std::endl;
    }

    const uint64_t max_bucket_size = buckets_stats.max_bucket_size();
    const uint64_t log2_max_bucket_size = std::ceil(std::log2(max_bucket_size));

    uint64_t num_partitions = constants::max_l - constants::min_l + 1;
    if (max_bucket_size < min_size) {
        num_partitions = 0;
    } else if (max_bucket_size < (1ULL << constants::max_l)) {
        num_partitions = log2_max_bucket_size - constants::min_l;
    }
    assert(num_partitions <= 8);

    if (build_config.verbose) {
        std::cout << "num_buckets_larger_than_1_not_in_skew_index "
                  << num_buckets_larger_than_1_not_in_skew_index << "/"
                  << buckets_stats.num_buckets() << " ("
                  << (num_buckets_larger_than_1_not_in_skew_index * 100.0) /
                         buckets_stats.num_buckets()
                  << "%)" << std::endl;
        std::cout << "num_buckets_in_skew_index " << num_buckets_in_skew_index << "/"
                  << buckets_stats.num_buckets() << " ("
                  << (num_buckets_in_skew_index * 100.0) / buckets_stats.num_buckets() << "%)"
                  << std::endl;
        std::cout << "max_bucket_size " << max_bucket_size << std::endl;
        std::cout << "log2_max_bucket_size " << log2_max_bucket_size << std::endl;
        std::cout << "num_partitions in skew index " << num_partitions << std::endl;
    }

    /* Materialize strings_offsets now (it's needed below to decode
       pos_in_seq into absolute offsets when emitting heavy-bucket kmer
       requests). `d.m_spss.strings` is materialized later (step 8) or
       stream-saved directly. */
    strings_offsets_builder.build(d.m_spss.strings_offsets);

    /* Precompute the layout of mid_load_buckets from the bucket-size
       histogram. begin_buckets_of_size[s] is the start offset (in
       positions, not bits) of size-s bucket positions in mid_load_buckets. */
    std::vector<uint32_t> begin_buckets_of_size(min_size + 1, 0);
    for (uint64_t s = 3; s <= min_size; ++s) {
        begin_buckets_of_size[s] = static_cast<uint32_t>(  //
            begin_buckets_of_size[s - 1] + buckets_stats.num_buckets_of_size(s - 1) * (s - 1));
    }
    d.m_ssi.begin_buckets_of_size = std::move(begin_buckets_of_size);

    /* All step-7.1 outputs are spilled to disk; the in-RAM dictionary
       fields stay empty (they're populated later either from disk for
       --check or substituted by the streaming saver). */
    const uint64_t step7_run_id = pthash::clock_type::now().time_since_epoch().count();
    auto step7_path = [&](std::string const& tag) {
        std::stringstream ss;
        ss << build_config.tmp_dirname << "/sshash.tmp.run_" << step7_run_id << "." << tag
           << ".bin";
        return ss.str();
    };

    spilled.control_codewords_path = step7_path("control_codewords");
    spilled.mid_load_buckets_path = step7_path("mid_load_buckets");
    spilled.heavy_load_buckets_path = step7_path("heavy_load_buckets");

    /* Streaming writers for the two compact_vectors that get strictly
       monotonic indices during the combined pass (control_codewords:
       indexed by bucket_id == mphf hash, monotonic across buckets in
       file order; heavy_load_buckets: indexed by a single monotone
       cursor advanced inside the heavy branch). */
    streaming_compact_vector_writer control_codewords_writer;
    control_codewords_writer.open(spilled.control_codewords_path, num_minimizers,
                                  num_bits_for_control);
    streaming_compact_vector_writer heavy_load_writer;
    heavy_load_writer.open(spilled.heavy_load_buckets_path,
                           num_minimizer_positions_of_buckets_in_skew_index, num_bits_per_offset);

    /* mid_load: per-size tmp files of raw uint64_t positions. Each file is
       written monotonically within its size class. After the combined
       pass we stream them in size order through a streaming
       compact_vector writer to assemble the final mid_load_buckets file. */
    auto mid_load_per_size_path = [&](uint64_t s) {
        std::stringstream ss;
        ss << build_config.tmp_dirname << "/sshash.tmp.run_" << step7_run_id << ".mid_load_size_"
           << s << ".bin";
        return ss.str();
    };
    std::vector<std::ofstream> mid_load_per_size(min_size + 1);
    for (uint64_t s = 2; s <= min_size; ++s) {
        if (buckets_stats.num_buckets_of_size(s) == 0) continue;
        mid_load_per_size[s].open(mid_load_per_size_path(s),
                                  std::ofstream::binary | std::ofstream::trunc);
        if (!mid_load_per_size[s].is_open()) {
            throw std::runtime_error("cannot open mid_load per-size tmp file");
        }
    }

    std::vector<uint64_t> list_id_per_size(min_size + 1, 0);
    uint64_t heavy_load_cursor = 0;

    /* Per-partition kmer counts; filled during the heavy branch of the
       combined pass below. */
    std::vector<uint64_t> num_kmers_in_partition(num_partitions, 0);

    /* Skew-index tmp file naming. */
    const uint64_t skew_run_id = pthash::clock_type::now().time_since_epoch().count();
    auto request_run_filename = [&](uint64_t id) {
        std::stringstream ss;
        ss << build_config.tmp_dirname << "/sshash.tmp.run_" << skew_run_id << ".kmer_requests."
           << id << ".bin";
        return ss.str();
    };
    auto skew_partition_filename = [&](uint64_t pid) {
        std::stringstream ss;
        ss << build_config.tmp_dirname << "/sshash.tmp.run_" << skew_run_id << ".skew_kmers." << pid
           << ".bin";
        return ss.str();
    };

    /* External-sort buffer for kmer-extraction requests (formerly step 7.2
       phase A; now folded into the combined pass). Capped at ram_limit/8
       so heap fragmentation across steps doesn't push peak RSS past the
       --ram-limit budget. */
    std::atomic<uint64_t> num_request_runs{0};
    const uint64_t request_buffer_capacity =
        std::max<uint64_t>(uint64_t(1) << 16, (build_config.ram_limit_in_GiB * essentials::GiB) /
                                                  (8 * sizeof(kmer_extraction_request)));
    std::vector<kmer_extraction_request> request_buffer;
    request_buffer.reserve(request_buffer_capacity);
    auto flush_request_buffer = [&]() {
        if (request_buffer.empty()) return;
        parallel_sort(request_buffer, build_config.num_threads,
                      [](kmer_extraction_request const& a, kmer_extraction_request const& b) {
                          return a.starting_pos < b.starting_pos;
                      });
        const uint64_t id = num_request_runs.fetch_add(1);
        const std::string fn = request_run_filename(id);
        std::ofstream out(fn, std::ofstream::binary);
        if (!out.is_open()) throw std::runtime_error("cannot open file");
        out.write(reinterpret_cast<char const*>(request_buffer.data()),
                  request_buffer.size() * sizeof(kmer_extraction_request));
        out.close();
        request_buffer.clear();
    };

    /* Map bucket size → partition_id for heavy buckets. num_partitions <= 8
       so this loop is constant time. */
    auto partition_for_size = [&](uint64_t bucket_size) -> uint64_t {
        assert(bucket_size > min_size);
        uint64_t pid = 0;
        uint64_t upper = 2 * min_size;
        while (bucket_size > upper && pid + 1 < num_partitions) {
            upper *= 2;
            ++pid;
        }
        return pid;
    };

    /*
        Combined pass: stream the merged minimizers file once and, per
        bucket, write the appropriate part of the sparse index DIRECTLY TO
        DISK via streaming compact_vector writers (control_codewords and
        heavy_load_buckets) or per-size raw-value tmp files (mid_load).
        For heavy buckets we also emit kmer-extraction requests in-line.

        Buckets are visited in mphf-hash (= bucket_id) order, so writes to
        control_codewords are strictly monotonic. heavy_load_cursor is also
        monotonic across the whole pass. mid_load per-size cursors are
        each monotonic within their size class.
    */
    {
        streaming_minimizer_bucket_reader reader;
        reader.open(minimizers_filename);
        std::vector<minimizer_tuple> bucket_buf;
        while (reader.has_next_bucket()) {
            const uint64_t bucket_id = reader.next_bucket(bucket_buf);
            uint64_t bucket_size = 0;
            {
                uint64_t prev = constants::invalid_uint64;
                for (auto const& mt : bucket_buf) {
                    if (mt.pos_in_seq != prev) {
                        ++bucket_size;
                        prev = mt.pos_in_seq;
                    }
                }
            }

            if (bucket_size == 1) {
                /* Singleton: code = |offset|0|, LSB = 0. */
                const uint64_t code = bucket_buf.front().pos_in_seq << 1;
                assert(code < (uint64_t(1) << num_bits_for_control));
                control_codewords_writer.set(bucket_id, code);
            } else if (bucket_size <= min_size) {
                /* Mid-load: write positions to per-size raw file at the
                   per-size cursor; assign the next list_id for this size. */
                const uint64_t list_id = list_id_per_size[bucket_size]++;
                const uint64_t code =
                    (((list_id << constants::min_l) | (bucket_size - 2)) << 2) | 1;
                assert(code < (uint64_t(1) << num_bits_for_control));
                control_codewords_writer.set(bucket_id, code);

                auto& out = mid_load_per_size[bucket_size];
                uint64_t prev_pos_in_seq = constants::invalid_uint64;
                for (auto const& mt : bucket_buf) {
                    if (mt.pos_in_seq != prev_pos_in_seq) {
                        const uint64_t v = mt.pos_in_seq;
                        out.write(reinterpret_cast<char const*>(&v), sizeof(uint64_t));
                        prev_pos_in_seq = mt.pos_in_seq;
                    }
                }
            } else {
                /* Heavy: write positions at the monotone heavy_load_cursor,
                   set the codeword (encodes the start offset and partition
                   id), and emit kmer-extraction requests for each
                   super-kmer in the bucket. */
                const uint64_t partition_id = partition_for_size(bucket_size);
                assert(partition_id < num_partitions);
                const uint64_t bucket_begin = heavy_load_cursor;
                const uint64_t code = (((bucket_begin << 3) | partition_id) << 2) | 3;
                assert(code < (uint64_t(1) << num_bits_for_control));
                control_codewords_writer.set(bucket_id, code);

                uint32_t pos_in_bucket = uint32_t(-1);
                uint64_t prev_pos_in_seq = constants::invalid_uint64;
                for (auto const& mt : bucket_buf) {
                    num_kmers_in_partition[partition_id] += mt.num_kmers_in_super_kmer;
                    if (mt.pos_in_seq != prev_pos_in_seq) {
                        heavy_load_writer.set(heavy_load_cursor++, mt.pos_in_seq);
                        prev_pos_in_seq = mt.pos_in_seq;
                        ++pos_in_bucket;
                    }
                    assert(mt.pos_in_seq >= mt.pos_in_kmer);
                    const uint64_t abs_offset =
                        d.m_spss.strings_offsets.decode(mt.pos_in_seq).absolute_offset;
                    const uint64_t starting_pos = abs_offset - mt.pos_in_kmer;
                    if (request_buffer.size() == request_buffer_capacity) flush_request_buffer();
                    request_buffer.emplace_back(starting_pos, uint32_t(partition_id), pos_in_bucket,
                                                uint32_t(mt.num_kmers_in_super_kmer));
                }
            }
        }
        reader.close();
        flush_request_buffer();
    }

    /* Finalize the directly-streamed compact_vector files. */
    control_codewords_writer.finalize();
    heavy_load_writer.finalize();

    /* Close per-size mid_load files. */
    for (uint64_t s = 2; s <= min_size; ++s) {
        if (mid_load_per_size[s].is_open()) mid_load_per_size[s].close();
    }

    /* Concatenate per-size mid_load files in size order into the final
       mid_load_buckets compact_vector file via the streaming writer.
       Each per-size file holds raw uint64_t values written monotonically
       within its size class; we just stream them through, packing into
       num_bits_per_offset-bit fields at the precomputed begin offset for
       each size. */
    {
        streaming_compact_vector_writer mid_load_writer;
        mid_load_writer.open(spilled.mid_load_buckets_path,
                             num_minimizer_positions_of_buckets_larger_than_1, num_bits_per_offset);
        uint64_t global_index = 0;
        for (uint64_t s = 2; s <= min_size; ++s) {
            const uint64_t expected = buckets_stats.num_buckets_of_size(s) * s;
            if (expected == 0) continue;
            std::ifstream in(mid_load_per_size_path(s), std::ifstream::binary);
            if (!in.is_open()) {
                throw std::runtime_error("cannot reopen mid_load per-size tmp file");
            }
            for (uint64_t i = 0; i != expected; ++i) {
                uint64_t v;
                in.read(reinterpret_cast<char*>(&v), sizeof(uint64_t));
                if (in.gcount() != static_cast<std::streamsize>(sizeof(uint64_t))) {
                    throw std::runtime_error("mid_load per-size tmp file truncated");
                }
                mid_load_writer.set(global_index++, v);
            }
            in.close();
            std::remove(mid_load_per_size_path(s).c_str());
        }
        mid_load_writer.finalize();
    }

    timer.stop();
    build_stats.add("step 7.1 (build sparse index)",
                    musec_as_seconds_str(uint64_t(timer.elapsed())).c_str());
    if (build_config.verbose) {
        print_time(uint64_t(timer.elapsed()), "step 7.1 (build sparse index)");
    }
    timer.reset();

    if (num_buckets_in_skew_index == 0) {
        if (build_config.verbose) buckets_stats.print_less();
        return;
    }

    /*
        step 2. build skew index

        Phases (B) and (C) below; phase (A) was folded into the combined
        sparse pass above. Phase (B) extracts k-mers from `strings` in a
        single forward sweep guided by the externally-sorted requests, and
        phase (C) builds the per-partition MPHF + cvb_positions on disk
        from the per-partition kmer files.
    */
    timer.start();

    if (build_config.verbose) {
        uint64_t total_kmers_in_skew = 0;
        for (uint64_t p = 0; p != num_partitions; ++p) {
            total_kmers_in_skew += num_kmers_in_partition[p];
            std::cout << "  partition = " << p
                      << ": num kmers in partition = " << num_kmers_in_partition[p] << std::endl;
        }
        std::cout << "num kmers in skew index = " << total_kmers_in_skew << " ("
                  << (total_kmers_in_skew * 100.0) / buckets_stats.num_kmers() << "%)" << std::endl;
    }

    /* (B) sequential extraction over `strings` -> per-partition kmer tmp files */
    {
        struct request_run_names_iterator {
            request_run_names_iterator(std::string const& tmp_dirname, uint64_t skew_run_id)
                : i(0), skew_run_id(skew_run_id), tmp_dirname(tmp_dirname) {}

            std::string operator*() const {
                std::stringstream ss;
                ss << tmp_dirname << "/sshash.tmp.run_" << skew_run_id << ".kmer_requests." << i
                   << ".bin";
                return ss.str();
            }
            void operator++() { ++i; }

            uint64_t i;
            uint64_t skew_run_id;
            std::string tmp_dirname;
        };

        request_run_names_iterator names_it(build_config.tmp_dirname, skew_run_id);
        file_merging_iterator<kmer_extraction_request> merger(names_it, num_request_runs.load());

        std::vector<std::ofstream> partition_writers(num_partitions);
        for (uint64_t p = 0; p != num_partitions; ++p) {
            if (num_kmers_in_partition[p] == 0) continue;
            partition_writers[p].open(skew_partition_filename(p),
                                      std::ofstream::binary | std::ofstream::trunc);
            if (!partition_writers[p].is_open()) {
                throw std::runtime_error("cannot open skew-partition tmp file");
            }
        }

        const uint64_t k = build_config.k;
        const bool canonical = build_config.canonical;
        auto strings_reader = strings_builder.make_reader();
        kmer_iterator<Kmer, disk_backed_strings::reader> kmer_it(strings_reader, k);

        while (merger.has_next())  //
        {
            const kmer_extraction_request req = *merger;
            kmer_it.at(Kmer::bits_per_char * req.starting_pos);
            for (uint32_t i = 0; i != req.num_kmers_in_super_kmer; ++i) {
                Kmer kmer = kmer_it.get();
                if (canonical) {
                    Kmer kmer_rc = kmer;
                    kmer_rc.reverse_complement_inplace(k);
                    kmer = std::min(kmer, kmer_rc);
                }
                auto& w = partition_writers[req.partition_id];
                skew_kmer_record_t<Kmer> rec{kmer.bits, req.pos_in_bucket};
                w.write(reinterpret_cast<char const*>(&rec), sizeof(rec));
                kmer_it.next();
            }
            merger.next();
        }
        merger.close();

        for (auto& w : partition_writers) {
            if (w.is_open()) w.close();
        }

        for (uint64_t i = 0; i != num_request_runs.load(); ++i) {
            std::remove(request_run_filename(i).c_str());
        }
    }

    /*
        (C) per-partition MPHF + cvb_positions build, both on disk.

        Per partition:
          (1) Build MPHF in external memory by streaming the partition's
              kmer file (pthash spills hashes to tmp_dirname under its own
              ram budget).
          (2) Stream-read the kmer file, compute F(kmer), emit
              (F(kmer), pos_in_bucket) tuples to disk; external-sort by
              F(kmer); stream sorted tuples through a
              streaming_compact_vector_writer to produce the partition's
              cvb_positions tmp file.

        Only the MPHF itself is held in RAM (pthash returns it as an
        in-memory struct); cvb_positions is fully spilled.
    */
    {
        spilled.skew_positions_paths.assign(num_partitions, std::string());
        spilled.skew_mphfs_paths.assign(num_partitions, std::string());
        std::vector<kmers_pthash_type<Kmer>> mphfs;
        mphfs.resize(num_partitions);

        pthash::build_configuration mphf_build_config;
        mphf_build_config.lambda = build_config.lambda + 2.0;
        mphf_build_config.alpha = 0.94;
        mphf_build_config.seed = util::get_seed_for_hash_function(build_config);
        mphf_build_config.verbose = false;
        util::configure_mphf_threads_and_partition(mphf_build_config, build_config.num_threads,
                                                   build_config.ram_limit_in_GiB,
                                                   build_config.verbose, "skew partition MPHF");
        mphf_build_config.ram = (build_config.ram_limit_in_GiB * essentials::GiB) / 2;
        mphf_build_config.tmp_dir = build_config.tmp_dirname;

        const uint64_t pos_run_basename_id = pthash::clock_type::now().time_since_epoch().count();
        auto pos_run_filename = [&](uint64_t partition_id, uint64_t id) {
            std::stringstream ss;
            ss << build_config.tmp_dirname << "/sshash.tmp.run_" << pos_run_basename_id
               << ".pos_runs.p" << partition_id << "." << id << ".bin";
            return ss.str();
        };
        auto skew_positions_filename = [&](uint64_t partition_id) {
            std::stringstream ss;
            ss << build_config.tmp_dirname << "/sshash.tmp.run_" << pos_run_basename_id
               << ".skew_positions.p" << partition_id << ".bin";
            return ss.str();
        };
        auto skew_mphf_filename = [&](uint64_t partition_id) {
            std::stringstream ss;
            ss << build_config.tmp_dirname << "/sshash.tmp.run_" << pos_run_basename_id
               << ".skew_mphf.p" << partition_id << ".bin";
            return ss.str();
        };

        /* Capped at ram_limit/8: this buffer is alive during phase C
           alongside pthash's parallel-build memory and the currently-
           building partition's MPHF, so it has to share the RAM budget. */
        const uint64_t pos_buffer_capacity = std::max<uint64_t>(
            uint64_t(1) << 16,
            (build_config.ram_limit_in_GiB * essentials::GiB) / (8 * sizeof(position_tuple)));

        uint64_t lower = min_size;
        uint64_t upper = 2 * lower;
        uint64_t num_bits_per_pos = constants::min_l + 1;
        if (num_partitions == 1) {
            upper = max_bucket_size;
            num_bits_per_pos = log2_max_bucket_size;
        }

        for (uint64_t partition_id = 0; partition_id != num_partitions; ++partition_id)  //
        {
            const uint64_t n = num_kmers_in_partition[partition_id];

            if (build_config.verbose) {
                std::cout << "  lower = " << lower << "; upper = " << upper
                          << "; num_bits_per_pos = " << num_bits_per_pos
                          << "; num_kmers_in_partition = " << n << std::endl;
            }

            if (n > 0)  //
            {
                const std::string kmer_fn = skew_partition_filename(partition_id);

                if (build_config.verbose) {
                    const uint64_t avg_partition_size =
                        pthash::compute_avg_partition_size(n, mphf_build_config);
                    const uint64_t pthash_num_partitions =
                        pthash::compute_num_partitions(n, avg_partition_size);
                    assert(pthash_num_partitions > 0);
                    std::cout << "    building MPHF (external memory) with "
                              << mphf_build_config.num_threads << " threads and "
                              << pthash_num_partitions
                              << " partitions (avg. partition size = " << avg_partition_size
                              << ")..." << std::endl;
                }

                /* (1) Build the MPHF by streaming kmers from the partition file. */
                auto& F = mphfs[partition_id];
                {
                    skew_partition_kmer_iterator<Kmer> iter;
                    iter.open(kmer_fn);
                    F.build_in_external_memory(iter, n, mphf_build_config);
                    iter.close();
                }

                if (build_config.verbose) {
                    std::cout << "    built mphs[" << partition_id << "] for " << F.num_keys()
                              << " kmers; bits/key = "
                              << static_cast<double>(F.num_bits()) / F.num_keys() << std::endl;
                }

                /* (2a) Stream-read kmer file, compute F(kmer), externally
                       sort (F(kmer), pos_in_bucket) tuples by F(kmer). */
                std::atomic<uint64_t> pos_num_runs{0};
                {
                    std::vector<position_tuple> pos_buffer;
                    pos_buffer.reserve(pos_buffer_capacity);
                    auto flush_pos_buffer = [&]() {
                        if (pos_buffer.empty()) return;
                        parallel_sort(pos_buffer, build_config.num_threads,
                                      [](position_tuple const& a, position_tuple const& b) {
                                          return a.mphf_pos < b.mphf_pos;
                                      });
                        const uint64_t id = pos_num_runs.fetch_add(1);
                        std::ofstream out(pos_run_filename(partition_id, id),
                                          std::ofstream::binary);
                        if (!out.is_open()) {
                            throw std::runtime_error("cannot open positions tuple run file");
                        }
                        out.write(reinterpret_cast<char const*>(pos_buffer.data()),
                                  pos_buffer.size() * sizeof(position_tuple));
                        out.close();
                        pos_buffer.clear();
                    };

                    buffered_record_stream<skew_kmer_record_t<Kmer>> rec_stream;
                    rec_stream.open(kmer_fn);
                    for (uint64_t i = 0; i != n; ++i) {
                        assert(!rec_stream.empty());
                        auto const& rec = rec_stream.current();
                        Kmer kmer;
                        kmer.bits = rec.kmer_bits;
                        const uint64_t pos = F(kmer);
                        if (pos_buffer.size() == pos_buffer_capacity) flush_pos_buffer();
                        pos_buffer.emplace_back(pos, rec.pib);
                        rec_stream.advance();
                    }
                    rec_stream.close();
                    std::remove(kmer_fn.c_str());
                    flush_pos_buffer();
                }

                /* (2b) Stream sorted tuples through the streaming
                       compact_vector writer to produce the partition's
                       cvb_positions tmp file. */
                {
                    spilled.skew_positions_paths[partition_id] =
                        skew_positions_filename(partition_id);
                    streaming_compact_vector_writer pos_writer;
                    pos_writer.open(spilled.skew_positions_paths[partition_id], n,
                                    num_bits_per_pos);

                    struct pos_run_names_iterator {
                        pos_run_names_iterator(uint64_t partition_id,
                                               std::function<std::string(uint64_t, uint64_t)> fn)
                            : i(0), partition_id(partition_id), fn(std::move(fn)) {}
                        std::string operator*() { return fn(partition_id, i); }
                        void operator++() { ++i; }
                        uint64_t i;
                        uint64_t partition_id;
                        std::function<std::string(uint64_t, uint64_t)> fn;
                    };
                    pos_run_names_iterator names_it(partition_id, pos_run_filename);
                    file_merging_iterator<position_tuple> merger(names_it, pos_num_runs.load());
                    while (merger.has_next()) {
                        position_tuple pt = *merger;
                        pos_writer.set(pt.mphf_pos, pt.pib);
                        merger.next();
                    }
                    merger.close();
                    pos_writer.finalize();
                }

                /* Cleanup the position-tuple run files. */
                for (uint64_t i = 0; i != pos_num_runs.load(); ++i) {
                    std::remove(pos_run_filename(partition_id, i).c_str());
                }

                if (build_config.verbose) {
                    std::cout << "    built positions[" << partition_id << "] for " << n
                              << " kmers; bits/key = " << num_bits_per_pos << std::endl;
                }

                /* Spill the partition's MPHF to disk (no longer needed
                   during build) and free its in-RAM footprint. The
                   accumulating skew MPHFs were the dominant resident
                   memory at the end of phase C. */
                spilled.skew_mphfs_paths[partition_id] = skew_mphf_filename(partition_id);
                essentials::save(F, spilled.skew_mphfs_paths[partition_id].c_str());
                F = kmers_pthash_type<Kmer>{};
            }

            /* advance partition state for the next iteration */
            lower = upper;
            upper = 2 * lower;
            num_bits_per_pos += 1;
            if (partition_id + 1 == num_partitions - 1) {
                upper = max_bucket_size;
                num_bits_per_pos = log2_max_bucket_size;
            }
        }

        d.m_ssi.ski.mphfs = std::move(mphfs);
        /* d.m_ssi.ski.positions stays empty here; it will be populated
           either by step 8 (materialize) or substituted at stream-save. */
    }

    timer.stop();

    build_stats.add("step 7.2 (build skew index)",
                    musec_as_seconds_str(uint64_t(timer.elapsed())).c_str());

    if (build_config.verbose) {
        print_time(uint64_t(timer.elapsed()), "step 7.2 (build skew index)");
        buckets_stats.print_less();
    }
}

}  // namespace sshash
