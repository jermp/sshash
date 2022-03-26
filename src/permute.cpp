#include <iostream>

#include "../include/gz/zip_stream.cpp"
#include "../external/pthash/external/essentials/include/essentials.hpp"
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/builder/cover.hpp"
#include "../include/builder/util.hpp"

using namespace sshash;

struct permute_data {
    permute_data() : num_runs_abundances(0), num_sequences(0) {}
    uint64_t num_runs_abundances;
    uint64_t num_sequences;
    std::vector<vertex> vertices;
};

void parse_file(std::istream& is, permute_data& data, build_configuration const& build_config) {
    std::string sequence;
    uint64_t k = build_config.k;
    uint64_t num_bases = 0;
    uint64_t num_kmers = 0;
    uint64_t seq_len = 0;

    uint64_t sum_of_abundances = 0;
    uint64_t num_sequences_diff_abs = 0;  // num sequences whose kmers have different abundances
    uint64_t num_sequences_all_mfa = 0;   // num sequences whose kmers have same abundance == mfa
    data.num_runs_abundances = 0;
    data.num_sequences = 0;

    auto parse_header = [&]() {
        if (sequence.empty()) return;

        /*
            Heder format:
            >[id] LN:i:[seq_len] ab:Z:[ab_seq]
            where [ab_seq] is a space-separated sequence of integer counters (the abundances),
            whose length is equal to [seq_len]-k+1
        */

        // example header: '>12 LN:i:41 ab:Z:2 2 2 2 2 2 2 2 2 2 2'

        expect(sequence[0], '>');
        uint64_t i = 0;
        i = sequence.find_first_of(' ', i);
        if (i == std::string::npos) throw parse_runtime_error();

        i += 1;
        expect(sequence[i + 0], 'L');
        expect(sequence[i + 1], 'N');
        expect(sequence[i + 2], ':');
        expect(sequence[i + 3], 'i');
        expect(sequence[i + 4], ':');
        i += 5;
        uint64_t j = sequence.find_first_of(' ', i);
        if (j == std::string::npos) throw parse_runtime_error();

        seq_len = std::strtoull(sequence.data() + i, nullptr, 10);
        i = j + 1;
        expect(sequence[i + 0], 'a');
        expect(sequence[i + 1], 'b');
        expect(sequence[i + 2], ':');
        expect(sequence[i + 3], 'Z');
        expect(sequence[i + 4], ':');
        i += 5;

        bool kmers_have_all_mfa = true;
        bool kmers_have_different_abundances = false;
        uint64_t front = constants::invalid;
        uint64_t back = constants::invalid;

        for (uint64_t j = 0, prev_abundance = constants::invalid; j != seq_len - k + 1; ++j) {
            uint64_t abundance = std::strtoull(sequence.data() + i, nullptr, 10);
            sum_of_abundances += abundance;
            i = sequence.find_first_of(' ', i) + 1;

            /* set front and back */
            if (j == 0) front = abundance;
            if (j == seq_len - k) back = abundance;

            /* accumulate statistics */
            if (abundance != constants::most_frequent_abundance) kmers_have_all_mfa = false;
            if (j > 0 and abundance != prev_abundance) kmers_have_different_abundances = true;

            /* count the number of runs */
            if (abundance != prev_abundance) data.num_runs_abundances += 1;

            prev_abundance = abundance;
        }

        num_sequences_diff_abs += kmers_have_different_abundances;
        num_sequences_all_mfa += kmers_have_all_mfa;

        constexpr bool sign = true;
        data.vertices.emplace_back(data.num_sequences, front, back, sign);
    };

    while (!is.eof()) {
        std::getline(is, sequence);  // header sequence
        parse_header();

        std::getline(is, sequence);  // DNA sequence
        if (sequence.size() < k) continue;

        if (++data.num_sequences % 100000 == 0) {
            std::cout << "read " << data.num_sequences << " sequences, " << num_bases << " bases, "
                      << num_kmers << " kmers" << std::endl;
        }

        num_bases += sequence.size();
        num_kmers += sequence.size() - k + 1;

        if (seq_len != sequence.size()) {
            std::cout << "ERROR: expected a sequence of length " << seq_len
                      << " but got one of length " << sequence.size() << std::endl;
            throw std::runtime_error("file is malformed");
        }
    }

    assert(data.vertices.size() == data.num_sequences);
    assert(data.num_runs_abundances >= data.num_sequences);

    std::cout << "read " << data.num_sequences << " sequences, " << num_bases << " bases, "
              << num_kmers << " kmers" << std::endl;
    std::cout << "sum_of_abundances " << sum_of_abundances << std::endl;
    std::cout << "num_sequences whose kmers have different abundances: " << num_sequences_diff_abs
              << "/" << data.num_sequences << " ("
              << (num_sequences_diff_abs * 100.0) / data.num_sequences << "%)" << std::endl;
    std::cout << "num_sequences whose kmers all have the same abundance != mfa: "
              << (data.num_sequences - num_sequences_diff_abs) << "/" << data.num_sequences << " ("
              << ((data.num_sequences - num_sequences_diff_abs) * 100.0) / data.num_sequences
              << "%)" << std::endl;
    std::cout << "num_sequences whose kmers all have the same abundance == mfa: "
              << num_sequences_all_mfa << "/" << (data.num_sequences - num_sequences_diff_abs)
              << " ("
              << (num_sequences_all_mfa * 100.0) / (data.num_sequences - num_sequences_diff_abs)
              << "%)" << std::endl;
}

permute_data parse_file(std::string const& filename, build_configuration const& build_config) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    std::cout << "reading file '" << filename << "'..." << std::endl;
    permute_data data;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        parse_file(zis, data, build_config);
    } else {
        parse_file(is, data, build_config);
    }
    is.close();
    return data;
}

struct lines_iterator {
    lines_iterator(uint8_t const* begin, uint8_t const* end) : m_begin(begin), m_end(end) {
        advance_to_next();
    }

    void advance_to_next() {
        m_header = read_line();
        m_sequence = read_line();
    }

    std::string const& header() const { return m_header; }
    std::string const& sequence() const { return m_sequence; }
    bool has_next() const { return !m_header.empty() and !m_sequence.empty(); }

private:
    uint8_t const* m_begin;
    uint8_t const* m_end;
    std::string m_header;
    std::string m_sequence;

    std::string read_line() {
        uint8_t const* begin = m_begin;
        while (m_begin != m_end and *m_begin++ != '\n')
            ;
        if (begin == m_begin) return std::string("");
        return std::string(reinterpret_cast<const char*>(begin), m_begin - begin - 1);
    }
};

void permute_and_write(std::istream& is, std::string const& output_filename,
                       std::string const& tmp_dirname, pthash::compact_vector const& permutation) {
    constexpr uint64_t limit = 1 * essentials::GB;
    std::vector<std::pair<std::string, std::string>> buffer;  // (header, dna)

    std::string run_identifier =
        std::to_string(pthash::clock_type::now().time_since_epoch().count());

    std::string header_sequence;
    std::string dna_sequence;
    uint64_t num_sequences = permutation.size();
    uint64_t num_bases = 0;
    uint64_t bytes = 0;
    uint64_t num_files_to_merge = 0;

    auto get_tmp_output_filename = [&](uint64_t id) {
        return tmp_dirname + "/tmp.run" + run_identifier + "." + std::to_string(id);
    };

    auto sort_and_flush = [&]() {
        if (buffer.empty()) return;

        std::cout << "sorting buffer..." << std::endl;
        std::sort(buffer.begin(), buffer.end(), [&](auto const& p_x, auto const& p_y) {
            assert(p_x.first.front() == '>');
            assert(p_y.first.front() == '>');
            uint64_t seq_id_x = std::strtoull(p_x.first.data() + 1, nullptr, 10);
            uint64_t seq_id_y = std::strtoull(p_y.first.data() + 1, nullptr, 10);
            return permutation[seq_id_x] < permutation[seq_id_y];
        });

        auto tmp_output_filename = get_tmp_output_filename(num_files_to_merge);
        std::cout << "saving to file '" << tmp_output_filename << "'..." << std::endl;
        std::ofstream out(tmp_output_filename.c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open file");
        for (auto const& seq : buffer) out << seq.first << '\n' << seq.second << '\n';
        out.close();

        buffer.clear();
        bytes = 0;
        num_files_to_merge += 1;
    };

    for (uint64_t i = 0; i != num_sequences; ++i) {
        std::getline(is, header_sequence);
        std::getline(is, dna_sequence);

        uint64_t seq_bytes = header_sequence.size() + dna_sequence.size() + 16;
        if (bytes + seq_bytes > limit) sort_and_flush();

        bytes += seq_bytes;
        buffer.emplace_back(header_sequence, dna_sequence);
        num_bases += dna_sequence.size();

        if (i != 0 and i % 100000 == 0) {
            std::cout << "read " << i << " sequences, " << num_bases << " bases" << std::endl;
        }
    }
    sort_and_flush();

    std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases" << std::endl;

    /* merge */
    {
        assert(num_files_to_merge > 0);
        std::cout << "files to merge = " << num_files_to_merge << std::endl;

        std::vector<lines_iterator> iterators;
        std::vector<uint32_t> idx_heap;
        iterators.reserve(num_files_to_merge);
        idx_heap.reserve(num_files_to_merge);
        std::vector<mm::file_source<uint8_t>> mm_files(num_files_to_merge);

        auto heap_idx_comparator = [&](uint32_t i, uint32_t j) {
            assert(iterators[i].header().front() == '>');
            assert(iterators[j].header().front() == '>');
            uint64_t seq_id_x = std::strtoull(iterators[i].header().data() + 1, nullptr, 10);
            uint64_t seq_id_y = std::strtoull(iterators[j].header().data() + 1, nullptr, 10);
            return permutation[seq_id_x] > permutation[seq_id_y];
        };

        auto advance_heap_head = [&]() {
            auto idx = idx_heap.front();
            iterators[idx].advance_to_next();
            if (iterators[idx].has_next()) {  // percolate down the head
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
        for (uint64_t i = 0; i != num_files_to_merge; ++i) {
            auto tmp_output_filename = get_tmp_output_filename(i);
            mm_files[i].open(tmp_output_filename, mm::advice::sequential);
            iterators.emplace_back(mm_files[i].data(), mm_files[i].data() + mm_files[i].size());
            idx_heap.push_back(i);
        }
        std::make_heap(idx_heap.begin(), idx_heap.end(), heap_idx_comparator);

        std::ofstream out(output_filename.c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open file");

        uint64_t num_written_sequences = 0;
        while (!idx_heap.empty()) {
            auto const& it = iterators[idx_heap.front()];
            out << it.header() << '\n' << it.sequence() << '\n';
            num_written_sequences += 1;
            if (num_written_sequences % 100000 == 0) {
                std::cout << "written sequences = " << num_written_sequences << "/" << num_sequences
                          << std::endl;
            }
            advance_heap_head();
        }
        std::cout << "written sequences = " << num_written_sequences << "/" << num_sequences
                  << std::endl;
        out.close();
        assert(num_written_sequences == num_sequences);

        /* remove tmp files */
        for (uint64_t i = 0; i != num_files_to_merge; ++i) {
            mm_files[i].close();
            auto tmp_output_filename = get_tmp_output_filename(i);
            std::remove(tmp_output_filename.c_str());
        }
    }
}

void permute_and_write(std::string const& input_filename, std::string const& output_filename,
                       std::string const& tmp_dirname, pthash::compact_vector const& permutation) {
    std::ifstream is(input_filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + input_filename + "'");
    std::cout << "reading file '" << input_filename << "'..." << std::endl;
    if (util::ends_with(input_filename, ".gz")) {
        zip_istream zis(is);
        permute_and_write(zis, output_filename, tmp_dirname, permutation);
    } else {
        permute_and_write(is, output_filename, tmp_dirname, permutation);
    }
    is.close();
}

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);

    /* mandatory arguments */
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:\n"
               "\t- without duplicate nor invalid kmers\n"
               "\t- one DNA sequence per line\n"
               "\t- with also kmers' abundances.\n"
               "\tFor example, it could be the de Bruijn graph topology output by BCALM.");

    parser.add("k", "K-mer length (must be <= " + std::to_string(constants::max_k) + ").");

    /* optional arguments */
    parser.add("output_filename", "Output file where the permuted collection will be written.",
               "-o", false);
    parser.add("tmp_dirname",
               "Temporary directory used for merging in external memory. Default is directory '" +
                   constants::default_tmp_dirname + "'.",
               "-d", false);

    if (!parser.parse()) return 1;

    auto input_filename = parser.get<std::string>("input_filename");
    auto k = parser.get<uint64_t>("k");

    build_configuration build_config;
    build_config.k = k;

    std::string output_filename = input_filename + ".permuted";
    if (parser.parsed("output_filename")) {
        output_filename = parser.get<std::string>("output_filename");
    }

    std::string tmp_dirname = constants::default_tmp_dirname;
    if (parser.parsed("tmp_dirname")) tmp_dirname = parser.get<std::string>("tmp_dirname");

    std::string permutation_filename = tmp_dirname + "/tmp.permutation";

    auto data = parse_file(input_filename, build_config);

    {
        uint64_t R_lower = data.num_runs_abundances - data.vertices.size() + 1;
        std::cout << "The trivial lower bound (too optimistic) assumes we are able to concatenate "
                     "all sequences : R_lo = "
                  << R_lower << std::endl;

        std::sort(data.vertices.begin(), data.vertices.end(), [](auto const& x, auto const& y) {
            /* sort on front */
            if (x.front != y.front) return x.front < y.front;
            if (x.back != y.back) return x.back < y.back;
            return x.id < y.id;
        });

        /* (abundance, num_seqs_with_front=abundance) */
        // We assume there are less than 2^32 sequences and that
        // the largest abundance fits into a 32-bit uint.
        std::unordered_map<uint32_t, uint32_t> abundance_map;

        uint64_t prev_front = data.vertices.front().front;
        uint64_t count = 0;
        for (auto const& vertex : data.vertices) {
            if (vertex.front != prev_front) {
                abundance_map[prev_front] = count;
                count = 0;
            }
            count += 1;
            prev_front = vertex.front;
        }
        abundance_map[prev_front] = count;

        uint64_t R = data.num_runs_abundances;
        for (auto const& vertex : data.vertices) {
            uint64_t back = vertex.back;
            auto it = abundance_map.find(back);
            if (it != abundance_map.cend()) {  // found
                if ((*it).second > 0) {        // if it is 0, we cannot find a match

                    /* We clearly cannot create more mergings than num_sequences - 1,
                       thus R cannot be lower than R_lower. */
                    if (R == R_lower) break;

                    (*it).second -= 1;
                    R -= 1;
                }
            }
        }
        std::cout << "Computed lower bound: R_hi = " << R << std::endl;
    }

    {
        /* compute cover */
        cover c(data.num_sequences, data.num_runs_abundances);
        assert(data.vertices.size() == data.num_sequences);
        c.compute(data.vertices);
        c.save(permutation_filename);
        std::vector<vertex>().swap(data.vertices);
    }

    /* permute */
    pthash::compact_vector permutation;
    {
        std::ifstream is(permutation_filename.c_str());
        if (!is.good()) {
            throw std::runtime_error("error in opening the file '" + permutation_filename + "'");
        }
        pthash::compact_vector::builder cv_builder(
            data.num_sequences,
            data.num_sequences == 1 ? 1 : std::ceil(std::log2(data.num_sequences)));
        for (uint64_t i = 0; i != data.num_sequences; ++i) {
            uint64_t position = 0;
            is >> position;
            cv_builder.set(position, i);
        }
        is.close();
        cv_builder.build(permutation);
    }

    permute_and_write(input_filename, output_filename, tmp_dirname, permutation);
    std::remove(permutation_filename.c_str());

    return 0;
}
