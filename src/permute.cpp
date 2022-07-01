#include <iostream>

#include "../include/gz/zip_stream.hpp"
#include "../external/pthash/external/essentials/include/essentials.hpp"
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/cover.hpp"
#include "../include/builder/util.hpp"

using namespace sshash;

struct permute_data {
    permute_data() : num_runs_weights(0), num_sequences(0) {}
    uint64_t num_runs_weights;
    uint64_t num_sequences;
    std::vector<node> nodes;
};

void parse_file(std::istream& is, permute_data& data, build_configuration const& build_config) {
    std::string sequence;
    uint64_t k = build_config.k;
    uint64_t num_bases = 0;
    uint64_t num_kmers = 0;
    uint64_t seq_len = 0;

    uint64_t sum_of_weights = 0;
    uint64_t num_sequences_diff_weights = 0;  // num sequences whose kmers have different weights
    uint64_t num_sequences_all_mfw = 0;       // num sequences whose kmers have same weight == mfw

    data.num_runs_weights = 0;
    data.num_sequences = 0;

    /* count number of distinct weights */
    std::unordered_set<uint32_t> distinct_weights;

    auto parse_header = [&]() {
        if (sequence.empty()) return;

        /*
            Heder format:
            >[id] LN:i:[seq_len] ab:Z:[ab_seq]
            where [ab_seq] is a space-separated sequence of integer counters (the weights),
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

        bool kmers_have_all_mfw = true;
        bool kmers_have_different_weights = false;
        uint64_t front = constants::invalid_uint64;
        uint64_t back = constants::invalid_uint64;

        for (uint64_t j = 0, prev_weight = constants::invalid_uint64; j != seq_len - k + 1; ++j) {
            uint64_t weight = std::strtoull(sequence.data() + i, nullptr, 10);
            sum_of_weights += weight;
            i = sequence.find_first_of(' ', i) + 1;

            distinct_weights.insert(weight);

            /* set front and back */
            if (j == 0) front = weight;
            if (j == seq_len - k) back = weight;

            /* accumulate statistics */
            if (weight != constants::most_frequent_weight) kmers_have_all_mfw = false;
            if (j > 0 and weight != prev_weight) kmers_have_different_weights = true;

            /* count the number of runs */
            if (weight != prev_weight) data.num_runs_weights += 1;

            prev_weight = weight;
        }

        num_sequences_diff_weights += kmers_have_different_weights;
        num_sequences_all_mfw += kmers_have_all_mfw;

        constexpr bool sign = true;
        data.nodes.emplace_back(data.num_sequences, front, back, sign);
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

    assert(data.nodes.size() == data.num_sequences);
    assert(data.num_runs_weights >= data.num_sequences);

    std::cout << "read " << data.num_sequences << " sequences, " << num_bases << " bases, "
              << num_kmers << " kmers" << std::endl;
    std::cout << "found " << distinct_weights.size() << " distinct weights" << std::endl;
    std::cout << "sum_of_weights " << sum_of_weights << std::endl;
    std::cout << "num_sequences whose kmers have different weights: " << num_sequences_diff_weights
              << "/" << data.num_sequences << " ("
              << (num_sequences_diff_weights * 100.0) / data.num_sequences << "%)" << std::endl;
    std::cout << "num_sequences whose kmers all have the same weight != mfw: "
              << (data.num_sequences - num_sequences_diff_weights) << "/" << data.num_sequences
              << " ("
              << ((data.num_sequences - num_sequences_diff_weights) * 100.0) / data.num_sequences
              << "%)" << std::endl;
    std::cout << "num_sequences whose kmers all have the same weight == mfw: "
              << num_sequences_all_mfw << "/" << (data.num_sequences - num_sequences_diff_weights)
              << " ("
              << (num_sequences_all_mfw * 100.0) / (data.num_sequences - num_sequences_diff_weights)
              << "%)" << std::endl;

    /* computation of lower bound to the number of final weight runs */
    {
        struct info {
            /* how many times an end-point weight appears */
            uint32_t freq;

            /* if this flag is true, then the weight w always appears in nodes of the form
            (w,w) and its freq is always even. */
            bool all_equal;
        };

        std::unordered_map<uint32_t, info> endpoint_weights;

        for (auto const& node : data.nodes) {
            uint64_t front = node.front;
            uint64_t back = node.back;
            {
                auto it = endpoint_weights.find(front);
                if (it == endpoint_weights.cend()) {
                    endpoint_weights[front] = {1, true};
                } else {
                    (*it).second.freq += 1;
                }
            }
            {
                auto it = endpoint_weights.find(back);
                if (it == endpoint_weights.cend()) {
                    endpoint_weights[back] = {1, true};
                } else {
                    (*it).second.freq += 1;
                }
            }
            if (front != back) {
                endpoint_weights[back].all_equal = false;
                endpoint_weights[front].all_equal = false;
            }
        }

        uint64_t sum_of_freqs = 0;
        uint64_t num_endpoints = 0;
        for (auto const& p : endpoint_weights) {
            uint64_t freq = p.second.freq;
            sum_of_freqs += freq;

            /* special case: weight will appear as singleton, so count twice */
            if (p.second.all_equal) {
                num_endpoints += 2;
                continue;
            }

            /* if the excess of frequency is odd, then the weight will appear as end-point */
            if (freq % 2 == 1) num_endpoints += 1;
        }

        assert(sum_of_freqs == 2 * data.num_sequences);
        assert(num_endpoints % 2 == 0);
        (void)sum_of_freqs;
        uint64_t num_distinct_weights = distinct_weights.size();
        uint64_t num_switch_points = data.num_runs_weights - data.nodes.size();
        uint64_t num_walks = num_endpoints / 2;
        std::cout << "(estimated) num_walks = " << num_walks << std::endl;
        std::cout << "num_runs_weights = " << data.num_runs_weights << std::endl;
        std::cout << "num_switch_points = " << num_switch_points << std::endl;

        uint64_t R_lo = std::max<uint64_t>(num_distinct_weights, num_switch_points + 1);
        uint64_t R_hi = num_switch_points + num_walks;
        std::cout << "R_lo = " << R_lo << std::endl;
        std::cout << "R_hi = " << R_hi << std::endl;
    }
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

void reverse_header(std::string const& input, std::string& output, uint64_t k) {
    // Example header:
    // >2 LN:i:61 ab:Z:4 4 4 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1
    // Expected output:
    // >2 LN:i:61 ab:Z:1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 4 4 4

    /* skip validation of input string */

    uint64_t i = 0;
    i = input.find_first_of(' ', i);
    if (i == std::string::npos) throw parse_runtime_error();

    i += 6;  // skip ' LN:i':
    uint64_t j = input.find_first_of(' ', i);
    if (j == std::string::npos) throw parse_runtime_error();

    uint64_t seq_len = std::strtoull(input.data() + i, nullptr, 10);
    i = j + 6;  // skip ' ab:Z:'
    output.append(input.data(), input.data() + i);

    std::vector<uint32_t> weights;
    weights.reserve(seq_len - k + 1);
    for (uint64_t j = 0; j != seq_len - k + 1; ++j) {
        uint64_t weight = std::strtoull(input.data() + i, nullptr, 10);
        weights.push_back(weight);
        i = input.find_first_of(' ', i) + 1;
    }

    std::reverse(weights.begin(), weights.end());
    for (auto weight : weights) output.append(std::to_string(weight) + " ");
}

void permute_and_write(std::istream& is, std::string const& output_filename,
                       std::string const& tmp_dirname, pthash::compact_vector const& permutation,
                       pthash::bit_vector const& signs, uint64_t k) {
    constexpr uint64_t limit = 1 * essentials::GB;
    std::vector<std::pair<std::string, std::string>> buffer;  // (header, dna)

    std::string run_identifier =
        std::to_string(pthash::clock_type::now().time_since_epoch().count());

    std::string header_sequence;
    std::string header_sequence_r;

    std::string dna_sequence;
    std::string dna_sequence_rc;

    uint64_t num_sequences = permutation.size();
    uint64_t num_bases = 0;
    uint64_t bytes = 0;
    uint64_t num_files_to_merge = 0;

    auto get_tmp_output_filename = [&](uint64_t id) {
        return tmp_dirname + "/sshash.tmp.run" + run_identifier + "." + std::to_string(id);
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

        if (!signs[i]) {
            /* compute reverse complement of dna_sequence
               and reverse the weights in header_sequence */
            dna_sequence_rc.resize(dna_sequence.size());
            header_sequence_r.clear();
            util::compute_reverse_complement(dna_sequence.data(), dna_sequence_rc.data(),
                                             dna_sequence.size());
            reverse_header(header_sequence, header_sequence_r, k);
            dna_sequence.swap(dna_sequence_rc);
            header_sequence.swap(header_sequence_r);
        }

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
                       std::string const& tmp_dirname, pthash::compact_vector const& permutation,
                       pthash::bit_vector const& signs, uint64_t k) {
    std::ifstream is(input_filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + input_filename + "'");
    std::cout << "reading file '" << input_filename << "'..." << std::endl;
    if (util::ends_with(input_filename, ".gz")) {
        zip_istream zis(is);
        permute_and_write(zis, output_filename, tmp_dirname, permutation, signs, k);
    } else {
        permute_and_write(is, output_filename, tmp_dirname, permutation, signs, k);
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
               "\t- with also kmers' weights.\n"
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
        /* compute cover */
        cover c(data.num_sequences, data.num_runs_weights);
        assert(data.nodes.size() == data.num_sequences);
        c.compute(data.nodes);
        c.save(permutation_filename);
        std::vector<node>().swap(data.nodes);
    }

    pthash::compact_vector permutation;
    pthash::bit_vector signs;
    {
        std::ifstream is(permutation_filename.c_str());
        if (!is.good()) {
            throw std::runtime_error("error in opening the file '" + permutation_filename + "'");
        }
        pthash::compact_vector::builder cv_builder(
            data.num_sequences,
            data.num_sequences == 1 ? 1 : std::ceil(std::log2(data.num_sequences)));
        pthash::bit_vector_builder bv_builder(data.num_sequences);
        for (uint64_t i = 0; i != data.num_sequences; ++i) {
            uint64_t position = 0;
            bool sign = 0;
            is >> position;
            is >> sign;
            cv_builder.set(position, i);
            bv_builder.set(position, sign);
        }
        is.close();
        cv_builder.build(permutation);
        signs.build(&bv_builder);
    }

    /* permute and save to output file */
    permute_and_write(input_filename, output_filename, tmp_dirname, permutation, signs, k);
    std::remove(permutation_filename.c_str());

    return 0;
}
