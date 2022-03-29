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
    uint64_t num_distinct_abundances;
    uint64_t num_distinct_links;
    std::vector<node> nodes;
};

struct pair_hash {
    inline uint64_t operator()(std::pair<uint32_t, uint32_t> const& v) const {
        return static_cast<uint64_t>(v.first) * 31 + static_cast<uint64_t>(v.second);
    }
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
    data.num_distinct_abundances = 0;

    std::unordered_set<uint32_t> distinct_abundances;  // count number of distinct abundances

    std::unordered_map<std::pair<uint32_t, uint32_t>, uint32_t, pair_hash> distinct_links_freq;
    std::unordered_map<uint32_t, uint32_t> indegree;
    std::unordered_map<uint32_t, uint32_t> outdegree;

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

            distinct_abundances.insert(abundance);

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
        data.nodes.emplace_back(data.num_sequences, front, back, sign);

        if (front != back) {
            {
                auto it = indegree.find(back);
                if (it != indegree.cend()) {
                    (*it).second += 1;
                } else {
                    indegree[back] = 1;
                }
            }
            {
                auto it = outdegree.find(front);
                if (it != outdegree.cend()) {
                    (*it).second += 1;
                } else {
                    outdegree[front] = 1;
                }
            }
            {
                auto it = distinct_links_freq.find({front, back});
                if (it != distinct_links_freq.cend()) {
                    (*it).second += 1;
                } else {
                    distinct_links_freq[{front, back}] = 1;
                }
            }
        }
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
    assert(data.num_runs_abundances >= data.num_sequences);
    data.num_distinct_abundances = distinct_abundances.size();
    data.num_distinct_links = distinct_links_freq.size();

    {
        /* (front_abundance, num_seqs_with_front=front_abundance) */
        std::unordered_map<uint32_t, uint32_t> front_abundance_freqs;

        /* (back_abundance, num_seqs_with_back=back_abundance) */
        std::unordered_map<uint32_t, uint32_t> back_abundance_freqs;

        for (auto const& node : data.nodes) {
            auto it_front = front_abundance_freqs.find(node.front);
            if (it_front == front_abundance_freqs.cend()) {
                front_abundance_freqs[node.front] = 1;
            } else {
                (*it_front).second += 1;
            }
            auto it_back = back_abundance_freqs.find(node.back);
            if (it_back == back_abundance_freqs.cend()) {
                back_abundance_freqs[node.back] = 1;
            } else {
                (*it_back).second += 1;
            }
        }

        std::vector<uint32_t> distinct_abundances_vec;
        distinct_abundances_vec.reserve(distinct_abundances.size());
        for (auto x : distinct_abundances) distinct_abundances_vec.push_back(x);
        std::sort(distinct_abundances_vec.begin(), distinct_abundances_vec.end());

        uint64_t will_appear = 0;
        for (auto ab : distinct_abundances_vec) {
            int64_t freq_front = 0;
            int64_t freq_back = 0;
            auto it_front = front_abundance_freqs.find(ab);
            auto it_back = back_abundance_freqs.find(ab);
            if (it_front != front_abundance_freqs.cend()) { freq_front = (*it_front).second; }
            if (it_back != back_abundance_freqs.cend()) { freq_back = (*it_back).second; }

            /* abundance only appears for internal kmers */
            if (freq_front == 0 and freq_back == 0) continue;

            // std::cout << "ab:" << ab << " - freq_front:" << freq_front
            //           << " - freq_back:" << freq_back;

            uint64_t in_count = 0;
            uint64_t out_count = 0;
            auto it_in_count = indegree.find(ab);
            auto it_out_count = outdegree.find(ab);
            if (it_in_count != indegree.cend()) { in_count = (*it_in_count).second; }
            if (it_out_count != outdegree.cend()) { out_count = (*it_out_count).second; }

            // std::cout << " - indegree:" << in_count << " - outdegree:" << out_count;

            /* special case */
            if (in_count == 0 and out_count == 0) {
                // std::cout << " - **2**\n";
                will_appear += 2;  // will appear as singleton, so count twice
                continue;
            }

            uint64_t diff = std::abs(freq_front - freq_back);
            if (diff % 2 == 0) {
                // std::cout << " - 0\n";
            } else {
                // std::cout << " - 1\n";
                will_appear += 1;
            }
        }

        uint64_t num_runs_abundances_internal = data.num_runs_abundances - data.nodes.size();
        std::cout << "will_appear = " << will_appear / 2 << std::endl;
        std::cout << "computed lower bound = " << (num_runs_abundances_internal + will_appear / 2)
                  << std::endl;

        // uint64_t cuts = 0;
        // for (uint64_t i = 0; i != distinct_abundances_vec.size(); ++i) {
        //     uint64_t front = distinct_abundances_vec[i];
        //     if (i + 1 == distinct_abundances_vec.size()) break;
        //     uint64_t back = distinct_abundances_vec[i + 1];

        //     std::cout << "trying to glue [" << front << "," << back << "]; ";

        //     int64_t freq_forward = 0;  // signed
        //     int64_t freq_bckward = 0;  // signed
        //     auto it_forward = distinct_links_freq.find({front, back});
        //     auto it_bckward = distinct_links_freq.find({back, front});
        //     if (it_forward != distinct_links_freq.cend()) { freq_forward = (*it_forward).second;
        //     } if (it_bckward != distinct_links_freq.cend()) { freq_bckward =
        //     (*it_bckward).second; }

        //     std::cout << "freq of link [" << front << "," << back << "] = " << freq_forward << ";
        //     "; std::cout << "freq of link [" << back << "," << front << "] = " << freq_bckward <<
        //     "; ";

        //     /* if the difference in abs value is even, then the link
        //     disappears and we cannot glue anymore, i.e., we originate a cut */
        //     uint64_t diff = std::abs(freq_forward - freq_bckward);
        //     if (diff % 2 == 0) {
        //         // std::cout << "NOT GLUEABLE";
        //         // cuts += 1;

        //         int64_t freq_front = 0;  // signed
        //         int64_t freq_back = 0;   // signed
        //         auto it_front = front_abundance_freqs.find(front);
        //         auto it_back = back_abundance_freqs.find(front);
        //         if (it_front != front_abundance_freqs.cend()) { freq_front = (*it_front).second;
        //         } if (it_back != back_abundance_freqs.cend()) { freq_back = (*it_back).second; }

        //         bool front_appears = (std::abs(freq_front - freq_back) % 2 != 0);
        //         // or (freq_front == 1 and freq_back == 1);

        //         freq_front = 0;  // signed
        //         freq_back = 0;   // signed
        //         it_front = front_abundance_freqs.find(back);
        //         it_back = back_abundance_freqs.find(back);
        //         if (it_front != front_abundance_freqs.cend()) { freq_front = (*it_front).second;
        //         } if (it_back != back_abundance_freqs.cend()) { freq_back = (*it_back).second; }

        //         bool back_appears = (std::abs(freq_front - freq_back) % 2 != 0);
        //         // or (freq_front == 1 and freq_back == 1);

        //         if (front_appears and back_appears) {  // and link is absent
        //             std::cout << "NOT GLUEABLE";
        //             cuts += 1;
        //         }

        //     } else {
        //         std::cout << "glueable";
        //     }

        //     std::cout << '\n';
        // }
        // std::cout << "cuts = " << cuts << std::endl;
    }

    std::cout << "read " << data.num_sequences << " sequences, " << num_bases << " bases, "
              << num_kmers << " kmers" << std::endl;
    std::cout << "found " << data.num_distinct_abundances << " distint abundances" << std::endl;
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

    std::vector<uint32_t> abundances;
    abundances.reserve(seq_len - k + 1);
    for (uint64_t j = 0; j != seq_len - k + 1; ++j) {
        uint64_t abundance = std::strtoull(input.data() + i, nullptr, 10);
        abundances.push_back(abundance);
        i = input.find_first_of(' ', i) + 1;
    }

    std::reverse(abundances.begin(), abundances.end());
    for (auto abundance : abundances) { output.append(std::to_string(abundance) + " "); }
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

        if (!signs[i]) {
            /* compute reverse complement of dna_sequence
               and reverse the abundances in header_sequence */
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
        /* R_lo = max(c,r-p+1)
           where c = num. distinct abundances
                 r = num. runs in the abundances
                 p = num. sequences */
        uint64_t num_runs_abundances_internal = data.num_runs_abundances - data.nodes.size();
        std::cout << "num_runs_abundances = " << data.num_runs_abundances << std::endl;
        std::cout << "num_runs_abundances_internal = " << num_runs_abundances_internal << std::endl;
        std::cout << "num_distinct_links = " << data.num_distinct_links << std::endl;

        uint64_t R_lo = data.num_runs_abundances - data.nodes.size() + 1;
        uint64_t R_hi = R_lo;
        if (R_lo < data.num_distinct_abundances) {
            /* clearly we cannot have less than data.num_distinct_abundances runs */
            R_lo = data.num_distinct_abundances;
            // R_hi = data.num_runs_abundances_internal;
            R_hi = data.num_distinct_abundances + num_runs_abundances_internal -
                   data.num_distinct_links;
        }
        // else {
        // We assume there are less than 2^32 sequences and that
        // the largest abundance fits into a 32-bit uint.

        /* (front_abundance, num_seqs_with_front=front_abundance) */
        std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> front_abundance_freqs;

        /* (back_abundance, num_seqs_with_back=back_abundance) */
        std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> back_abundance_freqs;

        uint64_t offset = 0;
        for (auto const& node : data.nodes) {
            auto it_front = front_abundance_freqs.find(node.front);
            if (it_front == front_abundance_freqs.cend()) {
                front_abundance_freqs[node.front] = {1, offset};
            } else {
                (*it_front).second.first += 1;
                (*it_front).second.second = -1;
            }
            auto it_back = back_abundance_freqs.find(node.back);
            if (it_back == back_abundance_freqs.cend()) {
                back_abundance_freqs[node.back] = {1, offset};
            } else {
                (*it_back).second.first += 1;
                (*it_back).second.second = -1;
            }
            ++offset;
        }

        // std::cout << "front_abundance_freqs.size() = " << front_abundance_freqs.size() <<
        // std::endl; std::cout << "back_abundance_freqs.size() = " << back_abundance_freqs.size()
        // << std::endl;

        // for (auto const& p : back_abundance_freqs) {
        //     std::cout << "back-ab:" << p.first << " - freq:" << p.second.first << std::endl;
        // }
        // for (auto const& p : front_abundance_freqs) {
        //     std::cout << "front-ab:" << p.first << " - freq:" << p.second.first << std::endl;
        // }

        for (auto const& p : back_abundance_freqs) {
            /* if there is only one occurrence of this back abundance
               and it does not appear in any front abundance, then it cannot be matched */
            if (p.second.first == 1) {
                auto it = front_abundance_freqs.find(p.first);
                if (it == front_abundance_freqs.cend()) { R_hi += 1; }
                // else {
                //     if ((*it).second.first == 1 and p.second.second == (*it).second.second) {
                //         // uint64_t offset = p.second.second;
                //         // auto node = data.nodes[offset];
                //         // std::cout << "node " << node.id << ":[" << node.front << "," <<
                //         // node.back
                //         //           << "]" << std::endl;
                //         R_hi += 1;
                //     }
                // }
            }
            // else {
            //     auto it = front_abundance_freqs.find(p.first);
            //     if (it != front_abundance_freqs.cend()) {
            //         // abundance is present both as back and front

            //         // and appears an even number of times, will originate a walk
            //         //   with both sides unmatchable
            //         if ((*it).second.first % 2 == 0 and p.second.first % 2 == 0) { R_hi += 1; }
            //     }
            // }
        }

        // for (auto const& p : front_abundance_freqs) {
        //     if (p.second.first == 1 and
        //         back_abundance_freqs.find(p.first) == back_abundance_freqs.cend()) {
        //         R_hi += 1;
        //     }
        // }

        /*
        note that we do NOT do the same for front_abundance_freqs,
        namely the following piece of code:

            for (auto const& p : front_abundance_freqs) {
                if (p.second == 1 and
                    back_abundance_freqs.find(p.first) == back_abundance_freqs.cend()) {
                    R_hi += 1;
                }
            }

        because otherwise we will end up in potentially counting twice the number of cuts.
        */
        // }
        std::cout << "R_lo = " << R_lo << std::endl;
        std::cout << "R_hi = " << R_hi << std::endl;
    }

    {
        /* compute cover */
        cover c(data.num_sequences, data.num_runs_abundances);
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
