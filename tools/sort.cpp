#pragma once

#include "external/pthash/include/pthash.hpp"
#include "external/pthash/external/bits/external/essentials/include/essentials.hpp"
// #include "external/gz/zip_stream.hpp"
// #include "external/gz/zip_stream.cpp"
#include "include/builder/util.hpp"

namespace sshash {

void sort(std::istream& is, std::string const& output_filename,
          std::string const& tmp_dirname)  //
{
    constexpr uint64_t limit_in_bytes = 4 * essentials::GB;
    constexpr uint64_t limit_in_num_sequences = 1000000;

    std::vector<std::string> buffer;
    buffer.reserve(limit_in_num_sequences);

    std::string run_identifier =
        std::to_string(pthash::clock_type::now().time_since_epoch().count());

    std::string sequence;
    uint64_t num_sequences = 0;
    uint64_t num_bases = 0;
    uint64_t bytes = 0;
    uint64_t num_files_to_merge = 0;

    auto get_tmp_output_filename = [&](uint64_t id) {
        return tmp_dirname + "/sshash.tmp.run" + run_identifier + "." + std::to_string(id);
    };

    auto sort_and_flush = [&]() {
        if (buffer.empty()) return;

        std::cout << "sorting buffer..." << std::endl;
        std::sort(buffer.begin(), buffer.end(),
                  [&](auto const& x, auto const& y) { return x.length() < y.length(); });

        auto tmp_output_filename = get_tmp_output_filename(num_files_to_merge);
        std::cout << "saving to file '" << tmp_output_filename << "'..." << std::endl;
        std::ofstream out(tmp_output_filename.c_str());
        if (!out.is_open()) throw std::runtime_error("cannot open file");
        for (auto const& s : buffer) out << s << '\n';
        out.close();

        buffer.clear();
        bytes = 0;
        num_files_to_merge += 1;
    };

    while (true) {
        std::getline(is, sequence);  // header sequence
        std::getline(is, sequence);  // DNA sequence
        if (is.eof()) break;

        if (++num_sequences % 100000 == 0) {
            std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases"
                      << std::endl;
        }

        const uint64_t n = sequence.size();
        num_bases += n;

        uint64_t seq_bytes = n + 16;  // overhead
        if (bytes + seq_bytes > limit_in_bytes or buffer.size() == limit_in_num_sequences) {
            sort_and_flush();
        }

        bytes += seq_bytes;
        buffer.emplace_back(sequence);
        num_bases += n;
    }
    sort_and_flush();

    std::cout << "read " << num_sequences << " sequences, " << num_bases << " bases" << std::endl;

    /* merge */
    {
        assert(num_files_to_merge > 0);
        std::cout << "files to merge = " << num_files_to_merge << std::endl;

        struct lines_iterator {
            lines_iterator(uint8_t const* begin, uint8_t const* end) : m_begin(begin), m_end(end) {
                advance_to_next();
            }

            void advance_to_next() { m_sequence = read_line(); }

            std::string const& sequence() const { return m_sequence; }
            bool has_next() const { return !m_sequence.empty(); }

        private:
            uint8_t const* m_begin;
            uint8_t const* m_end;
            std::string m_sequence;

            std::string read_line() {
                uint8_t const* begin = m_begin;
                while (m_begin != m_end and *m_begin++ != '\n')
                    ;
                if (begin == m_begin) return std::string("");
                return std::string(reinterpret_cast<const char*>(begin), m_begin - begin - 1);
            }
        };

        std::vector<lines_iterator> iterators;
        std::vector<uint32_t> idx_heap;
        iterators.reserve(num_files_to_merge);
        idx_heap.reserve(num_files_to_merge);
        std::vector<mm::file_source<uint8_t>> mm_files(num_files_to_merge);

        auto heap_idx_comparator = [&](uint32_t i, uint32_t j) {
            return iterators[i].sequence().length() > iterators[j].sequence().length();
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
            out << '>' << '\n' << it.sequence() << '\n';
            num_written_sequences += 1;
            if (num_written_sequences % 1000000 == 0) {
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

void sort(std::string const& input_filename, std::string const& output_filename,
          std::string const& tmp_dirname)  //
{
    std::ifstream is(input_filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + input_filename + "'");
    std::cout << "reading file '" << input_filename << "'..." << std::endl;
    if (util::ends_with(input_filename, ".gz")) {
        zip_istream zis(is);
        sort(zis, output_filename, tmp_dirname);
    } else {
        sort(is, output_filename, tmp_dirname);
    }
    is.close();
}

template <class kmer_t>
int sort(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);

    /* Required arguments. */
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not.",
               "-i", true);
    parser.add("output_filename", "Output filename.", "-o", true);

    /* Optional arguments. */
    parser.add("tmp_dirname",
               "Temporary directory used for merging in external memory. Default "
               "is directory '" +
                   constants::default_tmp_dirname + "'.",
               "-d", false);

    if (!parser.parse()) return 0;

    auto input_filename = parser.get<std::string>("input_filename");

    std::string output_filename = parser.get<std::string>("output_filename");
    std::string tmp_dirname = constants::default_tmp_dirname;
    if (parser.parsed("tmp_dirname")) {
        tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(tmp_dirname);
    }

    sort(input_filename, output_filename, tmp_dirname);

    return 0;
}

}  // namespace sshash