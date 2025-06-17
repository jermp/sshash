#include "external/gz/zip_stream.hpp"

#include "include/dictionary.hpp"
#include "include/util.hpp"
#include "include/streaming_query.hpp"

namespace sshash {

template <class kmer_t, typename Query>
streaming_query_report streaming_query_from_fasta_file_multiline(dictionary<kmer_t> const* dict,
                                                                 std::istream& is) {
    streaming_query_report report;
    buffered_lines_iterator it(is);
    std::string buffer;
    const uint64_t k = dict->k();
    Query query(dict);
    query.reset();
    while (!it.eof()) {
        bool empty_line_was_read = it.fill_buffer(buffer);
        const uint64_t num_kmers = buffer.size() - k + 1;
        report.num_kmers += num_kmers;
        for (uint64_t i = 0; i != num_kmers; ++i) {
            char const* kmer = buffer.data() + i;
            query.lookup_advanced(kmer);
        }
        if (empty_line_was_read) { /* re-start the kmers' buffer */
            buffer.clear();
            query.reset();
        } else {
            if (buffer.size() > k - 1) {
                std::copy(buffer.data() + buffer.size() - k + 1, buffer.data() + buffer.size(),
                          buffer.data());
                buffer.resize(k - 1);
            }
        }
    }
    report.num_searches = query.num_searches();
    report.num_extensions = query.num_extensions();
    report.num_positive_kmers = query.num_positive_lookups();
    report.num_negative_kmers = query.num_negative_lookups();
    report.num_invalid_kmers = query.num_invalid_lookups();
    assert(report.num_kmers ==
           report.num_positive_kmers + report.num_negative_kmers + report.num_invalid_kmers);
    return report;
}

template <class kmer_t, typename Query>
streaming_query_report streaming_query_from_fasta_file(dictionary<kmer_t> const* dict,
                                                       std::istream& is) {
    streaming_query_report report;
    std::string line;
    const uint64_t k = dict->k();
    Query query(dict);
    while (!is.eof()) {
        query.reset();
        std::getline(is, line);  // skip first header line
        std::getline(is, line);
        if (line.size() < k) continue;
        const uint64_t num_kmers = line.size() - k + 1;
        report.num_kmers += num_kmers;
        for (uint64_t i = 0; i != num_kmers; ++i) {
            char const* kmer = line.data() + i;
            query.lookup_advanced(kmer);
        }
    }
    report.num_searches = query.num_searches();
    report.num_extensions = query.num_extensions();
    report.num_positive_kmers = query.num_positive_lookups();
    report.num_negative_kmers = query.num_negative_lookups();
    report.num_invalid_kmers = query.num_invalid_lookups();
    assert(report.num_kmers ==
           report.num_positive_kmers + report.num_negative_kmers + report.num_invalid_kmers);
    return report;
}

template <class kmer_t, typename Query>
streaming_query_report streaming_query_from_fastq_file(dictionary<kmer_t> const* dict,
                                                       std::istream& is) {
    streaming_query_report report;
    std::string line;
    const uint64_t k = dict->k();
    Query query(dict);
    while (!is.eof()) {
        query.reset();
        /* We assume the file is well-formed, i.e., there are exactly 4 lines per read. */
        std::getline(is, line);  // skip first header line
        std::getline(is, line);
        if (line.size() >= k) {
            const uint64_t num_kmers = line.size() - k + 1;
            report.num_kmers += num_kmers;
            for (uint64_t i = 0; i != line.size() - k + 1; ++i) {
                char const* kmer = line.data() + i;
                query.lookup_advanced(kmer);
            }
        }
        std::getline(is, line);  // skip '+'
        std::getline(is, line);  // skip score
    }
    report.num_searches = query.num_searches();
    report.num_extensions = query.num_extensions();
    report.num_positive_kmers = query.num_positive_lookups();
    report.num_negative_kmers = query.num_negative_lookups();
    report.num_invalid_kmers = query.num_invalid_lookups();
    assert(report.num_kmers ==
           report.num_positive_kmers + report.num_negative_kmers + report.num_invalid_kmers);
    return report;
}

template <class kmer_t, typename Query>
streaming_query_report streaming_query_from_fasta_file(dictionary<kmer_t> const* dict,
                                                       std::istream& is, bool multiline) {
    if (multiline) return streaming_query_from_fasta_file_multiline<kmer_t, Query>(dict, is);
    return streaming_query_from_fasta_file<kmer_t, Query>(dict, is);
}

template <class kmer_t>
streaming_query_report dictionary<kmer_t>::streaming_query_from_file(std::string const& filename,
                                                                     bool multiline) const {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    streaming_query_report report;

    if (util::ends_with(filename, ".fa.gz") or util::ends_with(filename, ".fasta.gz")) {
        zip_istream zis(is);

        if (canonical()) {
            report = streaming_query_from_fasta_file<kmer_t, streaming_query<kmer_t, true>>(
                this, zis, multiline);
        } else {
            report = streaming_query_from_fasta_file<kmer_t, streaming_query<kmer_t, false>>(
                this, zis, multiline);
        }
    } else if (util::ends_with(filename, ".fq.gz") or util::ends_with(filename, ".fastq.gz")) {
        if (multiline) {
            std::cout << "==> Warning: option 'multiline' is only valid for FASTA files, not FASTQ."
                      << std::endl;
        }
        zip_istream zis(is);
        if (canonical()) {
            report =
                streaming_query_from_fastq_file<kmer_t, streaming_query<kmer_t, true>>(this, zis);
        } else {
            report =
                streaming_query_from_fastq_file<kmer_t, streaming_query<kmer_t, false>>(this, zis);
        }
    } else if (util::ends_with(filename, ".fa") or util::ends_with(filename, ".fasta")) {
        if (canonical()) {
            report = streaming_query_from_fasta_file<kmer_t, streaming_query<kmer_t, true>>(
                this, is, multiline);
        } else {
            report = streaming_query_from_fasta_file<kmer_t, streaming_query<kmer_t, false>>(
                this, is, multiline);
        }
    } else if (util::ends_with(filename, ".fq") or util::ends_with(filename, ".fastq")) {
        if (multiline) {
            std::cout << "==> Warning: option 'multiline' is only valid for FASTA files, not FASTQ."
                      << std::endl;
        }
        if (canonical()) {
            report =
                streaming_query_from_fastq_file<kmer_t, streaming_query<kmer_t, true>>(this, is);
        } else {
            report =
                streaming_query_from_fastq_file<kmer_t, streaming_query<kmer_t, false>>(this, is);
        }
    } else {
        std::cerr << "unsupported query file format" << std::endl;
    }

    is.close();
    return report;
}

}  // namespace sshash