#pragma once

#include "../dictionary.hpp"
#include "../util.hpp"

#include "../gz/zip_stream.hpp"
#include "streaming_query_canonical_parsing.hpp"
#include "streaming_query_regular_parsing.hpp"

namespace sshash {

template <typename Query>
streaming_query_report streaming_query_from_fasta_file_multiline(dictionary const* dict,
                                                                 std::istream& is) {
    streaming_query_report report;
    buffered_lines_iterator it(is);
    std::string buffer;
    uint64_t k = dict->k();
    Query query(dict);
    query.start();
    while (!it.eof()) {
        bool empty_line_was_read = it.fill_buffer(buffer);
        for (uint64_t i = 0; i != buffer.size() - k + 1; ++i) {
            char const* kmer = buffer.data() + i;
            auto answer = query.lookup_advanced(kmer);
            report.num_kmers += 1;
            report.num_positive_kmers += answer.kmer_id != constants::invalid_uint64;
        }
        if (empty_line_was_read) { /* re-start the kmers' buffer */
            buffer.clear();
            query.start();
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
    return report;
}

template <typename Query>
streaming_query_report streaming_query_from_fasta_file(dictionary const* dict, std::istream& is) {
    streaming_query_report report;
    std::string line;
    uint64_t k = dict->k();
    Query query(dict);
    while (!is.eof()) {
        query.start();
        std::getline(is, line);  // skip first header line
        std::getline(is, line);
        if (line.size() < k) continue;
        for (uint64_t i = 0; i != line.size() - k + 1; ++i) {
            char const* kmer = line.data() + i;
            auto answer = query.lookup_advanced(kmer);
            report.num_kmers += 1;
            report.num_positive_kmers += answer.kmer_id != constants::invalid_uint64;
        }
    }
    report.num_searches = query.num_searches();
    report.num_extensions = query.num_extensions();
    return report;
}

template <typename Query>
streaming_query_report streaming_query_from_fastq_file(dictionary const* dict, std::istream& is) {
    streaming_query_report report;
    std::string line;
    uint64_t k = dict->k();
    Query query(dict);
    while (!is.eof()) {
        query.start();
        /* We assume the file is well-formed, i.e., there are exactly 4 lines per read. */
        std::getline(is, line);  // skip first header line
        std::getline(is, line);
        if (line.size() >= k) {
            for (uint64_t i = 0; i != line.size() - k + 1; ++i) {
                char const* kmer = line.data() + i;
                auto answer = query.lookup_advanced(kmer);
                report.num_kmers += 1;
                report.num_positive_kmers += answer.kmer_id != constants::invalid_uint64;
            }
        }
        std::getline(is, line);  // skip '+'
        std::getline(is, line);  // skip score
    }
    report.num_searches = query.num_searches();
    report.num_extensions = query.num_extensions();
    return report;
}

template <typename Query>
streaming_query_report streaming_query_from_fasta_file(dictionary const* dict, std::istream& is,
                                                       bool multiline) {
    if (multiline) return streaming_query_from_fasta_file_multiline<Query>(dict, is);
    return streaming_query_from_fasta_file<Query>(dict, is);
}

streaming_query_report dictionary::streaming_query_from_file(std::string const& filename,
                                                             bool multiline) const {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    streaming_query_report report;

    if (util::ends_with(filename, ".fa.gz") or util::ends_with(filename, ".fasta.gz")) {
        zip_istream zis(is);

        if (canonicalized()) {
            report = streaming_query_from_fasta_file<streaming_query_canonical_parsing>(this, zis,
                                                                                        multiline);
        } else {
            report = streaming_query_from_fasta_file<streaming_query_regular_parsing>(this, zis,
                                                                                      multiline);
        }

    } else if (util::ends_with(filename, ".fq.gz") or util::ends_with(filename, ".fastq.gz")) {
        if (multiline) {
            std::cout << "==> Warning: option 'multiline' is only valid for FASTA files, not FASTQ."
                      << std::endl;
        }
        zip_istream zis(is);

        if (canonicalized()) {
            report = streaming_query_from_fastq_file<streaming_query_canonical_parsing>(this, zis);
        } else {
            report = streaming_query_from_fastq_file<streaming_query_regular_parsing>(this, zis);
        }

    } else if (util::ends_with(filename, ".fa") or util::ends_with(filename, ".fasta")) {
        if (canonicalized()) {
            report = streaming_query_from_fasta_file<streaming_query_canonical_parsing>(this, is,
                                                                                        multiline);
        } else {
            report = streaming_query_from_fasta_file<streaming_query_regular_parsing>(this, is,
                                                                                      multiline);
        }

    } else if (util::ends_with(filename, ".fq") or util::ends_with(filename, ".fastq")) {
        if (multiline) {
            std::cout << "==> Warning: option 'multiline' is only valid for FASTA files, not FASTQ."
                      << std::endl;
        }
        if (canonicalized()) {
            report = streaming_query_from_fastq_file<streaming_query_canonical_parsing>(this, is);
        } else {
            report = streaming_query_from_fastq_file<streaming_query_regular_parsing>(this, is);
        }

    } else {
        std::cerr << "unsupported query file format" << std::endl;
    }

    is.close();
    return report;
}

}  // namespace sshash