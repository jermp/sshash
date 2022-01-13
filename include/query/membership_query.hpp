#pragma once

#include "../dictionary.hpp"
#include "../util.hpp"
#include "../gz/zip_stream.cpp"
#include "membership_query_canonical_parsing.cpp"
#include "membership_query_regular_parsing.cpp"

namespace sshash {

template <typename Query>
dictionary::membership_query_result membership_query_from_fasta_file_multiline(
    dictionary const* dict, std::istream& is) {
    dictionary::membership_query_result result;
    buffered_lines_iterator it(is);
    std::string buffer;
    uint64_t k = dict->k();
    Query q(dict);
    q.start();
    while (!it.eof()) {
        bool empty_line_was_read = it.fill_buffer(buffer);
        for (uint64_t i = 0; i != buffer.size() - k + 1; ++i) {
            char const* kmer = buffer.data() + i;
            auto answer = q.is_member(kmer);
            result.num_kmers += 1;
            result.num_valid_kmers += answer.is_valid;
            result.num_positive_kmers += answer.is_member;
        }
        if (empty_line_was_read) { /* re-start the kmers' buffer */
            buffer.clear();
            q.start();
        } else {
            if (buffer.size() > k - 1) {
                std::copy(buffer.data() + buffer.size() - k + 1, buffer.data() + buffer.size(),
                          buffer.data());
                buffer.resize(k - 1);
            }
        }
    }
    result.num_searches = q.num_searches;
    result.num_extensions = q.num_extensions;
    return result;
}

template <typename Query>
dictionary::membership_query_result membership_query_from_fasta_file(dictionary const* dict,
                                                                     std::istream& is) {
    dictionary::membership_query_result result;
    std::string line;
    uint64_t k = dict->k();
    Query q(dict);
    while (!is.eof()) {
        q.start();
        std::getline(is, line);  // skip first header line
        std::getline(is, line);
        if (line.size() < k) continue;
        for (uint64_t i = 0; i != line.size() - k + 1; ++i) {
            char const* kmer = line.data() + i;
            auto answer = q.is_member(kmer);
            result.num_kmers += 1;
            result.num_valid_kmers += answer.is_valid;
            result.num_positive_kmers += answer.is_member;
        }
    }
    result.num_searches = q.num_searches;
    result.num_extensions = q.num_extensions;
    return result;
}

template <typename Query>
dictionary::membership_query_result membership_query_from_fastq_file(dictionary const* dict,
                                                                     std::istream& is) {
    dictionary::membership_query_result result;
    std::string line;
    uint64_t k = dict->k();
    Query q(dict);
    while (!is.eof()) {
        q.start();
        /* We assume the file is well-formed, i.e., there are exactly 4 lines per read. */
        std::getline(is, line);  // skip first header line
        std::getline(is, line);
        if (line.size() >= k) {
            for (uint64_t i = 0; i != line.size() - k + 1; ++i) {
                char const* kmer = line.data() + i;
                auto answer = q.is_member(kmer);
                result.num_kmers += 1;
                result.num_valid_kmers += answer.is_valid;
                result.num_positive_kmers += answer.is_member;
            }
        }
        std::getline(is, line);  // skip '+'
        std::getline(is, line);  // skip score
    }
    result.num_searches = q.num_searches;
    result.num_extensions = q.num_extensions;
    return result;
}

template <typename Query>
dictionary::membership_query_result membership_query_from_fasta_file(dictionary const* dict,
                                                                     std::istream& is,
                                                                     bool multiline) {
    if (multiline) return membership_query_from_fasta_file_multiline<Query>(dict, is);
    return membership_query_from_fasta_file<Query>(dict, is);
}

dictionary::membership_query_result dictionary::membership_query_from_file(
    std::string const& filename, bool multiline) const {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    dictionary::membership_query_result result;

    if (util::ends_with(filename, ".fa.gz") or util::ends_with(filename, ".fasta.gz")) {
        zip_istream zis(is);

        if (canonicalized()) {
            result = membership_query_from_fasta_file<membership_query_canonical_parsing>(
                this, zis, multiline);
        } else {
            result = membership_query_from_fasta_file<membership_query_regular_parsing>(this, zis,
                                                                                        multiline);
        }

    } else if (util::ends_with(filename, ".fq.gz") or util::ends_with(filename, ".fastq.gz")) {
        if (multiline) {
            std::cout << "==> Warning: option 'multiline' is only valid for FASTA files, not FASTQ."
                      << std::endl;
        }
        zip_istream zis(is);

        if (canonicalized()) {
            result =
                membership_query_from_fastq_file<membership_query_canonical_parsing>(this, zis);
        } else {
            result = membership_query_from_fastq_file<membership_query_regular_parsing>(this, zis);
        }

    } else if (util::ends_with(filename, ".fa") or util::ends_with(filename, ".fasta")) {
        if (canonicalized()) {
            result = membership_query_from_fasta_file<membership_query_canonical_parsing>(
                this, is, multiline);
        } else {
            result = membership_query_from_fasta_file<membership_query_regular_parsing>(this, is,
                                                                                        multiline);
        }

    } else if (util::ends_with(filename, ".fq") or util::ends_with(filename, ".fastq")) {
        if (multiline) {
            std::cout << "==> Warning: option 'multiline' is only valid for FASTA files, not FASTQ."
                      << std::endl;
        }
        if (canonicalized()) {
            result = membership_query_from_fastq_file<membership_query_canonical_parsing>(this, is);
        } else {
            result = membership_query_from_fastq_file<membership_query_regular_parsing>(this, is);
        }
    }

    is.close();
    return result;
}

}  // namespace sshash