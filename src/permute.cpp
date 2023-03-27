#include "../include/cover/parse_file.hpp"
#include "../include/cover/cover.hpp"

using namespace sshash;

int permute(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);

    /* Required arguments. */
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip "
               "(.gz) or not:\n"
               "\t- without duplicate nor invalid kmers\n"
               "\t- one DNA sequence per line\n"
               "\t- with also kmers' weights.\n"
               "\tFor example, it could be the de Bruijn graph topology output "
               "by BCALM.",
               "-i", true);
    parser.add("k", "K-mer length (must be <= " + std::to_string(constants::max_k) + ").", "-k",
               true);

    /* Optional arguments. */
    parser.add("output_filename", "Output file where the permuted collection will be written.",
               "-o", false);
    parser.add("tmp_dirname",
               "Temporary directory used for merging in external memory. Default "
               "is directory '" +
                   constants::default_tmp_dirname + "'.",
               "-d", false);

    if (!parser.parse()) return 1;

    auto input_filename = parser.get<std::string>("input_filename");
    auto k = parser.get<uint64_t>("k");

    /* validate k */
    if (k == 0) {
        std::cerr << "k must be > 0" << std::endl;
        return 1;
    }
    if (k > constants::max_k) {
        std::cerr << "k must be less <= " + std::to_string(constants::max_k) +
                         " but got k = " + std::to_string(k)
                  << std::endl;
        return 1;
    }

    build_configuration build_config;
    build_config.k = k;

    std::string output_filename = input_filename + ".permuted";
    if (parser.parsed("output_filename")) {
        output_filename = parser.get<std::string>("output_filename");
    }

    std::string tmp_dirname = constants::default_tmp_dirname;
    if (parser.parsed("tmp_dirname")) tmp_dirname = parser.get<std::string>("tmp_dirname");

    std::string permutation_filename = tmp_dirname + "/tmp.permutation";

    auto data = parse_weighted_file(input_filename, build_config);
    assert(data.nodes.size() == data.num_sequences);

    {
        /* compute cover */
        cover c(data.num_sequences, data.num_runs_weights);
        c.swap(data.nodes);
        c.compute();
        c.save(permutation_filename);
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