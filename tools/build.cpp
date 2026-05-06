using namespace sshash;

int build(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);

    /* Required arguments. */
    parser.add(
        "input_filename",
        "Must be a FASTA file (.fa/fasta extension) or a CUTTLEFISH file (.cfseg extension, "
        "produced by CUTTLEFISH with option -f 3) "
        "compressed with gzip (.gz) or not:\n"
        "\t- without duplicate nor invalid kmers\n"
        "\t- one DNA sequence per line.\n"
        "\tFor example, it could be the de Bruijn graph topology output by BCALM or CUTTLEFISH.",
        "-i", true);
    parser.add("k", "K-mer length (must be <= " + std::to_string(default_kmer_t::max_k) + ").",
               "-k", true);
    parser.add("m", "Minimizer length (must be < k).", "-m", true);

    /* Optional arguments. */
    parser.add("seed",
               "Seed for construction (default is " + std::to_string(constants::seed) + ").", "-s",
               false);
    parser.add("t", "Number of threads (default is 1).", "-t", false);
    parser.add("lambda",
               "A (floating point) constant that trades construction speed for space effectiveness "
               "of minimal perfect hashing. "
               "A reasonable value lies between 3.0 and 10.0 (default is " +
                   std::to_string(constants::lambda) + ").",
               "-a", false);
    parser.add("output_filename", "Output file name where the data structure will be serialized.",
               "-o", false);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory. Default is directory '" +
            constants::default_tmp_dirname + "'.",
        "-d", false);
    parser.add("RAM",
               "RAM limit in GiB. Default value is " +
                   std::to_string(constants::default_ram_limit_in_GiB) + " GiB.",
               "-g", false);
    parser.add("canonical",
               "This option results in a trade-off between index space and lookup time.",
               "--canonical", false, true);
    parser.add("weighted", "Also store the weights in compressed format.", "--weighted", false,
               true);
    parser.add("check", "Check correctness after construction.", "--check", false, true);
    parser.add("verbose", "Verbose output during construction.", "--verbose", false, true);

    if (!parser.parse()) return 0;

    auto input_filename = parser.get<std::string>("input_filename");
    auto k = parser.get<uint64_t>("k");
    auto m = parser.get<uint64_t>("m");

    build_configuration build_config;
    build_config.k = k;
    build_config.m = m;

    if (parser.parsed("seed")) build_config.seed = parser.get<uint64_t>("seed");
    if (parser.parsed("lambda")) build_config.lambda = parser.get<double>("lambda");
    build_config.canonical = parser.get<bool>("canonical");
    build_config.weighted = parser.get<bool>("weighted");
    build_config.verbose = parser.get<bool>("verbose");
    if (parser.parsed("tmp_dirname")) {
        build_config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(build_config.tmp_dirname);
    }
    if (parser.get<uint64_t>("RAM")) {
        build_config.ram_limit_in_GiB = parser.get<uint64_t>("RAM");
    }
    if (parser.parsed("t")) build_config.num_threads = parser.get<uint64_t>("t");

    // build_config.print();

    bool check = parser.get<bool>("check");
    bool has_output = parser.parsed("output_filename");

    /* Always build via the streaming-save path: peak RAM is bounded by
       the build phase only. If the caller didn't pass -o, write to a
       tmp file in `tmp_dirname` and delete it after the build (or after
       the --check verification). */
    std::string output_filename;
    if (has_output) {
        output_filename = parser.get<std::string>("output_filename");
    } else {
        std::stringstream ss;
        ss << build_config.tmp_dirname << "/sshash.tmp.run_"
           << pthash::clock_type::now().time_since_epoch().count() << ".index.bin";
        output_filename = ss.str();
    }

    {
        dictionary_type dict;
        essentials::logger("building data structure...");
        dict.build(input_filename, build_config, output_filename);
        essentials::logger("DONE");
    }

    if (check) {
        dictionary_type dict;
        open_dictionary(dict, output_filename, /*mmap=*/true, build_config.verbose);
        check_correctness_lookup_access(dict, input_filename);
        check_correctness_navigational_kmer_query(dict, input_filename);
        check_correctness_navigational_string_query(dict);
        if (build_config.weighted) check_correctness_weights(dict, input_filename);
        check_correctness_kmer_iterator(dict);
        check_correctness_string_iterator(dict);
    }

    if (!has_output) std::remove(output_filename.c_str());

    return 0;
}
