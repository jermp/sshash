[![Build](https://github.com/jermp/sshash/actions/workflows/build.yml/badge.svg)](https://github.com/jermp/sshash/actions/workflows/build.yml)
[![CodeQL](https://github.com/jermp/sshash/actions/workflows/codeql.yml/badge.svg)](https://github.com/jermp/sshash/actions/workflows/codeql.yml)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sshash/badges/platforms.svg)](https://anaconda.org/bioconda/sshash)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/sshash/README.html)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7772316.svg)](https://doi.org/10.5281/zenodo.7772316)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7239205.svg)](https://doi.org/10.5281/zenodo.7239205)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17582116.svg)](https://doi.org/10.5281/zenodo.17582116)

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="img/sshash_on_dark.png">
  <img src="img/sshash.png" width="350" alt="Logo">
</picture>

**SSHash** is a compressed dictionary data structure for k-mers
(strings of length k over the DNA alphabet {A,C,G,T}), based on **S**parse and **S**kew **Hash**ing.

**NEWS:** A Rust port of SSHash is available [here](https://github.com/COMBINE-lab/sshash-rs)!

The data structure is described in the following papers (most recent first):

* [Optimizing sparse and skew hashing: faster k-mer dictionaries](https://www.biorxiv.org/content/10.64898/2026.01.21.700884v1) [1]
* [On weighted k-mer dictionaries](https://almob.biomedcentral.com/articles/10.1186/s13015-023-00226-2) [2,3]
* [Sparse and skew hashing of k-mers](https://doi.org/10.1093/bioinformatics/btac245) [4]

**Please, cite these papers if you use SSHash.**

For a dictionary of n k-mers,
two basic queries are supported:

- i = **Lookup**(x), where i is in [0,n) if the k-mer x is found in the dictionary or i = -1 otherwise;
- x = **Access**(i), where x is the k-mer associated to the identifier i.

If also the weights of the k-mers (their frequency counts) are stored in the dictionary, then the dictionary is said to be *weighted* and it also supports:

- w = **Weight**(i), where i is a given k-mer identifier and w is the weight of the k-mer.

Other supported queries are:

- **Membership Queries**: determine if a given k-mer is present in the dictionary or not.
- **Streaming Queries**: stream through all k-mers of a given DNA file
(.fasta or .fastq formats) to determine their membership to the dictionary.
- **Navigational Queries**: given a k-mer x[1..k] determine if x[2..k]+c is present (forward neighbourhood) and if c+x[1..k-1] is present (backward neighbourhood), for c in {A,C,G,T} ('+' here means string concatenation).
SSHash internally stores a set of strings, each associated to a distinct identifier.
If a string identifier is specified for a navigational query (rather than a k-mer), then the backward neighbourhood of the first k-mer and the forward neighbourhood of the last k-mer in the string are returned.

If you are interested in a **membership-only** version of SSHash, have a look at [SSHash-Lite](https://github.com/jermp/sshash-lite). It also works for input files with duplicate k-mers (e.g., [matchtigs](https://github.com/algbio/matchtigs) [5]). For a query sequence S and a given coverage threshold E in [0,1], the sequence is considered to be present in the dictionary if at least E*(|S|-k+1) of the k-mers of S are positive.

**NOTE**: It is assumed that two k-mers being the *reverse complement* of each other are the same.

#### Table of contents
* [Compiling the Code](#compiling-the-code)
* [Dependencies](#dependencies)
* [Tools and Usage](#tools-and-usage)
* [Examples](#Examples)
* [Input Files](#input-files)
* [Create a New Release](#create-a-new-release)
* [Benchmarks](#benchmarks)
* [References](#references)

Compiling the Code
------------------

The code is tested on Linux with `gcc` and on Mac with `clang`.
To build the code, [`CMake`](https://cmake.org/) is required.

Clone the repository with

    git clone --recursive https://github.com/jermp/sshash.git

If you have cloned the repository **without** `--recursive`, be sure you pull the dependencies with the following command before
compiling:

    git submodule update --init --recursive

To compile the code for a release environment (see file `CMakeLists.txt` for the used compilation flags), it is sufficient to do the following:

    mkdir build
    cd build
    cmake ..
    make -j

**NOTE**: For best performance on `x86` architectures, the option `-D SSHASH_USE_ARCH_NATIVE` can be specified as well.

For a testing environment, use the following instead:

    mkdir debug_build
    cd debug_build
    cmake .. -D CMAKE_BUILD_TYPE=Debug -D SSHASH_USE_SANITIZERS=On
    make -j

### Encoding of Nucleotides

SSHash uses by default the following 2-bit encoding of nucleotides.

	 A     65     01000.00.1 -> 00
	 C     67     01000.01.1 -> 01
	 G     71     01000.11.1 -> 11
	 T     84     01010.10.0 -> 10

	 a     97     01100.00.1 -> 00
	 c     99     01100.01.1 -> 01
	 g    103     01100.11.1 -> 11
	 t    116     01110.10.0 -> 10

If you want to use the "traditional" encoding

	 A     65     01000001 -> 00
	 C     67     01000011 -> 01
	 G     71     01000111 -> 10
	 T     84     01010100 -> 11

	 a     97     01100001 -> 00
	 c     99     01100011 -> 01
	 g    103     01100111 -> 10
	 t    116     01110100 -> 11

for compatibility issues with other software, then
compile SSHash with the flag `-DSSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING=On`.

### K-mer Length

By default, SSHash uses a maximum k-mer length of 31.
If you want to support k-mer lengths up to (and including) 63,
compile the library with the flag `-DSSHASH_USE_MAX_KMER_LENGTH_63=On`.

Dependencies
------------

The repository has minimal dependencies: it only uses the [PTHash](https://github.com/jermp/pthash) library (for minimal perfect hashing), and `zlib` to read gzip-compressed streams.

To automatically pull the PTHash dependency, just clone the repo with
`--recursive` as explained in [Compiling the Code](#compiling-the-code).

If you do not have `zlib` installed, you can do

    sudo apt-get install zlib1g

if you are on Linux/Ubuntu, or

    brew install zlib

if you have a Mac.

Tools and Usage
---------------

There is one executable called `sshash` after the compilation, which can be used to run a tool.
Run `./sshash` as follows to see a list of available tools.

For large-scale indexing, it could be necessary to increase the number of file descriptors that can be opened simultaneously:

	ulimit -n 2048

Examples
--------

For the examples, we are going to use some collections
of *stitched unitigs* from the directory `data/unitigs_stitched`.

**Important note:** The value of k used during the formation of the unitigs
is indicated in the name of each file and the dictionaries
**must** be built with that value as well to ensure correctness.

For example, `data/unitigs_stitched/ecoli4_k31_ust.fa.gz` indicates the value k = 31, whereas `data/unitigs_stitched/se.ust.k63.fa.gz` indicates the value k = 63.

For all the examples below, we are going to use k = 31.

(The directory `data/unitigs_stitched/with_weights` contains some files with k-mers' weights too.)

In the section [Input Files](#input-files), we explain how
such collections of stitched unitigs can be obtained from raw FASTA files.

### Example 1

    ./sshash build -i ../data/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz -k 31 -m 13 --check -o salmonella_enterica.sshash

This example builds a dictionary for the k-mers read from the file `../data/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz`,
with k = 31 and m = 13. It also check the correctness of the dictionary (`--check` option) and serializes the index on disk to the file `salmonella_enterica.sshash`.

To run a performance benchmark after construction of the index,
use:

    ./sshash bench -i salmonella_enterica.sshash

To also store the weights, use the option `--weighted`:

    ./sshash build -i ../data/unitigs_stitched/with_weights/salmonella_enterica.ust.k31.fa.gz -k 31 -m 13 --weighted --check --verbose

### Example 2

    ./sshash build -i ../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz -k 31 -m 15 -o salmonella_100.sshash

This example builds a dictionary from the input file `../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz` (a pangenome consisting in 100 genomes of *Salmonella Enterica*), with k = 31, m = 15, and l = 2. It also serializes the index on disk to the file `salmonella_100.sshash`.

To perform some streaming membership queries, use:

    ./sshash query -i salmonella_100.sshash -q ../data/queries/SRR5833294.10K.fastq.gz

if your queries are meant to be read from a FASTQ file, or

    ./sshash query -i salmonella_100.sshash -q ../data/queries/salmonella_enterica.fasta.gz --multiline

if your queries are to be read from a (multi-line) FASTA file.

### Example 3

    ./sshash build -i ../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz -k 31 -m 13 --canonical -o salmonella_100.canon.sshash

This example builds a dictionary from the input file `../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz` (same used in Example 2), with k = 31, m = 13, and with the canonical parsing modality (option `--canonical`). The dictionary is serialized on disk to the file `salmonella_100.canon.sshash`.

The "canonical" version of the dictionary offers more speed for only a little space increase, especially under low-hit workloads -- when the majority of k-mers are not found in the dictionary. (For all details, refer to the paper.)

Below a comparison between the dictionary built in Example 2 (not canonical)
and the one just built (Example 3, canonical).

    ./sshash query -i salmonella_100.sshash -q ../data/queries/SRR5833294.10K.fastq.gz

    ./sshash query -i salmonella_100.canon.sshash -q ../data/queries/SRR5833294.10K.fastq.gz

Both queries should originate the following report (reported here for reference):

    ==== query report:
    num_kmers = 460000
    num_positive_kmers = 46 (0.01%)
    num_searches = 42/46 (91.3043%)
    num_extensions = 4/46 (8.69565%)

The canonical dictionary can be twice as fast as the regular dictionary
for low-hit workloads, even on this tiny example, for only +0.3 bits/k-mer.

### Example 4

    ./sshash permute -i ../data/unitigs_stitched/with_weights/ecoli_sakai.ust.k31.fa.gz -k 31 -o ecoli_sakai.permuted.fa

This command re-orders (and possibly reverse-complement) the strings in the collection as to *minimize* the number of runs in the weights and, hence, optimize the encoding of the weights.
The result is saved to the file `ecoli_sakai.permuted.fa`.

In this example for the E.Coli collection (Sakai strain) we reduce the number of runs in the weights from 5820 to 3723.

Then use the `build` command as usual to build the permuted collection:

    ./sshash build -i ecoli_sakai.permuted.fa -k 31 -m 13 --weighted --verbose

The index built on the permuted collection
optimizes the storage space for the weights which results in a 15.1X better space than the empirical entropy of the weights.

For reference, the index built on the original collection:

    ./sshash build -i ../data/unitigs_stitched/with_weights/ecoli_sakai.ust.k31.fa.gz -k 31 -m 13 --weighted --verbose

already achieves a 12.4X better space than the empirical entropy.

Input Files
-----------

SSHash is meant to index k-mers from collections that **do not contain duplicates
nor invalid k-mers** (strings containing symbols different from {A,C,G,T}).
These collections can be obtained, for example, by extracting the maximal unitigs of a de Bruijn graph, or eulertigs, using the [GGCAT](https://github.com/algbio/ggcat) algorithm.

**NOTE**: Input files are expected to have **one DNA sequence per line**. If a sequence spans multiple lines (e.g., multi-fasta), the lines should be concatenated before indexing.

#### Datasets

The script `scripts/download_and_preprocess_datasets.sh` of [this release](https://github.com/jermp/sshash/releases/tag/v3.0.0)
contains all the needed steps to download and pre-process
the datasets that we used in [1].

For the experiments in [2] and [3], we used the datasets available at [https://doi.org/10.5281/zenodo.7772316](https://doi.org/10.5281/zenodo.7772316).

For the latest benchmarks maintained in [this other repository](https://github.com/jermp/kmer_sets_benchmark)
we used the datasets described at [https://zenodo.org/records/17582116](https://zenodo.org/records/17582116).

#### Weights

Using the option `-all-abundance-counts` of [BCALM2](https://github.com/GATB/bcalm), it is possible to also include the abundance counts of the k-mers in the BCALM2 output. Then, use the option `-a 1` of [UST](https://github.com/jermp/UST) to include such counts in the stitched unitigs.

Create a New Release
--------------------

It is recommended to create a new release with the script `script/create_release.sh` which
**also includes the source code for the dependencies** in `external`
(this is not done by GitHub).

To create a new release, run the following command *from the parent directory*:

    bash script/create_release.sh --format zip [RELEASE-NAME]

for example

    bash script/create_release.sh --format zip v4.0.0.tar.gz

**Note 1**: The sha256 hash code printed at the end is needed for distribution via Bioconda.

**Note 2**: Avoid dashes in the name of the release because Bioconda does not like them.

Benchmarks
----------

The directory [`benchmarks`](/benchmarks) includes some performance benchmarks.

References
----------

* [1] Giulio Ermanno Pibiri and Rob Patro. [Optimizing sparse and skew hashing: faster k-mer dictionaries](https://www.biorxiv.org/content/10.64898/2026.01.21.700884v1). BioRxiv. 2026.
* [2] Giulio Ermanno Pibiri. [On weighted k-mer dictionaries](https://almob.biomedcentral.com/articles/10.1186/s13015-023-00226-2). Algorithms for Molecular Biology (ALGOMB). 2023.
* [3] Giulio Ermanno Pibiri. [On weighted k-mer dictionaries](https://drops.dagstuhl.de/opus/volltexte/2022/17043/). International Workshop on Algorithms in Bioinformatics (WABI). 2022.
* [4] Giulio Ermanno Pibiri. [Sparse and skew hashing of k-mers](https://doi.org/10.1093/bioinformatics/btac245). Bioinformatics. 2022.
* [5] Schmidt, S., Khan, S., Alanko, J., Pibiri, G. E., and Tomescu, A. I. [Matchtigs: minimum plain text representation of k-mer sets](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02968-z). Genome Biology 24, 136. 2023.
