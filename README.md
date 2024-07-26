[![CodeQL](https://github.com/jermp/sshash/actions/workflows/codeql.yml/badge.svg)](https://github.com/jermp/sshash/actions/workflows/codeql.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7772316.svg)](https://doi.org/10.5281/zenodo.7772316)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7239205.svg)](https://doi.org/10.5281/zenodo.7239205)

SSHash
======

This is a compressed dictionary data structure for k-mers
(strings of length k over the DNA alphabet {A,C,G,T}), based on **S**parse and **S**kew **Hash**ing.

The data structure is described in the following papers:

* [Sparse and Skew Hashing of K-Mers](https://doi.org/10.1093/bioinformatics/btac245) [1]
* [On Weighted K-Mer Dictionaries](https://almob.biomedcentral.com/articles/10.1186/s13015-023-00226-2) [2,3]

**Please, cite these papers if you use SSHash.**

For a dictionary of n k-mers,
two basic queries are supported:

- i = **Lookup**(g), where i is in [0,n) if the k-mer g is found in the dictionary or i = -1 otherwise;
- g = **Access**(i), where g is the k-mer associated to the identifier i.

If also the weights of the k-mers (their frequency counts) are stored in the dictionary, then the dictionary is said to be *weighted* and it also supports:

- w = **Weight**(i), where i is a given k-mer identifier and w is the weight of the k-mer.

Other supported queries are:

- **Membership Queries**: determine if a given k-mer is present in the dictionary or not.
- **Streaming Queries**: stream through all k-mers of a given DNA file
(.fasta or .fastq formats) to determine their membership to the dictionary.
- **Navigational Queries**: given a k-mer g[1..k] determine if g[2..k]+x is present (forward neighbourhood) and if x+g[1..k-1] is present (backward neighbourhood), for x = A, C, G, T ('+' here means string concatenation).
SSHash internally stores a set of strings, called *contigs* in the following, each associated to a distinct identifier.
If a contig identifier is specified for a navigational query (rather than a k-mer), then the backward neighbourhood of the first k-mer and the forward neighbourhood of the last k-mer in the contig are returned.

If you are interested in a **membership-only** version of SSHash, have a look at [SSHash-Lite](https://github.com/jermp/sshash-lite). It also works for input files with duplicate k-mers (e.g., [matchtigs](https://github.com/algbio/matchtigs) [4]). For a query sequence S and a given coverage threshold E in [0,1], the sequence is considered to be present in the dictionary if at least E*(|S|-k+1) of the k-mers of S are positive.

**NOTE**: It is assumed that two k-mers being the *reverse complement* of each other are the same.

#### Table of contents
* [Compiling the Code](#compiling-the-code)
* [Dependencies](#dependencies)
* [Tools and Usage](#tools-and-usage)
* [Examples](#Examples)
* [Input Files](#input-files)
* [Benchmarks](#benchmarks)
* [Author](#author)
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

    == SSHash: (S)parse and (S)kew (Hash)ing of k-mers =========================

    Usage: ./sshash <tool> ...

    Available tools:
      build                  build a dictionary
      query                  query a dictionary
      check                  check correctness of a dictionary
      bench                  run performance tests for a dictionary
      dump                   write super-k-mers of a dictionary to a fasta file
      permute                permute a weighted input file
      compute-statistics     compute index statistics

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

    ./sshash build -i ../data/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz -k 31 -m 13 --check --bench -o salmonella_enterica.index

This example builds a dictionary for the k-mers read from the file `../data/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz`,
with k = 31 and m = 13. It also check the correctness of the dictionary (`--check` option), run a performance benchmark (`--bench` option), and serializes the index on disk to the file `salmonella_enterica.index`.

To run a performance benchmark after construction of the index,
use:

    ./sshash bench -i salmonella_enterica.index

To also store the weights, use the option `--weighted`:

    ./sshash build -i ../data/unitigs_stitched/with_weights/salmonella_enterica.ust.k31.fa.gz -k 31 -m 13 --weighted --check --verbose

### Example 2

    ./sshash build -i ../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz -k 31 -m 15 -l 2 -o salmonella_100.index

This example builds a dictionary from the input file `../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz` (a pangenome consisting in 100 genomes of *Salmonella Enterica*), with k = 31, m = 15, and l = 2. It also serializes the index on disk to the file `salmonella_100.index`.

To perform some streaming membership queries, use:

    ./sshash query -i salmonella_100.index -q ../data/queries/SRR5833294.10K.fastq.gz

if your queries are meant to be read from a FASTQ file, or

    ./sshash query -i salmonella_100.index -q ../data/queries/salmonella_enterica.fasta.gz --multiline

if your queries are to be read from a (multi-line) FASTA file.

### Example 3

    ./sshash build -i ../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz -k 31 -m 13 -l 4 -s 347692 --canonical-parsing -o salmonella_100.canon.index

This example builds a dictionary from the input file `../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz` (same used in Example 2), with k = 31, m = 13, l = 4, using a seed 347692 for construction (`-s 347692`), and with the canonical parsing modality (option `--canonical-parsing`). The dictionary is serialized on disk to the file `salmonella_100.canon.index`.

The "canonical" version of the dictionary offers more speed for only a little space increase (for a suitable choice of parameters m and l), especially under low-hit workloads -- when the majority of k-mers are not found in the dictionary. (For all details, refer to the paper.)

Below a comparison between the dictionary built in Example 2 (not canonical)
and the one just built (Example 3, canonical).

    ./sshash query -i salmonella_100.index -q ../data/queries/SRR5833294.10K.fastq.gz

    ./sshash query -i salmonella_100.canon.index -q ../data/queries/SRR5833294.10K.fastq.gz

The canonical dictionary can be twice as fast as the regular dictionary
for low-hit workloads, even on this tiny example, for only +0.4 bits/k-mer.

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
These collections can be obtained, for example, by extracting the maximal unitigs of a de Bruijn graph.

To do so, we can use the tool [BCALM2](https://github.com/GATB/bcalm).
This tool builds a compacted de Bruijn graph and outputs its maximal unitigs.
From the output of BCALM2, we can then *stitch* (i.e., glue) some unitigs to reduce the number of nucleotides. The stitiching process is carried out using the [UST](https://github.com/jermp/UST) tool.

**NOTE**: Input files are expected to have **one DNA sequence per line**. If a sequence spans multiple lines (e.g., multi-fasta), the lines should be concatenated before indexing.

Below we provide a complete example (assuming both BCALM2 and UST are installed correctly) that downloads the Human (GRCh38) Chromosome 13 and extracts the maximal stitiched unitigs for k = 31.

    mkdir DNA_datasets
    wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz -O DNA_datasets/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
    ~/bcalm/build/bcalm -in ~/DNA_datasets/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz -kmer-size 31 -abundance-min 1 -nb-cores 8
    ~/UST/ust -k 31 -i ~/Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa
    gzip Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa.ust.fa
    rm ~/Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa

#### Datasets

The script `scripts/download_and_preprocess_datasets.sh`
contains all the needed steps to download and pre-process
the datasets that we used in [1].

For the experiments in [2] and [3], we used the datasets available on [Zenodo](https://doi.org/10.5281/zenodo.7772316).

#### Weights
Using the option `-all-abundance-counts` of BCALM2, it is possible to also include the abundance counts of the k-mers in the BCALM2 output. Then, use the option `-a 1` of UST to include such counts in the stitched unitigs.

Benchmarks
----------

For some example benchmarks, see the folder `/benchmarks`.

Some more large-scale benchmarks below.

*Pinus Taeda* ("pine", [GCA_000404065.3](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/404/065/GCA_000404065.3_Ptaeda2.0/GCA_000404065.3_Ptaeda2.0_genomic.fna.gz)) and *Ambystoma Mexicanum* ("axolotl", [GCA_002915635.2](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/915/635/GCA_002915635.3_AmbMex60DD/GCA_002915635.3_AmbMex60DD_genomic.fna.gz))
are some of the largest genome assemblies, respectively counting
10,508,232,575 and 17,987,935,180 distinct k-mers for k = 31.

After running BCALM2 and UST, we build the indexes as follows.

    ./sshash build -i ~/DNA_datasets.larger/GCA_000404065.3_Ptaeda2.0_genomic.ust_k31.fa.gz -k 31 -m 20 -l 6 -c 7 -o pinus.m20.index
    ./sshash build -i ~/DNA_datasets.larger/GCA_000404065.3_Ptaeda2.0_genomic.ust_k31.fa.gz -k 31 -m 19 -l 6 -c 7 --canonical-parsing -o pinus.m19.canon.index
    ./sshash build -i ~/DNA_datasets.larger/GCA_002915635.3_AmbMex60DD_genomic.ust_k31.fa.gz -k 31 -m 21 -l 6 -c 7 -o axolotl.m21.index
    ./sshash build -i ~/DNA_datasets.larger/GCA_002915635.3_AmbMex60DD_genomic.ust_k31.fa.gz -k 31 -m 20 -l 6 -c 7 --canonical-parsing -o axolotl.m20.canon.index

The following table summarizes the space of the dictionaries.

| Dictionary        |Pine       || Axolotl  ||
|:------------------|:---:|:----:|:---:|:---:|
|                   | GB     | bits/k-mer  | GB    | bits/k-mer |
| SSHash, regular   | 13.21  | 10.06       | 22.28 | 9.91       |
| SSHash, canonical | 14.94  | 11.37       | 25.03 | 11.13      |



To query the dictionaries, we use [SRR17023415](https://www.ebi.ac.uk/ena/browser/view/SRR17023415) fastq reads
(23,891,117 reads, each of 150 bases) for the pine,
and [GSM5747680](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5747680) multi-line fasta (15,548,160 lines) for the axolotl.

Timings have been collected on an Intel Xeon Platinum 8276L CPU @ 2.20GHz,
using a single thread.

| Dictionary        |Pine       || Axolotl  ||
|:------------------|:---:|:----:|:---:|:---:|
|                   |(>75% hits)||(>86% hits)|
|                   | tot (min) | avg (ns/k-mer) | tot (min) | avg (ns/k-mer) |
| SSHash, regular   | 19.2      | 400            | 4.2       | 269            |
| SSHash, canonical | 14.8      | 310            | 3.2       | 208            |

Below the complete query reports.

    ./sshash query -i pinus.m20.index -q ~/DNA_datasets.larger/queries/SRR17023415_1.fastq.gz
    ==== query report:
    num_kmers = 2866934040
    num_valid_kmers = 2866783488 (99.9947% of kmers)
    num_positive_kmers = 2151937575 (75.0645% of valid kmers)
    num_searches = 418897117/2151937575 (19.466%)
    num_extensions = 1733040458/2151937575 (80.534%)
    elapsed = 1146.58 sec / 19.1097 min / 399.933 ns/kmer

    ./sshash query -i pinus.m19.canon.index -q ~/DNA_datasets.larger/queries/SRR17023415_1.fastq.gz
    ==== query report:
    num_kmers = 2866934040
    num_valid_kmers = 2866783488 (99.9947% of kmers)
    num_positive_kmers = 2151937575 (75.0645% of valid kmers)
    num_searches = 359426304/2151937575 (16.7025%)
    num_extensions = 1792511271/2151937575 (83.2975%)
    elapsed = 889.779 sec / 14.8297 min / 310.359 ns/kmer

    ./sshash query -i axolotl.m21.index -q ~/DNA_datasets.larger/queries/Axolotl.Trinity.CellReports2017.fasta.gz --multiline
    ==== query report:
    num_kmers = 931366757
    num_valid_kmers = 748445346 (80.3599% of kmers)
    num_positive_kmers = 650467884 (86.9092% of valid kmers)
    num_searches = 124008258/650467884 (19.0645%)
    num_extensions = 526459626/650467884 (80.9355%)
    elapsed = 250.173 sec / 4.16955 min / 268.608 ns/kmer

    ./sshash query -i axolotl.m20.canon.index -q ~/DNA_datasets.larger/queries/Axolotl.Trinity.CellReports2017.fasta.gz --multiline
    ==== query report:
    num_kmers = 931366757
    num_valid_kmers = 748445346 (80.3599% of kmers)
    num_positive_kmers = 650467884 (86.9092% of valid kmers)
    num_searches = 106220473/650467884 (16.3299%)
    num_extensions = 544247411/650467884 (83.6701%)
    elapsed = 193.871 sec / 3.23119 min / 208.158 ns/kmer

Author
------

[Giulio Ermanno Pibiri](https://jermp.github.io) - <giulioermanno.pibiri@unive.it>

References
-----
* [1] Giulio Ermanno Pibiri. [Sparse and Skew Hashing of K-Mers](https://doi.org/10.1093/bioinformatics/btac245). Bioinformatics. 2022.
* [2] Giulio Ermanno Pibiri. [On Weighted K-Mer Dictionaries](https://drops.dagstuhl.de/opus/volltexte/2022/17043/). International Workshop on Algorithms in Bioinformatics (WABI). 2022.
* [3] Giulio Ermanno Pibiri. [On Weighted K-Mer Dictionaries](https://almob.biomedcentral.com/articles/10.1186/s13015-023-00226-2). Algorithms for Molecular Biology (ALGOMB). 2023.
* [4] Schmidt, S., Khan, S., Alanko, J., Pibiri, G. E., and Tomescu, A. I. [Matchtigs: minimum plain text representation of k-mer sets](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02968-z). Genome Biology 24, 136. 2023.
