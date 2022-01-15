[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/jermp/sshash.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/jermp/sshash/context:cpp)


SSHash
======

This is a compressed dictionary data structure for k-mers
(strings of length k over the DNA alphabet {A,C,G,T}), based on **S**parse and **S**kew **Hash**ing.

**A (pre-print) paper describing the data structure will appear soon.**

For a dictionary of n k-mers,
two basic queries are supported:

- i = Lookup(g), where i is in [0,n) if the k-mer g is found in the dictionary or i = -1 otherwise;
- g = Access(i), where g is the k-mer associated to the identifier i.

A membership query (determine if a given k-mer is present in the dictionary or not) is, therefore, supported by means of the lookup query.
The dictionary can also stream through all k-mers of a given DNA file
(.fasta or .fastq formats) to determine their membership to the dictionary.

**NOTE**: The Lookup query assumes that two k-mers being the *reverse complement* of each other are the same.

#### Table of contents
* [Compiling the Code](#compiling-the-code)
* [Dependencies](#dependencies)
* [Build a Dictionary](#build-a-dictionary)
* [Examples](#Examples)
* [Input Files](#input-files)

Compiling the Code
------------------

The code is tested on Linux with `gcc` and on Mac with `clang`.
To build the code, [`CMake`](https://cmake.org/) is required.

Clone the repository with

	git clone --recursive https://github.com/jermp/sshash.git

If you have cloned the repository without `--recursive`, you will need to perform the following commands before
compiling:

    git submodule init
    git submodule update

To compile the code for a release environment (see file `CMakeLists.txt` for the used compilation flags), it is sufficient to do the following:

    mkdir build
    cd build
    cmake ..
    make -j

For a testing environment, use the following instead:

    mkdir debug_build
    cd debug_build
    cmake .. -D CMAKE_BUILD_TYPE=Debug -D SSHASH_USE_SANITIZERS=On
    make

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

Build a Dictionary
------------------

The driver program
called `build` can be used to build a dictionary.

From within the directory
where the code was compiled (see the section [Compiling the Code](#compiling-the-code)), run the command:

	./build --help

to show the usage of the driver program (reported below for convenience).

	Usage: ./build [-h,--help] input_filename k m [-s seed] [-n max_num_kmers] [-l l] [-c c] [--canonical-parsing] [-o output_filename] [--check] [--bench] [--verbose]
	
	 input_filename
		Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:
		- without duplicate nor invalid kmers
		- one DNA sequence per line.
		For example, it could be the de Bruijn graph topology output by BCALM.
	
	 k
		K-mer length.
	
	 m
		Minimizer length (must be <= k).
	
	 [-s seed]
		Seed for construction (default is 1).
	
	 [-n max_num_kmers]
		Build the dictionary from at most this number of k-mers.
	
	 [-l l]
		A (integer) constant that controls the space/time trade-off of the dictionary. A reasonable values lies between 2 and 12 (default is 6).
	
	 [-c c]
		A (floating point) constant that trades construction speed for space effectiveness of minimal perfect hashing. A reasonable value lies between 3.0 and 10.0 (default is 3.000000).
	
	 [--canonical-parsing]
		Canonical parsing of k-mers. This option changes the parsing and results in a trade-off between index space and lookup time.
	
	 [-o output_filename]
		Output file name where the data structure will be serialized.
	
	 [--check]
		Check correctness after construction.
	
	 [--bench]
		Run benchmark after construction.
	
	 [--verbose]
		Verbose output during construction.
	
	 [-h,--help]
		Print this help text and silently exits.
		

Examples
--------

For the examples, we are going to use some collections
of *stitched unitigs* from the directory `../data/unitigs_stitched`.
These collections were built for k = 31, so dictionaries should be built with k = 31 as well to ensure correctness.

In the section [Input Files](#input-files), we explain how
such collections of stitched unitigs can be obtained from raw FASTA files.

### Example 1

	./build ../data/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz 31 13 --check --bench -o salmonella_enterica.index

This example builds a dictionary for the k-mers read from the file `../data/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz`,
with k = 31 and m = 13. It also check the correctness of the dictionary (`--check` option), run a performance benchmark (`--bench` option), and serializes the index on disk to the file `salmonella_enterica.index`.

To run a performance benchmark after construction of the index,
use:

	./bench salmonella_enterica.index
	
### Example 2

	./build ../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz 31 15 -l 2 -o salmonella_100.index

This example builds a dictionary from the input file `../data/unitigs_stitched/salmonella_100_k31_ust.fa.gz` (a pangenome consisting in 100 genomes of *Salmonella Enterica*), with k = 31, m = 15, and l = 2. It also serializes the index on disk to the file `salmonella_100.index`.

To perform some streaming membership queries, use:

	./query salmonella_100.index ../data/queries/SRR5833294.10K.fastq.gz

if your queries are meant to be read from a FASTQ file, or

	./query salmonella_100.index ../data/queries/salmonella_enterica.fasta.gz --multiline
	
if your queries are to be read from a (multi-line) FASTA file.


Input Files
-----------

SSHash is meant to index k-mers from collections that do not contain duplicates
nor invalid k-mers (strings containing symbols different from {A,C,G,T}).
These collections can be obtained, for example, by extracting the maximal unitigs of a de Bruijn graph.

Doing so is easy to do using the tool [BCALM2](https://github.com/GATB/bcalm).
This tool builds a compacted de Bruijn graph and outputs its maximal unitigs.
From the output of BCALM2, we can then *stitch* (i.e., glue) some unitigs to reduce the number of nucleotides. The stitiching process is carried out using the [UST](https://github.com/medvedevgroup/UST) tool.

Below we provide a complete example (assuming both BCALM2 and UST are installed correctly) that downloads the Human (GRCh38) Chromosome 13 and extracts the maximal stitiched unitigs for k = 31.

	mkdir DNA_datasets
	wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz -O DNA_datasets/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
	~/bcalm/build/bcalm -in ~/DNA_datasets/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz -kmer-size 31 -abundance-min 1 -nb-cores 8
	~/UST/ust -k 31 -i ~/Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa
	gzip Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa.ust.fa
	rm ~/Homo_sapiens.GRCh38.dna.chromosome.13.fa.unitigs.fa

See also the script `scripts/download_and_preprocess_datasets.sh`
for precise arguments.
