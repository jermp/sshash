Benchmarks
----------

For these benchmarks we used the whole genomes
of the following organisms:

- C. Elegans ("elegans")
- Gadus Morhua ("cod")
- Falco Tinnunculus ("kestrel")
- Homo Sapiens ("human")

for k = 31, 47, and 63.

The datasets can be downloaded from [zenodo](https://zenodo.org/record/7239205#.Y-K61OzMI-Q).

For streaming queries we used FASTQ files downloaded from [ENA](https://www.ebi.ac.uk/ena/browser/home), using accession numbers:

- elegans: SRR16288382
- cod: SRR12858649
- kestrel: SRR11449743
- human: SRR5833294

The query times are relative to the following configuration:

- Processor: Intel i9-9900K @ 3.60 GHz;
- Compiler: gcc 11.2.0;
- OS: GNU/Linux 5.13.0-52-generic x86_64.