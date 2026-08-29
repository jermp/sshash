[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17582116.svg)](https://doi.org/10.5281/zenodo.17582116)

Benchmarks
----------

For these benchmarks we used the datasets available here
[https://zenodo.org/records/17582116](https://zenodo.org/records/17582116).

To run the benchmarks, from within the `build` directory, run

    python3 ../script/build.py <log_label> <input_datasets_dir> <output_index_dir>
    python3 ../script/bench.py <log_label> <input_index_dir>
    python3 ../script/streaming-query.py <log_label> <input_index_dir> <input_queries_dir>

where `<log_label>` should be replaced by a suitable basename, e.g., the current date.

These are the results obtained on 28/08/26 (see logs [here](results-28-08-26))
on a machine equipped with an AMD Ryzen Threadripper PRO 7985WX processor clocked at 5.40GHz.
The code was compiled with `gcc` 13.3.0.
The indexes were build with a max RAM usage of 16 GB and 64 threads.
Queries were run using one thread, instead.

| k  | Collection | m  | Space (bits/kmer) | Space (total GB) | Building time (m:ss) | Positive random lookup (µs/kmer) | Negative random lookup (µs/kmer) | Random Access (µs/kmer) | Streaming Lookup high-hit (ns/kmer) |
|----|------------|:--:|:-----------------:|:----------------:|:--------------------:|:--------------------------------:|:--------------------------------:|:-----------------------:|:-----------------------------------:|
||
| 31 | Cod        | 20 | 7.94              | 0.50             | 0:27                 | 0.43                             | 0.36                             | 0.27                    | 24                                  |
|    | Kestrel    | 20 | 7.52              | 1.08             | 1:00                 | 0.44                             | 0.39                             | 0.28                    | 39                                  |
|    | Human      | 21 | 8.77              | 2.75             | 2:35                 | 0.56                             | 0.43                             | 0.35                    | 76                                  |
|    | NCBI-v     | 19 | 7.43              | 0.35             | 0:14                 | 0.41                             | 0.35                             | 0.26                    | 26                                  |
|    | SE         | 21 | 10.24             | 1.14             | 0:54                 | 0.62                             | 0.40                             | 0.36                    | 163                                 |
|    | HPRC       | 21 | 10.64             | 4.94             | 4:36                 | 0.70                             | 0.46                             | 0.53                    | 84                                  |
||
| 63 | Cod        | 24 | 4.62              | 0.32             | 0:24                 | 0.57                             | 0.45                             | 0.29                    | 37                                  |
|    | Kestrel    | 24 | 3.84              | 0.55             | 0:19                 | 0.54                             | 0.48                             | 0.32                    | 44                                  |
|    | Human      | 25 | 4.92              | 1.70             | 1:13                 | 0.66                             | 0.51                             | 0.36                    | 110                                 |
|    | NCBI-v     | 23 | 4.06              | 0.21             | 0:07                 | 0.53                             | 0.44                             | 0.28                    | 30                                  |
|    | SE         | 31 | 7.42              | 1.41             | 1:20                 | 0.97                             | 0.50                             | 0.40                    | 278                                 |
|    | HPRC       | 31 | 7.63              | 5.65             | 4:46                 | 0.94                             | 0.58                             | 0.62                    | 137                                 |

The results of a run, for both k=31 and k=63 at once, can be exported to CSV
format (every collected quantity) or to a markdown table (the main quantities)
with

    python3 ../benchmarks/print_table.py <results_dir>
    python3 ../benchmarks/print_table.py <results_dir> --md

where `<results_dir>` is the directory written by the scripts above (e.g.,
`results-21-01-26`), containing the `k31` and `k63` subdirectories.

Note that the scripts now produce a single set of result files per `k`
(`build.json`, `bench.json`, and `streaming-queries.json`), since the
regular/canonical distinction is gone: there is one indexing modality only.
Result directories archived before this change instead hold a `regular-` and a
`canon-` file for each of those; to re-read them, use the `print_csv.py`
script from the corresponding commit.
