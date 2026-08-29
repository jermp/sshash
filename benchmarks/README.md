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

These are the results obtained on 21/01/26 (see logs [here](results-21-01-26))
on a machine equipped with an AMD Ryzen Threadripper PRO 7985WX processor clocked at 5.40GHz.
The code was compiled with `gcc` 13.3.0.

The indexes were build with a max RAM usage of 16 GB and 64 threads.
Queries were run using one thread, instead.

![](results-21-01-26/results.png)

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
