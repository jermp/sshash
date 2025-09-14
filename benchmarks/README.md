[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7239205.svg)](https://doi.org/10.5281/zenodo.7239205)

Benchmarks
----------

For these benchmarks we used the whole genomes of the following organisms:

- Gadus Morhua ("Cod")
- Falco Tinnunculus ("Kestrel")
- Homo Sapiens ("Human")

for k = 31 and 63.

The datasets and queries used in these benchmarks can be downloaded
by running the script

```
bash download-datasets.sh
```

To run the benchmarks, from within the `build` directory, run

```
bash ../script/build.sh [prefix]
bash ../script/bench.sh [prefix]
bash ../script/streaming-query-high-hit.sh [prefix]
bash ../script/streaming-query-low-hit.sh [prefix]
```

where `[prefix]` should be replaced by a suitable basename, e.g., the current date.

These are the results obtained on 14/09/25 (see logs [here](results-14-09-25))
on a machine equipped with an Intel Xeon W-2245 CPU @ 3.90GHz, and running Ubuntu 18.04.6.
The code was compiled with `gcc` 10.3.

![](results-14-09-25/results.png)
