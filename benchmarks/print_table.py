#!/usr/bin/env python3

"""
Print the benchmark results of both k=31 and k=63 as one table.

The input directory must contain the two subdirectories `k31` and `k63`,
each with the files written by the benchmark scripts:
`build.json`, `bench.json`, and `streaming-queries.json`.

The default output is CSV, with every collected quantity; with `--md`,
a markdown table with the main quantities is printed instead.
"""

import argparse
import json
import os
import sys
import math
from statistics import mean, StatisticsError

K_VALUES = [31, 63]

# fixed dataset order and display names
DATASETS = ["cod", "kestrel", "human", "ncbi-virus", "se", "hprc"]
DISPLAY_NAME = {
    "cod": "Cod",
    "kestrel": "Kestrel",
    "human": "Human",
    "ncbi-virus": "NCBI-v",
    "se": "SE",
    "hprc": "HPRC",
}

def format_time(microseconds):
    seconds = microseconds / 1_000_000
    minutes = int(seconds // 60)
    seconds = int(seconds % 60)
    return f"{minutes}:{seconds:02d}"

def dataset_of(filename):
    return os.path.basename(filename).split(".")[0]

def parse_build_file(path):
    """Parse build JSONL file: one record per dataset."""
    results = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                print(f"Skipping invalid JSON line in {path}", file=sys.stderr)
                continue

            num_kmers = int(d["num_kmers"])
            index_bytes = int(d["index_size_in_bytes"])
            build_time_us = int(d["total_build_time_in_microsec"])

            results[dataset_of(d["input_filename"])] = {
                "m": d["m"],
                "bits_per_kmer": f"{(index_bytes * 8) / num_kmers:.2f}",
                "total_GB": f"{index_bytes / 1e9:.2f}",
                "build_time": format_time(build_time_us),
            }
    return results

def parse_bench_file(path):
    """Parse benchmark JSONL file and average the runs per dataset."""
    runs = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                print(f"Skipping invalid JSON line in {path}", file=sys.stderr)
                continue

            entry = runs.setdefault(dataset_of(d["index_filename"]),
                                    {"pos": [], "neg": [], "access": [], "iter": []})
            entry["pos"].append(float(d["positive lookup (avg_nanosec_per_kmer)"]))
            entry["neg"].append(float(d["negative lookup (avg_nanosec_per_kmer)"]))
            entry["access"].append(float(d["access (avg_nanosec_per_kmer)"]))
            entry["iter"].append(float(d["iterator (avg_nanosec_per_kmer)"]))

    results = {}
    for ds, v in runs.items():
        try:
            results[ds] = {
                "pos_us": f"{mean(v['pos']) / 1000:.2f}",
                "neg_us": f"{mean(v['neg']) / 1000:.2f}",
                "access_us": f"{mean(v['access']) / 1000:.2f}",
                "iter_ns": f"{mean(v['iter']):.2f}",
            }
        except StatisticsError:
            results[ds] = {"pos_us": "NA", "neg_us": "NA", "access_us": "NA", "iter_ns": "NA"}
    return results

def parse_streaming_file(path):
    """Parse streaming queries JSONL file: one record per dataset."""
    results = {}
    if not os.path.exists(path):
        return results

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                print(f"Skipping invalid JSON line in {path}", file=sys.stderr)
                continue

            num_kmers = int(d["num_kmers"])
            num_pos = int(d["num_positive_kmers"])
            num_ext = int(d["num_extensions"])
            elapsed_ms = int(d["elapsed_millisec"])

            results[dataset_of(d["index_filename"])] = {
                "ns_per_kmer": f"{int(math.ceil(elapsed_ms * 1e6 / num_kmers))}",
                "hit_rate": f"{(num_pos / num_kmers) * 100 if num_kmers else 0:.2f}",
                "extension_rate": f"{(num_ext / num_pos) * 100 if num_pos else 0:.2f}",
            }
    return results

def collect(input_dir):
    """One row per (k, dataset), in fixed order."""
    rows = []
    for k in K_VALUES:
        k_dir = os.path.join(input_dir, f"k{k}")
        if not os.path.isdir(k_dir):
            print(f"Error: expected subdirectory '{k_dir}'", file=sys.stderr)
            sys.exit(1)
        build = parse_build_file(os.path.join(k_dir, "build.json"))
        bench = parse_bench_file(os.path.join(k_dir, "bench.json"))
        stream = parse_streaming_file(os.path.join(k_dir, "streaming-queries.json"))
        first_row_of_k = len(rows)
        for ds in DATASETS:
            if ds not in build:
                continue
            row = {"k": str(k), "Collection": DISPLAY_NAME[ds], "first_of_k": False}
            row.update(build[ds])
            row.update(bench.get(ds, {"pos_us": "NA", "neg_us": "NA",
                                      "access_us": "NA", "iter_ns": "NA"}))
            row.update(stream.get(ds, {"ns_per_kmer": "NA", "hit_rate": "NA",
                                       "extension_rate": "NA"}))
            rows.append(row)
        if len(rows) > first_row_of_k:
            rows[first_row_of_k]["first_of_k"] = True
    return rows

def print_csv(rows):
    print("k,Collection,m,bits_per_kmer,total_GB,build_time,"
          "positive_lookup_us,negative_lookup_us,access_us,iteration_ns,"
          "streaming_ns_per_kmer,hit_rate,extension_rate")
    for r in rows:
        print(f"{r['k']},{r['Collection']},{r['m']},{r['bits_per_kmer']},{r['total_GB']},"
              f"{r['build_time']},{r['pos_us']},{r['neg_us']},{r['access_us']},{r['iter_ns']},"
              f"{r['ns_per_kmer']},{r['hit_rate']},{r['extension_rate']}")

MD_COLUMNS = [
    # (header, row key, alignment: 'left' or 'center')
    ("k", "k", "left"),
    ("Collection", "Collection", "left"),
    ("m", "m", "center"),
    ("Space (bits/kmer)", "bits_per_kmer", "center"),
    ("Space (total GB)", "total_GB", "center"),
    ("Building time (m:ss)", "build_time", "center"),
    ("Positive random lookup (µs/kmer)", "pos_us", "center"),
    ("Negative random lookup (µs/kmer)", "neg_us", "center"),
    ("Random Access (µs/kmer)", "access_us", "center"),
    ("Streaming Lookup high-hit (ns/kmer)", "ns_per_kmer", "center"),
]

def print_md(rows):
    """The k value is shown on the first row of its block only,
    and each block is preceded by an empty '||' row."""
    cells = [[str(r[key]) if key != "k" or r["first_of_k"] else ""
              for _, key, _ in MD_COLUMNS] for r in rows]
    widths = [max(len(header), *(len(row[i]) for row in cells))
              for i, (header, _, _) in enumerate(MD_COLUMNS)]

    def line(values):
        return "| " + " | ".join(v.ljust(w) for v, w in zip(values, widths)) + " |"

    print(line([header for header, _, _ in MD_COLUMNS]))
    print("|" + "|".join(f":{'-' * w}:" if align == "center" else "-" * (w + 2)
                         for w, (_, _, align) in zip(widths, MD_COLUMNS)) + "|")
    for r, row_cells in zip(rows, cells):
        if r["first_of_k"]:
            print("||")
        print(line(row_cells))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", help="directory containing the k31 and k63 subdirectories")
    parser.add_argument("--md", action="store_true",
                        help="print a markdown table instead of CSV")
    args = parser.parse_args()

    rows = collect(args.input_dir)
    if args.md:
        print_md(rows)
    else:
        print_csv(rows)

if __name__ == "__main__":
    main()
