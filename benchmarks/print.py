#!/usr/bin/env python3

import sys
import json
import os
from statistics import mean, StatisticsError

def format_time(microseconds):
    seconds = microseconds / 1_000_000
    minutes = int(seconds // 60)
    seconds = int(seconds % 60)
    return f"{minutes}:{seconds:02d}"

def parse_build_file(path, canonical_flag):
    """Parse build JSONL file."""
    results = []
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

            bits_per_kmer = (index_bytes * 8) / num_kmers
            gb = index_bytes / 1e9
            build_time_fmt = format_time(build_time_us)

            fname = os.path.basename(d["input_filename"])
            collection = fname.split(".")[0].capitalize()
            k = d["k"]

            results.append({
                "k": k,
                "Collection": collection,
                "m": d["m"],
                "canonical": "yes" if canonical_flag else "no",
                "bits_per_kmer": f"{bits_per_kmer:.2f}",
                "total_GB": f"{gb:.2f}",
                "build_time": build_time_fmt
            })
    return results

def parse_bench_file(path, canonical_flag):
    """Parse benchmark JSONL file and average per collection."""
    lookup_data = {}
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

            fname = os.path.basename(d["index_filename"])
            collection = fname.split(".")[0].capitalize()
            m = d["m"]
            k = d["k"]
            canonical = "yes" if canonical_flag else "no"

            key = (collection, m, canonical)
            entry = lookup_data.setdefault(key, {
                "k": k,
                "pos": [], "neg": [], "access": [], "iter": []
            })
            entry["pos"].append(float(d["positive lookup (avg_nanosec_per_kmer)"]))
            entry["neg"].append(float(d["negative lookup (avg_nanosec_per_kmer)"]))
            entry["access"].append(float(d["access (avg_nanosec_per_kmer)"]))
            entry["iter"].append(float(d["iterator (avg_nanosec_per_kmer)"]))

    # average the results
    for k, v in lookup_data.items():
        try:
            lookup_data[k] = {
                "k": v["k"],
                "pos": f"{mean(v['pos'])/1000:.2f}",
                "neg": f"{mean(v['neg'])/1000:.2f}",
                "access": f"{mean(v['access'])/1000:.2f}",
                "iter": f"{mean(v['iter']):.2f}",
            }
        except StatisticsError:
            lookup_data[k] = {"k": v["k"], "pos": "NA", "neg": "NA", "access": "NA", "iter": "NA"}
    return lookup_data

def main():
    if len(sys.argv) != 5:
        print("Usage: index_stats_csv.py regular-build.json canon-build.json regular-bench.json canon-bench.json", file=sys.stderr)
        sys.exit(1)

    reg_build_path, canon_build_path, reg_bench_path, canon_bench_path = sys.argv[1:]

    reg_build = parse_build_file(reg_build_path, False)
    canon_build = parse_build_file(canon_build_path, True)
    reg_bench = parse_bench_file(reg_bench_path, False)
    canon_bench = parse_bench_file(canon_bench_path, True)

    # merge everything
    all_builds = reg_build + canon_build
    lookup_all = {**reg_bench, **canon_bench}

    # CSV header
    print("k,Collection,m,canonical,bits_per_kmer,total_GB,build_time,positive_lookup_ns,negative_lookup_ns,access_ns,iteration_ns")

    for r in sorted(all_builds, key=lambda x: (int(x["k"]), x["Collection"], x["canonical"])):
        key = (r["Collection"], r["m"], r["canonical"])
        lookup = lookup_all.get(key, {"pos": "NA", "neg": "NA", "access": "NA", "iter": "NA", "k": r["k"]})

        print(f"{r['k']},{r['Collection']},{r['m']},{r['canonical']},{r['bits_per_kmer']},{r['total_GB']},{r['build_time']},{lookup['pos']},{lookup['neg']},{lookup['access']},{lookup['iter']}")

if __name__ == "__main__":
    main()

