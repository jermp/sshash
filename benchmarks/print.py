#!/usr/bin/env python3

import sys
import json
import os
from statistics import mean, StatisticsError
import math

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


def parse_streaming_file(path, canonical_flag):
    """Parse streaming queries JSON file."""
    stream_data = {}
    if not os.path.exists(path):
        return stream_data

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
            canonical = "yes" if canonical_flag else "no"

            key = (collection, canonical)
            num_kmers = int(d["num_kmers"])
            num_pos = int(d["num_positive_kmers"])
            num_ext = int(d["num_extensions"])
            elapsed_ms = int(d["elapsed_millisec"])

            ns_per_kmer = int(math.ceil(elapsed_ms * 1e6 / num_kmers))
            hit_rate = (num_pos / num_kmers) * 100 if num_kmers else 0
            extension_rate = (num_ext / num_pos) * 100 if num_pos else 0

            stream_data[key] = {
                "ns_per_kmer": f"{ns_per_kmer}",
                "hit_rate": f"{hit_rate:.2f}",
                "extension_rate": f"{extension_rate:.2f}"
            }
    return stream_data


def main():
    if len(sys.argv) != 2:
        print("Usage: print.py input_dir", file=sys.stderr)
        sys.exit(1)

    input_dir = sys.argv[1]
    reg_build_path = input_dir + "/regular-build.json"
    canon_build_path = input_dir + "/canon-build.json"
    reg_bench_path = input_dir + "/regular-bench.json"
    canon_bench_path = input_dir + "/canon-bench.json"
    reg_stream_path = input_dir + "/regular-streaming-queries-high-hit.json"
    canon_stream_path = input_dir + "/canon-streaming-queries-high-hit.json"

    reg_build = parse_build_file(reg_build_path, False)
    canon_build = parse_build_file(canon_build_path, True)
    reg_bench = parse_bench_file(reg_bench_path, False)
    canon_bench = parse_bench_file(canon_bench_path, True)
    reg_stream = parse_streaming_file(reg_stream_path, False)
    canon_stream = parse_streaming_file(canon_stream_path, True)

    # merge everything
    all_builds = reg_build + canon_build
    lookup_all = {**reg_bench, **canon_bench}
    stream_all = {**reg_stream, **canon_stream}

    # CSV header
    print("k,Collection,m,canonical,bits_per_kmer,total_GB,build_time,positive_lookup_ns,negative_lookup_ns,access_ns,iteration_ns,ns_per_kmer,hit_rate,extension_rate")

    for r in sorted(all_builds, key=lambda x: (int(x["k"]), x["Collection"], x["canonical"])):
        lookup = lookup_all.get(
            (r["Collection"], r["m"], r["canonical"]), # key
            {"pos": "NA", "neg": "NA", "access": "NA", "iter": "NA", "k": r["k"]})
        stream = stream_all.get(
            (r["Collection"], r["canonical"]), # key
            {"ns_per_kmer": "NA", "hit_rate": "NA", "extension_rate": "NA"})

        print(f"{r['k']},{r['Collection']},{r['m']},{r['canonical']},{r['bits_per_kmer']},{r['total_GB']},{r['build_time']},{lookup['pos']},{lookup['neg']},{lookup['access']},{lookup['iter']},{stream['ns_per_kmer']},{stream['hit_rate']},{stream['extension_rate']}")

if __name__ == "__main__":
    main()
