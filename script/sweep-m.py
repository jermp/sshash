#!/usr/bin/env python3

import os
import subprocess
import sys
from pathlib import Path

if len(sys.argv) < 4:
    print("Usage: python3 sweep-m.py <log_label> <input_datasets_dir> <output_index_dir>")
    sys.exit(1)

log_label = sys.argv[1]
datasets_dir = Path(sys.argv[2]).resolve()
index_dir = Path(sys.argv[3]).resolve()
tmp_dir = datasets_dir / "tmp_dir"
results_dir = Path(f"results-{log_label}")
threads = 16
g = 16

# 1. Target datasets
datasets = ["human", "se"]

# 2. Define the sweeps for m (minimizer length)
m_sweeps_k31 = {
    "human": [17, 19, 21, 23, 25],
    "se":    [17, 19, 21, 23, 25]
}

m_sweeps_k63 = {
    "human": [21, 23, 25, 27, 29],
    "se":    [23, 25, 27, 29, 31]
}

# --- Utilities ---
def run_cmd(cmd, cwd=None, append_to=None):
    print(f"[RUN] {' '.join(cmd)}")
    if append_to:
        with open(append_to, "a") as f:
            subprocess.run(cmd, cwd=cwd, stdout=f, stderr=f, check=True)
    else:
        subprocess.run(cmd, cwd=cwd, check=True)

def build_project(max_k63: bool):
    flag = "On" if max_k63 else "Off"
    print(f"\n=== Building SSHASH (MAX_KMER_LENGTH_63={flag}) ===\n")
    run_cmd([
        "cmake", "..",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_CXX_COMPILER=/usr/bin/g++",
        "-DSSHASH_USE_ARCH_NATIVE=On",
        "-DSSHASH_USE_SANITIZERS=Off",
        f"-DSSHASH_USE_MAX_KMER_LENGTH_63={flag}"
    ])
    run_cmd(["make", "-j"])

def build_sshash(k, dataset, m_val):
    # Differentiate results dir by m_val
    out_dir = results_dir / f"m{m_val}" / f"k{k}"
    out_dir.mkdir(parents=True, exist_ok=True)

    log_file = out_dir / "build.log"
    json_file = out_dir / "build.json"
    time_file = out_dir / "build.time.log"

    input_file = datasets_dir / f"{dataset}.k{k}.eulertigs.fa.gz"

    # Append m_val to the output filename
    output_file = index_dir / f"{dataset}.k{k}.m{m_val}"

    print(f"\n>>> Building {dataset} (k={k}, m={m_val})\n")

    # Clean tmp directory (should be empty after each build anyway)
    subprocess.run(f"rm -rf {tmp_dir}/*", shell=True, check=True)

    cmd = [
        "/usr/bin/time", "-v", "-a", "-o", str(time_file),
        "./sshash", "build",
        "-i", str(input_file),
        "-k", str(k),
        "-m", str(m_val),
        "-g", str(g),
        "-t", str(threads),
        "--verbose",
        "-d", str(tmp_dir),
        "-o", f"{output_file}.sshash"
    ]

    # Append stdout to .log, stderr to .json
    with open(log_file, "a") as log, open(json_file, "a") as js:
        subprocess.run(cmd, stdout=log, stderr=js, check=True)

def run_bench(k, dataset, m_val, runs=3):
    """Run SSHASH benchmark for a specific dataset and m."""
    # Store results in the specific m_val / k folder
    out_dir = results_dir / f"m{m_val}" / f"k{k}"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_file = out_dir / "bench.log"
    json_file = out_dir / "bench.json"

    # Match the index naming scheme that includes m_val
    index_path = index_dir / f"{dataset}.k{k}.m{m_val}.sshash"

    print(f"\n>>> Benchmarking {dataset} (k={k}, m={m_val})\n")
    for i in range(runs):
        print(f"  ==> run {i+1}/{runs}")
        cmd = ["./sshash", "bench", "-i", str(index_path)]
        # Append stdout to .log, stderr to .json
        with open(log_file, "a") as log, open(json_file, "a") as js:
            subprocess.run(cmd, stdout=log, stderr=js, check=True)


# --- Main pipeline ---
index_dir.mkdir(parents=True, exist_ok=True)
results_dir.mkdir(parents=True, exist_ok=True)
tmp_dir.mkdir(parents=True, exist_ok=True)

print("\n=======================================================")
print("   STARTING SWEEP FOR m (MINIMIZER LENGTH) ")
print("=======================================================\n")

# --- k = 31 Sweep ---
build_project(max_k63=False)

for dataset in datasets:
    for current_m in m_sweeps_k31[dataset]:
        build_sshash(31, dataset, current_m)
        run_bench(31, dataset, current_m)

# --- k = 63 Sweep ---
build_project(max_k63=True)

for dataset in datasets:
    for current_m in m_sweeps_k63[dataset]:
        build_sshash(63, dataset, current_m)
        run_bench(63, dataset, current_m)

# Restore default compilation at the end
print("\nRestoring default compilation (max_k63=False)...")
build_project(max_k63=False)

print("\n All SSHash indexes built and benchmarked successfully across varying 'm' values. \n")
