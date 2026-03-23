#!/usr/bin/env python3

import os
import re
import subprocess
import sys
from pathlib import Path

if len(sys.argv) < 4:
    print("Usage: python3 sweep-min-l.py <log_label> <input_datasets_dir> <output_index_dir>")
    sys.exit(1)

log_label = sys.argv[1]
datasets_dir = Path(sys.argv[2]).resolve()
index_dir = Path(sys.argv[3]).resolve()
tmp_dir = datasets_dir / "tmp_dir"
results_dir = Path(f"results-{log_label}")
threads = 16
g = 16

# 1. Target only human and se datasets
datasets = [
    "human", "se"
]

m_values_k31 = {
    "human": 21, "se": 21
}

m_values_k63 = {
    "human": 25, "se": 31
}

# Values of min_l to benchmark
l_values_to_test = [4, 5, 6, 7, 8]

# Assuming script is run from the `build` directory
constants_hpp_path = Path("../include/constants.hpp").resolve()

# --- Utilities ---
def run_cmd(cmd, cwd=None, append_to=None):
    print(f"[RUN] {' '.join(cmd)}")
    if append_to:
        with open(append_to, "a") as f:
            subprocess.run(cmd, cwd=cwd, stdout=f, stderr=f, check=True)
    else:
        subprocess.run(cmd, cwd=cwd, check=True)

def update_constants_hpp(min_l):
    """
    Updates the min_l and max_l values in include/constants.hpp.
    max_l is set to min_l + 7 to satisfy the static_assert constraints.
    """
    if not constants_hpp_path.exists():
        print(f"Error: Could not find {constants_hpp_path}")
        sys.exit(1)
        
    with open(constants_hpp_path, 'r') as f:
        content = f.read()

    # Update min_l
    content = re.sub(r'constexpr uint64_t min_l = \d+;', 
                     f'constexpr uint64_t min_l = {min_l};', 
                     content)
    
    # Update max_l safely
    max_l = min_l + 7
    content = re.sub(r'constexpr uint64_t max_l = \d+;', 
                     f'constexpr uint64_t max_l = {max_l};', 
                     content)

    with open(constants_hpp_path, 'w') as f:
        f.write(content)
        
    print(f"\n[CONFIG] Updated {constants_hpp_path.name} -> min_l = {min_l}, max_l = {max_l}")

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

def build_sshash(k, canonical, m_values, l_val):
    # Differentiate results dir by l_val
    mode_dir = results_dir / f"l{l_val}" / f"k{k}"
    mode_dir.mkdir(parents=True, exist_ok=True)

    mode = "canon" if canonical else "regular"
    log_file = mode_dir / f"{mode}-build.log"
    json_file = mode_dir / f"{mode}-build.json"
    time_file = mode_dir / f"{mode}-build.time.log"

    for dataset in datasets:
        m_val = m_values[dataset]
        input_file = datasets_dir / f"{dataset}.k{k}.eulertigs.fa.gz"
        
        # 2. Append l_val to the output filename
        output_file = index_dir / f"{dataset}.k{k}.l{l_val}"
        if canonical:
            output_file = str(output_file) + ".canon"

        print(f"\n>>> Building {dataset} (k={k}, m={m_val}, l={l_val}, mode={mode})\n")

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
        if canonical:
            cmd.append("--canonical")

        # Append stdout to .log, stderr to .json
        with open(log_file, "a") as log, open(json_file, "a") as js:
            subprocess.run(cmd, stdout=log, stderr=js, check=True)

def run_bench(k, canonical, l_val, runs=3):
    """Run SSHASH benchmark for all datasets."""
    mode = "canon" if canonical else "regular"
    
    # Store results in the specific l_val / k folder
    out_dir = results_dir / f"l{l_val}" / f"k{k}"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_file = out_dir / f"{mode}-bench.log"
    json_file = out_dir / f"{mode}-bench.json"

    for dataset in datasets:
        # Match the new index naming scheme that includes l_val
        suffix = f".k{k}.l{l_val}.canon.sshash" if canonical else f".k{k}.l{l_val}.sshash"
        index_path = index_dir / f"{dataset}{suffix}"

        print(f"\n>>> Benchmarking {dataset} (k={k}, l={l_val}, mode={mode})\n")
        for i in range(runs):
            print(f"  ==> run {i+1}/{runs}")
            cmd = ["./sshash", "bench", "-i", str(index_path)]
            # Append stdout to .log, stderr to .json
            with open(log_file, "a") as log, open(json_file, "a") as js:
                subprocess.run(cmd, stdout=log, stderr=js, check=True)


# --- Main pipeline ---
index_dir.mkdir(parents=True, exist_ok=True)
results_dir.mkdir(parents=True, exist_ok=True)

for current_l in l_values_to_test:
    print(f"\n=======================================================")
    print(f"   STARTING SWEEP FOR min_l = {current_l}")
    print(f"=======================================================\n")
    
    # Update the header file
    update_constants_hpp(current_l)
    
    # Build and benchmark for k = 31 (Regular)
    build_project(max_k63=False)
    build_sshash(31, False, m_values_k31, current_l)
    run_bench(31, False, current_l)
    
    # Build and benchmark for k = 31 (Canonical)
    build_sshash(31, True, m_values_k31, current_l)
    run_bench(31, True, current_l)

    # Build and benchmark for k = 63 (Regular)
    build_project(max_k63=True)
    build_sshash(63, False, m_values_k63, current_l)
    run_bench(63, False, current_l)
    
    # Build and benchmark for k = 63 (Canonical)
    build_sshash(63, True, m_values_k63, current_l)
    run_bench(63, True, current_l)

# Restore default constants file at the end
print("\nRestoring default constants.hpp (min_l = 6)...")
update_constants_hpp(6)
build_project(max_k63=False)

print("\n All SSHash indexes built and benchmarked successfully across varying 'min_l' values. \n")