#!/usr/bin/env python3
import os
import subprocess
import sys
from pathlib import Path

# ------------------------------
#   Argument parsing
# ------------------------------
if len(sys.argv) < 3:
    print("Usage: python3 bench.py <log_label> <input_index_dir>")
    sys.exit(1)

log_label = sys.argv[1]
index_dir = Path(sys.argv[2]).resolve()

# ------------------------------
#   Global configuration
# ------------------------------
results_dir = Path(f"results-{log_label}")

datasets = [
    "cod", "kestrel", "human", "ncbi-virus", "se", "hprc"
]

# ------------------------------
#   Utility functions
# ------------------------------
def run_cmd(cmd, cwd=None):
    """Run a shell command and print it."""
    print(f"[RUN] {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)


def build_project(max_k63: bool):
    """Run cmake + make with max_k63 = True/False."""
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


def run_bench(k, runs = 3):
    """Run SSHASH benchmark for all datasets."""
    out_dir = results_dir / f"k{k}"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_file = out_dir / "bench.log"
    json_file = out_dir / "bench.json"

    for dataset in datasets:
        index_path = index_dir / f"{dataset}.k{k}.sshash"

        print(f"\n>>> Benchmarking {dataset} (k={k})\n")
        for i in range(runs):
            print(f"  ==> run {i+1}/{runs}")
            cmd = ["./sshash", "bench", "-i", str(index_path)]
            # Append stdout to .log, stderr to .json
            with open(log_file, "a") as log, open(json_file, "a") as js:
                subprocess.run(cmd, stdout=log, stderr=js, check=True)

# ------------------------------
#   Prepare directories
# ------------------------------
results_dir.mkdir(parents=True, exist_ok=True)
(index_dir).mkdir(parents=True, exist_ok=True)

(results_dir / "k31").mkdir(exist_ok=True)
(results_dir / "k63").mkdir(exist_ok=True)

# ------------------------------
#   Run benchmarks
# ------------------------------

# --- Build for k=31 ---
build_project(max_k63=False)
run_bench(31)

# --- Build for k=63 ---
build_project(max_k63=True)
run_bench(63)

# --- Restore to default ---
build_project(max_k63=False)

print("\n All SSHash benchmark runs completed successfully. \n")
