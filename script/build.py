#!/usr/bin/env python3

import os
import subprocess
import sys
from pathlib import Path

if len(sys.argv) < 4:
    print("Usage: python3 build.py <log_label> <input_datasets_dir> <output_index_dir>")
    sys.exit(1)

log_label = sys.argv[1]
datasets_dir = Path(sys.argv[2]).resolve()
index_dir = Path(sys.argv[3]).resolve()
tmp_dir = datasets_dir / "tmp_dir"
results_dir = Path(f"results-{log_label}")
threads = 64
g = 16

datasets = [
    "cod", "kestrel", "human", "ncbi-virus", "se", "hprc"
]

m_values_k31 = {
    "cod": 20, "kestrel": 20, "human": 21, "ncbi-virus": 19, "se": 21, "hprc": 21
}

m_values_k63 = {
    "cod": 24, "kestrel": 24, "human": 25, "ncbi-virus": 23, "se": 31, "hprc": 31
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


def build_sshash(k, canonical, m_values):
    mode_dir = results_dir / f"k{k}"
    mode_dir.mkdir(parents=True, exist_ok=True)

    mode = "canon" if canonical else "regular"
    log_file = mode_dir / f"{mode}-build.log"
    json_file = mode_dir / f"{mode}-build.json"
    time_file = mode_dir / f"{mode}-build.time.log"

    for dataset in datasets:
        m_val = m_values[dataset]
        input_file = datasets_dir / f"{dataset}.k{k}.eulertigs.fa.gz"
        output_file = index_dir / f"{dataset}.k{k}"
        if canonical:
            output_file = str(output_file) + ".canon"

        print(f"\n>>> Building {dataset} (k={k}, m={m_val}, mode={mode})\n")

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


# --- Main pipeline ---
index_dir.mkdir(parents=True, exist_ok=True)
results_dir.mkdir(parents=True, exist_ok=True)

# k = 31
build_project(max_k63=False)
build_sshash(31, False, m_values_k31)
build_sshash(31, True, m_values_k31)

# k = 63
build_project(max_k63=True)
build_sshash(63, False, m_values_k63)
build_sshash(63, True, m_values_k63)

# rebuild back to default
build_project(max_k63=False)

print("\n All SSHash indexes built successfully. \n")
