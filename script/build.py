#!/usr/bin/env python3

import os
import subprocess
import sys
from pathlib import Path

if len(sys.argv) < 3:
    print("Usage: run_sshash_builds.py <log_label> <output_index_dir>")
    sys.exit(1)

log_label = sys.argv[1]
index_dir = Path(sys.argv[2]).resolve()

prefix = Path("/mnt/hd2/pibiri/DNA")
tmp_dir = prefix / "tmp_dir"
results_dir = Path(f"results-{log_label}")
threads = 64
g = 16

datasets = [
    "cod", "kestrel", "human", "axolotl", "hprc",
    "ec", "se", "ncbi-virus", "jgi_fungi.batch-0"
]

m_values_k31 = {
    "cod": 20, "kestrel": 20, "human": 21, "axolotl": 21,
    "hprc": 21, "ec": 21, "se": 21, "ncbi-virus": 19, "jgi_fungi.batch-0": 21
}

m_values_k63 = {
    "cod": 24, "kestrel": 24, "human": 25, "axolotl": 25,
    "hprc": 31, "ec": 31, "se": 31, "ncbi-virus": 23, "jgi_fungi.batch-0": 25
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


def build_sshash(k, canonical_flag, mode, m_values):
    """Perform sshash builds for a given k, mode, and canonical flag."""
    mode_dir = results_dir / f"k{k}"
    mode_dir.mkdir(parents=True, exist_ok=True)

    log_file = mode_dir / f"{mode}-build.log"
    json_file = mode_dir / f"{mode}-build.json"
    time_file = mode_dir / f"{mode}-build.time.log"

    for dataset in datasets:
        m_val = m_values[dataset]
        input_file = prefix / f"eulertigs/{dataset}.k{k}.eulertigs.fa.gz"
        output_file = index_dir / f"{dataset}.k{k}"
        if canonical_flag:
            output_file = output_file.with_suffix(".canon")

        print(f"\n>>> Building {dataset} (k={k}, m={m_val}, mode={mode})\n")

        # Clean tmp directory
        subprocess.run(["rm", "-rf", str(tmp_dir / "*")], shell=True)

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
        if canonical_flag:
            cmd.append("--canonical")

        # Append stdout to .log, stderr to .json
        with open(log_file, "a") as log, open(json_file, "a") as js:
            subprocess.run(cmd, stdout=log, stderr=js, check=True)


# --- Main pipeline ---
index_dir.mkdir(parents=True, exist_ok=True)
results_dir.mkdir(parents=True, exist_ok=True)

# k = 31
build_project(max_k63=False)
build_sshash(31, canonical_flag=False, mode="regular", m_values=m_values_k31)
build_sshash(31, canonical_flag=True, mode="canon", m_values=m_values_k31)

# k = 63
build_project(max_k63=True)
build_sshash(63, canonical_flag=False, mode="regular", m_values=m_values_k63)
build_sshash(63, canonical_flag=True, mode="canon", m_values=m_values_k63)

# rebuild back to default
build_project(max_k63=False)

print("\n All SSHash indexes built successfully.\n")
