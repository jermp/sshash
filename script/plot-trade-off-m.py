#!/usr/bin/env python3

import os
import sys
import json
import re
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def parse_results(results_dir):
    """
    Traverses the results directory and extracts space and time metrics.
    Calculates bits/k-mer from bytes and num_kmers.
    Extracts dataset name from filenames.
    """
    results_dir = Path(results_dir)
    data = []

    # Regex to extract m and k from folder names
    m_pattern = re.compile(r'm(\d+)')
    k_pattern = re.compile(r'k(\d+)')

    for m_dir in results_dir.glob('m*'):
        m_match = m_pattern.search(m_dir.name)
        if not m_match: continue
        m_val = int(m_match.group(1))

        for k_dir in m_dir.glob('k*'):
            k_match = k_pattern.search(k_dir.name)
            if not k_match: continue
            k_val = int(k_match.group(1))

            for mode in ['regular', 'canon']:
                build_json = k_dir / f"{mode}-build.json"
                bench_json = k_dir / f"{mode}-bench.json"

                if not build_json.exists() or not bench_json.exists():
                    continue

                # 1. Parse Build JSON for Space (bits/k-mer)
                space_dict = {}
                with open(build_json, 'r') as f:
                    for line in f:
                        try:
                            j = json.loads(line)

                            # Extract dataset name from "/mnt/.../human.k31.eulertigs.fa.gz"
                            filename = os.path.basename(j.get("input_filename", ""))
                            ds = filename.split('.')[0] if filename else "unknown"

                            if "index_size_in_bytes" in j and "num_kmers" in j:
                                bytes_size = float(j["index_size_in_bytes"])
                                num_kmers = float(j["num_kmers"])
                                # Calculate bits per k-mer
                                bits_per_kmer = (bytes_size * 8.0) / num_kmers
                                space_dict[ds] = bits_per_kmer
                        except json.JSONDecodeError:
                            continue

                # 2. Parse Bench JSON for Query Time (ns/kmer)
                time_dict = {}
                count_dict = {}
                with open(bench_json, 'r') as f:
                    for line in f:
                        try:
                            j = json.loads(line)

                            # Extract dataset name from ".../human.k31.m17.sshash"
                            filename = os.path.basename(j.get("index_filename", ""))
                            ds = filename.split('.')[0] if filename else "unknown"

                            # Use positive lookup time
                            t_str = j.get("positive lookup (avg_nanosec_per_kmer)")
                            if t_str is not None:
                                t = float(t_str)
                                time_dict[ds] = time_dict.get(ds, 0.0) + t
                                count_dict[ds] = count_dict.get(ds, 0) + 1
                        except json.JSONDecodeError:
                            continue

                # 3. Combine and store
                for ds in space_dict.keys():
                    if ds in time_dict and count_dict[ds] > 0:
                        # Average the 3 benchmark runs
                        avg_time = time_dict[ds] / count_dict[ds]
                        data.append({
                            'Dataset': ds,
                            'k': k_val,
                            'm': m_val,
                            'Mode': mode,
                            'Space (bits/k-mer)': space_dict[ds],
                            'Query Time (ns/k-mer)': avg_time
                        })

    return pd.DataFrame(data)

def plot_tradeoff(df, output_img="tradeoff_plot.png"):
    """
    Generates a space-time trade-off plot.
    Different lines for datasets/k/modes, points vary by 'm'.
    """
    if df.empty:
        print("No data parsed! Please check the JSON keys in the script.")
        return

    # Enforce categorical order so the legend is populated exactly how we want:
    # Human before SE, and regular before canon
    df['Dataset'] = pd.Categorical(df['Dataset'], categories=['human', 'se'], ordered=True)
    df['Mode'] = pd.Categorical(df['Mode'], categories=['regular', 'canon'], ordered=True)
    df = df.sort_values(by=['Dataset', 'k', 'Mode'])

    plt.figure(figsize=(10, 8))
    plt.style.use('seaborn-v0_8-whitegrid')

    # Group by Dataset, k, and Mode with sort=False to preserve our categorical ordering
    groups = df.groupby(['Dataset', 'k', 'Mode'], sort=False)

    for (dataset, k, mode), group in groups:
        # Sort by m to make the line connect logically
        group = group.sort_values(by='m')

        label = f"{'Human' if dataset == 'human' else 'SE'} (k={k}, {mode})"

        # Color logic: Red for Human, Blue for SE. Darker if canonical.
        if dataset == 'human':
            color = 'firebrick' if mode == 'canon' else 'lightcoral'
        else: # se
            color = 'royalblue' if mode == 'canon' else 'lightskyblue'

        # Marker logic: Circle for k=31, Square for k=63
        marker = 'o' if k == 31 else 's'

        # Plot line and scatter (added linewidth=2.5 for thicker lines)
        plt.plot(group['Space (bits/k-mer)'], group['Query Time (ns/k-mer)'],
                 linestyle='-', color=color, alpha=0.7, linewidth=2.5)
        plt.scatter(group['Space (bits/k-mer)'], group['Query Time (ns/k-mer)'],
                    label=label, color=color, marker=marker, s=80, edgecolor='k', zorder=5)

        # Find min and max m for this group to filter annotations
        min_m = group['m'].min()
        max_m = group['m'].max()

        # Annotate points with 'm' values
        for _, row in group.iterrows():
            current_m = row['m']

            # If k=63, ONLY annotate if m is the smallest or largest in this sweep
            if row['k'] == 63 and current_m not in (min_m, max_m):
                continue

            plt.annotate(f"{int(current_m)}",
                         (row['Space (bits/k-mer)'], row['Query Time (ns/k-mer)']),
                         textcoords="offset points",
                         xytext=(5,5),
                         ha='left', fontsize=10)

    # Force x-axis ticks to be 1 unit (1 bit) apart
    plt.gca().xaxis.set_major_locator(MultipleLocator(1))

    # Make axis labels bold
    plt.xlabel("Index Space (bits/k-mer)", fontsize=12, fontweight='bold')
    plt.ylabel("Positive Lookup Time (ns/k-mer)", fontsize=12, fontweight='bold')
    plt.legend(title="Configuration", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    print(f"Plot saved successfully as '{output_img}'!")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 plot_tradeoff.py <path_to_results_dir>")
        sys.exit(1)

    results_directory = sys.argv[1]
    df = parse_results(results_directory)
    
    print("Extracted Data:")
    print(df.to_string(index=False))
    
    plot_tradeoff(df, output_img="sshash_tradeoff_m.png")
