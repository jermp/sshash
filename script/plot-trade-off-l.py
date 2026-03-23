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

    # Regex to extract l and k from folder names
    l_pattern = re.compile(r'l(\d+)')
    k_pattern = re.compile(r'k(\d+)')

    for l_dir in results_dir.glob('l*'):
        l_match = l_pattern.search(l_dir.name)
        if not l_match: continue
        l_val = int(l_match.group(1))

        for k_dir in l_dir.glob('k*'):
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

                            # Extract dataset name from ".../human.k31.l4.sshash"
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
                        # Average the benchmark runs
                        avg_time = time_dict[ds] / count_dict[ds]
                        data.append({
                            'Dataset': ds,
                            'k': k_val,
                            'l': l_val,
                            'Mode': mode,
                            'Space (bits/k-mer)': space_dict[ds],
                            'Query Time (ns/k-mer)': avg_time
                        })

    return pd.DataFrame(data)

def plot_tradeoff(df, output_img="tradeoff_plot_l.pdf"):
    """
    Generates a space-time trade-off plot.
    Different lines for datasets/k/modes, points vary by 'l'.
    """
    if df.empty:
        print("No data parsed! Please check the JSON keys in the script.")
        return

    # Enforce categorical order so the legend is populated exactly how we want:
    # Human before SE, and regular before canon
    df['Dataset'] = pd.Categorical(df['Dataset'], categories=['human', 'se'], ordered=True)
    df['Mode'] = pd.Categorical(df['Mode'], categories=['regular', 'canon'], ordered=True)
    df = df.sort_values(by=['Dataset', 'k', 'Mode'])

    plt.figure(figsize=(5, 10))
    plt.style.use('seaborn-v0_8-whitegrid')

    # Group by Dataset, k, and Mode with sort=False to preserve our categorical ordering
    groups = df.groupby(['Dataset', 'k', 'Mode'], sort=False)

    for (dataset, k, mode), group in groups:
        # Sort by l to make the line connect logically
        group = group.sort_values(by='l')

        label = f"{'Human' if dataset == 'human' else 'SE'} (k={k}, {mode})"
        
        # Color logic: Red for Human, Blue for SE. Darker if canonical.
        if dataset == 'human':
            color = 'firebrick' if mode == 'canon' else 'lightcoral'
        else: # se
            color = 'royalblue' if mode == 'canon' else 'lightskyblue'

        # Marker logic: Circle for k=31, Square for k=63
        marker = 'o' if k == 31 else 's'

        # Plot line and scatter (linewidth=2.5 for thicker lines)
        plt.plot(group['Space (bits/k-mer)'], group['Query Time (ns/k-mer)'],
                 linestyle='-', color=color, alpha=0.7, linewidth=2.5)
        plt.scatter(group['Space (bits/k-mer)'], group['Query Time (ns/k-mer)'],
                    label=label, color=color, marker=marker, s=80, edgecolor='k', zorder=5)

        # Annotate ALL points with 'l' values
        for _, row in group.iterrows():
            current_l = row['l']
                
            plt.annotate(f"{int(current_l)}",
                         (row['Space (bits/k-mer)'], row['Query Time (ns/k-mer)']),
                         textcoords="offset points",
                         xytext=(5,5),
                         ha='left', fontsize=10)

    # Force x-axis ticks to be 1 unit (1 bit) apart
    plt.gca().xaxis.set_major_locator(MultipleLocator(1))

    # Make axis labels bold
    plt.xlabel("Index Space (bits/k-mer)", fontsize=10
            # , fontweight='bold'
            )
    plt.ylabel("Positive Lookup Time (ns/k-mer)", fontsize=10
            # , fontweight='bold'
            )

    # Legend modifications: Move above plot, centered, multi-column
    plt.legend(title="Configuration",
               bbox_to_anchor=(0.5, 1.02),
               loc='lower center',
               ncol=2,
               borderaxespad=0.)

    plt.tight_layout()

    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    print(f"Plot saved successfully as '{output_img}'!")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 plot_tradeoff_l.py <path_to_results_dir>")
        sys.exit(1)

    results_directory = sys.argv[1]
    df = parse_results(results_directory)

    print("Extracted Data:")
    print(df.to_string(index=False))

    plot_tradeoff(df, output_img="sshash_tradeoff_l.pdf")
