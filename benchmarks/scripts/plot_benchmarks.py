#!/usr/bin/env python3
"""
Plot benchmark results from CSV files.

Usage:
    python plot_benchmarks.py results/bench_10M_5M_*.csv
    python plot_benchmarks.py results/bench_10M_5M_*.csv --output plots/

Requirements:
    pip install matplotlib pandas
"""

import argparse
import csv
import sys
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import pandas as pd
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False


def load_data(csv_files: list[str]) -> 'pd.DataFrame':
    """Load and combine CSV files into DataFrame."""
    dfs = []
    for f in csv_files:
        df = pd.read_csv(f)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def plot_speedup(df: 'pd.DataFrame', output_dir: Path):
    """Plot speedup comparison."""
    # Get GRIT results only
    grit_df = df[df['tool'] == 'grit'].copy()

    for dist in grit_df['distribution'].unique():
        dist_df = grit_df[grit_df['distribution'] == dist]

        fig, ax = plt.subplots(figsize=(10, 6))

        commands = dist_df['command'].tolist()
        speedups = dist_df['speedup_vs_bedtools'].astype(float).tolist()

        colors = ['#2ecc71' if s >= 5 else '#f39c12' if s >= 2 else '#e74c3c' for s in speedups]

        bars = ax.barh(commands, speedups, color=colors)
        ax.axvline(x=1, color='red', linestyle='--', label='bedtools baseline')

        ax.set_xlabel('Speedup (x faster than bedtools)')
        ax.set_title(f'GRIT vs bedtools Speedup ({dist.title()} Distribution)')
        ax.legend()

        # Add value labels
        for bar, val in zip(bars, speedups):
            ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{val:.1f}x', va='center')

        plt.tight_layout()
        output_path = output_dir / f'speedup_{dist}.png'
        plt.savefig(output_path, dpi=150)
        print(f"Saved: {output_path}")
        plt.close()


def plot_memory(df: 'pd.DataFrame', output_dir: Path):
    """Plot memory comparison."""
    for dist in df['distribution'].unique():
        dist_df = df[df['distribution'] == dist]

        fig, ax = plt.subplots(figsize=(12, 6))

        commands = dist_df[dist_df['tool'] == 'bedtools']['command'].unique()
        x = range(len(commands))
        width = 0.35

        bt_mem = []
        grit_mem = []
        for cmd in commands:
            bt = dist_df[(dist_df['command'] == cmd) & (dist_df['tool'] == 'bedtools')]['rss_mb']
            grit = dist_df[(dist_df['command'] == cmd) & (dist_df['tool'] == 'grit')]['rss_mb']
            bt_mem.append(float(bt.iloc[0]) if len(bt) > 0 else 0)
            grit_mem.append(float(grit.iloc[0]) if len(grit) > 0 else 0)

        ax.bar([i - width/2 for i in x], bt_mem, width, label='bedtools', color='#e74c3c')
        ax.bar([i + width/2 for i in x], grit_mem, width, label='GRIT', color='#2ecc71')

        ax.set_xlabel('Command')
        ax.set_ylabel('Peak Memory (MB)')
        ax.set_title(f'Memory Usage Comparison ({dist.title()} Distribution)')
        ax.set_xticks(x)
        ax.set_xticklabels(commands, rotation=45, ha='right')
        ax.legend()
        ax.set_yscale('log')

        plt.tight_layout()
        output_path = output_dir / f'memory_{dist}.png'
        plt.savefig(output_path, dpi=150)
        print(f"Saved: {output_path}")
        plt.close()


def main():
    parser = argparse.ArgumentParser(description="Plot benchmark results")
    parser.add_argument("csv_files", nargs="+", help="CSV files to process")
    parser.add_argument("--output", "-o", default="plots", help="Output directory")
    args = parser.parse_args()

    if not HAS_PLOTTING:
        print("Error: matplotlib and pandas required")
        print("Install with: pip install matplotlib pandas")
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_data(args.csv_files)

    plot_speedup(df, output_dir)
    plot_memory(df, output_dir)

    print(f"\nPlots saved to {output_dir}/")


if __name__ == "__main__":
    main()
