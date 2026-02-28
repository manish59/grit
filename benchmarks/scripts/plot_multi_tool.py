#!/usr/bin/env python3
"""
Multi-Tool Benchmark Visualization

Generates comparison charts for GRIT vs other genomic tools.

Usage:
    python plot_multi_tool.py results.csv -o output_dir/
    python plot_multi_tool.py results.csv --format svg
"""

import argparse
import csv
import sys
from pathlib import Path
from collections import defaultdict

# Try to import matplotlib, but provide useful output without it
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# Tool colors - consistent across all charts
TOOL_COLORS = {
    "grit": "#2ecc71",       # Green - our tool
    "bedtools": "#3498db",   # Blue
    "bedops": "#9b59b6",     # Purple
    "granges": "#e74c3c",    # Red
    "pyranges": "#f39c12",   # Orange
    "polars_bio": "#1abc9c", # Teal
}

TOOL_ORDER = ["grit", "bedtools", "bedops", "granges", "pyranges", "polars_bio"]


def load_results(csv_path: str) -> dict:
    """Load benchmark results from CSV."""
    results = defaultdict(dict)

    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            op = row["operation"]
            tool = row["tool"]
            time_s = row["time_s"]
            mem_mb = row["memory_mb"]
            speedup = row.get("speedup_vs_bedtools", "NA")

            if time_s != "NA":
                results[op][tool] = {
                    "time": float(time_s),
                    "memory": float(mem_mb) if mem_mb != "NA" else None,
                    "speedup": float(speedup) if speedup != "NA" else None,
                }

    return dict(results)


def print_text_table(results: dict):
    """Print a text-based comparison table."""
    # Collect all tools that have results
    all_tools = set()
    for op_results in results.values():
        all_tools.update(op_results.keys())

    # Order tools
    tools = [t for t in TOOL_ORDER if t in all_tools]

    # Print header
    print("\n" + "=" * 80)
    print("MULTI-TOOL BENCHMARK RESULTS")
    print("=" * 80)

    # Time comparison
    print("\n### Execution Time (seconds) ###\n")
    header = f"{'Operation':<12}"
    for tool in tools:
        header += f" | {tool:>10}"
    print(header)
    print("-" * len(header))

    for op, op_results in sorted(results.items()):
        row = f"{op:<12}"
        for tool in tools:
            if tool in op_results:
                time = op_results[tool]["time"]
                row += f" | {time:>10.3f}"
            else:
                row += f" | {'-':>10}"
        print(row)

    # Speedup comparison
    print("\n### Speedup vs bedtools ###\n")
    header = f"{'Operation':<12}"
    for tool in tools:
        if tool != "bedtools":
            header += f" | {tool:>10}"
    print(header)
    print("-" * len(header))

    for op, op_results in sorted(results.items()):
        row = f"{op:<12}"
        bedtools_time = op_results.get("bedtools", {}).get("time")

        for tool in tools:
            if tool == "bedtools":
                continue
            if tool in op_results and bedtools_time:
                speedup = bedtools_time / op_results[tool]["time"]
                row += f" | {speedup:>9.2f}x"
            else:
                row += f" | {'-':>10}"
        print(row)

    # Memory comparison
    print("\n### Peak Memory (MB) ###\n")
    header = f"{'Operation':<12}"
    for tool in tools:
        header += f" | {tool:>10}"
    print(header)
    print("-" * len(header))

    for op, op_results in sorted(results.items()):
        row = f"{op:<12}"
        for tool in tools:
            if tool in op_results and op_results[tool]["memory"]:
                mem = op_results[tool]["memory"]
                row += f" | {mem:>10.1f}"
            else:
                row += f" | {'-':>10}"
        print(row)

    print("\n" + "=" * 80)


def generate_markdown_table(results: dict, output_path: str):
    """Generate a markdown comparison table."""
    all_tools = set()
    for op_results in results.values():
        all_tools.update(op_results.keys())
    tools = [t for t in TOOL_ORDER if t in all_tools]

    with open(output_path, "w") as f:
        f.write("# Multi-Tool Benchmark Results\n\n")

        # Time table
        f.write("## Execution Time (seconds)\n\n")
        f.write("| Operation |")
        for tool in tools:
            f.write(f" {tool} |")
        f.write("\n")
        f.write("|-----------|")
        for _ in tools:
            f.write("-------:|")
        f.write("\n")

        for op, op_results in sorted(results.items()):
            f.write(f"| {op} |")
            for tool in tools:
                if tool in op_results:
                    time = op_results[tool]["time"]
                    f.write(f" {time:.3f} |")
                else:
                    f.write(" - |")
            f.write("\n")

        # Speedup table
        f.write("\n## Speedup vs bedtools\n\n")
        f.write("| Operation |")
        for tool in tools:
            if tool != "bedtools":
                f.write(f" {tool} |")
        f.write("\n")
        f.write("|-----------|")
        for tool in tools:
            if tool != "bedtools":
                f.write("-------:|")
        f.write("\n")

        for op, op_results in sorted(results.items()):
            f.write(f"| {op} |")
            bedtools_time = op_results.get("bedtools", {}).get("time")
            for tool in tools:
                if tool == "bedtools":
                    continue
                if tool in op_results and bedtools_time:
                    speedup = bedtools_time / op_results[tool]["time"]
                    # Highlight GRIT speedups
                    if tool == "grit" and speedup > 1:
                        f.write(f" **{speedup:.2f}x** |")
                    else:
                        f.write(f" {speedup:.2f}x |")
                else:
                    f.write(" - |")
            f.write("\n")

        # Memory table
        f.write("\n## Peak Memory (MB)\n\n")
        f.write("| Operation |")
        for tool in tools:
            f.write(f" {tool} |")
        f.write("\n")
        f.write("|-----------|")
        for _ in tools:
            f.write("-------:|")
        f.write("\n")

        for op, op_results in sorted(results.items()):
            f.write(f"| {op} |")
            for tool in tools:
                if tool in op_results and op_results[tool]["memory"]:
                    mem = op_results[tool]["memory"]
                    f.write(f" {mem:.1f} |")
                else:
                    f.write(" - |")
            f.write("\n")

    print(f"Markdown table saved to: {output_path}")


def plot_time_comparison(results: dict, output_dir: Path, fmt: str = "png"):
    """Generate time comparison bar chart."""
    if not HAS_MATPLOTLIB:
        print("matplotlib not installed, skipping chart generation")
        return

    operations = sorted(results.keys())
    all_tools = set()
    for op_results in results.values():
        all_tools.update(op_results.keys())
    tools = [t for t in TOOL_ORDER if t in all_tools]

    fig, ax = plt.subplots(figsize=(12, 6))

    x = range(len(operations))
    width = 0.12
    multiplier = 0

    for tool in tools:
        times = []
        for op in operations:
            if tool in results[op]:
                times.append(results[op][tool]["time"])
            else:
                times.append(0)

        offset = width * multiplier
        bars = ax.bar(
            [i + offset for i in x],
            times,
            width,
            label=tool,
            color=TOOL_COLORS.get(tool, "#888888"),
        )
        multiplier += 1

    ax.set_xlabel("Operation")
    ax.set_ylabel("Time (seconds)")
    ax.set_title("Execution Time Comparison")
    ax.set_xticks([i + width * (len(tools) - 1) / 2 for i in x])
    ax.set_xticklabels(operations)
    ax.legend(loc="upper right")
    ax.set_yscale("log")

    plt.tight_layout()
    output_path = output_dir / f"time_comparison.{fmt}"
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Time comparison chart saved to: {output_path}")


def plot_speedup_chart(results: dict, output_dir: Path, fmt: str = "png"):
    """Generate speedup comparison chart."""
    if not HAS_MATPLOTLIB:
        return

    operations = sorted(results.keys())
    all_tools = set()
    for op_results in results.values():
        all_tools.update(op_results.keys())
    tools = [t for t in TOOL_ORDER if t in all_tools and t != "bedtools"]

    fig, ax = plt.subplots(figsize=(12, 6))

    x = range(len(operations))
    width = 0.15
    multiplier = 0

    for tool in tools:
        speedups = []
        for op in operations:
            bedtools_time = results[op].get("bedtools", {}).get("time")
            tool_time = results[op].get(tool, {}).get("time")
            if bedtools_time and tool_time:
                speedups.append(bedtools_time / tool_time)
            else:
                speedups.append(0)

        offset = width * multiplier
        bars = ax.bar(
            [i + offset for i in x],
            speedups,
            width,
            label=tool,
            color=TOOL_COLORS.get(tool, "#888888"),
        )
        multiplier += 1

    ax.axhline(y=1, color="gray", linestyle="--", linewidth=1, label="bedtools baseline")
    ax.set_xlabel("Operation")
    ax.set_ylabel("Speedup (higher is better)")
    ax.set_title("Speedup vs bedtools")
    ax.set_xticks([i + width * (len(tools) - 1) / 2 for i in x])
    ax.set_xticklabels(operations)
    ax.legend(loc="upper right")

    plt.tight_layout()
    output_path = output_dir / f"speedup_comparison.{fmt}"
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Speedup chart saved to: {output_path}")


def plot_memory_comparison(results: dict, output_dir: Path, fmt: str = "png"):
    """Generate memory usage comparison chart."""
    if not HAS_MATPLOTLIB:
        return

    operations = sorted(results.keys())
    all_tools = set()
    for op_results in results.values():
        all_tools.update(op_results.keys())
    tools = [t for t in TOOL_ORDER if t in all_tools]

    fig, ax = plt.subplots(figsize=(12, 6))

    x = range(len(operations))
    width = 0.12
    multiplier = 0

    for tool in tools:
        memories = []
        for op in operations:
            mem = results[op].get(tool, {}).get("memory")
            memories.append(mem if mem else 0)

        offset = width * multiplier
        bars = ax.bar(
            [i + offset for i in x],
            memories,
            width,
            label=tool,
            color=TOOL_COLORS.get(tool, "#888888"),
        )
        multiplier += 1

    ax.set_xlabel("Operation")
    ax.set_ylabel("Peak Memory (MB)")
    ax.set_title("Memory Usage Comparison")
    ax.set_xticks([i + width * (len(tools) - 1) / 2 for i in x])
    ax.set_xticklabels(operations)
    ax.legend(loc="upper right")
    ax.set_yscale("log")

    plt.tight_layout()
    output_path = output_dir / f"memory_comparison.{fmt}"
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Memory comparison chart saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate multi-tool benchmark visualizations")
    parser.add_argument("csv_file", help="Path to benchmark results CSV")
    parser.add_argument("-o", "--output", default=".", help="Output directory for charts")
    parser.add_argument("--format", choices=["png", "svg", "pdf"], default="png", help="Output format")
    parser.add_argument("--text", action="store_true", help="Print text table only")
    parser.add_argument("--markdown", action="store_true", help="Generate markdown table")

    args = parser.parse_args()

    if not Path(args.csv_file).exists():
        print(f"Error: CSV file not found: {args.csv_file}", file=sys.stderr)
        sys.exit(1)

    results = load_results(args.csv_file)

    if not results:
        print("Error: No results found in CSV file", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Always print text table
    print_text_table(results)

    if args.markdown:
        generate_markdown_table(results, output_dir / "results.md")

    if not args.text:
        if HAS_MATPLOTLIB:
            plot_time_comparison(results, output_dir, args.format)
            plot_speedup_chart(results, output_dir, args.format)
            plot_memory_comparison(results, output_dir, args.format)
        else:
            print("\nNote: Install matplotlib for chart generation: pip install matplotlib")


if __name__ == "__main__":
    main()
