#!/usr/bin/env python3
"""
Generate Markdown performance tables from benchmark CSV results.

Usage:
    python generate_markdown_table.py results/bench_10M_5M_*.csv
    python generate_markdown_table.py results/bench_10M_5M_*.csv --output README_PERF.md
"""

import argparse
import csv
import sys
from pathlib import Path
from collections import defaultdict


def parse_csv(filepath: str) -> list[dict]:
    """Parse benchmark CSV file."""
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        return list(reader)


def generate_markdown_table(rows: list[dict], distribution: str = "uniform") -> str:
    """Generate markdown table from benchmark results."""
    # Filter by distribution
    filtered = [r for r in rows if r.get('distribution', 'uniform') == distribution]

    # Group by command
    commands = defaultdict(dict)
    for row in filtered:
        cmd = row['command']
        tool = row['tool']
        commands[cmd][tool] = row

    # Build table
    lines = []
    lines.append(f"### {distribution.title()} Distribution\n")
    lines.append("| Command | bedtools | GRIT | Speedup | BT Memory | GRIT Memory |")
    lines.append("|---------|----------|------|---------|-----------|-------------|")

    for cmd in sorted(commands.keys()):
        tools = commands[cmd]
        bt = tools.get('bedtools', {})
        grit = tools.get('grit', {})

        bt_time = bt.get('time_s', 'N/A')
        grit_time = grit.get('time_s', 'N/A')
        speedup = grit.get('speedup_vs_bedtools', 'N/A')
        bt_mem = bt.get('rss_mb', 'N/A')
        grit_mem = grit.get('rss_mb', 'N/A')
        correctness = grit.get('correctness', 'N/A')

        # Format values
        if bt_time != 'N/A':
            bt_time = f"{float(bt_time):.2f}s"
        if grit_time != 'N/A':
            grit_time = f"{float(grit_time):.2f}s"
        if speedup != 'N/A':
            speedup = f"**{float(speedup):.1f}x**"
        if bt_mem != 'N/A':
            bt_mem = format_memory(float(bt_mem))
        if grit_mem != 'N/A':
            grit_mem = format_memory(float(grit_mem))

        status = "PASS" if correctness == "PASS" else "FAIL"
        lines.append(f"| {cmd} | {bt_time} | {grit_time} | {speedup} | {bt_mem} | {grit_mem} |")

    return "\n".join(lines)


def format_memory(mb: float) -> str:
    """Format memory in human-readable format."""
    if mb >= 1024:
        return f"{mb/1024:.1f} GB"
    return f"{mb:.0f} MB"


def main():
    parser = argparse.ArgumentParser(description="Generate markdown tables from benchmark CSV")
    parser.add_argument("csv_files", nargs="+", help="CSV files to process")
    parser.add_argument("--output", "-o", help="Output markdown file")
    args = parser.parse_args()

    all_rows = []
    for filepath in args.csv_files:
        all_rows.extend(parse_csv(filepath))

    output_lines = ["# GRIT Benchmark Results\n"]
    output_lines.append("Generated from benchmark CSV files.\n")

    # Generate tables for each distribution
    for dist in ["uniform", "clustered"]:
        table = generate_markdown_table(all_rows, dist)
        if table.count('\n') > 3:  # Has data
            output_lines.append(table)
            output_lines.append("")

    output = "\n".join(output_lines)

    if args.output:
        Path(args.output).write_text(output)
        print(f"Written to {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
