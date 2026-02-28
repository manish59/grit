#!/usr/bin/env python3
"""
Performance benchmarks comparing pygrit Python API with grit CLI.

This script measures:
1. Execution time
2. Memory usage
3. Throughput (intervals/second)

Run with: python bench_performance.py
"""

import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pygrit


def generate_bed_file(path: Path, num_intervals: int, chrom_count: int = 5) -> None:
    """Generate a sorted BED file with the specified number of intervals."""
    with open(path, "w") as f:
        intervals_per_chrom = num_intervals // chrom_count
        for chrom_idx in range(1, chrom_count + 1):
            chrom = f"chr{chrom_idx}"
            for i in range(intervals_per_chrom):
                start = i * 1000
                end = start + 500
                f.write(f"{chrom}\t{start}\t{end}\n")


def run_cli_benchmark(cmd: list[str]) -> tuple[float, int]:
    """Run CLI command and return (time_seconds, memory_kb)."""
    # Use /usr/bin/time for memory measurement on macOS/Linux
    if sys.platform == "darwin":
        time_cmd = ["/usr/bin/time", "-l"] + cmd
        result = subprocess.run(time_cmd, capture_output=True, text=True)
        # Parse macOS time output for memory
        stderr = result.stderr
        time_line = stderr.split("\n")[0]
        elapsed = float(time_line.split()[0])
        # Find "maximum resident set size" line
        for line in stderr.split("\n"):
            if "maximum resident set size" in line:
                memory_kb = int(line.strip().split()[0]) // 1024
                break
        else:
            memory_kb = 0
    else:
        # Linux
        time_cmd = ["/usr/bin/time", "-v"] + cmd
        result = subprocess.run(time_cmd, capture_output=True, text=True)
        stderr = result.stderr
        elapsed = 0.0
        memory_kb = 0
        for line in stderr.split("\n"):
            if "Elapsed" in line:
                # Parse mm:ss.ss format
                time_str = line.split()[-1]
                parts = time_str.split(":")
                if len(parts) == 2:
                    elapsed = float(parts[0]) * 60 + float(parts[1])
                else:
                    elapsed = float(time_str)
            if "Maximum resident set size" in line:
                memory_kb = int(line.split()[-1])

    return elapsed, memory_kb


def run_python_benchmark(func, *args, **kwargs) -> tuple[float, list]:
    """Run Python function and return (time_seconds, result)."""
    import tracemalloc

    tracemalloc.start()
    start = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = time.perf_counter() - start
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return elapsed, peak // 1024, result


def format_size(kb: int) -> str:
    """Format size in KB to human readable."""
    if kb >= 1024:
        return f"{kb / 1024:.1f} MB"
    return f"{kb} KB"


def format_rate(intervals: int, seconds: float) -> str:
    """Format throughput as intervals/second."""
    if seconds == 0:
        return "N/A"
    rate = intervals / seconds
    if rate >= 1_000_000:
        return f"{rate / 1_000_000:.2f}M/s"
    if rate >= 1_000:
        return f"{rate / 1_000:.2f}K/s"
    return f"{rate:.2f}/s"


def check_grit_available() -> bool:
    """Check if grit CLI is available."""
    try:
        subprocess.run(["grit", "--version"], capture_output=True)
        return True
    except FileNotFoundError:
        return False


def main():
    print("=" * 70)
    print("pygrit Performance Benchmark")
    print("=" * 70)
    print()

    if not check_grit_available():
        print("WARNING: grit CLI not found. Skipping CLI comparisons.")
        cli_available = False
    else:
        cli_available = True

    # Test sizes
    sizes = [10_000, 100_000, 1_000_000]

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        for size in sizes:
            print(f"\n--- Benchmark: {size:,} intervals ---\n")

            # Generate test files
            a_file = tmpdir / f"a_{size}.bed"
            b_file = tmpdir / f"b_{size}.bed"
            generate_bed_file(a_file, size)
            generate_bed_file(b_file, size // 2)

            # Count overlaps for throughput calculation
            result = pygrit.intersect(str(a_file), str(b_file))
            overlap_count = len(result) if result else 0

            print(f"Input A: {size:,} intervals")
            print(f"Input B: {size // 2:,} intervals")
            print(f"Overlaps: {overlap_count:,}")
            print()

            # Benchmark: intersect
            print("INTERSECT:")

            # Python API
            py_time, py_mem, _ = run_python_benchmark(
                pygrit.intersect,
                str(a_file),
                str(b_file),
                output=str(tmpdir / "py_out.bed"),
            )
            print(f"  Python: {py_time:.3f}s, {format_size(py_mem)}, {format_rate(overlap_count, py_time)}")

            if cli_available:
                cli_time, cli_mem = run_cli_benchmark([
                    "grit", "intersect",
                    "--streaming", "--assume-sorted",
                    "-a", str(a_file),
                    "-b", str(b_file),
                ])
                print(f"  CLI:    {cli_time:.3f}s, {format_size(cli_mem)}, {format_rate(overlap_count, cli_time)}")
                if py_time > 0 and cli_time > 0:
                    ratio = py_time / cli_time
                    print(f"  Ratio:  {ratio:.2f}x (Python/CLI)")

            # Benchmark: merge
            print("\nMERGE:")

            py_time, py_mem, _ = run_python_benchmark(
                pygrit.merge,
                str(a_file),
                output=str(tmpdir / "py_merge.bed"),
            )
            print(f"  Python: {py_time:.3f}s, {format_size(py_mem)}")

            if cli_available:
                cli_time, cli_mem = run_cli_benchmark([
                    "grit", "merge",
                    "--assume-sorted",
                    "-i", str(a_file),
                ])
                print(f"  CLI:    {cli_time:.3f}s, {format_size(cli_mem)}")
                if py_time > 0 and cli_time > 0:
                    ratio = py_time / cli_time
                    print(f"  Ratio:  {ratio:.2f}x (Python/CLI)")

            # Benchmark: subtract
            print("\nSUBTRACT:")

            py_time, py_mem, _ = run_python_benchmark(
                pygrit.subtract,
                str(a_file),
                str(b_file),
                output=str(tmpdir / "py_sub.bed"),
            )
            print(f"  Python: {py_time:.3f}s, {format_size(py_mem)}")

            if cli_available:
                cli_time, cli_mem = run_cli_benchmark([
                    "grit", "subtract",
                    "--streaming", "--assume-sorted",
                    "-a", str(a_file),
                    "-b", str(b_file),
                ])
                print(f"  CLI:    {cli_time:.3f}s, {format_size(cli_mem)}")
                if py_time > 0 and cli_time > 0:
                    ratio = py_time / cli_time
                    print(f"  Ratio:  {ratio:.2f}x (Python/CLI)")

    print("\n" + "=" * 70)
    print("Benchmark complete")
    print("=" * 70)


if __name__ == "__main__":
    main()
