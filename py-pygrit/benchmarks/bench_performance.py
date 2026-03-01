#!/usr/bin/env python3
"""
Performance benchmarks comparing pygrit Python API with grit CLI.

This script measures:
1. Execution time
2. Memory usage
3. Throughput (intervals/second)

Run with: python bench_performance.py
"""

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


def generate_unsorted_bed_file(path: Path, num_intervals: int, chrom_count: int = 5) -> None:
    """Generate an unsorted BED file for sort benchmarks."""
    import random
    random.seed(42)
    intervals = []
    intervals_per_chrom = num_intervals // chrom_count
    for chrom_idx in range(1, chrom_count + 1):
        chrom = f"chr{chrom_idx}"
        for i in range(intervals_per_chrom):
            start = random.randint(0, 10_000_000)
            end = start + random.randint(100, 1000)
            intervals.append((chrom, start, end))
    random.shuffle(intervals)
    with open(path, "w") as f:
        for chrom, start, end in intervals:
            f.write(f"{chrom}\t{start}\t{end}\n")


def generate_genome_file(path: Path, chrom_count: int = 5) -> None:
    """Generate a genome file with chromosome sizes."""
    with open(path, "w") as f:
        for chrom_idx in range(1, chrom_count + 1):
            f.write(f"chr{chrom_idx}\t100000000\n")


def run_cli_benchmark(cmd: list[str]) -> tuple[float, int] | None:
    """Run CLI command and return (time_seconds, memory_kb), or None on error."""
    try:
        # Use /usr/bin/time for memory measurement on macOS/Linux
        if sys.platform == "darwin":
            time_cmd = ["/usr/bin/time", "-l"] + cmd
            result = subprocess.run(time_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                return None
            # Parse macOS time output for memory
            stderr = result.stderr
            time_line = stderr.split("\n")[0]
            try:
                elapsed = float(time_line.split()[0])
            except (ValueError, IndexError):
                return None
            # Find "maximum resident set size" line
            memory_kb = 0
            for line in stderr.split("\n"):
                if "maximum resident set size" in line:
                    try:
                        memory_kb = int(line.strip().split()[0]) // 1024
                    except (ValueError, IndexError):
                        pass
                    break
        else:
            # Linux
            time_cmd = ["/usr/bin/time", "-v"] + cmd
            result = subprocess.run(time_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                return None
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
    except Exception:
        return None


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


def benchmark_command(name: str, py_func, py_args, py_kwargs, cli_cmd: list[str] | None,
                       cli_available: bool, count: int = 0):
    """Run a benchmark for a single command."""
    print(f"\n{name}:")

    # Python API
    py_time, py_mem, _ = run_python_benchmark(py_func, *py_args, **py_kwargs)
    rate_str = f", {format_rate(count, py_time)}" if count > 0 else ""
    print(f"  Python: {py_time:.3f}s, {format_size(py_mem)}{rate_str}")

    if cli_available and cli_cmd:
        result = run_cli_benchmark(cli_cmd)
        if result is not None:
            cli_time, cli_mem = result
            rate_str = f", {format_rate(count, cli_time)}" if count > 0 else ""
            print(f"  CLI:    {cli_time:.3f}s, {format_size(cli_mem)}{rate_str}")
            if py_time > 0 and cli_time > 0:
                ratio = py_time / cli_time
                print(f"  Ratio:  {ratio:.2f}x (Python/CLI)")
        else:
            print("  CLI:    (command failed)")

    return py_time


def main():
    print("=" * 70)
    print("pygrit Performance Benchmark - All Commands")
    print("=" * 70)
    print()

    if not check_grit_available():
        print("WARNING: grit CLI not found. Skipping CLI comparisons.")
        cli_available = False
    else:
        cli_available = True

    # Test sizes - use smaller size for comprehensive test
    sizes = [50_000, 200_000]

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        for size in sizes:
            print(f"\n{'=' * 70}")
            print(f"Benchmark: {size:,} intervals")
            print("=" * 70)

            # Generate test files
            a_file = tmpdir / f"a_{size}.bed"
            b_file = tmpdir / f"b_{size}.bed"
            c_file = tmpdir / f"c_{size}.bed"
            unsorted_file = tmpdir / f"unsorted_{size}.bed"
            genome_file = tmpdir / "genome.txt"

            generate_bed_file(a_file, size)
            generate_bed_file(b_file, size // 2)
            generate_bed_file(c_file, size // 3)
            generate_unsorted_bed_file(unsorted_file, size)
            generate_genome_file(genome_file)

            # Count overlaps for throughput calculation
            result = pygrit.intersect(str(a_file), str(b_file))
            overlap_count = len(result) if result else 0

            print(f"\nInput A: {size:,} intervals")
            print(f"Input B: {size // 2:,} intervals")
            print(f"Input C: {size // 3:,} intervals")
            print(f"Overlaps (A∩B): {overlap_count:,}")

            # ===== Original Commands =====
            print("\n--- Original Commands ---")

            # INTERSECT
            benchmark_command(
                "INTERSECT",
                pygrit.intersect,
                [str(a_file), str(b_file)],
                {"output": str(tmpdir / "py_out.bed")},
                ["grit", "intersect", "--streaming", "--assume-sorted",
                 "-a", str(a_file), "-b", str(b_file)],
                cli_available,
                overlap_count,
            )

            # MERGE
            benchmark_command(
                "MERGE",
                pygrit.merge,
                [str(a_file)],
                {"output": str(tmpdir / "py_merge.bed")},
                ["grit", "merge", "--assume-sorted", "-i", str(a_file)],
                cli_available,
            )

            # SUBTRACT
            benchmark_command(
                "SUBTRACT",
                pygrit.subtract,
                [str(a_file), str(b_file)],
                {"output": str(tmpdir / "py_sub.bed")},
                ["grit", "subtract", "--streaming", "--assume-sorted",
                 "-a", str(a_file), "-b", str(b_file)],
                cli_available,
            )

            # COVERAGE
            benchmark_command(
                "COVERAGE",
                pygrit.coverage,
                [str(a_file), str(b_file)],
                {"output": str(tmpdir / "py_cov.bed")},
                ["grit", "coverage", "--assume-sorted",
                 "-a", str(a_file), "-b", str(b_file)],
                cli_available,
            )

            # CLOSEST
            benchmark_command(
                "CLOSEST",
                pygrit.closest,
                [str(a_file), str(b_file)],
                {"output": str(tmpdir / "py_closest.bed")},
                ["grit", "closest", "--assume-sorted",
                 "-a", str(a_file), "-b", str(b_file)],
                cli_available,
            )

            # WINDOW
            benchmark_command(
                "WINDOW",
                pygrit.window,
                [str(a_file), str(b_file)],
                {"output": str(tmpdir / "py_window.bed"), "window": 1000},
                ["grit", "window", "--assume-sorted",
                 "-a", str(a_file), "-b", str(b_file), "-w", "1000"],
                cli_available,
            )

            # ===== New Commands =====
            print("\n--- New Commands ---")

            # SORT
            benchmark_command(
                "SORT",
                pygrit.sort,
                [str(unsorted_file)],
                {"output": str(tmpdir / "py_sorted.bed")},
                ["grit", "sort", "-i", str(unsorted_file)],
                cli_available,
                size,
            )

            # SLOP
            benchmark_command(
                "SLOP",
                pygrit.slop,
                [str(a_file), str(genome_file)],
                {"both": 100.0, "output": str(tmpdir / "py_slop.bed")},
                ["grit", "slop", "-i", str(a_file), "-g", str(genome_file), "-b", "100"],
                cli_available,
                size,
            )

            # COMPLEMENT
            benchmark_command(
                "COMPLEMENT",
                pygrit.complement,
                [str(a_file), str(genome_file)],
                {"output": str(tmpdir / "py_complement.bed")},
                ["grit", "complement", "-i", str(a_file), "-g", str(genome_file)],
                cli_available,
            )

            # GENOMECOV
            benchmark_command(
                "GENOMECOV",
                pygrit.genomecov,
                [str(a_file), str(genome_file)],
                {"output": str(tmpdir / "py_genomecov.txt"), "bg": True},
                ["grit", "genomecov", "-i", str(a_file), "-g", str(genome_file), "-bg"],
                cli_available,
            )

            # JACCARD
            benchmark_command(
                "JACCARD",
                pygrit.jaccard,
                [str(a_file), str(b_file)],
                {"output": str(tmpdir / "py_jaccard.txt")},
                ["grit", "jaccard", "-a", str(a_file), "-b", str(b_file)],
                cli_available,
            )

            # MULTIINTER
            benchmark_command(
                "MULTIINTER",
                pygrit.multiinter,
                [[str(a_file), str(b_file), str(c_file)]],
                {"output": str(tmpdir / "py_multiinter.bed")},
                ["grit", "multiinter", "-i", str(a_file), str(b_file), str(c_file)],
                cli_available,
            )

    print("\n" + "=" * 70)
    print("Benchmark complete")
    print("=" * 70)


if __name__ == "__main__":
    main()
