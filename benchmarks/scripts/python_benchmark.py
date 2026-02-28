#!/usr/bin/env python3
"""
Python Benchmark Wrapper for pyranges and polars-bio

Provides a unified interface for benchmarking genomic interval operations
using Python libraries, comparable to CLI tool benchmarks.

Usage:
    python python_benchmark.py <library> <operation> <file_a> [file_b] [options]

Examples:
    python python_benchmark.py pyranges intersect a.bed b.bed -o output.bed
    python python_benchmark.py polars-bio merge a.bed -o output.bed
    python python_benchmark.py pyranges closest a.bed b.bed -o output.bed
"""

import argparse
import sys
import time
import tracemalloc
from pathlib import Path


def load_bed_pyranges(filepath: str):
    """Load a BED file into a PyRanges object."""
    import pyranges as pr

    # PyRanges expects specific column names
    return pr.read_bed(filepath)


def load_bed_polars_bio(filepath: str):
    """Load a BED file for polars-bio operations."""
    import polars as pl

    # Read as standard BED format
    df = pl.read_csv(
        filepath,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start", "end"] + [f"col{i}" for i in range(3, 12)],
        truncate_ragged_lines=True,
    )
    # Keep only first 3 columns for basic operations
    return df.select(["chrom", "start", "end"])


def save_bed_pyranges(pr_obj, filepath: str):
    """Save PyRanges object to BED file."""
    pr_obj.to_bed(filepath)


def save_bed_polars(df, filepath: str):
    """Save Polars DataFrame to BED file."""
    df.write_csv(filepath, separator="\t", include_header=False)


# =============================================================================
# PyRanges Operations
# =============================================================================

def pyranges_intersect(file_a: str, file_b: str, output: str, **kwargs):
    """Find overlapping intervals between two BED files."""
    import pyranges as pr

    a = pr.read_bed(file_a)
    b = pr.read_bed(file_b)
    result = a.join(b)
    result.to_bed(output)
    return len(result)


def pyranges_merge(file_a: str, output: str, **kwargs):
    """Merge overlapping intervals."""
    import pyranges as pr

    a = pr.read_bed(file_a)
    result = a.merge()
    result.to_bed(output)
    return len(result)


def pyranges_subtract(file_a: str, file_b: str, output: str, **kwargs):
    """Subtract intervals in B from A."""
    import pyranges as pr

    a = pr.read_bed(file_a)
    b = pr.read_bed(file_b)
    result = a.subtract(b)
    result.to_bed(output)
    return len(result)


def pyranges_closest(file_a: str, file_b: str, output: str, **kwargs):
    """Find closest intervals."""
    import pyranges as pr

    a = pr.read_bed(file_a)
    b = pr.read_bed(file_b)
    result = a.nearest(b)
    result.to_bed(output)
    return len(result)


def pyranges_coverage(file_a: str, file_b: str, output: str, **kwargs):
    """Calculate coverage of A over B."""
    import pyranges as pr

    a = pr.read_bed(file_a)
    b = pr.read_bed(file_b)
    result = b.coverage(a)
    result.to_bed(output)
    return len(result)


def pyranges_sort(file_a: str, output: str, **kwargs):
    """Sort intervals."""
    import pyranges as pr

    a = pr.read_bed(file_a)
    result = a.sort()
    result.to_bed(output)
    return len(result)


# =============================================================================
# Polars-Bio Operations
# =============================================================================

def polars_bio_intersect(file_a: str, file_b: str, output: str, **kwargs):
    """Find overlapping intervals using polars-bio."""
    import polars_bio as pb

    result = pb.overlap(file_a, file_b, output_type="inner")
    result.write_csv(output, separator="\t", include_header=False)
    return len(result)


def polars_bio_merge(file_a: str, output: str, **kwargs):
    """Merge overlapping intervals using polars-bio."""
    import polars_bio as pb

    result = pb.merge(file_a)
    result.write_csv(output, separator="\t", include_header=False)
    return len(result)


def polars_bio_closest(file_a: str, file_b: str, output: str, **kwargs):
    """Find closest intervals using polars-bio."""
    import polars_bio as pb

    result = pb.nearest(file_a, file_b)
    result.write_csv(output, separator="\t", include_header=False)
    return len(result)


def polars_bio_coverage(file_a: str, file_b: str, output: str, **kwargs):
    """Calculate coverage using polars-bio."""
    import polars_bio as pb

    result = pb.coverage(file_a, file_b)
    result.write_csv(output, separator="\t", include_header=False)
    return len(result)


# =============================================================================
# Operation Registry
# =============================================================================

OPERATIONS = {
    "pyranges": {
        "intersect": pyranges_intersect,
        "merge": pyranges_merge,
        "subtract": pyranges_subtract,
        "closest": pyranges_closest,
        "coverage": pyranges_coverage,
        "sort": pyranges_sort,
    },
    "polars-bio": {
        "intersect": polars_bio_intersect,
        "merge": polars_bio_merge,
        "closest": polars_bio_closest,
        "coverage": polars_bio_coverage,
    },
}


def list_operations(library: str = None):
    """List available operations for a library or all libraries."""
    if library:
        if library in OPERATIONS:
            print(f"{library}: {', '.join(OPERATIONS[library].keys())}")
        else:
            print(f"Unknown library: {library}")
    else:
        for lib, ops in OPERATIONS.items():
            print(f"{lib}: {', '.join(ops.keys())}")


def run_benchmark(library: str, operation: str, file_a: str, file_b: str = None,
                  output: str = None, **kwargs):
    """Run a benchmark and return timing/memory metrics."""

    if library not in OPERATIONS:
        print(f"Unknown library: {library}", file=sys.stderr)
        print(f"Available: {', '.join(OPERATIONS.keys())}", file=sys.stderr)
        sys.exit(1)

    if operation not in OPERATIONS[library]:
        print(f"Unknown operation '{operation}' for {library}", file=sys.stderr)
        print(f"Available: {', '.join(OPERATIONS[library].keys())}", file=sys.stderr)
        sys.exit(1)

    func = OPERATIONS[library][operation]

    # Determine if operation needs one or two files
    needs_two_files = operation in ["intersect", "subtract", "closest", "coverage"]

    if needs_two_files and not file_b:
        print(f"Operation '{operation}' requires two input files", file=sys.stderr)
        sys.exit(1)

    # Default output
    if not output:
        output = "/dev/null"

    # Start memory tracking
    tracemalloc.start()

    # Run benchmark
    start_time = time.perf_counter()

    try:
        if needs_two_files:
            result_count = func(file_a, file_b, output, **kwargs)
        else:
            result_count = func(file_a, output, **kwargs)
    except Exception as e:
        print(f"Error running {library} {operation}: {e}", file=sys.stderr)
        tracemalloc.stop()
        sys.exit(1)

    end_time = time.perf_counter()

    # Get memory stats
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    elapsed_time = end_time - start_time
    peak_mb = peak / (1024 * 1024)

    return {
        "library": library,
        "operation": operation,
        "time_seconds": elapsed_time,
        "peak_memory_mb": peak_mb,
        "result_count": result_count,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Python benchmark wrapper for pyranges and polars-bio"
    )
    parser.add_argument("library", nargs="?", help="Library to use (pyranges, polars-bio)")
    parser.add_argument("operation", nargs="?", help="Operation to perform")
    parser.add_argument("file_a", nargs="?", help="First input BED file")
    parser.add_argument("file_b", nargs="?", help="Second input BED file (for two-file operations)")
    parser.add_argument("-o", "--output", help="Output file path")
    parser.add_argument("--list", action="store_true", help="List available operations")
    parser.add_argument("--json", action="store_true", help="Output results as JSON")
    parser.add_argument("--csv", action="store_true", help="Output results as CSV line")

    args = parser.parse_args()

    if args.list:
        list_operations(args.library)
        sys.exit(0)

    if not args.library or not args.operation or not args.file_a:
        parser.print_help()
        sys.exit(1)

    results = run_benchmark(
        library=args.library,
        operation=args.operation,
        file_a=args.file_a,
        file_b=args.file_b,
        output=args.output,
    )

    if args.json:
        import json
        print(json.dumps(results, indent=2))
    elif args.csv:
        print(f"{results['library']},{results['operation']},{results['time_seconds']:.4f},{results['peak_memory_mb']:.2f},{results['result_count']}")
    else:
        print(f"Library:     {results['library']}")
        print(f"Operation:   {results['operation']}")
        print(f"Time:        {results['time_seconds']:.4f}s")
        print(f"Peak Memory: {results['peak_memory_mb']:.2f} MB")
        print(f"Results:     {results['result_count']} intervals")


if __name__ == "__main__":
    main()
