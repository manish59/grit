# GRIT v0.1.0 - Genomic Range Interval Toolkit

High-performance genomic interval operations implemented in Rust.

## Installation

```bash
# From crates.io
cargo install grit-genomics

# The CLI command is 'grit'
grit --help
```

## Features

- **10 BED file commands**: intersect, subtract, merge, coverage, closest, window, complement, jaccard, slop, sort
- **Streaming algorithms**: O(k) memory where k = max overlapping intervals
- **3-13x faster** than bedtools with **up to 1000x less memory**
- **Byte-for-byte correctness**: SHA256-verified output matching bedtools

## Benchmarks (10M × 5M intervals)

| Command | Speedup | Memory Reduction |
|---------|---------|------------------|
| intersect | 4.1x | 19x less |
| subtract | 5.8x | 19x less |
| coverage | 8.1x | 135x less |
| closest | 4.7x | 59x less |
| window | 13.5x | 138x less |
| jaccard | 3.1x | 1160x less |

## Usage

```bash
# Intersect two BED files (streaming mode)
grit intersect -a regions.bed -b reads.bed --streaming --assume-sorted

# Merge overlapping intervals
grit merge -i input.bed --assume-sorted

# Calculate coverage
grit coverage -a regions.bed -b reads.bed --assume-sorted
```

## Full Changelog

Initial public release.
