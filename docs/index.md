---
layout: default
title: GRIT - Genomic Range Interval Toolkit
---

# GRIT Documentation

**GRIT** (Genomic Range Interval Toolkit) is a high-performance Rust implementation of common BED file operations. It provides streaming algorithms with O(k) memory complexity, achieving **3-13x speedup** over bedtools with **up to 1000x less memory**.

## Installation

```bash
# From crates.io (recommended)
cargo install grit-genomics

# From Homebrew (macOS/Linux)
brew tap manish59/grit && brew install grit

# From Bioconda (coming soon)
conda install -c bioconda grit-genomics
```

## Quick Start

```bash
# Sort a BED file
grit sort -i unsorted.bed > sorted.bed

# Find overlapping intervals (streaming mode for large files)
grit intersect -a regions.bed -b reads.bed --streaming --assume-sorted

# Merge overlapping intervals
grit merge -i intervals.bed --assume-sorted

# Calculate coverage
grit coverage -a genes.bed -b reads.bed --assume-sorted
```

## Commands

| Command | Description |
|---------|-------------|
| [sort](commands/sort.html) | Sort a BED file by chromosome and position |
| [merge](commands/merge.html) | Merge overlapping intervals |
| [intersect](commands/intersect.html) | Find overlapping intervals between two BED files |
| [subtract](commands/subtract.html) | Remove intervals in A that overlap with B |
| [closest](commands/closest.html) | Find the closest interval in B for each interval in A |
| [window](commands/window.html) | Find intervals in B within a window of A |
| [coverage](commands/coverage.html) | Calculate coverage of A intervals by B intervals |
| [slop](commands/slop.html) | Extend intervals by a given number of bases |
| [complement](commands/complement.html) | Return intervals NOT covered by the input |
| [genomecov](commands/genomecov.html) | Compute genome-wide coverage |
| [jaccard](commands/jaccard.html) | Calculate Jaccard similarity between two BED files |
| [multiinter](commands/multiinter.html) | Identify common intervals across multiple files |

## Global Options

All commands support these global options:

| Option | Description |
|--------|-------------|
| `-t, --threads <N>` | Number of threads (default: number of CPUs) |
| `--bedtools-compatible` | Match bedtools behavior for zero-length intervals |
| `-h, --help` | Print help |
| `-V, --version` | Print version |

## Performance Tips

For maximum performance with sorted input files:

```bash
# Use --streaming for constant memory usage
grit intersect -a large_a.bed -b large_b.bed --streaming --assume-sorted

# Skip sort validation if you know input is sorted
grit merge -i sorted.bed --assume-sorted
```

## Benchmarks

Tested on 10M × 5M intervals:

| Command | Speedup | Memory Reduction |
|---------|---------|------------------|
| intersect | 4.1x | 19x less |
| subtract | 5.8x | 19x less |
| coverage | 8.1x | 135x less |
| closest | 4.7x | 59x less |
| window | 13.5x | 138x less |
| jaccard | 3.1x | 1160x less |

## Links

- [GitHub Repository](https://github.com/manish59/grit)
- [crates.io](https://crates.io/crates/grit-genomics)
- [Report Issues](https://github.com/manish59/grit/issues)
