---
layout: default
title: Commands
nav_order: 7
has_children: true
---

# GRIT Command Reference

GRIT: Genomic Range Interval Toolkit - High-performance genomic interval operations.

## Global Options

| Flag | Description |
|------|-------------|
| `-t, --threads` | Number of threads to use (default: number of CPUs) |
| `--bedtools-compatible` | Normalize zero-length intervals to 1bp for bedtools compatibility |

## Command Index

| Command | Description | Streaming | Example |
|---------|-------------|:---------:|---------|
| [sort](EXAMPLES/sort.md) | Sort BED file by chromosome and position | - | `grit sort -i input.bed` |
| [merge](EXAMPLES/merge.md) | Merge overlapping intervals | Yes | `grit merge -i input.bed` |
| [intersect](EXAMPLES/intersect.md) | Find overlapping intervals between two files | Yes | `grit intersect -a a.bed -b b.bed` |
| [subtract](EXAMPLES/subtract.md) | Remove intervals in A that overlap with B | Yes | `grit subtract -a a.bed -b b.bed` |
| [closest](EXAMPLES/closest.md) | Find closest interval in B for each A | Yes | `grit closest -a a.bed -b b.bed` |
| [window](EXAMPLES/window.md) | Find intervals within a window | Yes | `grit window -a a.bed -b b.bed` |
| [coverage](EXAMPLES/coverage.md) | Calculate coverage of A by B | Yes | `grit coverage -a a.bed -b b.bed` |
| [slop](EXAMPLES/slop.md) | Extend intervals by given bases | - | `grit slop -i input.bed -g genome.txt` |
| [complement](EXAMPLES/complement.md) | Return uncovered regions | Yes | `grit complement -i input.bed -g genome.txt` |
| [genomecov](EXAMPLES/genomecov.md) | Compute genome-wide coverage | - | `grit genomecov -i input.bed -g genome.txt` |
| [jaccard](EXAMPLES/jaccard.md) | Calculate Jaccard similarity | - | `grit jaccard -a a.bed -b b.bed` |
| [multiinter](EXAMPLES/multiinter.md) | Find common intervals across files | - | `grit multiinter -i a.bed b.bed c.bed` |
| [generate](EXAMPLES/generate.md) | Generate synthetic datasets | - | `grit generate --sizes 1M` |

## Streaming Support

Commands marked with "Yes" in the Streaming column support streaming mode with O(k) memory usage, where k is the maximum number of overlapping intervals at any position.

### Enabling Streaming Mode

```bash
# Requires sorted input
grit intersect -a a.bed -b b.bed --streaming --assume-sorted
```

### Performance Flags

| Flag | Description |
|------|-------------|
| `--streaming` | Use streaming algorithm (constant memory) |
| `--assume-sorted` | Skip sorted validation (faster for pre-sorted input) |
| `--stats` | Print execution statistics to stderr |

## Bedtools Compatibility Mode

GRIT uses strict half-open interval semantics by default. Zero-length intervals (where start == end) represent points and do not overlap with themselves.

For bedtools compatibility, use the `--bedtools-compatible` flag:

```bash
grit --bedtools-compatible intersect -a snps.bed -b regions.bed
```

This normalizes zero-length intervals to 1bp intervals during parsing, matching bedtools behavior.

## Input/Output

### Standard Input

Most commands support reading from stdin using `-` as the filename:

```bash
cat input.bed | grit sort -i - | grit merge --assume-sorted
```

### Piping Commands

```bash
grit sort -i unsorted.bed | grit merge --assume-sorted | grit intersect -a - -b features.bed --streaming
```

## Memory Usage

| Mode | Memory | Use Case |
|------|--------|----------|
| Streaming | O(k) | Large sorted files |
| Parallel | O(n) | Unsorted files, complex operations |
| In-memory | O(n) | Full file manipulation |

Where:
- k = maximum overlapping intervals at any position
- n = total number of intervals

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Error (invalid input, I/O error, etc.) |

## See Also

- [Tutorial](TUTORIAL.md) - Getting started guide
- [README](../README.md) - Project overview
