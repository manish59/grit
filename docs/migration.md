---
layout: default
title: Migrating from bedtools
---

# Migrating from bedtools

GRIT is designed as a drop-in replacement for bedtools. This guide covers command mappings, GRIT-specific optimizations, and common workflows.

## Command Comparison

| bedtools | GRIT (basic) | GRIT (optimized) |
|----------|--------------|------------------|
| `bedtools intersect -a A.bed -b B.bed` | `grit intersect -a A.bed -b B.bed` | `grit intersect -a A.bed -b B.bed --streaming --assume-sorted` |
| `bedtools intersect -a A.bed -b B.bed -sorted` | `grit intersect -a A.bed -b B.bed` | `grit intersect -a A.bed -b B.bed --streaming --assume-sorted` |
| `bedtools subtract -a A.bed -b B.bed` | `grit subtract -a A.bed -b B.bed` | `grit subtract -a A.bed -b B.bed --streaming --assume-sorted` |
| `bedtools merge -i A.bed` | `grit merge -i A.bed` | `grit merge -i A.bed --assume-sorted` |
| `bedtools closest -a A.bed -b B.bed` | `grit closest -a A.bed -b B.bed` | `grit closest -a A.bed -b B.bed --streaming --assume-sorted` |
| `bedtools coverage -a A.bed -b B.bed -sorted` | `grit coverage -a A.bed -b B.bed` | `grit coverage -a A.bed -b B.bed --assume-sorted` |
| `bedtools window -a A.bed -b B.bed -w 1000` | `grit window -a A.bed -b B.bed -w 1000` | `grit window -a A.bed -b B.bed -w 1000 --assume-sorted` |
| `bedtools sort -i A.bed` | `grit sort -i A.bed` | `grit sort -i A.bed` |
| `bedtools slop -i A.bed -g genome.txt -b 100` | `grit slop -i A.bed -g genome.txt -b 100` | Same |
| `bedtools complement -i A.bed -g genome.txt` | `grit complement -i A.bed -g genome.txt` | `grit complement -i A.bed -g genome.txt --assume-sorted` |
| `bedtools jaccard -a A.bed -b B.bed` | `grit jaccard -a A.bed -b B.bed` | Same |

## Key GRIT Flags

| Flag | Description | When to Use |
|------|-------------|-------------|
| `--streaming` | O(k) memory mode | Large files (>1GB), memory-constrained systems |
| `--assume-sorted` | Skip sort validation | Pre-sorted files for faster startup |
| `--allow-unsorted` | Auto-sort in memory | Unsorted input (uses more memory) |
| `-g, --genome` | Validate chromosome order | Ensure genome-specific ordering |
| `--bedtools-compatible` | Match bedtools behavior | Zero-length interval handling |

## Performance Modes

```bash
# Basic (validates input, loads into memory)
grit intersect -a A.bed -b B.bed

# Streaming (constant memory, requires sorted input)
grit intersect -a A.bed -b B.bed --streaming

# Maximum performance (skip validation, streaming)
grit intersect -a A.bed -b B.bed --streaming --assume-sorted

# Handle unsorted input (auto-sorts in memory)
grit intersect -a unsorted.bed -b B.bed --allow-unsorted
```

## Common Workflow Migration

### bedtools workflow
```bash
bedtools sort -i raw.bed > sorted.bed
bedtools merge -i sorted.bed > merged.bed
bedtools intersect -a merged.bed -b features.bed -sorted > result.bed
```

### GRIT equivalent (faster)
```bash
grit sort -i raw.bed > sorted.bed
grit merge -i sorted.bed --assume-sorted > merged.bed
grit intersect -a merged.bed -b features.bed --streaming --assume-sorted > result.bed
```

### GRIT pipeline (fastest - no intermediate files)
```bash
grit sort -i raw.bed | grit merge -i - --assume-sorted | grit intersect -a - -b features.bed --streaming --assume-sorted > result.bed
```

## Global Options

All commands support these options:

| Option | Description |
|--------|-------------|
| `-t, --threads <N>` | Number of threads (default: all CPUs) |
| `--bedtools-compatible` | Normalize zero-length intervals to 1bp for bedtools parity |
| `-h, --help` | Show help for any command |
| `-V, --version` | Show version |

```bash
# Run with 4 threads
grit -t 4 intersect -a file1.bed -b file2.bed

# Enable bedtools-compatible mode for zero-length intervals
grit --bedtools-compatible intersect -a snps.bed -b features.bed

# Get help for a specific command
grit intersect --help
```

## Zero-Length Interval Differences

GRIT uses strict half-open interval semantics by default. Zero-length intervals (`start == end`) contain no bases and don't overlap with anything.

bedtools treats zero-length intervals as 1bp intervals. To match this behavior:

```bash
grit --bedtools-compatible intersect -a snps.bed -b features.bed
```

See [Input Validation](VALIDATION.html) for more details.

## Streaming Mode

For very large files, streaming mode processes data with constant O(k) memory:

```bash
# Intersect
grit intersect -a a.bed -b b.bed --streaming > result.bed

# Subtract
grit subtract -a a.bed -b b.bed --streaming > result.bed

# Closest
grit closest -a a.bed -b b.bed --streaming > result.bed
```

**Memory comparison:**

| Mode | Memory Usage | Best For |
|------|--------------|----------|
| Default (parallel) | O(n + m) | Maximum speed |
| Streaming | O(k) ≈ 2 MB | Large files, low RAM |

## Validation Differences

| Scenario | bedtools | GRIT |
|----------|----------|------|
| Unsorted input | Silent wrong results | Error with fix suggestion |
| Wrong chromosome order | Depends on command | Error with `-g` flag |
| Zero-length intervals | Treats as 1bp | Strict (or `--bedtools-compatible`) |

GRIT validates input by default to prevent silent failures. See [Input Validation](VALIDATION.html) for details.
