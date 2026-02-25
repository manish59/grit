---
layout: default
title: genomecov - GRIT Documentation
---

# grit genomecov

Compute genome-wide coverage.

## Usage

```bash
grit genomecov [OPTIONS] -i <INPUT> -g <GENOME>
```

## Options

| Option | Description |
|--------|-------------|
| `-i, --input <FILE>` | Input BED file |
| `-g, --genome <FILE>` | Genome file (chromosome sizes) |
| `-d, --per-base` | Report depth at each position (1-based) |
| `--bg` | Report BedGraph format (non-zero regions only) |
| `--bga` | Report BedGraph format (including zero coverage) |
| `--scale <F>` | Scale depth by factor (default: 1.0) |
| `--streaming` | Use streaming mode (O(k) memory) |
| `--assume-sorted` | Skip sorted validation |

## Examples

### Default histogram

```bash
# Generate coverage histogram
grit genomecov -i reads.bed -g genome.txt > histogram.txt
```

### BedGraph output

```bash
# BedGraph format (for visualization)
grit genomecov -i reads.bed -g genome.txt --bg > coverage.bedgraph

# Include zero-coverage regions
grit genomecov -i reads.bed -g genome.txt --bga > coverage_all.bedgraph
```

### Per-base depth

```bash
# Report depth at every position
grit genomecov -i reads.bed -g genome.txt -d > per_base.txt
```

### Scaled coverage

```bash
# Scale by RPM factor
grit genomecov -i reads.bed -g genome.txt --bg --scale 0.001 > rpm.bedgraph
```

### Streaming mode

```bash
# For large files, use streaming
grit genomecov -i sorted_reads.bed -g genome.txt --bg --streaming --assume-sorted
```

## Output Formats

**Default (histogram):**
```
chr1    0    1000000    0.80
chr1    1    200000     0.16
chr1    2    50000      0.04
genome  0    5000000    0.75
genome  1    1000000    0.20
genome  2    250000     0.05
```

| Column | Description |
|--------|-------------|
| 1 | Chromosome (or "genome" for total) |
| 2 | Coverage depth |
| 3 | Number of bases at this depth |
| 4 | Fraction of chromosome/genome |

**BedGraph (--bg):**
```
chr1    0      1000   5
chr1    1000   2000   3
chr1    3000   4000   7
```

**Per-base (-d):**
```
chr1    1    2
chr1    2    2
chr1    3    3
...
```

## Use Cases

### RNA-seq coverage track
```bash
grit genomecov -i aligned.bed -g genome.txt --bg > coverage.bedgraph
```

### Normalized coverage (RPM)
```bash
# Calculate scaling factor
total_reads=$(wc -l < reads.bed)
scale=$(echo "scale=10; 1000000 / $total_reads" | bc)

grit genomecov -i reads.bed -g genome.txt --bg --scale $scale > rpm.bedgraph
```

### Coverage statistics
```bash
grit genomecov -i reads.bed -g genome.txt > stats.txt
# Parse the "genome" lines for overall statistics
```

## Performance

```bash
# Fastest with streaming mode
grit genomecov -i sorted.bed -g genome.txt --bg --streaming --assume-sorted
```

[← Back to Commands](../index.html)
