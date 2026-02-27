---
layout: default
title: coverage
parent: Commands
nav_order: 7
---

# grit coverage

Calculate coverage of A intervals by B intervals.

## Usage

```bash
grit coverage [OPTIONS] -a <FILE_A> -b <FILE_B>
```

## Options

| Option | Description |
|--------|-------------|
| `-a, --file-a <FILE>` | Input BED file A (regions) |
| `-b, --file-b <FILE>` | Input BED file B (reads/features) |
| `--hist` | Report a histogram of coverage |
| `-d, --per-base` | Report depth at each position |
| `--mean` | Report mean depth |
| `--assume-sorted` | Skip sorted validation |
| `-g, --genome <FILE>` | Genome file for validation |

## Examples

### Basic coverage

```bash
# Calculate coverage of genes by reads
grit coverage -a genes.bed -b reads.bed > coverage.bed

# With sorted input (faster)
grit coverage -a genes.bed -b reads.bed --assume-sorted
```

### Coverage histogram

```bash
# Generate coverage histogram
grit coverage -a regions.bed -b reads.bed --hist > histogram.txt
```

### Per-base depth

```bash
# Report depth at each position
grit coverage -a genes.bed -b reads.bed -d > per_base.bed
```

### Mean coverage

```bash
# Report mean coverage depth
grit coverage -a genes.bed -b reads.bed --mean > mean_coverage.bed
```

## Output

**Default output** (7 columns):
```
chr1    100    200    3    75    100    0.7500000
```

| Column | Description |
|--------|-------------|
| 1-3 | Chromosome, start, end (from A) |
| 4 | Number of B intervals overlapping |
| 5 | Number of bases covered |
| 6 | Length of A interval |
| 7 | Fraction of A covered |

**With --hist:**
```
chr1    100    200    0    25
chr1    100    200    1    50
chr1    100    200    2    25
all     0      100    25
all     1      100    50
all     2      100    25
```

**With -d (per-base):**
```
chr1    100    2
chr1    101    2
chr1    102    3
...
```

## Use Cases

### RNA-seq gene coverage
```bash
grit coverage -a genes.bed -b aligned_reads.bed --assume-sorted > gene_coverage.bed
```

### ChIP-seq peak coverage
```bash
grit coverage -a peaks.bed -b reads.bed --mean > peak_depths.bed
```

### Exome coverage analysis
```bash
grit coverage -a exome_targets.bed -b wes_reads.bed --hist > coverage_hist.txt
```

## Performance

```bash
# Fastest with sorted input
grit coverage -a sorted_regions.bed -b sorted_reads.bed --assume-sorted
```

[← Back to Commands](../index.html)
