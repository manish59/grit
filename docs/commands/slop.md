---
layout: default
title: slop
parent: Commands
nav_order: 8
---

# grit slop

Extend intervals by a given number of bases.

## Usage

```bash
grit slop [OPTIONS] -i <INPUT> -g <GENOME>
```

## Options

| Option | Description |
|--------|-------------|
| `-i, --input <FILE>` | Input BED file |
| `-g, --genome <FILE>` | Genome file (chromosome sizes) |
| `-b, --both <N>` | Extend both sides by N bases |
| `-l, --left <N>` | Extend left side by N bases |
| `-r, --right <N>` | Extend right side by N bases |
| `-s, --strand` | Use strand info (left=upstream, right=downstream) |
| `--pct` | Interpret values as fraction of interval size |

## Examples

### Extend both sides

```bash
# Extend each interval by 100bp on both sides
grit slop -i regions.bed -g genome.txt -b 100 > extended.bed

# Extend by 1kb
grit slop -i peaks.bed -g genome.txt -b 1000 > peaks_1kb.bed
```

### Asymmetric extension

```bash
# Extend 500bp left, 200bp right
grit slop -i genes.bed -g genome.txt -l 500 -r 200 > extended.bed
```

### Strand-aware extension

```bash
# Extend 1kb upstream, 500bp downstream (strand-aware)
grit slop -i genes.bed -g genome.txt -l 1000 -r 500 -s > extended.bed
```

### Percentage-based extension

```bash
# Extend by 10% of interval size on each side
grit slop -i regions.bed -g genome.txt -b 0.1 --pct > extended.bed

# Double the interval size (50% each side)
grit slop -i regions.bed -g genome.txt -b 0.5 --pct > doubled.bed
```

## Genome File Format

The genome file specifies chromosome sizes:

```
chr1    248956422
chr2    242193529
chr3    198295559
chrX    156040895
chrY    57227415
```

Create from a FASTA index:
```bash
cut -f1,2 genome.fa.fai > genome.txt
```

## Output

**Input:**
```
chr1    100    200
```

**With -b 50:**
```
chr1    50    250
```

**Note:** Coordinates are clamped to chromosome boundaries (0 to chrom_size).

## Visual Example

```
Original:           |--------|
With -b 100:  |-----|--------|-----|
              ^100bp         ^100bp

With -l 200 -r 50:
             |------|--------|--|
             ^200bp          ^50bp
```

## Use Cases

### Promoter regions
```bash
# Get 2kb upstream of TSS
grit slop -i tss.bed -g genome.txt -l 2000 -r 0 -s > promoters.bed
```

### Peak extension
```bash
# Extend ChIP-seq peaks by 250bp
grit slop -i peaks.bed -g genome.txt -b 250 > extended_peaks.bed
```

[← Back to Commands](../index.html)
