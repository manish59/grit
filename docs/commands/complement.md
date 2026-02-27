---
layout: default
title: complement
parent: Commands
nav_order: 9
---

# grit complement

Return intervals NOT covered by the input BED file.

## Usage

```bash
grit complement [OPTIONS] -i <INPUT> -g <GENOME>
```

## Options

| Option | Description |
|--------|-------------|
| `-i, --input <FILE>` | Input BED file |
| `-g, --genome <FILE>` | Genome file (chromosome sizes) |
| `--assume-sorted` | Assume input is sorted (O(1) memory streaming) |

## Examples

### Basic complement

```bash
# Get regions NOT covered by input
grit complement -i covered.bed -g genome.txt > gaps.bed

# With sorted input (faster, less memory)
grit complement -i sorted.bed -g genome.txt --assume-sorted > gaps.bed
```

### Find intergenic regions

```bash
# Get regions between genes
grit complement -i genes.bed -g genome.txt > intergenic.bed
```

### Find uncovered regions

```bash
# Find gaps in sequencing coverage
grit complement -i covered_regions.bed -g genome.txt > uncovered.bed
```

## Genome File Format

```
chr1    248956422
chr2    242193529
chr3    198295559
```

## Output

**Input (covered.bed):**
```
chr1    100    200
chr1    300    400
```

**Genome (genome.txt):**
```
chr1    500
```

**Output (gaps):**
```
chr1    0      100
chr1    200    300
chr1    400    500
```

## Visual Example

```
Chromosome:  |----------------------------------------|
             0                                      500

Input:            |-----|     |-----|
                  100-200     300-400

Complement:  |----|     |-----|     |-----------------|
             0-100      200-300     400-500
```

## Use Cases

### Intergenic regions
```bash
grit complement -i genes.bed -g genome.txt > intergenic.bed
```

### Coverage gaps
```bash
grit complement -i aligned_regions.bed -g genome.txt > gaps.bed
```

### Accessible chromatin (inverse of closed)
```bash
grit complement -i closed_chromatin.bed -g genome.txt > open.bed
```

## Performance

```bash
# Streaming mode with sorted input
grit complement -i sorted.bed -g genome.txt --assume-sorted
```

[← Back to Commands](../index.html)
