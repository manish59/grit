---
layout: default
title: sort - GRIT Documentation
---

# grit sort

Sort a BED file by chromosome and position.

## Usage

```bash
grit sort [OPTIONS] -i <INPUT>
```

## Options

| Option | Description |
|--------|-------------|
| `-i, --input <FILE>` | Input BED file (use `-` for stdin) |
| `-g, --genome <FILE>` | Genome file for chromosome ordering |
| `--sizeA` | Sort by interval size (ascending) |
| `--sizeD` | Sort by interval size (descending) |
| `-r, --reverse` | Reverse the sort order |
| `--chrThenSizeA` | Sort by chromosome name only |
| `--stats` | Print sorting statistics to stderr |

## Examples

### Basic sorting

```bash
# Sort by chromosome and position
grit sort -i unsorted.bed > sorted.bed

# Sort from stdin
cat unsorted.bed | grit sort -i - > sorted.bed
```

### Sort by interval size

```bash
# Sort by size, smallest first
grit sort -i regions.bed --sizeA > by_size_asc.bed

# Sort by size, largest first
grit sort -i regions.bed --sizeD > by_size_desc.bed
```

### Custom chromosome order

```bash
# Use genome file for chromosome ordering
grit sort -i regions.bed -g genome.txt > sorted.bed
```

**genome.txt format:**
```
chr1    248956422
chr2    242193529
chr3    198295559
...
```

### Reverse sort

```bash
# Reverse the default sort order
grit sort -i regions.bed -r > reversed.bed
```

## Output

Output is a sorted BED file with the same columns as input.

**Input:**
```
chr2    100    200
chr1    500    600
chr1    100    200
```

**Output:**
```
chr1    100    200
chr1    500    600
chr2    100    200
```

## Performance

- Uses parallel sorting for large files
- Memory usage: O(n) where n = number of intervals

[← Back to Commands](../index.html)
