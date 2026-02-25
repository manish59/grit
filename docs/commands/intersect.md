---
layout: default
title: intersect - GRIT Documentation
---

# grit intersect

Find overlapping intervals between two BED files.

## Usage

```bash
grit intersect [OPTIONS] -a <FILE_A> -b <FILE_B>
```

## Options

| Option | Description |
|--------|-------------|
| `-a, --file-a <FILE>` | Input BED file A |
| `-b, --file-b <FILE>` | Input BED file B |
| `--wa` | Write original A entry |
| `--wb` | Write original B entry |
| `-u, --unique` | Only report unique A intervals |
| `-v, --no-overlap` | Only report A intervals with NO overlap |
| `-f, --fraction <F>` | Minimum overlap fraction for A (0.0-1.0) |
| `-r, --reciprocal` | Require reciprocal fraction overlap |
| `-c, --count` | Report the number of overlaps |
| `--streaming` | Use streaming mode (constant memory) |
| `--assume-sorted` | Skip sorted validation |
| `--allow-unsorted` | Allow unsorted input (uses O(n) memory) |
| `-g, --genome <FILE>` | Genome file for chromosome order validation |

## Examples

### Basic intersection

```bash
# Find overlapping regions
grit intersect -a regions.bed -b reads.bed > overlaps.bed

# Streaming mode for large files
grit intersect -a regions.bed -b reads.bed --streaming --assume-sorted
```

### Report original entries

```bash
# Report original A entry for each overlap
grit intersect -a genes.bed -b peaks.bed --wa > results.bed

# Report both A and B entries
grit intersect -a genes.bed -b peaks.bed --wa --wb > results.bed
```

### Unique intervals

```bash
# Only report each A interval once
grit intersect -a genes.bed -b reads.bed -u > unique_genes.bed
```

### Find non-overlapping intervals

```bash
# Report A intervals that DON'T overlap with B
grit intersect -a all_regions.bed -b excluded.bed -v > remaining.bed
```

### Minimum overlap fraction

```bash
# Require at least 50% overlap
grit intersect -a genes.bed -b reads.bed -f 0.5 > significant_overlaps.bed

# Require reciprocal 50% overlap (both A and B must overlap 50%)
grit intersect -a set1.bed -b set2.bed -f 0.5 -r > reciprocal.bed
```

### Count overlaps

```bash
# Count how many B intervals overlap each A interval
grit intersect -a genes.bed -b reads.bed -c > counts.bed
```

## Output

**Default output** (intersection coordinates):
```
chr1    150    200
```

**With --wa --wb** (both entries):
```
chr1    100    200    chr1    150    250
```

**With -c** (counts):
```
chr1    100    200    5
```

## Performance

For large files, use streaming mode:

```bash
# Fastest: streaming with pre-sorted input
grit intersect -a sorted_a.bed -b sorted_b.bed --streaming --assume-sorted
```

| Mode | Memory | Requirements |
|------|--------|--------------|
| Streaming | O(k) | Sorted input |
| Default | O(n) | Any input |

[← Back to Commands](../index.html)
