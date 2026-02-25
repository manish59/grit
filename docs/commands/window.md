---
layout: default
title: window - GRIT Documentation
---

# grit window

Find intervals in B that are within a window of A.

## Usage

```bash
grit window [OPTIONS] -a <FILE_A> -b <FILE_B>
```

## Options

| Option | Description |
|--------|-------------|
| `-a, --file-a <FILE>` | Input BED file A |
| `-b, --file-b <FILE>` | Input BED file B |
| `-w, --window <N>` | Window size on both sides (default: 1000) |
| `-l, --left <N>` | Left window size |
| `-r, --right <N>` | Right window size |
| `-c, --count` | Report number of matches |
| `-v, --no-overlap` | Only report A intervals with no matches |
| `--assume-sorted` | Skip sorted validation |
| `-g, --genome <FILE>` | Genome file for validation |

## Examples

### Basic window search

```bash
# Find B intervals within 1kb of A intervals
grit window -a genes.bed -b enhancers.bed > nearby.bed

# Custom window size (5kb)
grit window -a genes.bed -b enhancers.bed -w 5000 > nearby_5kb.bed
```

### Asymmetric windows

```bash
# 2kb upstream, 500bp downstream
grit window -a tss.bed -b peaks.bed -l 2000 -r 500 > asymmetric.bed

# Only look upstream (left)
grit window -a genes.bed -b promoters.bed -l 5000 -r 0 > upstream.bed
```

### Count nearby intervals

```bash
# Count how many B intervals are near each A
grit window -a genes.bed -b snps.bed -w 10000 -c > snp_counts.bed
```

### Find isolated intervals

```bash
# Find A intervals with NO nearby B intervals
grit window -a genes.bed -b other_genes.bed -w 50000 -v > isolated.bed
```

## Output

**Default output:**
```
chr1    100    200    chr1    150    250
chr1    100    200    chr1    300    400
```

**With -c (counts):**
```
chr1    100    200    2
chr1    500    600    0
```

## Visual Example

```
Window: 100bp on each side

A interval:        |--------|
Window:      |-----|--------|-----|
B intervals:   x        x       x
               ^        ^       ^
            (within window, all reported)
```

## Performance

```bash
# Fastest with sorted input
grit window -a sorted_a.bed -b sorted_b.bed --assume-sorted
```

[← Back to Commands](../index.html)
