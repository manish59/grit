---
layout: default
title: subtract
parent: Commands
nav_order: 4
---

# grit subtract

Remove intervals in A that overlap with B.

## Usage

```bash
grit subtract [OPTIONS] -a <FILE_A> -b <FILE_B>
```

## Options

| Option | Description |
|--------|-------------|
| `-a, --file-a <FILE>` | Input BED file A |
| `-b, --file-b <FILE>` | Input BED file B |
| `-A, --remove-entire` | Remove entire A feature if any overlap |
| `-f, --fraction <F>` | Minimum overlap fraction required |
| `-r, --reciprocal` | Require reciprocal fraction overlap |
| `--streaming` | Use streaming mode (O(k) memory) |
| `--assume-sorted` | Skip sorted validation |
| `--allow-unsorted` | Allow unsorted input |
| `-g, --genome <FILE>` | Genome file for chromosome order validation |

## Examples

### Basic subtraction

```bash
# Remove overlapping portions
grit subtract -a regions.bed -b blacklist.bed > clean.bed

# Streaming mode for large files
grit subtract -a regions.bed -b blacklist.bed --streaming --assume-sorted
```

### Remove entire intervals

```bash
# Remove entire A interval if ANY overlap with B
grit subtract -a genes.bed -b repeats.bed -A > non_repeat_genes.bed
```

### Minimum overlap fraction

```bash
# Only subtract if overlap is at least 50%
grit subtract -a regions.bed -b mask.bed -f 0.5 > filtered.bed
```

## Output

**Input A:**
```
chr1    100    300
chr1    400    500
```

**Input B:**
```
chr1    150    200
```

**Output (default):**
```
chr1    100    150
chr1    200    300
chr1    400    500
```

**Output (with -A):**
```
chr1    400    500
```

## Visual Example

```
A:      |------------|        |------|
B:           |--|

Result: |----| |----|        |------|
        ^kept^ ^kept^        ^kept^
             ^removed^
```

## Performance

```bash
# Fastest for sorted input
grit subtract -a sorted_a.bed -b sorted_b.bed --streaming --assume-sorted
```

[← Back to Commands](../index.html)
