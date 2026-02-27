---
layout: default
title: closest
parent: Commands
nav_order: 5
---

# grit closest

Find the closest interval in B for each interval in A.

## Usage

```bash
grit closest [OPTIONS] -a <FILE_A> -b <FILE_B>
```

## Options

| Option | Description |
|--------|-------------|
| `-a, --file-a <FILE>` | Input BED file A |
| `-b, --file-b <FILE>` | Input BED file B |
| `-d, --distance` | Report distance in output |
| `-t, --tie <MODE>` | How to handle ties: `all`, `first`, `last` |
| `--io` | Ignore overlapping intervals |
| `--iu` | Ignore upstream intervals |
| `--id` | Ignore downstream intervals |
| `-D, --max-distance <N>` | Maximum distance to report |
| `--streaming` | Use streaming mode (O(k) memory) |
| `--assume-sorted` | Skip sorted validation |
| `--allow-unsorted` | Allow unsorted input |
| `-g, --genome <FILE>` | Genome file for validation |

## Examples

### Basic closest

```bash
# Find closest B interval for each A interval
grit closest -a genes.bed -b enhancers.bed > closest.bed

# Streaming mode for large files
grit closest -a genes.bed -b enhancers.bed --streaming --assume-sorted
```

### Report distance

```bash
# Include distance to closest interval
grit closest -a tss.bed -b peaks.bed -d > with_distance.bed
```

### Handle ties

```bash
# Report all equidistant intervals
grit closest -a genes.bed -b peaks.bed -t all > all_closest.bed

# Report only the first tie
grit closest -a genes.bed -b peaks.bed -t first > first_closest.bed
```

### Direction filtering

```bash
# Only find downstream intervals
grit closest -a tss.bed -b enhancers.bed --iu > downstream_only.bed

# Only find upstream intervals
grit closest -a tss.bed -b enhancers.bed --id > upstream_only.bed

# Ignore overlapping, find truly closest
grit closest -a genes.bed -b peaks.bed --io > non_overlapping.bed
```

### Maximum distance

```bash
# Only report if within 10kb
grit closest -a tss.bed -b enhancers.bed -D 10000 > within_10kb.bed
```

## Output

**Default output:**
```
chr1    100    200    chr1    250    300
```

**With -d (distance):**
```
chr1    100    200    chr1    250    300    50
```

**No match found:**
```
chr1    100    200    .    -1    -1
```

## Performance

```bash
# Fastest for sorted input
grit closest -a sorted_a.bed -b sorted_b.bed --streaming --assume-sorted
```

[← Back to Commands](../index.html)
