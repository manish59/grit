---
layout: default
title: merge
parent: Commands
nav_order: 2
---

# grit merge

Merge overlapping intervals into single intervals.

## Usage

```bash
grit merge [OPTIONS] -i <INPUT>
```

## Options

| Option | Description |
|--------|-------------|
| `-i, --input <FILE>` | Input BED file (use `-` for stdin) |
| `-d, --distance <N>` | Maximum distance between intervals to merge (default: 0) |
| `-s, --strand` | Require strand to match for merging |
| `-c, --count` | Report count of merged intervals |
| `--in-memory` | Use in-memory mode (handles unsorted input) |
| `--assume-sorted` | Skip sorted validation (faster) |
| `-g, --genome <FILE>` | Genome file for chromosome order validation |
| `--stats` | Print streaming statistics to stderr |

## Examples

### Basic merge

```bash
# Merge overlapping intervals
grit merge -i regions.bed > merged.bed

# With pre-sorted input (faster)
grit merge -i sorted.bed --assume-sorted > merged.bed
```

### Merge nearby intervals

```bash
# Merge intervals within 100bp of each other
grit merge -i regions.bed -d 100 > merged.bed

# Merge intervals within 1kb
grit merge -i regions.bed -d 1000 > merged.bed
```

### Strand-specific merge

```bash
# Only merge intervals on the same strand
grit merge -i stranded.bed -s > merged.bed
```

### Count merged intervals

```bash
# Report how many intervals were merged
grit merge -i regions.bed -c > merged_with_counts.bed
```

### Handle unsorted input

```bash
# Automatically sort and merge (uses more memory)
grit merge -i unsorted.bed --in-memory > merged.bed
```

## Output

**Input:**
```
chr1    100    200
chr1    150    250
chr1    300    400
```

**Output (basic):**
```
chr1    100    250
chr1    300    400
```

**Output (with -c):**
```
chr1    100    250    2
chr1    300    400    1
```

## Performance

- **Streaming mode** (default): O(k) memory, requires sorted input
- **In-memory mode**: O(n) memory, handles unsorted input

```bash
# Fastest for sorted input
grit merge -i sorted.bed --assume-sorted
```

[← Back to Commands](../index.html)
