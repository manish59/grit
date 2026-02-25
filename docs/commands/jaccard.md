---
layout: default
title: jaccard - GRIT Documentation
---

# grit jaccard

Calculate Jaccard similarity between two BED files.

## Usage

```bash
grit jaccard [OPTIONS] -a <FILE_A> -b <FILE_B>
```

## Options

| Option | Description |
|--------|-------------|
| `-a, --file-a <FILE>` | Input BED file A |
| `-b, --file-b <FILE>` | Input BED file B |

## Examples

### Basic Jaccard similarity

```bash
# Calculate Jaccard index between two BED files
grit jaccard -a peaks1.bed -b peaks2.bed
```

### Compare replicates

```bash
# Compare ChIP-seq replicates
grit jaccard -a rep1_peaks.bed -b rep2_peaks.bed
```

### Compare conditions

```bash
# Compare peaks between conditions
grit jaccard -a control_peaks.bed -b treatment_peaks.bed
```

## Output

```
intersection    union    jaccard    n_intersections
500000    2000000    0.25    150
```

| Column | Description |
|--------|-------------|
| intersection | Total bases in intersection |
| union | Total bases in union |
| jaccard | Jaccard index (intersection/union) |
| n_intersections | Number of intersecting interval pairs |

## Understanding Jaccard Index

The Jaccard index measures similarity between two sets:

```
J(A,B) = |A ∩ B| / |A ∪ B|
```

| Value | Interpretation |
|-------|----------------|
| 0.0 | No overlap |
| 0.5 | 50% similarity |
| 1.0 | Identical |

## Visual Example

```
A:    |-------|    |-------|
B:        |-------|    |-----|

Intersection: |--|         |-|
Union:    |-----------|---------|

Jaccard = intersection_bp / union_bp
```

## Use Cases

### Replicate concordance
```bash
# High Jaccard = good reproducibility
grit jaccard -a rep1.bed -b rep2.bed
# Output: jaccard = 0.85 (good concordance)
```

### Method comparison
```bash
# Compare peak callers
grit jaccard -a macs2_peaks.bed -b homer_peaks.bed
```

### Condition comparison
```bash
# Compare differential binding
grit jaccard -a wt_peaks.bed -b ko_peaks.bed
```

### Quality control
```bash
# Compare to gold standard
grit jaccard -a my_peaks.bed -b encode_peaks.bed
```

## Performance

GRIT uses an efficient streaming algorithm, using significantly less memory than bedtools for large files.

[← Back to Commands](../index.html)
