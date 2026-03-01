# Quick Start

This guide will get you up and running with pygrit.

## Installation

```bash
pip install grit-genomics
```

## Basic Concepts

pygrit works with **genomic intervals** - regions on chromosomes defined by:

- **Chromosome**: The chromosome name (e.g., "chr1", "chrX")
- **Start**: The 0-based start position (inclusive)
- **End**: The end position (exclusive)

This follows the BED format convention (0-based, half-open intervals).

## Your First Operations

### Creating Intervals

```python
import pygrit

# Create a single interval
interval = pygrit.Interval("chr1", 100, 200)
print(f"Interval: {interval}")        # chr1    100    200
print(f"Length: {len(interval)} bp")  # Length: 100 bp
```

### Reading BED Files

```python
# Read intervals from a BED file
intervals = pygrit.read_bed("regions.bed")
print(f"Loaded {len(intervals)} intervals")

# Access individual intervals
first = intervals[0]
print(f"First interval: {first.chrom}:{first.start}-{first.end}")
```

### Finding Overlaps

```python
# Find overlapping intervals between two files
overlaps = pygrit.intersect("features.bed", "regions.bed")
print(f"Found {len(overlaps)} overlapping intervals")

# Write results to a file
pygrit.intersect(
    "features.bed",
    "regions.bed",
    output="overlaps.bed"
)
```

### Merging Intervals

```python
# Merge overlapping intervals
merged = pygrit.merge("input.bed")

# Merge intervals within 1000bp of each other
merged = pygrit.merge("input.bed", distance=1000)
```

### Sorting a BED File

```python
# Sort an unsorted BED file
pygrit.sort("unsorted.bed", output="sorted.bed")
```

## Input Requirements

!!! warning "Sorted Input Required"
    Most file-based functions require **sorted BED files** (by chromosome, then start position).

    **Unsorted input will produce incorrect results without warning.**

Sort your files first:

=== "Python"
    ```python
    pygrit.sort("unsorted.bed", output="sorted.bed")
    ```

=== "grit CLI"
    ```bash
    grit sort -i unsorted.bed > sorted.bed
    ```

=== "Unix"
    ```bash
    sort -k1,1 -k2,2n unsorted.bed > sorted.bed
    ```

## File-Based vs In-Memory Operations

### File-Based Operations

Process files directly - suitable for larger files:

```python
# Returns list of Interval objects
overlaps = pygrit.intersect("a.bed", "b.bed")

# Or write directly to file
pygrit.intersect("a.bed", "b.bed", output="result.bed")
```

### In-Memory Operations

Load intervals into memory for manipulation:

```python
# Load into memory
intervals = pygrit.read_bed("regions.bed")
merged = intervals.merge(distance=100)

# Convert to list for iteration
for iv in merged.to_list():
    print(iv)
```

## Working with Results

### As Python Lists

```python
results = pygrit.intersect("a.bed", "b.bed")

for interval in results:
    print(f"{interval.chrom}:{interval.start}-{interval.end}")
```

### As NumPy Arrays

```python
import numpy as np

intervals = pygrit.read_bed("regions.bed")
arr = intervals.to_numpy()  # Shape: (n, 2) with [start, end]

# Calculate interval lengths
lengths = arr[:, 1] - arr[:, 0]
print(f"Mean length: {np.mean(lengths):.1f} bp")
```

## Common Options

Most file-based operations support these options:

| Option | Description |
|--------|-------------|
| `output` | Write to file instead of returning list |
| `fraction` | Minimum overlap fraction (0.0-1.0) |
| `reciprocal` | Require reciprocal overlap fraction |

### Example: Reciprocal Overlap

```python
# Require at least 50% overlap in both directions
pygrit.intersect(
    "a.bed",
    "b.bed",
    fraction=0.5,
    reciprocal=True,
    output="significant_overlaps.bed"
)
```

## Next Steps

- [File Operations Guide](../guide/file-operations.md) - Detailed guide on all commands
- [In-Memory Guide](../guide/in-memory.md) - Working with intervals in memory
- [API Reference](../api/index.md) - Complete API documentation
