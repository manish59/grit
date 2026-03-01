# pygrit

Python bindings for GRIT (Genomic Range Interval Toolkit) - high-performance genomic interval operations powered by Rust.

## Installation

```bash
pip install pygrit
```

Or build from source:

```bash
cd py-pygrit
maturin develop --release
```

## Quick Start

```python
import pygrit

# File-based streaming (recommended for large files)
# Uses O(k) memory where k = max overlapping intervals
results = pygrit.intersect("a.bed", "b.bed")
pygrit.merge("input.bed", output="merged.bed", distance=100)
pygrit.coverage("regions.bed", "reads.bed", output="coverage.bed")

# In-memory operations
intervals = pygrit.read_bed("regions.bed")
merged = intervals.merge(distance=100)

# Work with individual intervals
iv = pygrit.Interval("chr1", 100, 200)
print(len(iv))  # 100
```

## File-Based Streaming API

For large BED files (millions of intervals), use the file-based functions that maintain streaming behavior:

```python
import pygrit

# Intersect two BED files
results = pygrit.intersect("a.bed", "b.bed")  # Returns list of Interval
pygrit.intersect("a.bed", "b.bed", output="out.bed")  # Writes to file

# With options (matching bedtools flags)
pygrit.intersect(
    "a.bed", "b.bed",
    write_a=True,      # -wa
    write_b=True,      # -wb
    fraction=0.5,      # -f 0.5
    reciprocal=True,   # -r
    no_overlap=True,   # -v
)

# Merge overlapping intervals
merged = pygrit.merge("input.bed", distance=100)

# Subtract intervals
remaining = pygrit.subtract("a.bed", "b.bed")

# Calculate coverage
coverage = pygrit.coverage("regions.bed", "reads.bed")

# Find closest intervals
closest = pygrit.closest("a.bed", "b.bed")

# Find intervals within window
nearby = pygrit.window("a.bed", "b.bed", window=1000)

# Sort a BED file
pygrit.sort("unsorted.bed", output="sorted.bed")

# Extend intervals (slop)
pygrit.slop("regions.bed", "genome.txt", both=100.0, output="extended.bed")

# Find gaps between intervals
pygrit.complement("features.bed", "genome.txt", output="gaps.bed")

# Genome-wide coverage
pygrit.genomecov("reads.bed", "genome.txt", bg=True, output="coverage.bg")

# Jaccard similarity
pygrit.jaccard("set_a.bed", "set_b.bed", output="similarity.txt")

# Multi-file intersection
pygrit.multiinter(["a.bed", "b.bed", "c.bed"], output="overlap.bed")

# Generate test data
stats = pygrit.generate("./test_data", num_intervals=10000, seed=42)
```

## In-Memory API

For smaller datasets or programmatic interval manipulation:

```python
import pygrit
import numpy as np

# Read a BED file
intervals = pygrit.read_bed("regions.bed")

# Create intervals from a string
intervals = pygrit.parse_bed("""
chr1	100	200
chr1	150	250
chr1	300	400
""")

# Merge overlapping intervals
merged = intervals.merge(distance=100)

# Find intersections
other = pygrit.read_bed("other.bed")
overlaps = intervals.intersect(other, fraction=0.5)

# Find non-overlapping intervals
non_overlapping = intervals.non_overlapping(other)

# Work with NumPy arrays
arr = np.array([[0, 100], [150, 200]], dtype=np.int64)
intervals = pygrit.from_numpy("chr1", arr)

# Convert back to NumPy
output = intervals.to_numpy()
```

## API Reference

### File-Based Functions

| Function | Description | bedtools equivalent |
|----------|-------------|---------------------|
| `intersect(a, b, ...)` | Find overlapping intervals | `bedtools intersect` |
| `merge(input, ...)` | Merge overlapping intervals | `bedtools merge` |
| `subtract(a, b, ...)` | Remove overlapping regions | `bedtools subtract` |
| `coverage(a, b, ...)` | Calculate coverage depth | `bedtools coverage` |
| `closest(a, b, ...)` | Find nearest intervals | `bedtools closest` |
| `window(a, b, ...)` | Find intervals within distance | `bedtools window` |
| `sort(input, ...)` | Sort BED file | `bedtools sort` |
| `slop(input, genome, ...)` | Extend interval boundaries | `bedtools slop` |
| `complement(input, genome)` | Find gaps between intervals | `bedtools complement` |
| `genomecov(input, genome, ...)` | Genome-wide coverage | `bedtools genomecov` |
| `jaccard(a, b)` | Jaccard similarity coefficient | `bedtools jaccard` |
| `multiinter(inputs, ...)` | Multi-file intersection | `bedtools multiinter` |
| `generate(output_dir, ...)` | Generate synthetic BED files | `bedtools random` |

### Classes

#### `Interval`
A single genomic interval (0-based, half-open).

```python
iv = pygrit.Interval("chr1", 100, 200)
iv.chrom          # "chr1"
iv.start          # 100
iv.end            # 200
len(iv)           # 100
iv.overlaps(iv2)  # True/False
iv.overlap_length(iv2)  # overlap bases
iv.distance_to(iv2)     # distance or None
```

#### `IntervalSet`
A collection of intervals with bulk operations.

```python
intervals = pygrit.read_bed("file.bed")
len(intervals)                    # count
intervals[0]                      # first interval
intervals.merge(distance=0)       # merge overlapping
intervals.intersect(other)        # find overlaps
intervals.sort()                  # sort in place
intervals.to_list()               # list of Interval
intervals.to_numpy()              # numpy array
```

### I/O Functions

- `read_bed(path)` - Read intervals from a BED file
- `parse_bed(content)` - Parse intervals from a string
- `from_numpy(chrom, arr)` - Create IntervalSet from NumPy array (n, 2)

## Performance

pygrit maintains the streaming architecture of GRIT:
- O(k) memory complexity where k = max overlapping intervals
- GIL released during computation for full parallelism
- Zero-copy file I/O when writing to output files
- Rayon-based parallel processing

## Requirements

- Python >= 3.9
- NumPy >= 1.20
- Input files should be sorted for streaming functions (use `grit sort` or `sort -k1,1 -k2,2n`)

## License

MIT
