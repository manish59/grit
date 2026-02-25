# rbedtools Python Bindings

Fast, parallel genomic interval operations powered by Rust.

## Installation

```bash
pip install rbedtools
```

Or build from source:

```bash
cd py-rbedtools
maturin develop --release
```

## Usage

```python
import rbedtools as bt
import numpy as np

# Read a BED file
intervals = bt.read_bed("regions.bed")

# Create intervals from a string
intervals = bt.parse_bed("""
chr1	100	200
chr1	150	250
chr1	300	400
""")

# Merge overlapping intervals
merged = intervals.merge(distance=100)

# Find intersections
other = bt.read_bed("other.bed")
overlaps = intervals.intersect(other, fraction=0.5)

# Find non-overlapping intervals
non_overlapping = intervals.non_overlapping(other)

# Work with NumPy arrays
arr = np.array([[0, 100], [150, 200]], dtype=np.int64)
intervals = bt.from_numpy("chr1", arr)

# Convert back to NumPy
output = intervals.to_numpy()

# Access individual intervals
for interval in intervals.to_list():
    print(f"{interval.chrom}:{interval.start}-{interval.end}")
```

## API Reference

### Classes

#### `Interval`
A single genomic interval.

- `Interval(chrom: str, start: int, end: int)` - Create a new interval
- `.chrom` - Chromosome name
- `.start` - Start position (0-based)
- `.end` - End position (exclusive)
- `.overlaps(other)` - Check if intervals overlap
- `.overlap_length(other)` - Get overlap length
- `.distance_to(other)` - Get distance to another interval

#### `IntervalSet`
A collection of intervals.

- `.add(interval)` - Add an interval
- `.merge(distance=0)` - Merge overlapping intervals
- `.intersect(other, fraction=None, reciprocal=False)` - Find intersections
- `.non_overlapping(other)` - Find intervals with no overlap
- `.sort()` - Sort intervals in place
- `.to_list()` - Convert to list of Interval objects
- `.to_numpy()` - Convert to NumPy array

### Functions

- `read_bed(path)` - Read intervals from a BED file
- `parse_bed(content)` - Parse intervals from a string
- `from_numpy(chrom, arr)` - Create IntervalSet from NumPy array
- `merge(intervals, distance=0)` - Merge overlapping intervals
- `intersect(a, b, fraction=None, reciprocal=False)` - Find intersections

## Performance

rbedtools uses parallel processing via Rayon to achieve significant speedups
on multi-core systems. Operations are typically 3-8x faster than the original
bedtools for large datasets.

## License

MIT
