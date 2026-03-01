# In-Memory Operations

For smaller datasets or when you need programmatic access to interval data, pygrit provides in-memory classes: `Interval` and `IntervalSet`.

## Interval Class

The `Interval` class represents a single genomic interval.

### Creating Intervals

```python
import pygrit

# Create an interval
iv = pygrit.Interval("chr1", 100, 200)

# Access properties
print(iv.chrom)  # "chr1"
print(iv.start)  # 100
print(iv.end)    # 200
print(len(iv))   # 100 (length in base pairs)
```

### Interval Methods

#### Overlap Detection

```python
iv1 = pygrit.Interval("chr1", 100, 200)
iv2 = pygrit.Interval("chr1", 150, 250)
iv3 = pygrit.Interval("chr1", 300, 400)

# Check if intervals overlap
iv1.overlaps(iv2)  # True
iv1.overlaps(iv3)  # False

# Get overlap length
iv1.overlap_length(iv2)  # 50
iv1.overlap_length(iv3)  # 0
```

#### Distance Calculation

```python
iv1 = pygrit.Interval("chr1", 100, 200)
iv2 = pygrit.Interval("chr1", 300, 400)
iv3 = pygrit.Interval("chr2", 100, 200)

# Distance between intervals
iv1.distance_to(iv2)  # 100 (gap between them)
iv1.distance_to(iv3)  # None (different chromosomes)

# Overlapping intervals have distance 0
iv1.distance_to(pygrit.Interval("chr1", 150, 250))  # 0
```

### Conversion Methods

```python
iv = pygrit.Interval("chr1", 100, 200)

# To tuple
iv.to_tuple()  # ("chr1", 100, 200)

# String representations
str(iv)   # "chr1\t100\t200" (BED format)
repr(iv)  # "Interval('chr1', 100, 200)"
```

### Comparison and Hashing

Intervals are hashable and can be used in sets and dictionaries:

```python
iv1 = pygrit.Interval("chr1", 100, 200)
iv2 = pygrit.Interval("chr1", 100, 200)
iv3 = pygrit.Interval("chr1", 150, 250)

# Equality
iv1 == iv2  # True
iv1 == iv3  # False

# Use in sets
unique_intervals = {iv1, iv2, iv3}
len(unique_intervals)  # 2

# Use as dict keys
coverage = {iv1: 10.5, iv3: 8.2}
```

---

## IntervalSet Class

The `IntervalSet` class represents a collection of intervals with bulk operations.

### Creating IntervalSets

```python
# Empty set
intervals = pygrit.IntervalSet()

# From list of intervals
intervals = pygrit.IntervalSet.from_intervals([
    pygrit.Interval("chr1", 100, 200),
    pygrit.Interval("chr1", 150, 250),
    pygrit.Interval("chr1", 300, 400),
])

# From BED file
intervals = pygrit.read_bed("regions.bed")

# From BED string
content = "chr1\t100\t200\nchr1\t300\t400\n"
intervals = pygrit.parse_bed(content)
```

### Basic Operations

```python
intervals = pygrit.read_bed("regions.bed")

# Length
len(intervals)  # Number of intervals

# Indexing
first = intervals[0]
last = intervals[-1]

# Iteration
for iv in intervals.to_list():
    print(iv)
```

### Adding Intervals

```python
intervals = pygrit.IntervalSet()
intervals.add(pygrit.Interval("chr1", 100, 200))
intervals.add(pygrit.Interval("chr1", 300, 400))
```

### Merge

Merge overlapping or nearby intervals:

```python
intervals = pygrit.IntervalSet.from_intervals([
    pygrit.Interval("chr1", 100, 200),
    pygrit.Interval("chr1", 150, 250),  # Overlaps previous
    pygrit.Interval("chr1", 500, 600),  # Separate
])

# Merge overlapping
merged = intervals.merge()
# Result: [chr1:100-250, chr1:500-600]

# Merge with distance tolerance
merged = intervals.merge(distance=300)
# Result: [chr1:100-600] (all within 300bp)
```

### Intersect

Find intervals that overlap with another set:

```python
set_a = pygrit.IntervalSet.from_intervals([
    pygrit.Interval("chr1", 100, 200),
    pygrit.Interval("chr1", 300, 400),
])

set_b = pygrit.IntervalSet.from_intervals([
    pygrit.Interval("chr1", 150, 350),
])

# Find overlapping intervals from A
overlapping = set_a.intersect(set_b)

# With fraction requirement
overlapping = set_a.intersect(set_b, fraction=0.5)

# Reciprocal fraction
overlapping = set_a.intersect(set_b, fraction=0.5, reciprocal=True)
```

### Non-Overlapping

Find intervals that don't overlap with another set:

```python
non_overlapping = set_a.non_overlapping(set_b)
```

### Sort

Sort intervals in place by chromosome and start position:

```python
intervals = pygrit.IntervalSet.from_intervals([
    pygrit.Interval("chr2", 100, 200),
    pygrit.Interval("chr1", 300, 400),
    pygrit.Interval("chr1", 100, 200),
])

intervals.sort()
# Now: [chr1:100-200, chr1:300-400, chr2:100-200]
```

---

## NumPy Integration

### From NumPy to Intervals

```python
import numpy as np

# Create intervals from NumPy array
arr = np.array([
    [100, 200],
    [300, 400],
    [500, 600],
], dtype=np.int64)

intervals = pygrit.from_numpy("chr1", arr)
```

### From Intervals to NumPy

```python
intervals = pygrit.read_bed("regions.bed")

# Convert to NumPy array (shape: n x 2)
arr = intervals.to_numpy()

# Calculate statistics
lengths = arr[:, 1] - arr[:, 0]
print(f"Mean length: {np.mean(lengths):.1f}")
print(f"Total coverage: {np.sum(lengths)}")
```

---

## When to Use In-Memory vs File-Based

| Use Case | Recommended Approach |
|----------|---------------------|
| Large files (>1M intervals) | File-based operations |
| Small datasets (<100K intervals) | Either works well |
| Complex manipulation | In-memory (IntervalSet) |
| Pipeline processing | File-based (streaming) |
| Interactive analysis | In-memory (IntervalSet) |
| NumPy integration | In-memory with to_numpy() |
