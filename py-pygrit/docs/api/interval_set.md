# IntervalSet

A collection of genomic intervals with bulk operations.

## Overview

The `IntervalSet` class provides efficient storage and operations for multiple intervals.

```python
import pygrit

# Create from list
intervals = pygrit.IntervalSet.from_intervals([
    pygrit.Interval("chr1", 100, 200),
    pygrit.Interval("chr1", 300, 400),
])
```

## Constructors

### `IntervalSet()`

Create an empty IntervalSet.

**Returns:** `IntervalSet`

```python
intervals = pygrit.IntervalSet()
len(intervals)  # 0
```

---

### `IntervalSet.from_intervals(intervals)`

Create an IntervalSet from a list of Interval objects.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `intervals` | `list[Interval]` | List of Interval objects |

**Returns:** `IntervalSet`

```python
intervals = pygrit.IntervalSet.from_intervals([
    pygrit.Interval("chr1", 100, 200),
    pygrit.Interval("chr1", 300, 400),
])
```

---

## Methods

### `add(interval)`

Add an interval to the set.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `interval` | `Interval` | The interval to add |

**Returns:** `None`

```python
intervals = pygrit.IntervalSet()
intervals.add(pygrit.Interval("chr1", 100, 200))
intervals.add(pygrit.Interval("chr1", 300, 400))
len(intervals)  # 2
```

---

### `to_list()`

Convert to a list of Interval objects.

**Returns:** `list[Interval]`

```python
intervals = pygrit.read_bed("regions.bed")
for iv in intervals.to_list():
    print(f"{iv.chrom}:{iv.start}-{iv.end}")
```

---

### `merge(distance=0)`

Merge overlapping or nearby intervals.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `distance` | `int` | `0` | Maximum gap to bridge when merging |

**Returns:** `IntervalSet` - New IntervalSet with merged intervals

**Example:**

```python
intervals = pygrit.IntervalSet.from_intervals([
    pygrit.Interval("chr1", 100, 200),
    pygrit.Interval("chr1", 150, 250),  # Overlaps
    pygrit.Interval("chr1", 500, 600),  # Separate
])

# Merge overlapping only
merged = intervals.merge()
# Result: [chr1:100-250, chr1:500-600]

# Merge with distance tolerance
merged = intervals.merge(distance=300)
# Result: [chr1:100-600]
```

---

### `intersect(other, fraction=None, reciprocal=False)`

Find intervals that overlap with another IntervalSet.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `other` | `IntervalSet` | - | IntervalSet to intersect with |
| `fraction` | `float \| None` | `None` | Minimum overlap fraction (0.0-1.0) |
| `reciprocal` | `bool` | `False` | Require reciprocal overlap fraction |

**Returns:** `IntervalSet` - Intervals from self that overlap with other

**Example:**

```python
set_a = pygrit.IntervalSet.from_intervals([
    pygrit.Interval("chr1", 100, 200),  # 100bp
    pygrit.Interval("chr1", 300, 400),  # 100bp
])

set_b = pygrit.IntervalSet.from_intervals([
    pygrit.Interval("chr1", 150, 350),  # Overlaps both
])

# Basic intersection
overlapping = set_a.intersect(set_b)
# Both intervals overlap

# With fraction requirement (50% of A must overlap)
overlapping = set_a.intersect(set_b, fraction=0.5)
# First: 50bp overlap / 100bp = 50% ✓
# Second: 50bp overlap / 100bp = 50% ✓

# With reciprocal fraction
overlapping = set_a.intersect(set_b, fraction=0.3, reciprocal=True)
# Also requires 30% of B to overlap
```

---

### `non_overlapping(other)`

Find intervals that don't overlap with another IntervalSet.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `other` | `IntervalSet` | IntervalSet to check against |

**Returns:** `IntervalSet` - Intervals from self that don't overlap with other

**Example:**

```python
set_a = pygrit.IntervalSet.from_intervals([
    pygrit.Interval("chr1", 100, 200),
    pygrit.Interval("chr1", 300, 400),
    pygrit.Interval("chr1", 500, 600),
])

set_b = pygrit.IntervalSet.from_intervals([
    pygrit.Interval("chr1", 150, 250),  # Overlaps first
])

non_overlapping = set_a.non_overlapping(set_b)
# Result: [chr1:300-400, chr1:500-600]
```

---

### `sort()`

Sort intervals in place by chromosome and start position.

**Returns:** `None`

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

### `to_numpy()`

Convert to a NumPy array.

**Returns:** `np.ndarray` - Shape `(n, 2)` with columns `[start, end]`

!!! note
    This only returns start/end coordinates. Chromosome information is not included in the array.

```python
import numpy as np

intervals = pygrit.IntervalSet.from_intervals([
    pygrit.Interval("chr1", 100, 200),
    pygrit.Interval("chr1", 300, 400),
])

arr = intervals.to_numpy()
# array([[100, 200],
#        [300, 400]])

# Calculate lengths
lengths = arr[:, 1] - arr[:, 0]
```

---

## Special Methods

### `__len__()`

Return the number of intervals.

**Returns:** `int`

```python
intervals = pygrit.read_bed("regions.bed")
len(intervals)  # Number of intervals
```

### `__getitem__(index)`

Access interval by index.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `index` | `int` | Index (supports negative indexing) |

**Returns:** `Interval`

**Raises:** `IndexError` if index is out of range

```python
intervals = pygrit.read_bed("regions.bed")
first = intervals[0]
last = intervals[-1]
```

---

## Usage Patterns

### Iteration

```python
intervals = pygrit.read_bed("regions.bed")

# Using to_list()
for iv in intervals.to_list():
    print(iv)

# Using indexing
for i in range(len(intervals)):
    print(intervals[i])
```

### Filtering

```python
intervals = pygrit.read_bed("regions.bed")

# Filter by chromosome
chr1_intervals = [
    iv for iv in intervals.to_list()
    if iv.chrom == "chr1"
]

# Filter by length
long_intervals = [
    iv for iv in intervals.to_list()
    if len(iv) > 1000
]

# Create new IntervalSet from filtered
filtered = pygrit.IntervalSet.from_intervals(long_intervals)
```

### Chaining Operations

```python
intervals = pygrit.read_bed("regions.bed")

# Sort, merge, then filter
intervals.sort()
merged = intervals.merge(distance=100)

# Further operations
final = merged.non_overlapping(exclusions)
```
