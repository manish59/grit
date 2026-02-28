# pygrit API Reference

## Module: pygrit

Python bindings for GRIT (Genomic Range Interval Toolkit).

---

## Classes

### `Interval`

A genomic interval with chromosome, start, and end coordinates.

Coordinates are 0-based, half-open (BED format).

```python
class Interval:
    chrom: str  # Chromosome name
    start: int  # Start position (0-based, inclusive)
    end: int    # End position (exclusive)
```

#### Constructor

```python
Interval(chrom: str, start: int, end: int) -> Interval
```

Creates a new interval. Raises `ValueError` if start > end.

#### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `__len__()` | `int` | Length of interval (end - start) |
| `__repr__()` | `str` | `"Interval('chr1', 100, 200)"` |
| `__str__()` | `str` | BED format: `"chr1\t100\t200"` |
| `__eq__(other)` | `bool` | Equality comparison |
| `__hash__()` | `int` | Hash for use in sets/dicts |
| `overlaps(other)` | `bool` | True if intervals overlap |
| `overlap_length(other)` | `int` | Number of overlapping bases |
| `distance_to(other)` | `int | None` | Distance (0 if overlapping, None if different chromosomes) |
| `to_tuple()` | `tuple[str, int, int]` | `(chrom, start, end)` |

#### Example

```python
iv1 = pygrit.Interval("chr1", 100, 200)
iv2 = pygrit.Interval("chr1", 150, 250)

len(iv1)              # 100
iv1.overlaps(iv2)     # True
iv1.overlap_length(iv2)  # 50
iv1.distance_to(iv2)  # 0 (overlapping)
```

---

### `IntervalSet`

A collection of genomic intervals with bulk operations.

#### Constructor

```python
IntervalSet() -> IntervalSet
IntervalSet.from_intervals(intervals: list[Interval]) -> IntervalSet
```

#### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `__len__()` | `int` | Number of intervals |
| `__getitem__(idx)` | `Interval` | Access by index |
| `add(interval)` | `None` | Add an interval |
| `to_list()` | `list[Interval]` | Convert to list |
| `merge(distance=0)` | `IntervalSet` | Merge overlapping intervals |
| `intersect(other, fraction=None, reciprocal=False)` | `IntervalSet` | Find intersecting intervals |
| `non_overlapping(other)` | `IntervalSet` | Find intervals with no overlap |
| `sort()` | `None` | Sort in place by chrom, start |
| `to_numpy()` | `np.ndarray` | Shape (n, 2) array of [start, end] |

#### Example

```python
# Create from list
intervals = IntervalSet.from_intervals([
    Interval("chr1", 100, 200),
    Interval("chr1", 150, 250),
])

# Merge overlapping
merged = intervals.merge()  # 1 interval: 100-250

# With distance
merged = intervals.merge(distance=100)  # Merge if gap <= 100bp
```

---

## File-Based Streaming Functions

These functions operate on BED files using streaming algorithms with O(k) memory complexity.

### `intersect`

```python
def intersect(
    a: str,                    # Path to file A
    b: str,                    # Path to file B
    output: str | None = None, # Output file path (returns list if None)
    write_a: bool = False,     # Include original A record (-wa)
    write_b: bool = False,     # Include original B record (-wb)
    fraction: float | None = None,  # Minimum overlap fraction (-f)
    reciprocal: bool = False,  # Require reciprocal fraction (-r)
    count: bool = False,       # Report count (-c)
    unique: bool = False,      # Report each A once (-u)
    no_overlap: bool = False,  # Report non-overlapping A (-v)
) -> list[Interval] | None
```

Find overlapping intervals between two BED files.

**Returns:** List of Interval objects if `output` is None, otherwise None.

**Example:**

```python
# Get overlapping intervals
results = pygrit.intersect("a.bed", "b.bed")

# Write to file
pygrit.intersect("a.bed", "b.bed", output="out.bed")

# With 50% overlap requirement
pygrit.intersect("a.bed", "b.bed", fraction=0.5, reciprocal=True)

# Get non-overlapping A intervals
pygrit.intersect("a.bed", "b.bed", no_overlap=True)
```

---

### `merge`

```python
def merge(
    input: str,                # Path to input BED file
    output: str | None = None, # Output file path
    distance: int = 0,         # Max distance to merge (-d)
    strand: bool = False,      # Strand-specific merge (-s)
) -> list[Interval] | None
```

Merge overlapping intervals.

**Example:**

```python
# Merge overlapping intervals
merged = pygrit.merge("input.bed")

# Merge intervals within 100bp
merged = pygrit.merge("input.bed", distance=100)

# Strand-specific merge
merged = pygrit.merge("input.bed", strand=True)
```

---

### `subtract`

```python
def subtract(
    a: str,                    # Path to file A
    b: str,                    # Path to file B
    output: str | None = None, # Output file path
    remove_entire: bool = False,  # Remove entire A if any overlap (-A)
    fraction: float | None = None,  # Minimum overlap fraction
    reciprocal: bool = False,  # Require reciprocal fraction
) -> list[Interval] | None
```

Subtract B intervals from A intervals.

**Example:**

```python
# Subtract B from A
result = pygrit.subtract("a.bed", "b.bed")

# Remove entire A interval if any overlap
result = pygrit.subtract("a.bed", "b.bed", remove_entire=True)
```

---

### `coverage`

```python
def coverage(
    a: str,                    # Path to regions file
    b: str,                    # Path to reads/features file
    output: str | None = None, # Output file path
    histogram: bool = False,   # Report depth histogram
    mean: bool = False,        # Report mean depth
) -> str | None
```

Calculate coverage of A regions by B features.

**Example:**

```python
# Get coverage output
result = pygrit.coverage("regions.bed", "reads.bed")

# With mean depth
result = pygrit.coverage("regions.bed", "reads.bed", mean=True)
```

---

### `closest`

```python
def closest(
    a: str,                    # Path to file A
    b: str,                    # Path to file B
    output: str | None = None, # Output file path
    ignore_overlaps: bool = False,    # Ignore overlapping (-io)
    ignore_upstream: bool = False,    # Ignore upstream (-iu)
    ignore_downstream: bool = False,  # Ignore downstream (-id)
) -> str | None
```

Find closest B interval for each A interval.

**Example:**

```python
# Find closest
result = pygrit.closest("a.bed", "b.bed")

# Only downstream
result = pygrit.closest("a.bed", "b.bed", ignore_upstream=True)
```

---

### `window`

```python
def window(
    a: str,                    # Path to file A
    b: str,                    # Path to file B
    output: str | None = None, # Output file path
    window: int = 1000,        # Window size (bp)
    left: int | None = None,   # Left window (overrides window)
    right: int | None = None,  # Right window (overrides window)
    count: bool = False,       # Report count
    no_overlap: bool = False,  # Only non-overlapping
) -> str | None
```

Find B intervals within window distance of A intervals.

**Example:**

```python
# Find within 1kb window
result = pygrit.window("a.bed", "b.bed", window=1000)

# Asymmetric window
result = pygrit.window("a.bed", "b.bed", left=500, right=0)
```

---

## I/O Functions

### `read_bed`

```python
def read_bed(path: str) -> IntervalSet
```

Read intervals from a BED file.

### `parse_bed`

```python
def parse_bed(content: str) -> IntervalSet
```

Parse intervals from a BED-formatted string.

### `from_numpy`

```python
def from_numpy(chrom: str, arr: np.ndarray) -> IntervalSet
```

Create IntervalSet from NumPy array. Array shape must be (n, 2) with [start, end] columns.

---

## Constants

### `__version__`

Package version string.

```python
>>> pygrit.__version__
'0.1.0'
```

---

## Performance Notes

1. **Memory Complexity**: O(k) where k = maximum overlapping intervals at any point
2. **GIL Release**: All heavy computations release the GIL for threading
3. **Sorted Input**: File-based functions require sorted input for streaming
4. **File Output**: Using `output=` parameter is more memory-efficient than returning lists

## Error Handling

- `ValueError`: Invalid parameters (e.g., start > end)
- `IOError`: File not found or I/O errors
- `RuntimeError`: Processing errors (e.g., unsorted input)
