# Interval

A genomic interval with chromosome, start, and end coordinates.

## Overview

The `Interval` class represents a single genomic region using 0-based, half-open coordinates (BED format).

```python
import pygrit

iv = pygrit.Interval("chr1", 100, 200)
```

## Constructor

### `Interval(chrom, start, end)`

Create a new interval.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `chrom` | `str` | Chromosome name (e.g., "chr1", "chrX") |
| `start` | `int` | Start position (0-based, inclusive) |
| `end` | `int` | End position (exclusive) |

**Returns:** `Interval`

**Raises:** `ValueError` if `start > end`

**Example:**

```python
# Create an interval covering positions 100-199
iv = pygrit.Interval("chr1", 100, 200)

# Invalid: start > end
try:
    pygrit.Interval("chr1", 200, 100)
except ValueError as e:
    print(e)  # "start (200) cannot be greater than end (100)"
```

---

## Properties

### `chrom`

The chromosome name.

**Type:** `str`

```python
iv = pygrit.Interval("chr1", 100, 200)
iv.chrom  # "chr1"
```

### `start`

The start position (0-based, inclusive).

**Type:** `int`

```python
iv = pygrit.Interval("chr1", 100, 200)
iv.start  # 100
```

### `end`

The end position (exclusive).

**Type:** `int`

```python
iv = pygrit.Interval("chr1", 100, 200)
iv.end  # 200
```

---

## Methods

### `overlaps(other)`

Check if this interval overlaps another.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `other` | `Interval` | The interval to check against |

**Returns:** `bool` - `True` if intervals overlap, `False` otherwise

**Example:**

```python
iv1 = pygrit.Interval("chr1", 100, 200)
iv2 = pygrit.Interval("chr1", 150, 250)
iv3 = pygrit.Interval("chr1", 300, 400)
iv4 = pygrit.Interval("chr2", 100, 200)

iv1.overlaps(iv2)  # True (overlap at 150-200)
iv1.overlaps(iv3)  # False (no overlap)
iv1.overlaps(iv4)  # False (different chromosomes)
```

---

### `overlap_length(other)`

Calculate the number of overlapping bases.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `other` | `Interval` | The interval to check against |

**Returns:** `int` - Number of overlapping bases (0 if no overlap)

**Example:**

```python
iv1 = pygrit.Interval("chr1", 100, 200)
iv2 = pygrit.Interval("chr1", 150, 250)

iv1.overlap_length(iv2)  # 50 (positions 150-199)
```

---

### `distance_to(other)`

Calculate the distance to another interval.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `other` | `Interval` | The interval to measure distance to |

**Returns:** `int | None`

- `0` if intervals overlap
- Positive integer for the gap between non-overlapping intervals
- `None` if intervals are on different chromosomes

**Example:**

```python
iv1 = pygrit.Interval("chr1", 100, 200)
iv2 = pygrit.Interval("chr1", 150, 250)
iv3 = pygrit.Interval("chr1", 300, 400)
iv4 = pygrit.Interval("chr2", 100, 200)

iv1.distance_to(iv2)  # 0 (overlapping)
iv1.distance_to(iv3)  # 100 (gap from 200 to 300)
iv1.distance_to(iv4)  # None (different chromosomes)
```

---

### `to_tuple()`

Convert to a tuple.

**Returns:** `tuple[str, int, int]` - `(chrom, start, end)`

**Example:**

```python
iv = pygrit.Interval("chr1", 100, 200)
iv.to_tuple()  # ("chr1", 100, 200)

# Useful for unpacking
chrom, start, end = iv.to_tuple()
```

---

## Special Methods

### `__len__()`

Return the length of the interval.

**Returns:** `int` - `end - start`

```python
iv = pygrit.Interval("chr1", 100, 200)
len(iv)  # 100
```

### `__str__()`

Return BED format string.

**Returns:** `str` - Tab-separated `chrom\tstart\tend`

```python
iv = pygrit.Interval("chr1", 100, 200)
str(iv)  # "chr1\t100\t200"
```

### `__repr__()`

Return Python representation.

**Returns:** `str` - `Interval('chrom', start, end)`

```python
iv = pygrit.Interval("chr1", 100, 200)
repr(iv)  # "Interval('chr1', 100, 200)"
```

### `__eq__(other)`

Check equality with another interval.

**Returns:** `bool`

```python
iv1 = pygrit.Interval("chr1", 100, 200)
iv2 = pygrit.Interval("chr1", 100, 200)
iv3 = pygrit.Interval("chr1", 100, 201)

iv1 == iv2  # True
iv1 == iv3  # False
```

### `__hash__()`

Return hash value for use in sets and dicts.

**Returns:** `int`

```python
iv1 = pygrit.Interval("chr1", 100, 200)
iv2 = pygrit.Interval("chr1", 100, 200)

# Use in sets
intervals = {iv1, iv2}
len(intervals)  # 1 (duplicates removed)

# Use as dict keys
data = {iv1: "exon1"}
data[iv2]  # "exon1"
```

---

## Usage Patterns

### Working with Multiple Intervals

```python
intervals = [
    pygrit.Interval("chr1", 100, 200),
    pygrit.Interval("chr1", 150, 250),
    pygrit.Interval("chr1", 300, 400),
]

# Find all overlapping pairs
for i, iv1 in enumerate(intervals):
    for iv2 in intervals[i+1:]:
        if iv1.overlaps(iv2):
            print(f"{iv1} overlaps {iv2}")
```

### Converting to Other Formats

```python
iv = pygrit.Interval("chr1", 100, 200)

# To BED line
bed_line = str(iv)

# To tuple
t = iv.to_tuple()

# To dict
d = {"chrom": iv.chrom, "start": iv.start, "end": iv.end}
```
