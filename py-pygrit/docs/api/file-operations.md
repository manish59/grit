# File Operations

Streaming functions for processing BED files with O(k) memory complexity.

!!! important "Sorted Input Required"
    All file operations require sorted BED files (by chromosome, then start position).

---

## intersect

```python
def intersect(
    a: str,
    b: str,
    output: str | None = None,
    write_a: bool = False,
    write_b: bool = False,
    fraction: float | None = None,
    reciprocal: bool = False,
    count: bool = False,
    unique: bool = False,
    no_overlap: bool = False,
) -> list[Interval] | None
```

Find overlapping intervals between two BED files.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `a` | `str` | - | Path to file A |
| `b` | `str` | - | Path to file B |
| `output` | `str \| None` | `None` | Output file path. If None, returns list |
| `write_a` | `bool` | `False` | Include original A record in output |
| `write_b` | `bool` | `False` | Include original B record in output |
| `fraction` | `float \| None` | `None` | Minimum overlap as fraction of A (0.0-1.0) |
| `reciprocal` | `bool` | `False` | Require fraction overlap in both A and B |
| `count` | `bool` | `False` | Report count of overlaps per A interval |
| `unique` | `bool` | `False` | Report each A interval only once |
| `no_overlap` | `bool` | `False` | Report A intervals with no B overlap |

### Returns

- `list[Interval]` if `output` is None
- `None` if `output` is specified (results written to file)

### Raises

- `IOError`: File not found or I/O error
- `RuntimeError`: Processing error (e.g., unsorted input)
- `ValueError`: Invalid parameter values

### Examples

```python
import pygrit

# Basic intersection
overlaps = pygrit.intersect("features.bed", "regions.bed")

# Write to file
pygrit.intersect("a.bed", "b.bed", output="overlaps.bed")

# Include original records
pygrit.intersect("a.bed", "b.bed", write_a=True, write_b=True, output="out.bed")

# Require 50% overlap
pygrit.intersect("a.bed", "b.bed", fraction=0.5, output="out.bed")

# Reciprocal 50% overlap
pygrit.intersect("a.bed", "b.bed", fraction=0.5, reciprocal=True, output="out.bed")

# Report each A once
pygrit.intersect("a.bed", "b.bed", unique=True, output="out.bed")

# Count overlaps
pygrit.intersect("a.bed", "b.bed", count=True, output="out.bed")

# Non-overlapping A intervals
non_overlapping = pygrit.intersect("a.bed", "b.bed", no_overlap=True)
```

---

## merge

```python
def merge(
    input: str,
    output: str | None = None,
    distance: int = 0,
    strand: bool = False,
) -> list[Interval] | None
```

Merge overlapping or nearby intervals.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input` | `str` | - | Path to input BED file |
| `output` | `str \| None` | `None` | Output file path. If None, returns list |
| `distance` | `int` | `0` | Maximum gap to bridge when merging |
| `strand` | `bool` | `False` | Only merge intervals on same strand |

### Returns

- `list[Interval]` if `output` is None
- `None` if `output` is specified

### Examples

```python
# Basic merge
merged = pygrit.merge("input.bed")

# Merge intervals within 100bp
merged = pygrit.merge("input.bed", distance=100)

# Strand-specific merge
merged = pygrit.merge("input.bed", strand=True)

# Write to file
pygrit.merge("input.bed", output="merged.bed")
```

---

## subtract

```python
def subtract(
    a: str,
    b: str,
    output: str | None = None,
    remove_entire: bool = False,
    fraction: float | None = None,
    reciprocal: bool = False,
) -> list[Interval] | None
```

Subtract B intervals from A intervals.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `a` | `str` | - | Path to file A |
| `b` | `str` | - | Path to file B |
| `output` | `str \| None` | `None` | Output file path. If None, returns list |
| `remove_entire` | `bool` | `False` | Remove entire A interval if any overlap |
| `fraction` | `float \| None` | `None` | Minimum overlap fraction to subtract |
| `reciprocal` | `bool` | `False` | Require reciprocal fraction |

### Returns

- `list[Interval]` if `output` is None
- `None` if `output` is specified

### Examples

```python
# Basic subtraction
result = pygrit.subtract("features.bed", "exclude.bed")

# Remove entire interval on any overlap
result = pygrit.subtract("a.bed", "b.bed", remove_entire=True)

# Only subtract significant overlaps
result = pygrit.subtract("a.bed", "b.bed", fraction=0.5)
```

---

## coverage

```python
def coverage(
    a: str,
    b: str,
    output: str | None = None,
    histogram: bool = False,
    mean: bool = False,
) -> str | None
```

Calculate coverage of A regions by B features.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `a` | `str` | - | Path to regions file |
| `b` | `str` | - | Path to reads/features file |
| `output` | `str \| None` | `None` | Output file path. If None, returns string |
| `histogram` | `bool` | `False` | Report depth histogram |
| `mean` | `bool` | `False` | Report mean depth per region |

### Returns

- `str` if `output` is None (coverage output as string)
- `None` if `output` is specified

### Examples

```python
# Basic coverage
result = pygrit.coverage("regions.bed", "reads.bed")

# With mean depth
result = pygrit.coverage("regions.bed", "reads.bed", mean=True)

# Histogram output
result = pygrit.coverage("regions.bed", "reads.bed", histogram=True)
```

---

## closest

```python
def closest(
    a: str,
    b: str,
    output: str | None = None,
    ignore_overlaps: bool = False,
    ignore_upstream: bool = False,
    ignore_downstream: bool = False,
) -> str | None
```

Find closest B interval for each A interval.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `a` | `str` | - | Path to file A |
| `b` | `str` | - | Path to file B |
| `output` | `str \| None` | `None` | Output file path. If None, returns string |
| `ignore_overlaps` | `bool` | `False` | Skip overlapping intervals |
| `ignore_upstream` | `bool` | `False` | Only look downstream (3') |
| `ignore_downstream` | `bool` | `False` | Only look upstream (5') |

### Returns

- `str` if `output` is None (closest output as string)
- `None` if `output` is specified

### Examples

```python
# Find closest
result = pygrit.closest("queries.bed", "targets.bed")

# Ignore overlapping
result = pygrit.closest("a.bed", "b.bed", ignore_overlaps=True)

# Only downstream
result = pygrit.closest("a.bed", "b.bed", ignore_upstream=True)
```

---

## window

```python
def window(
    a: str,
    b: str,
    output: str | None = None,
    window: int = 1000,
    left: int | None = None,
    right: int | None = None,
    count: bool = False,
    no_overlap: bool = False,
) -> str | None
```

Find B intervals within window distance of A intervals.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `a` | `str` | - | Path to file A |
| `b` | `str` | - | Path to file B |
| `output` | `str \| None` | `None` | Output file path. If None, returns string |
| `window` | `int` | `1000` | Window size in base pairs |
| `left` | `int \| None` | `None` | Left window (overrides `window`) |
| `right` | `int \| None` | `None` | Right window (overrides `window`) |
| `count` | `bool` | `False` | Report count of B in window |
| `no_overlap` | `bool` | `False` | Only report non-overlapping |

### Returns

- `str` if `output` is None (window output as string)
- `None` if `output` is specified

### Examples

```python
# Find within 1kb
result = pygrit.window("genes.bed", "enhancers.bed", window=1000)

# Asymmetric window (5kb upstream, 1kb downstream)
result = pygrit.window("genes.bed", "enhancers.bed", left=5000, right=1000)

# Count mode
result = pygrit.window("a.bed", "b.bed", count=True)
```
