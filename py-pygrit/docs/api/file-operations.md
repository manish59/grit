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

---

## sort

```python
def sort(
    input: str,
    output: str | None = None,
    genome: str | None = None,
    reverse: bool = False,
) -> str | None
```

Sort a BED file by chromosome and start position.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input` | `str` | - | Path to input BED file |
| `output` | `str \| None` | `None` | Output file path. If None, returns string |
| `genome` | `str \| None` | `None` | Genome file for chromosome ordering |
| `reverse` | `bool` | `False` | Sort in reverse order |

### Returns

- `str` if `output` is None (sorted output as string)
- `None` if `output` is specified

### Examples

```python
# Basic sort
sorted_bed = pygrit.sort("unsorted.bed")

# Write to file
pygrit.sort("unsorted.bed", output="sorted.bed")

# Reverse sort
pygrit.sort("input.bed", reverse=True, output="reversed.bed")

# Use genome file for chromosome ordering
pygrit.sort("input.bed", genome="genome.txt", output="sorted.bed")
```

---

## slop

```python
def slop(
    input: str,
    genome: str,
    output: str | None = None,
    both: float | None = None,
    left: float | None = None,
    right: float | None = None,
    pct: bool = False,
) -> str | None
```

Extend intervals by a specified amount on each side.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input` | `str` | - | Path to input BED file |
| `genome` | `str` | - | Path to genome file (chromosome sizes) |
| `output` | `str \| None` | `None` | Output file path. If None, returns string |
| `both` | `float \| None` | `None` | Extend both sides by this amount |
| `left` | `float \| None` | `None` | Extend left (5') side |
| `right` | `float \| None` | `None` | Extend right (3') side |
| `pct` | `bool` | `False` | Interpret values as percentage of interval length |

### Returns

- `str` if `output` is None (extended intervals as string)
- `None` if `output` is specified

### Examples

```python
# Extend both sides by 100bp
result = pygrit.slop("regions.bed", "genome.txt", both=100.0)

# Asymmetric extension
result = pygrit.slop("regions.bed", "genome.txt", left=50.0, right=200.0)

# Percentage-based extension (50% of interval length)
result = pygrit.slop("regions.bed", "genome.txt", both=0.5, pct=True)

# Write to file
pygrit.slop("regions.bed", "genome.txt", both=100.0, output="extended.bed")
```

---

## complement

```python
def complement(
    input: str,
    genome: str,
    output: str | None = None,
) -> str | None
```

Calculate the complement of intervals (gaps between intervals).

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input` | `str` | - | Path to input BED file |
| `genome` | `str` | - | Path to genome file (chromosome sizes) |
| `output` | `str \| None` | `None` | Output file path. If None, returns string |

### Returns

- `str` if `output` is None (complement intervals as string)
- `None` if `output` is specified

### Examples

```python
# Get gaps between intervals
gaps = pygrit.complement("features.bed", "genome.txt")

# Write to file
pygrit.complement("features.bed", "genome.txt", output="gaps.bed")
```

---

## genomecov

```python
def genomecov(
    input: str,
    genome: str,
    output: str | None = None,
    bg: bool = False,
    bga: bool = False,
    scale: float = 1.0,
) -> str | None
```

Calculate genome-wide coverage.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input` | `str` | - | Path to input BED file |
| `genome` | `str` | - | Path to genome file (chromosome sizes) |
| `output` | `str \| None` | `None` | Output file path. If None, returns string |
| `bg` | `bool` | `False` | Output BedGraph format (non-zero regions only) |
| `bga` | `bool` | `False` | Output BedGraph format (all regions, including zero) |
| `scale` | `float` | `1.0` | Scale coverage by this factor |

### Returns

- `str` if `output` is None (coverage output as string)
- `None` if `output` is specified

### Examples

```python
# Histogram output (default)
result = pygrit.genomecov("reads.bed", "genome.txt")

# BedGraph format (non-zero only)
result = pygrit.genomecov("reads.bed", "genome.txt", bg=True)

# BedGraph with zero coverage regions
result = pygrit.genomecov("reads.bed", "genome.txt", bga=True)

# Scaled coverage
result = pygrit.genomecov("reads.bed", "genome.txt", bg=True, scale=0.5)
```

---

## jaccard

```python
def jaccard(
    a: str,
    b: str,
    output: str | None = None,
) -> str | None
```

Calculate Jaccard similarity between two BED files.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `a` | `str` | - | Path to file A |
| `b` | `str` | - | Path to file B |
| `output` | `str \| None` | `None` | Output file path. If None, returns string |

### Returns

- `str` if `output` is None (Jaccard statistics as string)
- `None` if `output` is specified

### Examples

```python
# Calculate Jaccard similarity
result = pygrit.jaccard("set_a.bed", "set_b.bed")
print(result)  # intersection, union, jaccard, n_intersections

# Write to file
pygrit.jaccard("set_a.bed", "set_b.bed", output="similarity.txt")
```

---

## multiinter

```python
def multiinter(
    inputs: list[str],
    output: str | None = None,
    cluster: bool = False,
) -> str | None
```

Find intervals that overlap across multiple BED files.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `inputs` | `list[str]` | - | List of paths to input BED files (minimum 2) |
| `output` | `str \| None` | `None` | Output file path. If None, returns string |
| `cluster` | `bool` | `False` | Cluster overlapping intervals |

### Returns

- `str` if `output` is None (multiinter output as string)
- `None` if `output` is specified

### Raises

- `ValueError`: If fewer than 2 input files provided

### Examples

```python
# Find regions covered by multiple files
result = pygrit.multiinter(["a.bed", "b.bed", "c.bed"])

# With clustering
result = pygrit.multiinter(["a.bed", "b.bed"], cluster=True)

# Write to file
pygrit.multiinter(["a.bed", "b.bed", "c.bed"], output="overlap.bed")
```

---

## generate

```python
def generate(
    output_dir: str,
    num_intervals: int = 10000,
    num_chroms: int = 5,
    chrom_size: int = 100000000,
    len_min: int = 50,
    len_max: int = 5000,
    mode: str = "uniform",
    num_files: int = 2,
    sorted: bool = True,
    seed: int | None = None,
) -> dict
```

Generate synthetic BED files for testing and benchmarking.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_dir` | `str` | - | Output directory path |
| `num_intervals` | `int` | `10000` | Number of intervals per file |
| `num_chroms` | `int` | `5` | Number of chromosomes |
| `chrom_size` | `int` | `100000000` | Size of each chromosome |
| `len_min` | `int` | `50` | Minimum interval length |
| `len_max` | `int` | `5000` | Maximum interval length |
| `mode` | `str` | `"uniform"` | Distribution mode: "uniform", "balanced", "clustered" |
| `num_files` | `int` | `2` | Number of files to generate |
| `sorted` | `bool` | `True` | Generate sorted output |
| `seed` | `int \| None` | `None` | Random seed for reproducibility |

### Returns

- `dict` with generation statistics:
    - `total_files`: Number of files generated
    - `total_intervals`: Total intervals across all files
    - `elapsed_secs`: Time taken in seconds

### Raises

- `ValueError`: If invalid mode specified

### Examples

```python
# Generate test files
stats = pygrit.generate("./test_data", num_intervals=10000)
print(f"Generated {stats['total_files']} files")

# Clustered distribution
stats = pygrit.generate("./data", mode="clustered", seed=42)

# Reproducible generation
stats = pygrit.generate("./data", num_intervals=5000, seed=12345)
```
