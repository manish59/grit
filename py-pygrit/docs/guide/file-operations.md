# File Operations

File-based operations are the core strength of pygrit. They use streaming algorithms that process data with **O(k) memory complexity**, where k is the maximum number of overlapping intervals at any genomic position.

## Requirements

!!! important "Sorted Input Required"
    All file-based operations require **sorted BED files**. Sort your files by chromosome and start position:
    ```bash
    sort -k1,1 -k2,2n input.bed > sorted.bed
    ```

## intersect

Find overlapping intervals between two files.

### Basic Usage

```python
import pygrit

# Return as list
overlaps = pygrit.intersect("features.bed", "regions.bed")

# Write to file (more memory efficient)
pygrit.intersect("features.bed", "regions.bed", output="overlaps.bed")
```

### Options

#### Write Original Records

```python
# Include original A record (-wa equivalent)
pygrit.intersect("a.bed", "b.bed", write_a=True, output="out.bed")

# Include both A and B records (-wa -wb equivalent)
pygrit.intersect("a.bed", "b.bed", write_a=True, write_b=True, output="out.bed")
```

#### Overlap Fraction

```python
# Require at least 50% of A to be overlapped
pygrit.intersect("a.bed", "b.bed", fraction=0.5, output="out.bed")

# Require reciprocal 50% overlap
pygrit.intersect("a.bed", "b.bed", fraction=0.5, reciprocal=True, output="out.bed")
```

#### Unique and Count

```python
# Report each A interval only once (-u equivalent)
pygrit.intersect("a.bed", "b.bed", unique=True, output="out.bed")

# Report overlap count for each A (-c equivalent)
pygrit.intersect("a.bed", "b.bed", count=True, output="out.bed")
```

#### Non-Overlapping

```python
# Report A intervals that DON'T overlap B (-v equivalent)
non_overlapping = pygrit.intersect("a.bed", "b.bed", no_overlap=True)
```

---

## merge

Merge overlapping or nearby intervals.

### Basic Usage

```python
# Merge overlapping intervals
merged = pygrit.merge("input.bed")

# Write to file
pygrit.merge("input.bed", output="merged.bed")
```

### Options

#### Distance

```python
# Merge intervals within 1000bp of each other
pygrit.merge("input.bed", distance=1000, output="merged.bed")
```

#### Strand-Specific

```python
# Merge only intervals on the same strand
pygrit.merge("input.bed", strand=True, output="merged.bed")
```

---

## subtract

Remove portions of A that overlap with B.

### Basic Usage

```python
# Subtract B from A
result = pygrit.subtract("features.bed", "exclude.bed")
```

### Options

#### Remove Entire Interval

```python
# Remove entire A interval if any overlap with B (-A equivalent)
pygrit.subtract("a.bed", "b.bed", remove_entire=True, output="out.bed")
```

#### Fraction Requirements

```python
# Only subtract if overlap >= 50% of A
pygrit.subtract("a.bed", "b.bed", fraction=0.5, output="out.bed")
```

---

## coverage

Calculate coverage of A by B.

### Basic Usage

```python
# Get coverage information
result = pygrit.coverage("regions.bed", "reads.bed")
```

### Options

#### Histogram

```python
# Report depth histogram
result = pygrit.coverage("regions.bed", "reads.bed", histogram=True)
```

#### Mean Depth

```python
# Report mean depth per region
result = pygrit.coverage("regions.bed", "reads.bed", mean=True)
```

---

## closest

Find the closest B interval for each A interval.

### Basic Usage

```python
result = pygrit.closest("queries.bed", "targets.bed")
```

### Options

#### Direction Control

```python
# Ignore overlapping intervals
result = pygrit.closest("a.bed", "b.bed", ignore_overlaps=True)

# Only find downstream (3') intervals
result = pygrit.closest("a.bed", "b.bed", ignore_upstream=True)

# Only find upstream (5') intervals
result = pygrit.closest("a.bed", "b.bed", ignore_downstream=True)
```

---

## window

Find B intervals within a window distance of A intervals.

### Basic Usage

```python
# Find B intervals within 1000bp of A intervals
result = pygrit.window("a.bed", "b.bed", window=1000)
```

### Options

#### Asymmetric Window

```python
# Different left and right windows
result = pygrit.window("a.bed", "b.bed", left=5000, right=1000)
```

#### Count Mode

```python
# Report count of B intervals in window
result = pygrit.window("a.bed", "b.bed", window=1000, count=True)
```

---

## Output Modes

### Return as List

When `output` is not specified, operations return a list of `Interval` objects:

```python
results = pygrit.intersect("a.bed", "b.bed")
for iv in results:
    print(f"{iv.chrom}:{iv.start}-{iv.end}")
```

### Write to File

When `output` is specified, results are written directly to the file:

```python
pygrit.intersect("a.bed", "b.bed", output="results.bed")
# Returns None, but results.bed contains the output
```

!!! tip "Memory Efficiency"
    Use `output=` for large result sets. Writing directly to a file avoids loading all results into memory.

---

## Error Handling

```python
try:
    results = pygrit.intersect("a.bed", "missing.bed")
except IOError as e:
    print(f"File error: {e}")
except RuntimeError as e:
    print(f"Processing error: {e}")  # e.g., unsorted input
except ValueError as e:
    print(f"Invalid parameter: {e}")
```
