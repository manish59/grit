# NumPy Integration

pygrit provides seamless integration with NumPy for numerical analysis of genomic intervals.

## Converting Between Formats

### IntervalSet to NumPy

```python
import numpy as np
import pygrit

# Read intervals
intervals = pygrit.read_bed("regions.bed")

# Convert to NumPy array
arr = intervals.to_numpy()
print(arr.shape)  # (n_intervals, 2)
print(arr.dtype)  # int64

# Columns are [start, end]
starts = arr[:, 0]
ends = arr[:, 1]
```

### NumPy to IntervalSet

```python
import numpy as np
import pygrit

# Create array
arr = np.array([
    [100, 200],
    [300, 400],
    [500, 600],
], dtype=np.int64)

# Convert to IntervalSet
intervals = pygrit.from_numpy("chr1", arr)
```

---

## Statistical Analysis

### Basic Statistics

```python
import numpy as np
import pygrit

intervals = pygrit.read_bed("peaks.bed")
arr = intervals.to_numpy()

# Calculate lengths
lengths = arr[:, 1] - arr[:, 0]

print(f"Count: {len(lengths)}")
print(f"Mean length: {np.mean(lengths):.1f} bp")
print(f"Median length: {np.median(lengths):.1f} bp")
print(f"Std dev: {np.std(lengths):.1f} bp")
print(f"Min: {np.min(lengths)} bp")
print(f"Max: {np.max(lengths)} bp")
print(f"Total coverage: {np.sum(lengths):,} bp")
```

### Percentiles

```python
import numpy as np
import pygrit

intervals = pygrit.read_bed("peaks.bed")
arr = intervals.to_numpy()
lengths = arr[:, 1] - arr[:, 0]

# Percentiles
percentiles = [10, 25, 50, 75, 90]
values = np.percentile(lengths, percentiles)

for p, v in zip(percentiles, values):
    print(f"P{p}: {v:.0f} bp")
```

### Histogram

```python
import numpy as np
import pygrit

intervals = pygrit.read_bed("peaks.bed")
arr = intervals.to_numpy()
lengths = arr[:, 1] - arr[:, 0]

# Create histogram
hist, bin_edges = np.histogram(lengths, bins=50)

print("Length Distribution:")
for i, count in enumerate(hist):
    if count > 0:
        print(f"  {bin_edges[i]:.0f}-{bin_edges[i+1]:.0f}: {count}")
```

---

## Filtering and Selection

### Filter by Length

```python
import numpy as np
import pygrit

intervals = pygrit.read_bed("regions.bed")
arr = intervals.to_numpy()

# Calculate lengths
lengths = arr[:, 1] - arr[:, 0]

# Filter: keep intervals > 1000 bp
mask = lengths > 1000
filtered_arr = arr[mask]

# Convert back
long_intervals = pygrit.from_numpy("chr1", filtered_arr)
print(f"Kept {len(long_intervals)} / {len(intervals)} intervals")
```

### Filter by Position

```python
import numpy as np
import pygrit

intervals = pygrit.read_bed("regions.bed")
arr = intervals.to_numpy()

# Keep intervals in a specific region
region_start = 1_000_000
region_end = 2_000_000

mask = (arr[:, 0] >= region_start) & (arr[:, 1] <= region_end)
region_intervals = pygrit.from_numpy("chr1", arr[mask])
```

### Top N by Length

```python
import numpy as np
import pygrit

intervals = pygrit.read_bed("peaks.bed")
arr = intervals.to_numpy()

# Get top 100 longest intervals
lengths = arr[:, 1] - arr[:, 0]
top_indices = np.argsort(lengths)[-100:][::-1]  # Top 100, descending

top_intervals = pygrit.from_numpy("chr1", arr[top_indices])
```

---

## Transformations

### Extend Intervals

```python
import numpy as np
import pygrit

intervals = pygrit.read_bed("peaks.bed")
arr = intervals.to_numpy()

# Extend by 500bp on each side
extension = 500
arr[:, 0] = np.maximum(0, arr[:, 0] - extension)  # Don't go below 0
arr[:, 1] = arr[:, 1] + extension

extended = pygrit.from_numpy("chr1", arr)
# Merge overlapping after extension
merged = extended.merge()
```

### Resize to Fixed Width

```python
import numpy as np
import pygrit

intervals = pygrit.read_bed("peaks.bed")
arr = intervals.to_numpy()

# Resize all intervals to 1000bp centered on midpoint
target_width = 1000
midpoints = (arr[:, 0] + arr[:, 1]) // 2
half_width = target_width // 2

new_arr = np.column_stack([
    np.maximum(0, midpoints - half_width),
    midpoints + half_width,
])

resized = pygrit.from_numpy("chr1", new_arr)
```

### Shift Intervals

```python
import numpy as np
import pygrit

intervals = pygrit.read_bed("regions.bed")
arr = intervals.to_numpy()

# Shift all intervals by 1000bp
shift = 1000
arr[:, 0] += shift
arr[:, 1] += shift

shifted = pygrit.from_numpy("chr1", arr)
```

---

## Generating Intervals

### Sliding Windows

```python
import numpy as np
import pygrit

def create_sliding_windows(chrom: str, chrom_length: int,
                           window_size: int, step: int) -> pygrit.IntervalSet:
    """Create sliding windows across a chromosome."""
    starts = np.arange(0, chrom_length - window_size + 1, step)
    ends = starts + window_size

    arr = np.column_stack([starts, ends])
    return pygrit.from_numpy(chrom, arr)

# Create 1kb windows with 500bp step
windows = create_sliding_windows("chr1", 10_000_000, 1000, 500)
print(f"Created {len(windows)} windows")
```

### Random Intervals

```python
import numpy as np
import pygrit

def create_random_intervals(chrom: str, chrom_length: int,
                            n_intervals: int, min_length: int,
                            max_length: int, seed: int = None) -> pygrit.IntervalSet:
    """Create random non-overlapping intervals."""
    if seed is not None:
        np.random.seed(seed)

    # Generate random starts
    starts = np.sort(np.random.randint(0, chrom_length - max_length, n_intervals))

    # Generate random lengths
    lengths = np.random.randint(min_length, max_length + 1, n_intervals)
    ends = starts + lengths

    # Ensure ends don't exceed chromosome length
    ends = np.minimum(ends, chrom_length)

    arr = np.column_stack([starts, ends])
    intervals = pygrit.from_numpy(chrom, arr)

    # Merge overlapping
    return intervals.merge()

# Create 1000 random intervals
random_intervals = create_random_intervals(
    "chr1", 100_000_000, 1000, 100, 1000, seed=42
)
```

### Tiled Regions

```python
import numpy as np
import pygrit

def tile_region(chrom: str, start: int, end: int, tile_size: int) -> pygrit.IntervalSet:
    """Tile a region with non-overlapping intervals."""
    starts = np.arange(start, end, tile_size)
    ends = np.minimum(starts + tile_size, end)

    arr = np.column_stack([starts, ends])
    return pygrit.from_numpy(chrom, arr)

# Tile a 1Mb region with 10kb tiles
tiles = tile_region("chr1", 1_000_000, 2_000_000, 10_000)
print(f"Created {len(tiles)} tiles")
```

---

## Integration with Pandas

### To DataFrame

```python
import numpy as np
import pandas as pd
import pygrit

intervals = pygrit.read_bed("regions.bed")

# Convert to DataFrame
data = []
for iv in intervals.to_list():
    data.append({
        "chrom": iv.chrom,
        "start": iv.start,
        "end": iv.end,
        "length": len(iv),
    })

df = pd.DataFrame(data)
print(df.head())
```

### From DataFrame

```python
import numpy as np
import pandas as pd
import pygrit

# DataFrame with interval data
df = pd.DataFrame({
    "chrom": ["chr1", "chr1", "chr1"],
    "start": [100, 300, 500],
    "end": [200, 400, 600],
})

# Group by chromosome and convert
all_intervals = []
for chrom, group in df.groupby("chrom"):
    arr = group[["start", "end"]].to_numpy()
    intervals = pygrit.from_numpy(chrom, arr)
    all_intervals.extend(intervals.to_list())

# Create combined IntervalSet
combined = pygrit.IntervalSet.from_intervals(all_intervals)
```

---

## Performance Tips

### Minimize Conversions

```python
# Bad: Converting repeatedly
for _ in range(1000):
    arr = intervals.to_numpy()  # Converts each time
    # ... do something ...

# Good: Convert once
arr = intervals.to_numpy()
for _ in range(1000):
    # ... do something with arr ...
```

### Use Vectorized Operations

```python
import numpy as np

# Bad: Python loop
lengths = []
for iv in intervals.to_list():
    lengths.append(len(iv))

# Good: Vectorized NumPy
arr = intervals.to_numpy()
lengths = arr[:, 1] - arr[:, 0]
```

### Pre-allocate Arrays

```python
import numpy as np

# When building intervals programmatically
n = 10000
arr = np.empty((n, 2), dtype=np.int64)

for i in range(n):
    arr[i, 0] = i * 1000
    arr[i, 1] = i * 1000 + 500

intervals = pygrit.from_numpy("chr1", arr)
```
