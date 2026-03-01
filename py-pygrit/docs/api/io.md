# I/O Functions

Functions for reading and parsing genomic interval data.

---

## read_bed

```python
def read_bed(path: str) -> IntervalSet
```

Read intervals from a BED file.

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `path` | `str` | Path to the BED file |

### Returns

`IntervalSet` containing all intervals from the file.

### Raises

- `IOError`: File not found or I/O error

### Example

```python
import pygrit

# Read BED file
intervals = pygrit.read_bed("regions.bed")

print(f"Loaded {len(intervals)} intervals")

# Access intervals
first = intervals[0]
print(f"First: {first.chrom}:{first.start}-{first.end}")
```

### Notes

- Reads the first 3 columns (chrom, start, end)
- Additional columns are ignored
- Empty lines and comment lines (starting with #) are skipped

---

## parse_bed

```python
def parse_bed(content: str) -> IntervalSet
```

Parse intervals from a BED-formatted string.

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `content` | `str` | BED-formatted string content |

### Returns

`IntervalSet` containing parsed intervals.

### Example

```python
import pygrit

# Parse BED content from string
content = """chr1\t100\t200
chr1\t300\t400
chr2\t500\t600
"""

intervals = pygrit.parse_bed(content)
print(f"Parsed {len(intervals)} intervals")
```

### Use Cases

```python
# Parse from subprocess output
import subprocess

result = subprocess.run(
    ["some_tool", "output.bed"],
    capture_output=True,
    text=True
)
intervals = pygrit.parse_bed(result.stdout)

# Parse from HTTP response
import requests

response = requests.get("https://example.com/data.bed")
intervals = pygrit.parse_bed(response.text)

# Parse from multi-line string
content = "\n".join([
    "chr1\t100\t200",
    "chr1\t300\t400",
])
intervals = pygrit.parse_bed(content)
```

---

## from_numpy

```python
def from_numpy(chrom: str, arr: np.ndarray) -> IntervalSet
```

Create an IntervalSet from a NumPy array.

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `chrom` | `str` | Chromosome name for all intervals |
| `arr` | `np.ndarray` | Array with shape `(n, 2)` containing `[start, end]` pairs |

### Returns

`IntervalSet` with `n` intervals, all on the specified chromosome.

### Raises

- `ValueError`: If array shape is not `(n, 2)`

### Example

```python
import numpy as np
import pygrit

# Create from NumPy array
arr = np.array([
    [100, 200],
    [300, 400],
    [500, 600],
], dtype=np.int64)

intervals = pygrit.from_numpy("chr1", arr)
print(f"Created {len(intervals)} intervals")

# Verify
for iv in intervals.to_list():
    print(iv)
# chr1    100    200
# chr1    300    400
# chr1    500    600
```

### Use Cases

```python
# Generate intervals programmatically
import numpy as np

# Create sliding windows
window_size = 1000
step = 500
starts = np.arange(0, 100000, step)
ends = starts + window_size

arr = np.column_stack([starts, ends])
windows = pygrit.from_numpy("chr1", arr)

# From pandas DataFrame
import pandas as pd

df = pd.DataFrame({
    "start": [100, 300, 500],
    "end": [200, 400, 600],
})

arr = df[["start", "end"]].to_numpy()
intervals = pygrit.from_numpy("chr1", arr)

# Random intervals
np.random.seed(42)
n_intervals = 1000
starts = np.sort(np.random.randint(0, 1000000, n_intervals))
lengths = np.random.randint(100, 1000, n_intervals)
ends = starts + lengths

arr = np.column_stack([starts, ends])
random_intervals = pygrit.from_numpy("chr1", arr)
```

### Notes

- All intervals will have the same chromosome
- For multiple chromosomes, create separate IntervalSets and process them
- Array dtype should be integer (int32, int64, uint32, etc.)

---

## Round-Trip Example

```python
import numpy as np
import pygrit

# Read from file
intervals = pygrit.read_bed("regions.bed")

# Convert to NumPy for analysis
arr = intervals.to_numpy()

# Analyze
lengths = arr[:, 1] - arr[:, 0]
print(f"Mean length: {np.mean(lengths):.1f}")
print(f"Total coverage: {np.sum(lengths)}")

# Modify (extend by 100bp on each side)
arr[:, 0] = np.maximum(0, arr[:, 0] - 100)
arr[:, 1] = arr[:, 1] + 100

# Convert back (assuming single chromosome)
extended = pygrit.from_numpy("chr1", arr)

# Merge overlapping extended intervals
merged = extended.merge()
```
