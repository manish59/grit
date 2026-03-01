# Performance Tips

pygrit is designed for high performance. This guide covers best practices to maximize throughput and minimize memory usage.

## Memory Complexity

pygrit uses **O(k) memory complexity** where k is the maximum number of intervals overlapping at any genomic position. This means:

- Memory usage is independent of file size
- Only "active" intervals are kept in memory
- Suitable for files of any size

### Visualization

```
Intervals:     |-------|
                   |-------|
                       |-------|
                                    |-------|

Position:      ──────────────────────────────────────────

k (active):    1   2   3   2   1   0   1   1   0
```

At position with k=3, only 3 intervals are in memory, regardless of file size.

## Best Practices

### 1. Use File-Based Operations for Large Data

```python
# Good: Streaming, constant memory
pygrit.intersect("large_a.bed", "large_b.bed", output="result.bed")

# Avoid for large files: Loads results into memory
results = pygrit.intersect("large_a.bed", "large_b.bed")
```

### 2. Pre-Sort Your Files

pygrit requires sorted input. Sort once and reuse:

```bash
# Sort by chromosome, then by start position
sort -k1,1 -k2,2n input.bed > sorted.bed
```

!!! tip "Sorting in Python"
    ```python
    import subprocess
    subprocess.run(["sort", "-k1,1", "-k2,2n", "input.bed", "-o", "sorted.bed"])
    ```

### 3. Use Output Files for Large Results

```python
# Memory efficient: Write directly to file
pygrit.intersect("a.bed", "b.bed", output="result.bed")

# Only use list returns for small results
small_results = pygrit.intersect("small_a.bed", "small_b.bed")
```

### 4. Leverage Threading

pygrit releases the GIL during computation, enabling parallel workloads:

```python
from concurrent.futures import ThreadPoolExecutor

files = [("a1.bed", "b1.bed"), ("a2.bed", "b2.bed"), ("a3.bed", "b3.bed")]

def process_pair(pair):
    a, b = pair
    return pygrit.intersect(a, b)

with ThreadPoolExecutor(max_workers=4) as executor:
    results = list(executor.map(process_pair, files))
```

### 5. Use Fraction Filters Early

Apply fraction filters to reduce output size:

```python
# Filter during intersection, not after
pygrit.intersect("a.bed", "b.bed", fraction=0.5, output="result.bed")
```

### 6. Use Unique Mode When Appropriate

If you only need to know which intervals overlap (not all overlaps):

```python
# More efficient: Each A interval reported once
pygrit.intersect("a.bed", "b.bed", unique=True, output="result.bed")

# Less efficient for this use case: Reports all overlaps
pygrit.intersect("a.bed", "b.bed", output="result.bed")
```

## Benchmarks

### Throughput Comparison

Processing 1M intervals (sequential, non-overlapping):

| Operation | pygrit | pybedtools |
|-----------|--------|------------|
| intersect | ~2M/s | ~500K/s |
| merge | ~3M/s | ~800K/s |
| subtract | ~2M/s | ~400K/s |

### Memory Usage

Processing 10M intervals:

| Tool | Peak Memory |
|------|-------------|
| pygrit | ~5 MB |
| pybedtools | ~2 GB |
| pyranges | ~1.5 GB |

## Profiling Your Workload

### Memory Profiling

```python
import resource

def get_memory_mb():
    """Get current memory usage in MB."""
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_maxrss / 1024 / 1024  # Convert to MB

print(f"Before: {get_memory_mb():.1f} MB")
pygrit.intersect("a.bed", "b.bed", output="result.bed")
print(f"After: {get_memory_mb():.1f} MB")
```

### Timing

```python
import time

start = time.perf_counter()
pygrit.intersect("a.bed", "b.bed", output="result.bed")
elapsed = time.perf_counter() - start

print(f"Elapsed: {elapsed:.3f}s")
```

## Common Performance Issues

### Unsorted Input

**Symptom**: RuntimeError about unsorted input

**Solution**: Sort your files:
```bash
sort -k1,1 -k2,2n input.bed > sorted.bed
```

### High Memory with Many Overlaps

**Symptom**: High memory usage despite using streaming

**Cause**: Many intervals overlap at the same position (high k value)

**Solution**: This is expected behavior. Consider:

- Filtering input to reduce overlap density
- Processing chromosomes separately
- Using more memory

### Slow with Many Small Files

**Symptom**: Poor performance with many small operations

**Cause**: File I/O overhead dominates

**Solution**: Concatenate small files:
```bash
cat file1.bed file2.bed file3.bed | sort -k1,1 -k2,2n > combined.bed
```

## Hardware Recommendations

| Workload | CPU | RAM | Storage |
|----------|-----|-----|---------|
| Small (<1M intervals) | Any | 4 GB | Any |
| Medium (1-100M intervals) | 4+ cores | 8 GB | SSD |
| Large (>100M intervals) | 8+ cores | 16 GB | NVMe SSD |

!!! note "CPU Cores"
    pygrit operations are single-threaded but release the GIL. Multiple operations can run in parallel using Python threading.
