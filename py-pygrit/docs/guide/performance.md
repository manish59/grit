# Performance Tips

This guide covers best practices for using pygrit effectively.

## Best Practices

### 1. Sort Your Input Files

!!! warning "Critical: Sorted Input Required"
    Most pygrit functions require **sorted BED files** (sorted by chromosome, then by start position).

    **Unsorted input will produce incorrect results without warning.**

    pygrit uses streaming algorithms that process intervals in a single pass. This approach
    requires sorted input to function correctly. Unlike the CLI which validates sort order
    by default, the Python API assumes your input is already sorted for maximum performance.

Sort your files first:

=== "Python"
    ```python
    pygrit.sort("unsorted.bed", output="sorted.bed")
    ```

=== "grit CLI"
    ```bash
    grit sort -i unsorted.bed > sorted.bed
    ```

=== "Unix"
    ```bash
    sort -k1,1 -k2,2n input.bed > sorted.bed
    ```

### 2. Use Output Files for Large Results

When processing large files, write directly to output files instead of returning results to Python:

```python
# Recommended for large files: Write directly to file
pygrit.intersect("a.bed", "b.bed", output="result.bed")

# For small results: Return to Python
small_results = pygrit.intersect("small_a.bed", "small_b.bed")
```

### 3. Apply Filters Early

Use filters during operations to reduce output size:

```python
# Filter during intersection
pygrit.intersect("a.bed", "b.bed", fraction=0.5, output="result.bed")

# Use unique mode when you only need to know which intervals overlap
pygrit.intersect("a.bed", "b.bed", unique=True, output="result.bed")
```

### 4. Use the Right Tool for the Job

| Data Size | Recommendation |
|-----------|----------------|
| Small (<100K intervals) | In-memory operations work well |
| Medium (100K - 10M) | File-based operations |
| Large (>10M) | File-based with output files |

## Profiling Your Workload

### Memory Profiling

```python
import tracemalloc

tracemalloc.start()
pygrit.intersect("a.bed", "b.bed", output="result.bed")
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

print(f"Peak memory: {peak / 1024 / 1024:.1f} MB")
```

### Timing

```python
import time

start = time.perf_counter()
pygrit.intersect("a.bed", "b.bed", output="result.bed")
elapsed = time.perf_counter() - start

print(f"Elapsed: {elapsed:.3f}s")
```

## Common Issues

### Unsorted Input

**Symptom**: RuntimeError about unsorted input

**Solution**: Sort your files first:

```python
pygrit.sort("unsorted.bed", output="sorted.bed")
```

### Missing Genome File

**Symptom**: Error about missing genome file

**Solution**: Create a genome file (tab-separated chromosome name and size):

```
chr1	248956422
chr2	242193529
chr3	198295559
```

Functions that require genome files: `slop`, `complement`, `genomecov`

### High Memory Usage

**Symptom**: High memory usage when processing large files

**Solutions**:

1. Use `output=` parameter to write directly to files
2. Process chromosomes separately
3. Use smaller input files

## Benchmark Your Own Data

Performance varies significantly based on:

- Data size and distribution
- Overlap density
- Hardware (CPU, disk speed, RAM)

We recommend benchmarking on your own datasets:

```python
import time
import tempfile

def benchmark(func, *args, **kwargs):
    start = time.perf_counter()
    func(*args, **kwargs)
    return time.perf_counter() - start

with tempfile.NamedTemporaryFile(suffix=".bed") as tmp:
    elapsed = benchmark(
        pygrit.intersect,
        "a.bed", "b.bed",
        output=tmp.name
    )
    print(f"intersect: {elapsed:.3f}s")
```
