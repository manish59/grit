---
layout: default
title: Performance Tuning
nav_order: 9
---

# Performance Tuning Guide

This guide helps you get the best performance from GRIT for your specific workload.

## Quick Decision Tree

```
Is your input sorted?
├── No → Use default mode (auto-sorts in memory)
└── Yes → Use --streaming --assume-sorted (fastest)

Does your data fit in memory?
├── No → Use --streaming mode
└── Yes → Use default parallel mode (faster)

Processing a pipeline?
├── Yes → Use --streaming with stdin/stdout
└── No → Use file paths for best I/O
```

## Mode Selection

### Parallel Mode (Default)

Best for:
- Files that fit in memory (< available RAM)
- Maximum speed on multi-core systems
- Unsorted input (will be sorted automatically)

```bash
grit intersect -a regions.bed -b features.bed > output.bed
```

Memory: O(n + m) where n, m are interval counts

### Streaming Mode

Best for:
- Large files exceeding available RAM
- Memory-constrained environments
- Pipeline processing (stdin/stdout)
- Already-sorted input

```bash
grit intersect -a large.bed -b large.bed --streaming --assume-sorted > output.bed
```

Memory: O(k) where k = max overlapping intervals (typically < 100)

## Benchmarks by Data Size

| File Size | Intervals | Recommended Mode | Expected Memory |
|-----------|-----------|------------------|-----------------|
| < 100 MB | < 1M | Parallel | < 500 MB |
| 100 MB - 1 GB | 1M - 10M | Either | 500 MB - 4 GB |
| > 1 GB | > 10M | Streaming | < 50 MB |

## Command-Specific Tips

### intersect

```bash
# Fastest for pre-sorted files
grit intersect -a sorted_a.bed -b sorted_b.bed --streaming --assume-sorted

# Use -u when you only need unique A intervals
grit intersect -a a.bed -b b.bed -u > unique_overlaps.bed

# Use -c for counting (faster than full output)
grit intersect -a a.bed -b b.bed -c > counts.bed
```

### merge

```bash
# Always streaming - very fast
grit merge -i input.bed > merged.bed

# Skip sort validation if pre-sorted
grit merge -i sorted.bed --assume-sorted
```

### sort

```bash
# Default: Fast parallel sort
grit sort -i input.bed > sorted.bed

# With genome file for chromosome ordering
grit sort -i input.bed -g genome.txt > sorted.bed
```

### coverage

```bash
# Streaming mode for large files
grit coverage -a regions.bed -b reads.bed --streaming > coverage.bed

# Mean coverage is faster than per-base
grit coverage -a regions.bed -b reads.bed --mean > mean_coverage.bed
```

## I/O Optimization

### Use Files Instead of Pipes for Large Data

```bash
# Slower: Multiple pipe operations
cat large.bed | grit sort | grit merge > output.bed

# Faster: Intermediate files (better I/O patterns)
grit sort -i large.bed > sorted.bed
grit merge -i sorted.bed > output.bed
```

### SSD vs HDD

- **SSD**: Use default settings
- **HDD**: Streaming mode benefits from sequential access patterns

### Compression

GRIT reads uncompressed BED files. Decompress first:

```bash
# Decompress then process
zcat input.bed.gz | grit intersect -a - -b features.bed > output.bed

# Or decompress to file (faster for repeated use)
gunzip -k input.bed.gz
grit intersect -a input.bed -b features.bed > output.bed
```

## Memory Management

### Estimating Memory Usage

**Parallel mode**:
```
Memory ≈ (n + m) × 100 bytes
```

Where n, m are interval counts.

Example: 10M intervals ≈ 1 GB RAM

**Streaming mode**:
```
Memory ≈ k × 100 bytes + buffer (< 10 MB)
```

Where k = max overlapping intervals.

Example: k=1000 ≈ 100 KB + 10 MB buffer

### Reducing Memory Usage

1. **Use streaming mode**:
   ```bash
   grit intersect -a large.bed -b large.bed --streaming --assume-sorted
   ```

2. **Process chromosomes separately**:
   ```bash
   for chr in chr1 chr2 chr3; do
       grep "^$chr\t" input.bed | grit merge > ${chr}_merged.bed
   done
   ```

3. **Filter early**:
   ```bash
   # Filter before expensive operations
   grep "^chr1\t" large.bed | grit intersect -a - -b features.bed
   ```

## CPU Optimization

### Thread Count

GRIT automatically uses all available cores. To limit:

```bash
# Limit to 4 threads (using RAYON_NUM_THREADS)
RAYON_NUM_THREADS=4 grit intersect -a a.bed -b b.bed
```

### When to Limit Threads

- Shared systems (be considerate)
- Memory-constrained (fewer threads = less memory)
- I/O-bound workloads (more threads won't help)

## Profiling Your Workload

### Timing

```bash
# Basic timing
time grit intersect -a a.bed -b b.bed > output.bed

# Detailed timing
/usr/bin/time -v grit intersect -a a.bed -b b.bed > output.bed 2>&1
```

### Memory Profiling (Linux)

```bash
# Peak memory usage
/usr/bin/time -v grit intersect -a a.bed -b b.bed 2>&1 | grep "Maximum resident"

# Detailed memory over time
valgrind --tool=massif grit intersect -a a.bed -b b.bed
```

### Memory Profiling (macOS)

```bash
# Using /usr/bin/time
/usr/bin/time -l grit intersect -a a.bed -b b.bed 2>&1 | grep "maximum resident"
```

## Common Performance Issues

### Slow: Unsorted Input with Streaming

**Symptom**: Streaming mode is slow or fails

**Cause**: Streaming requires sorted input; unsorted causes errors or re-sorting

**Fix**: Sort first or use parallel mode
```bash
grit sort -i unsorted.bed > sorted.bed
grit intersect -a sorted.bed -b sorted.bed --streaming --assume-sorted
```

### Slow: High Memory with Parallel Mode

**Symptom**: System swapping, slow performance

**Cause**: Files too large for available RAM

**Fix**: Use streaming mode
```bash
grit intersect -a large.bed -b large.bed --streaming --assume-sorted
```

### Slow: Many Small Files

**Symptom**: Slow when processing many files

**Cause**: I/O overhead dominates

**Fix**: Concatenate inputs first
```bash
cat file1.bed file2.bed file3.bed | grit sort > combined_sorted.bed
```

### Slow: Highly Overlapping Data

**Symptom**: Slow intersect operations

**Cause**: High overlap density increases output size

**Fix**: Use filtering options
```bash
# Only report unique A intervals
grit intersect -a a.bed -b b.bed -u > unique.bed

# Require minimum overlap
grit intersect -a a.bed -b b.bed -f 0.5 > filtered.bed
```

## Comparison with bedtools

GRIT is typically 3-15x faster than bedtools:

| Operation | GRIT | bedtools | Speedup |
|-----------|------|----------|---------|
| intersect (10M × 10M) | 12s | 45s | 3.8x |
| merge (10M intervals) | 3s | 18s | 6x |
| sort (10M intervals) | 5s | 25s | 5x |
| coverage (1M × 10M) | 8s | 95s | 12x |

*Benchmarks on AMD Ryzen 9 5900X, 32GB RAM, NVMe SSD*

### Why GRIT is Faster

1. **Rust vs C**: Modern optimizations, no runtime overhead
2. **Parallel by default**: Uses all CPU cores
3. **Cache-efficient**: Data structures optimized for CPU cache
4. **Zero-copy parsing**: Minimal memory allocations

## Python-Specific Performance

### Use Output Files for Large Results

```python
# Slower: Returns to Python (memory copy)
results = pygrit.intersect("a.bed", "b.bed")

# Faster: Writes directly to file
pygrit.intersect("a.bed", "b.bed", output="result.bed")
```

### GIL Release

GRIT releases Python's GIL during computation, enabling true parallelism:

```python
import threading
import pygrit

# These run in parallel (GIL released)
t1 = threading.Thread(target=pygrit.intersect, args=("a1.bed", "b.bed"),
                      kwargs={"output": "r1.bed"})
t2 = threading.Thread(target=pygrit.intersect, args=("a2.bed", "b.bed"),
                      kwargs={"output": "r2.bed"})
t1.start()
t2.start()
t1.join()
t2.join()
```

## Summary Checklist

- [ ] Input files are sorted (use `grit sort` if not)
- [ ] Using `--streaming --assume-sorted` for large files
- [ ] Using parallel mode for files that fit in RAM
- [ ] Filtering early to reduce data size
- [ ] Using appropriate output options (`-u`, `-c`, etc.)
- [ ] Using output files instead of returning to Python for large results
