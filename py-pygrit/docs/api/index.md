# API Reference

Complete API documentation for pygrit.

## Overview

pygrit provides two main ways to work with genomic intervals:

### Core Types

- **[`Interval`](interval.md)**: A single genomic interval with chromosome, start, and end
- **[`IntervalSet`](interval_set.md)**: A collection of intervals with bulk operations

### File-Based Operations

Streaming functions that process BED files with O(k) memory:

- **[`intersect`](file-operations.md#intersect)**: Find overlapping intervals
- **[`merge`](file-operations.md#merge)**: Merge overlapping intervals
- **[`subtract`](file-operations.md#subtract)**: Subtract intervals
- **[`coverage`](file-operations.md#coverage)**: Calculate coverage
- **[`closest`](file-operations.md#closest)**: Find closest intervals
- **[`window`](file-operations.md#window)**: Find intervals within a window

### I/O Functions

- **[`read_bed`](io.md#read_bed)**: Read intervals from a BED file
- **[`parse_bed`](io.md#parse_bed)**: Parse intervals from a string
- **[`from_numpy`](io.md#from_numpy)**: Create intervals from NumPy array

## Quick Reference

### Coordinate System

pygrit uses **0-based, half-open** coordinates (BED format):

- **Start**: 0-based, inclusive
- **End**: 0-based, exclusive
- **Length**: `end - start`

```python
# An interval covering bases 100-199 (100 bases)
iv = pygrit.Interval("chr1", 100, 200)
len(iv)  # 100
```

### Import

```python
import pygrit

# All public API
pygrit.Interval
pygrit.IntervalSet
pygrit.intersect
pygrit.merge
pygrit.subtract
pygrit.coverage
pygrit.closest
pygrit.window
pygrit.read_bed
pygrit.parse_bed
pygrit.from_numpy
pygrit.__version__
```

## Type Summary

| Type | Description |
|------|-------------|
| `Interval` | Single genomic interval |
| `IntervalSet` | Collection of intervals |
| `str` | File paths, chromosome names |
| `int` | Positions, distances |
| `float` | Fractions (0.0-1.0) |
| `bool` | Option flags |
| `np.ndarray` | NumPy arrays for bulk data |

## Exceptions

| Exception | When Raised |
|-----------|-------------|
| `ValueError` | Invalid parameters (e.g., start > end) |
| `IOError` | File not found or I/O errors |
| `RuntimeError` | Processing errors (e.g., unsorted input) |

## Thread Safety

All file-based operations release the GIL during computation, making them safe for use with Python threading:

```python
from concurrent.futures import ThreadPoolExecutor

with ThreadPoolExecutor(max_workers=4) as executor:
    futures = [
        executor.submit(pygrit.intersect, f"a{i}.bed", f"b{i}.bed")
        for i in range(10)
    ]
    results = [f.result() for f in futures]
```
