# API Reference

Complete API documentation for pygrit.

## Overview

pygrit provides two main ways to work with genomic intervals:

### Core Types

- **[`Interval`](interval.md)**: A single genomic interval with chromosome, start, and end
- **[`IntervalSet`](interval_set.md)**: A collection of intervals with bulk operations

### File-Based Operations

Functions that process BED files:

- **[`intersect`](file-operations.md#intersect)**: Find overlapping intervals
- **[`merge`](file-operations.md#merge)**: Merge overlapping intervals
- **[`subtract`](file-operations.md#subtract)**: Subtract intervals
- **[`coverage`](file-operations.md#coverage)**: Calculate coverage
- **[`closest`](file-operations.md#closest)**: Find closest intervals
- **[`window`](file-operations.md#window)**: Find intervals within a window
- **[`sort`](file-operations.md#sort)**: Sort BED files
- **[`slop`](file-operations.md#slop)**: Extend interval boundaries
- **[`complement`](file-operations.md#complement)**: Find gaps between intervals
- **[`genomecov`](file-operations.md#genomecov)**: Genome-wide coverage
- **[`jaccard`](file-operations.md#jaccard)**: Jaccard similarity
- **[`multiinter`](file-operations.md#multiinter)**: Multi-file intersection
- **[`generate`](file-operations.md#generate)**: Generate synthetic data

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

# Core types
pygrit.Interval
pygrit.IntervalSet

# File operations
pygrit.intersect
pygrit.merge
pygrit.subtract
pygrit.coverage
pygrit.closest
pygrit.window
pygrit.sort
pygrit.slop
pygrit.complement
pygrit.genomecov
pygrit.jaccard
pygrit.multiinter
pygrit.generate

# I/O
pygrit.read_bed
pygrit.parse_bed
pygrit.from_numpy

# Metadata
pygrit.__version__
```

## Input Requirements

!!! warning "Sorted Input Required"
    Most file-based functions require **sorted BED files** (by chromosome, then start position).

    **Unsorted input will produce incorrect results without warning.**

    Always sort first:
    ```python
    pygrit.sort("unsorted.bed", output="sorted.bed")
    ```

Functions requiring a **genome file**: `slop`, `complement`, `genomecov`

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
