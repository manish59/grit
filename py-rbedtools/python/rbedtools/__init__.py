"""
rbedtools: Fast, parallel genomic interval operations.

This module provides Python bindings for the rbedtools Rust library,
offering high-performance interval operations for genomic data analysis.

Example usage:

    >>> import rbedtools as bt
    >>>
    >>> # Read a BED file
    >>> intervals = bt.read_bed("regions.bed")
    >>>
    >>> # Merge overlapping intervals
    >>> merged = intervals.merge(distance=100)
    >>>
    >>> # Find intersections
    >>> other = bt.read_bed("other.bed")
    >>> overlaps = intervals.intersect(other, fraction=0.5)
    >>>
    >>> # Work with NumPy arrays
    >>> import numpy as np
    >>> arr = np.array([[0, 100], [150, 200]], dtype=np.int64)
    >>> intervals = bt.from_numpy("chr1", arr)
"""

from .rbedtools import (
    Interval,
    IntervalSet,
    read_bed,
    parse_bed,
    from_numpy,
    merge,
    intersect,
    __version__,
)

__all__ = [
    "Interval",
    "IntervalSet",
    "read_bed",
    "parse_bed",
    "from_numpy",
    "merge",
    "intersect",
    "__version__",
]
