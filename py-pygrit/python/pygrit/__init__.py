"""
pygrit: High-performance genomic interval operations.

Python bindings for GRIT (Genomic Range Interval Toolkit), providing
streaming algorithms for large-scale genomic data analysis.

Example usage:

    >>> import pygrit
    >>>
    >>> # File-based streaming (recommended for large files)
    >>> results = pygrit.intersect("a.bed", "b.bed")
    >>> pygrit.merge("input.bed", output="merged.bed", distance=100)
    >>>
    >>> # In-memory operations
    >>> intervals = pygrit.read_bed("regions.bed")
    >>> merged = intervals.merge(distance=100)
    >>>
    >>> # Work with individual intervals
    >>> iv = pygrit.Interval("chr1", 100, 200)
    >>> print(len(iv))  # 100
    >>>
    >>> # NumPy integration
    >>> import numpy as np
    >>> arr = np.array([[0, 100], [150, 200]], dtype=np.int64)
    >>> intervals = pygrit.from_numpy("chr1", arr)
"""

from pygrit.pygrit import (
    # Core types
    Interval,
    IntervalSet,
    # File-based streaming functions
    intersect,
    merge,
    subtract,
    coverage,
    closest,
    window,
    sort,
    slop,
    complement,
    genomecov,
    jaccard,
    multiinter,
    generate,
    # I/O utilities
    read_bed,
    parse_bed,
    from_numpy,
    # Metadata
    __version__,
)

__all__ = [
    # Core types
    "Interval",
    "IntervalSet",
    # File-based streaming functions
    "intersect",
    "merge",
    "subtract",
    "coverage",
    "closest",
    "window",
    "sort",
    "slop",
    "complement",
    "genomecov",
    "jaccard",
    "multiinter",
    "generate",
    # I/O utilities
    "read_bed",
    "parse_bed",
    "from_numpy",
    # Metadata
    "__version__",
]
