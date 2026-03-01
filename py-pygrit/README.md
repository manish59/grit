# grit-genomics

Python bindings for [GRIT](https://github.com/manish59/grit) (Genomic Range Interval Toolkit).

[![PyPI](https://img.shields.io/pypi/v/grit-genomics.svg)](https://pypi.org/project/grit-genomics/)
[![Python](https://img.shields.io/pypi/pyversions/grit-genomics.svg)](https://pypi.org/project/grit-genomics/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Installation

```bash
pip install grit-genomics
```

The package is imported as `pygrit`:

```python
import pygrit
print(pygrit.__version__)
```

### Building from Source

Requires [Rust](https://rustup.rs/) and [maturin](https://github.com/PyO3/maturin):

```bash
cd py-pygrit
pip install maturin
maturin develop --release
```

## Quick Start

```python
import pygrit

# Find overlapping intervals between two BED files
overlaps = pygrit.intersect("a.bed", "b.bed")
print(f"Found {len(overlaps)} overlapping intervals")

# Merge overlapping intervals
merged = pygrit.merge("input.bed")

# Write results to file
pygrit.intersect("a.bed", "b.bed", output="overlaps.bed")
```

## Available Functions

### File Operations

All functions accept file paths and optionally write to an output file.

| Function | Description | Requires |
|----------|-------------|----------|
| `intersect(a, b)` | Find overlapping intervals | Sorted input |
| `merge(input)` | Merge overlapping intervals | Sorted input |
| `subtract(a, b)` | Remove overlapping regions | Sorted input |
| `coverage(a, b)` | Calculate coverage depth | Sorted input |
| `closest(a, b)` | Find nearest intervals | Sorted input |
| `window(a, b)` | Find intervals within distance | Sorted input |
| `sort(input)` | Sort a BED file | - |
| `slop(input, genome)` | Extend interval boundaries | Genome file |
| `complement(input, genome)` | Find gaps between intervals | Genome file, sorted input |
| `genomecov(input, genome)` | Genome-wide coverage | Genome file, sorted input |
| `jaccard(a, b)` | Jaccard similarity coefficient | Sorted input |
| `multiinter(inputs)` | Multi-file intersection | Sorted input |
| `generate(output_dir)` | Generate synthetic BED files | - |

### Example Usage

```python
import pygrit

# intersect - find overlapping intervals
overlaps = pygrit.intersect("a.bed", "b.bed")
pygrit.intersect("a.bed", "b.bed", output="out.bed")
pygrit.intersect("a.bed", "b.bed", fraction=0.5, reciprocal=True)

# merge - combine overlapping intervals
merged = pygrit.merge("input.bed")
pygrit.merge("input.bed", distance=100)  # merge within 100bp

# subtract - remove overlapping regions
remaining = pygrit.subtract("a.bed", "b.bed")

# sort - sort a BED file
pygrit.sort("unsorted.bed", output="sorted.bed")

# slop - extend intervals (requires genome file)
pygrit.slop("regions.bed", "genome.txt", both=100.0)
pygrit.slop("regions.bed", "genome.txt", left=50.0, right=100.0)

# complement - find gaps
gaps = pygrit.complement("features.bed", "genome.txt")

# genomecov - coverage histogram or bedgraph
pygrit.genomecov("reads.bed", "genome.txt", bg=True)

# jaccard - similarity between two files
result = pygrit.jaccard("a.bed", "b.bed")

# multiinter - intersection across multiple files
pygrit.multiinter(["a.bed", "b.bed", "c.bed"])

# generate - create synthetic test data
stats = pygrit.generate("./test_data", num_intervals=10000, seed=42)
```

### In-Memory Operations

For programmatic interval manipulation:

```python
import pygrit

# Read intervals from file
intervals = pygrit.read_bed("regions.bed")
print(f"Loaded {len(intervals)} intervals")

# Parse from string
intervals = pygrit.parse_bed("chr1\t100\t200\nchr1\t150\t250\n")

# Create individual intervals
iv = pygrit.Interval("chr1", 100, 200)
print(iv.chrom, iv.start, iv.end)  # chr1 100 200
print(len(iv))  # 100

# Check overlap
iv2 = pygrit.Interval("chr1", 150, 250)
print(iv.overlaps(iv2))  # True
print(iv.overlap_length(iv2))  # 50

# Merge intervals
merged = intervals.merge(distance=0)

# Find intersections
other = pygrit.read_bed("other.bed")
overlaps = intervals.intersect(other)
```

### NumPy Integration

```python
import pygrit
import numpy as np

# Create from NumPy array (shape: n x 2)
arr = np.array([[100, 200], [300, 400]], dtype=np.int64)
intervals = pygrit.from_numpy("chr1", arr)

# Convert back to NumPy
output = intervals.to_numpy()
print(output.shape)  # (2, 2)
```

## Input Requirements

- **Sorted input**: Most functions require BED files sorted by chromosome and start position
- **Genome file**: Functions like `slop`, `complement`, `genomecov` require a genome file (tab-separated: chromosome name and size)

Example genome file:
```
chr1	248956422
chr2	242193529
chr3	198295559
```

Sort files with:
```bash
sort -k1,1 -k2,2n input.bed > sorted.bed
# or
grit sort -i input.bed > sorted.bed
```

## Requirements

- Python >= 3.9
- NumPy >= 1.20

## Platform Support

Pre-built wheels available for:
- macOS (ARM64)

For other platforms, pip will build from source (requires Rust toolchain).

## Testing

```bash
pip install pytest
pytest tests/
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- [GRIT](https://github.com/manish59/grit) - The underlying Rust implementation
- [bedtools](https://bedtools.readthedocs.io/) by Aaron Quinlan - The original genomic interval toolkit
- [PyO3](https://pyo3.rs/) - Rust bindings for Python
- [maturin](https://github.com/PyO3/maturin) - Build system for Rust Python extensions

## Links

- [PyPI Package](https://pypi.org/project/grit-genomics/)
- [GitHub Repository](https://github.com/manish59/grit)
- [GRIT Documentation](https://manish59.github.io/grit/)
