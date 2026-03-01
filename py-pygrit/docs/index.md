# grit-genomics

**Python bindings for GRIT (Genomic Range Interval Toolkit).**

[![PyPI](https://img.shields.io/pypi/v/grit-genomics.svg)](https://pypi.org/project/grit-genomics/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

grit-genomics provides Python bindings for [GRIT](https://github.com/manish59/grit), a genomic interval toolkit written in Rust. The package is imported as `pygrit`.

## Features

- **13 Commands**: intersect, merge, subtract, closest, window, coverage, sort, slop, complement, genomecov, jaccard, multiinter, generate
- **File Operations**: Process BED files with streaming algorithms
- **In-Memory Operations**: Work with intervals programmatically
- **NumPy Integration**: Convert between intervals and NumPy arrays

## Quick Example

```python
import pygrit

# Find overlapping intervals between two BED files
results = pygrit.intersect("genes.bed", "peaks.bed")

# Merge overlapping intervals
merged = pygrit.merge("regions.bed", distance=100)

# Sort a BED file
pygrit.sort("unsorted.bed", output="sorted.bed")

# Work with individual intervals
iv = pygrit.Interval("chr1", 1000, 2000)
print(f"Length: {len(iv)} bp")  # Length: 1000 bp
```

## Installation

```bash
pip install grit-genomics
```

```python
import pygrit
print(pygrit.__version__)
```

See the [Installation Guide](getting-started/installation.md) for more options.

## Documentation

- **[Getting Started](getting-started/quickstart.md)**: Installation and first steps
- **[User Guide](guide/file-operations.md)**: Tutorials and guides
- **[API Reference](api/index.md)**: Complete API documentation
- **[Examples](examples/workflows.md)**: Usage examples

## Requirements

- Python >= 3.9
- NumPy >= 1.20
- Most functions require sorted BED input files

## Acknowledgments

- [GRIT](https://github.com/manish59/grit) - The underlying Rust implementation
- [bedtools](https://bedtools.readthedocs.io/) by Aaron Quinlan - The original genomic interval toolkit that inspired this project

## License

MIT License. See [LICENSE](https://github.com/manish59/grit/blob/main/LICENSE) for details.
