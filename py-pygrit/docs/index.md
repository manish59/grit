# pygrit

**High-performance genomic interval operations for Python.**

[![PyPI version](https://badge.fury.io/py/pygrit.svg)](https://badge.fury.io/py/pygrit)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

pygrit provides Python bindings for [GRIT](https://github.com/manish59/grit) (Genomic Range Interval Toolkit), offering streaming algorithms for large-scale genomic data analysis with **O(k) memory complexity**.

## Features

- **Streaming Operations**: Process files of any size with constant memory usage
- **High Performance**: Rust-powered core with zero-copy data handling
- **BED Format Support**: Full compatibility with standard BED files
- **NumPy Integration**: Seamless interoperability with NumPy arrays
- **Thread-Safe**: GIL-releasing operations for multi-threaded workloads

## Quick Example

```python
import pygrit

# Find overlapping intervals between two BED files
results = pygrit.intersect("genes.bed", "peaks.bed")

# Merge overlapping intervals
merged = pygrit.merge("regions.bed", distance=100)

# Subtract one set from another
remaining = pygrit.subtract("features.bed", "exclude.bed")

# Work with individual intervals
iv = pygrit.Interval("chr1", 1000, 2000)
print(f"Length: {len(iv)} bp")  # Length: 1000 bp
```

## Why pygrit?

| Feature | pygrit | pybedtools | pyranges |
|---------|--------|------------|----------|
| Memory complexity | O(k) | O(n) | O(n) |
| Large file support | Streaming | Load all | Load all |
| Implementation | Rust | C/Python | Python/Cython |
| Thread safety | Yes | Limited | Limited |

## Installation

```bash
pip install pygrit
```

See the [Installation Guide](getting-started/installation.md) for more options.

## Documentation

- **[Getting Started](getting-started/quickstart.md)**: Installation and first steps
- **[User Guide](guide/file-operations.md)**: In-depth tutorials and guides
- **[API Reference](api/index.md)**: Complete API documentation
- **[Examples](examples/workflows.md)**: Real-world usage examples

## License

pygrit is released under the MIT License. See [LICENSE](https://github.com/manish59/grit/blob/main/LICENSE) for details.
