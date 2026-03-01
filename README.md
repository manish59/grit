# GRIT: Genomic Range Interval Toolkit

A genomic interval toolkit written in Rust, inspired by [bedtools](https://bedtools.readthedocs.io/).

[![CI](https://github.com/manish59/grit/actions/workflows/ci.yml/badge.svg)](https://github.com/manish59/grit/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/grit-genomics.svg)](https://crates.io/crates/grit-genomics)
[![PyPI](https://img.shields.io/pypi/v/grit-genomics.svg)](https://pypi.org/project/grit-genomics/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Features

- **Streaming algorithms**: Process large files with constant memory usage
- **13 commands**: intersect, merge, subtract, closest, window, coverage, sort, slop, complement, genomecov, jaccard, multiinter, generate
- **Python bindings**: Available via `pip install grit-genomics`
- **Compatible output**: Produces output matching bedtools format

---

## Installation

### Cargo (Rust)

```bash
cargo install grit-genomics
```

### PyPI (Python)

```bash
pip install grit-genomics
```

```python
import pygrit
result = pygrit.intersect("a.bed", "b.bed")
```

### Pre-built Binaries

Download from [GitHub Releases](https://github.com/manish59/grit/releases) for:
- Linux x86_64
- Linux ARM64
- macOS x86_64
- macOS ARM64

### From Source

```bash
git clone https://github.com/manish59/grit
cd grit
cargo install --path .
```

Verify installation: `grit --version`

---

## Quick Start

```bash
# Find overlapping intervals
grit intersect -a regions.bed -b features.bed > overlaps.bed

# Merge overlapping intervals
grit merge -i intervals.bed > merged.bed

# Sort a BED file
grit sort -i unsorted.bed > sorted.bed

# Streaming mode for large files (requires sorted input)
grit intersect -a large_a.bed -b large_b.bed --streaming --assume-sorted > result.bed
```

---

## Commands

| Command | Description | bedtools equivalent |
|---------|-------------|---------------------|
| `intersect` | Find overlapping intervals | `bedtools intersect` |
| `subtract` | Remove overlapping regions | `bedtools subtract` |
| `merge` | Combine overlapping intervals | `bedtools merge` |
| `sort` | Sort BED files | `bedtools sort` |
| `closest` | Find nearest intervals | `bedtools closest` |
| `window` | Find intervals within a window | `bedtools window` |
| `coverage` | Calculate interval coverage | `bedtools coverage` |
| `slop` | Extend intervals | `bedtools slop` |
| `complement` | Find gaps between intervals | `bedtools complement` |
| `genomecov` | Genome-wide coverage | `bedtools genomecov` |
| `jaccard` | Similarity coefficient | `bedtools jaccard` |
| `multiinter` | Multi-file intersection | `bedtools multiinter` |
| `generate` | Generate synthetic datasets | - |

Run `grit <command> --help` for detailed usage.

---

## Streaming Mode

For large files, use `--streaming --assume-sorted` flags to process data with constant memory:

```bash
# Requires pre-sorted input (sort by chromosome, then start position)
grit intersect -a sorted_a.bed -b sorted_b.bed --streaming --assume-sorted

# Sort first if needed
grit sort -i unsorted.bed > sorted.bed
```

**Note**: Streaming mode requires sorted input files. Use `grit sort` or standard Unix sort (`sort -k1,1 -k2,2n`) to prepare files.

---

## Benchmarks

Performance comparisons were run on synthetic datasets. Results vary based on data distribution, overlap patterns, and hardware.

### Test Conditions
- **Hardware**: Apple M1 Pro, 16GB RAM
- **Data**: Synthetic BED files with uniform and clustered distributions
- **Methodology**: Median of 3 runs, hyperfine for timing

### Results (10M intervals, uniform distribution)

| Command | bedtools | GRIT | Notes |
|---------|----------|------|-------|
| intersect | 6.8s | 1.5s | With `--streaming --assume-sorted` |
| merge | 3.7s | 0.3s | With `--assume-sorted` |
| subtract | 9.5s | 1.5s | With `--streaming --assume-sorted` |
| closest | 9.7s | 2.0s | With `--streaming --assume-sorted` |

Memory usage in streaming mode stays constant regardless of file size (tested up to 100M intervals).

<details>
<summary>Benchmark commands</summary>

```bash
# GRIT
grit intersect -a A.bed -b B.bed --streaming --assume-sorted
grit merge -i A.bed --assume-sorted

# bedtools
bedtools intersect -a A.bed -b B.bed -sorted
bedtools merge -i A.bed
```

</details>

**Disclaimer**: These benchmarks are from our testing. Your results may vary depending on data characteristics and hardware. We encourage users to benchmark on their own datasets.

---

## Python API

Install: `pip install grit-genomics`

```python
import pygrit

# File-based operations
overlaps = pygrit.intersect("a.bed", "b.bed")
pygrit.merge("input.bed", output="merged.bed")

# In-memory operations
intervals = pygrit.read_bed("regions.bed")
merged = intervals.merge(distance=100)

# Create intervals programmatically
iv = pygrit.Interval("chr1", 100, 200)
print(len(iv))  # 100
```

See [py-pygrit/README.md](py-pygrit/README.md) for full Python API documentation.

---

## Documentation

- [Command Reference](https://manish59.github.io/grit/) - All commands with examples
- [Python API](py-pygrit/README.md) - Python bindings documentation

---

## Contributing

Contributions welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-feature`)
3. Run tests (`cargo test`)
4. Commit changes (`git commit -m 'feat: add new feature'`)
5. Push to branch (`git push origin feature/new-feature`)
6. Open a Pull Request

---

## Testing

```bash
# Run all tests
cargo test

# Run with release optimizations
cargo test --release
```

---

## License

MIT License - see [LICENSE](LICENSE) for details.

---

## Acknowledgments

- [bedtools](https://bedtools.readthedocs.io/) by Aaron Quinlan - The original genomic interval toolkit that inspired this project
- [Rust](https://www.rust-lang.org/) - Systems programming language
- [PyO3](https://pyo3.rs/) - Rust bindings for Python

## Citation

If you use GRIT in your research, please cite:

```
GRIT: Genomic Range Interval Toolkit
https://github.com/manish59/grit
```

For bedtools methodology, please also cite:
> Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842.
