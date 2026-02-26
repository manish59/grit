# GRIT: Genomic Range Interval Toolkit

A high-performance genomic interval toolkit written in Rust. Drop-in replacement for bedtools with **3-15x faster** performance.

[![CI](https://github.com/manish59/grit/actions/workflows/ci.yml/badge.svg)](https://github.com/manish59/grit/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/manish59/grit/branch/main/graph/badge.svg)](https://codecov.io/gh/manish59/grit)
[![Crates.io](https://img.shields.io/crates/v/grit-genomics.svg)](https://crates.io/crates/grit-genomics)
[![docs.rs](https://docs.rs/grit-genomics/badge.svg)](https://docs.rs/grit-genomics)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/rust-1.85%2B-blue.svg)](https://www.rust-lang.org)

---

## Why GRIT?

| Feature | bedtools | GRIT |
|---------|----------|------|
| Speed | Baseline | **3-15x faster** |
| Memory (streaming) | N/A | **O(k) constant** |
| Large file support | Limited by RAM | Process 50GB+ on 4GB RAM |

---

## Benchmarks

Tested on 10M × 5M intervals ([full methodology](https://manish59.github.io/grit/benchmarks.html)):

### Uniform Distribution

| Command | bedtools | GRIT | Speedup | BT Memory | GRIT Memory | Reduction |
|---------|----------|------|---------|-----------|-------------|-----------|
| window | 32.18s | 2.10s | **15.3x** | 1.5 GB | 11 MB | 137x less |
| merge | 3.68s | 0.34s | **10.8x** | 2.6 MB | 2.8 MB | ~same |
| coverage | 16.53s | 1.84s | **9.0x** | 1.4 GB | 11 MB | 134x less |
| subtract | 9.49s | 1.47s | **6.5x** | 208 MB | 11 MB | 19x less |
| closest | 9.70s | 1.95s | **5.0x** | 670 MB | 11 MB | 59x less |
| intersect | 6.77s | 1.54s | **4.4x** | 208 MB | 11 MB | 19x less |
| jaccard | 4.98s | 1.59s | **3.1x** | 3.4 GB | 2.8 MB | 1230x less |

### Clustered Distribution (Real-world hotspots)

| Command | bedtools | GRIT | Speedup | BT Memory | GRIT Memory |
|---------|----------|------|---------|-----------|-------------|
| window | 28.80s | 1.97s | **14.6x** | 1.4 GB | 12 MB |
| subtract | 14.72s | 1.22s | **12.1x** | 1.3 GB | 11 MB |
| coverage | 14.59s | 1.50s | **9.7x** | 1.4 GB | 11 MB |
| merge | 2.17s | 0.31s | **7.0x** | 55 MB | 2.8 MB |
| closest | 9.51s | 1.80s | **5.3x** | 583 MB | 12 MB |
| intersect | 6.27s | 1.44s | **4.4x** | 207 MB | 11 MB |
| jaccard | 4.51s | 1.95s | **2.3x** | 3.4 GB | 3.4 MB |

---

## Installation

### Bioconda (Recommended for Python users)

```bash
conda install -c bioconda grit-genomics
```

### Homebrew (macOS/Linux)

```bash
brew install manish59/grit/grit
```

### Cargo (Rust users)

```bash
cargo install grit-genomics
```

### Pre-built Binaries

Download from [GitHub Releases](https://github.com/manish59/grit/releases) for Linux (x86_64, ARM64) and macOS (x86_64, ARM64).

### From Source

```bash
git clone https://github.com/manish59/grit && cd grit && cargo install --path .
```

Verify: `grit --version`

---

## Quick Start

```bash
# Find overlapping intervals
grit intersect -a regions.bed -b features.bed > overlaps.bed

# Merge overlapping intervals
grit merge -i intervals.bed > merged.bed

# Sort a BED file
grit sort -i unsorted.bed > sorted.bed

# Streaming mode for large files (minimal memory)
grit intersect -a large_a.bed -b large_b.bed --streaming --assume-sorted > result.bed
```

---

## Documentation

**Full documentation: [https://manish59.github.io/grit/](https://manish59.github.io/grit/)**

- [Command Reference](https://manish59.github.io/grit/) - All commands with examples
- [Migration from bedtools](https://manish59.github.io/grit/migration.html) - Drop-in replacement guide
- [Benchmarks](https://manish59.github.io/grit/benchmarks.html) - Performance methodology
- [Input Validation](https://manish59.github.io/grit/VALIDATION.html) - Sort order & genome validation

---

## Commands

| Command | Description |
|---------|-------------|
| `intersect` | Find overlapping intervals |
| `subtract` | Remove overlapping regions |
| `merge` | Combine overlapping intervals |
| `sort` | Sort BED files |
| `closest` | Find nearest intervals |
| `window` | Find intervals within a window |
| `coverage` | Calculate interval coverage |
| `slop` | Extend intervals |
| `complement` | Find gaps between intervals |
| `genomecov` | Genome-wide coverage |
| `jaccard` | Similarity coefficient |
| `multiinter` | Multi-file intersection |
| `generate` | Generate synthetic datasets |

Run `grit <command> --help` for usage details.

---

## Contributing

Contributions welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-feature`)
3. Commit changes (`git commit -m 'feat: add new feature'`)
4. Push to branch (`git push origin feature/new-feature`)
5. Open a Pull Request

---

## License

MIT License - see [LICENSE](LICENSE) for details.

---

## Acknowledgments

- [bedtools](https://bedtools.readthedocs.io/) by Aaron Quinlan - the inspiration for this project
