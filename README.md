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

## Installation

```bash
# From crates.io (recommended)
cargo install grit-genomics

# From Homebrew (macOS/Linux)
brew install manish59/grit/grit

# From source
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
