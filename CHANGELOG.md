# Changelog

All notable changes to GRIT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Cross-platform SHA256 support in benchmark scripts (Ubuntu compatibility)
- Centralized streaming utilities module (`src/streaming/`)
- `--assume-sorted` flag for faster processing of pre-sorted inputs
- Unified benchmark script (`bench.sh`) for comparing GRIT vs bedtools

### Changed
- Refactored streaming commands to use shared active set management
- Improved code organization with shared parsing and validation utilities

### Fixed
- Jaccard `n_intersections` overcounting issue
- Streaming window now preserves full columns in output
- Streaming closest now correctly handles downstream ties
- Various clippy warnings and code style improvements

## [0.1.0] - 2024-12-01

### Added
- Initial release of GRIT (Genomic Range Interval Toolkit)
- Core commands: `intersect`, `subtract`, `merge`, `sort`, `closest`, `window`, `coverage`
- Additional commands: `slop`, `complement`, `genomecov`, `jaccard`, `multiinter`
- Streaming mode with O(k) memory for large file processing
- Parallel processing with Rayon for multi-core utilization
- Fast radix sort implementation
- Zero-allocation parsing for performance-critical paths

### Performance
- 2.8-8.3x faster than bedtools on standard benchmarks
- O(k) memory streaming algorithms (k = max overlapping intervals)
- Memory-mapped I/O for large files
- Parallel chromosome processing

### Compatibility
- SHA256 parity with bedtools output
- Support for BED3, BED4, BED6, and extended BED formats
- Cross-platform: Linux, macOS, Windows
