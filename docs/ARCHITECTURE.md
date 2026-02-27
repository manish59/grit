---
layout: default
title: Architecture
nav_order: 6
---

# GRIT Architecture

This document describes the high-level architecture and design decisions in GRIT.

## Overview

GRIT (Genomic Range Interval Toolkit) is a high-performance Rust implementation of genomic interval operations. It prioritizes:

1. **Performance**: Zero-allocation hot paths, SIMD-accelerated parsing
2. **Memory efficiency**: O(k) streaming algorithms where k = max overlapping intervals
3. **Correctness**: SHA256 parity with bedtools for all commands
4. **Usability**: Drop-in replacement for bedtools CLI

## Module Structure

```
src/
├── main.rs              # CLI entry point and argument parsing
├── lib.rs               # Library exports
├── bed/
│   ├── mod.rs           # BED record types and parsing
│   ├── record.rs        # Interval record representation
│   └── parser.rs        # Zero-copy BED parsing
├── commands/
│   ├── mod.rs           # Command registry
│   ├── intersect.rs     # Intersection operations
│   ├── subtract.rs      # Subtraction operations
│   ├── merge.rs         # Merge operations
│   ├── sort.rs          # Sorting operations
│   ├── closest.rs       # Nearest neighbor search
│   ├── window.rs        # Window-based operations
│   ├── coverage.rs      # Coverage computation
│   ├── slop.rs          # Interval extension
│   ├── complement.rs    # Gap finding
│   ├── genomecov.rs     # Genome-wide coverage
│   ├── jaccard.rs       # Similarity coefficient
│   └── multiinter.rs    # Multi-file intersection
├── streaming/
│   ├── mod.rs           # Streaming infrastructure
│   ├── merge.rs         # Streaming merge
│   ├── intersect.rs     # Streaming intersection
│   └── sweep.rs         # Sweep-line algorithm
├── sort/
│   ├── mod.rs           # Sort implementations
│   └── radix.rs         # Radix sort for intervals
└── utils/
    ├── mod.rs           # Utility functions
    ├── genome.rs        # Genome file handling
    └── validation.rs    # Input validation
```

## Core Components

### BED Record Representation

```rust
pub struct BedRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub rest: Option<String>,  // Optional additional fields
}
```

**Design decisions:**
- Use `u64` for coordinates (supports chromosomes > 4GB)
- Store raw rest-of-line to avoid parsing unused fields
- Zero-copy parsing where possible

### Streaming Infrastructure

The streaming system enables O(k) memory processing:

```rust
pub trait StreamProcessor {
    fn process_interval(&mut self, record: &BedRecord) -> Result<Vec<BedRecord>>;
    fn flush(&mut self) -> Result<Vec<BedRecord>>;
}
```

See [STREAMING_MODEL.md](STREAMING_MODEL.md) for details.

### Sort Order Validation

GRIT validates input sort order to prevent silent failures:

```rust
pub struct SortValidator {
    last_chrom: Option<String>,
    last_start: u64,
    genome_order: Option<Vec<String>>,
}
```

See [VALIDATION.md](VALIDATION.md) for details.

## Command Architecture

Each command follows a consistent pattern:

1. **Argument parsing** via `clap` derive macros
2. **Input validation** (sort order, chromosome order)
3. **Core algorithm** (streaming or parallel)
4. **Output formatting** (bedtools-compatible)

### Example: Intersect

```rust
pub fn run(args: IntersectArgs) -> Result<()> {
    // 1. Open input files
    let a_reader = BedReader::open(&args.file_a)?;
    let b_reader = BedReader::open(&args.file_b)?;

    // 2. Validate sort order (unless --assume-sorted)
    if !args.assume_sorted {
        validate_sorted(&a_reader)?;
        validate_sorted(&b_reader)?;
    }

    // 3. Choose algorithm based on flags
    if args.streaming {
        streaming_intersect(a_reader, b_reader, &args)
    } else {
        parallel_intersect(a_reader, b_reader, &args)
    }
}
```

## Performance Optimizations

### Zero-Allocation Parsing

```rust
// Use memchr for fast field detection
fn parse_bed_line(line: &[u8]) -> Result<BedRecord> {
    let tab1 = memchr::memchr(b'\t', line)?;
    let tab2 = memchr::memchr(b'\t', &line[tab1+1..])?;
    // ...
}
```

### SIMD Integer Parsing

```rust
// Fast integer parsing without allocation
fn parse_u64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0u64, |acc, &b| acc * 10 + (b - b'0') as u64)
}
```

### Memory-Mapped I/O

```rust
// For large files, use mmap for efficient random access
let mmap = unsafe { Mmap::map(&file)? };
```

## Thread Model

GRIT uses Rayon for parallel processing:

- **Parallel mode**: Load intervals into memory, process in parallel
- **Streaming mode**: Single-threaded sweep-line for O(k) memory

```rust
// Parallel intersection
intervals_a.par_iter()
    .flat_map(|a| find_overlaps(a, &tree_b))
    .for_each(|result| output(result));
```

## Error Handling

GRIT uses `thiserror` for structured errors:

```rust
#[derive(Error, Debug)]
pub enum GritError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("file not sorted: position {0} after {1} on {2}")]
    UnsortedInput(u64, u64, String),

    #[error("invalid BED format at line {0}: {1}")]
    ParseError(usize, String),
}
```

## Testing Strategy

1. **Unit tests**: Each module has inline tests
2. **Integration tests**: End-to-end command tests
3. **Parity tests**: SHA256 comparison with bedtools
4. **Property tests**: Invariant verification

See the `tests/` directory for examples.
