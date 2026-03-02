---
layout: default
title: Rust Library Cookbook
nav_order: 8
---

# Rust Library Cookbook

This guide provides practical examples for using `grit_genomics` as a Rust library in your own projects.

## Getting Started

Add to your `Cargo.toml`:

```toml
[dependencies]
grit-genomics = "0.1"
```

## Basic Usage

### Reading BED Files

```rust
use grit_genomics::prelude::*;
use grit_genomics::bed::{read_intervals, parse_intervals, BedReader};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Read all intervals from a file
    let intervals = read_intervals("regions.bed")?;
    println!("Loaded {} intervals", intervals.len());

    // Parse from a string
    let content = "chr1\t100\t200\nchr1\t300\t400\n";
    let intervals = parse_intervals(content)?;

    // Streaming reader for large files
    let mut reader = BedReader::from_path("large.bed")?;
    while let Some(record) = reader.read_record()? {
        println!("{}\t{}\t{}", record.chrom, record.start, record.end);
    }

    Ok(())
}
```

### Working with Intervals

```rust
use grit_genomics::interval::{Interval, Strand};

fn main() {
    // Create an interval
    let iv = Interval::new("chr1", 100, 200);

    // With strand information
    let iv_stranded = Interval::with_strand("chr1", 100, 200, Strand::Forward);

    // Check overlap
    let other = Interval::new("chr1", 150, 250);
    if iv.overlaps(&other) {
        println!("Intervals overlap by {} bp", iv.overlap_length(&other));
    }

    // Calculate distance (returns None if overlapping)
    let distant = Interval::new("chr1", 300, 400);
    if let Some(dist) = iv.distance_to(&distant) {
        println!("Distance: {} bp", dist);
    }
}
```

## Command Examples

### Intersect: Find Overlapping Intervals

```rust
use grit_genomics::prelude::*;
use grit_genomics::bed::parse_intervals;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a = parse_intervals("chr1\t100\t200\nchr1\t300\t400\n")?;
    let b = parse_intervals("chr1\t150\t250\nchr1\t350\t450\n")?;

    // Build index for efficient lookups
    let b_index = IntervalIndex::from_intervals(b);

    // Find intersections
    let cmd = IntersectCommand::new();
    let results = cmd.find_intersections(&a, &b_index);

    for interval in results {
        println!("{}\t{}\t{}", interval.chrom, interval.start, interval.end);
    }

    Ok(())
}
```

### Intersect with Options

```rust
use grit_genomics::commands::IntersectCommand;

let mut cmd = IntersectCommand::new();
cmd.write_a = true;           // Include original A interval
cmd.write_b = true;           // Include original B interval
cmd.fraction = Some(0.5);     // Require 50% overlap
cmd.reciprocal = true;        // Require reciprocal overlap
cmd.unique = true;            // Report each A only once
```

### Merge: Combine Overlapping Intervals

```rust
use grit_genomics::prelude::*;
use grit_genomics::bed::parse_intervals;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let intervals = parse_intervals(
        "chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n"
    )?;

    let cmd = MergeCommand::new();
    let merged = cmd.merge(intervals);

    // Result: chr1:100-250 and chr1:300-400
    assert_eq!(merged.len(), 2);

    Ok(())
}
```

### Merge with Distance

```rust
use grit_genomics::commands::MergeCommand;

let mut cmd = MergeCommand::new();
cmd.distance = 100;  // Merge intervals within 100bp of each other
```

### Subtract: Remove Overlapping Regions

```rust
use grit_genomics::prelude::*;
use grit_genomics::bed::parse_intervals;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a = parse_intervals("chr1\t100\t300\n")?;
    let b = parse_intervals("chr1\t150\t200\n")?;

    let b_index = IntervalIndex::from_intervals(b);

    let cmd = SubtractCommand::new();
    let results = cmd.subtract(&a, &b_index);

    // Result: chr1:100-150 and chr1:200-300
    assert_eq!(results.len(), 2);

    Ok(())
}
```

### Closest: Find Nearest Intervals

```rust
use grit_genomics::prelude::*;
use grit_genomics::bed::parse_intervals;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a = parse_intervals("chr1\t100\t200\n")?;
    let b = parse_intervals("chr1\t300\t400\nchr1\t500\t600\n")?;

    let b_index = IntervalIndex::from_intervals(b);

    let cmd = ClosestCommand::new();
    let results = cmd.find_closest(&a, &b_index);

    // Each result is (a_interval, closest_b_interval, distance)
    for (a_iv, b_iv, dist) in results {
        println!("{} closest to {} (distance: {})",
            a_iv.start, b_iv.start, dist);
    }

    Ok(())
}
```

### Sort: Order Intervals

```rust
use grit_genomics::prelude::*;
use grit_genomics::bed::parse_intervals;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut intervals = parse_intervals(
        "chr2\t100\t200\nchr1\t300\t400\nchr1\t100\t200\n"
    )?;

    let cmd = SortCommand::new();
    cmd.sort(&mut intervals);

    // Now sorted: chr1:100-200, chr1:300-400, chr2:100-200

    Ok(())
}
```

## Streaming Operations

For large files that don't fit in memory, use streaming commands:

### Streaming Intersect

```rust
use grit_genomics::commands::StreamingIntersectCommand;
use std::io::{BufWriter, stdout};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = StreamingIntersectCommand::new();
    cmd.assume_sorted = true;  // Skip validation for speed

    let mut output = BufWriter::new(stdout());
    cmd.run("a.bed", "b.bed", &mut output)?;

    Ok(())
}
```

### Streaming Merge

```rust
use grit_genomics::commands::StreamingMergeCommand;
use std::io::{BufWriter, stdout};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = StreamingMergeCommand::new();
    cmd.distance = 100;

    let mut output = BufWriter::new(stdout());
    let stats = cmd.run("input.bed", &mut output)?;

    eprintln!("Merged {} intervals into {}",
        stats.input_count, stats.output_count);

    Ok(())
}
```

## Building an Index

For repeated queries against the same dataset:

```rust
use grit_genomics::prelude::*;
use grit_genomics::bed::read_intervals;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Load and index once
    let features = read_intervals("features.bed")?;
    let index = IntervalIndex::from_intervals(features);

    // Query multiple times
    let query = Interval::new("chr1", 1000, 2000);
    let overlapping: Vec<_> = index.find_overlapping(&query).collect();

    println!("Found {} overlapping features", overlapping.len());

    Ok(())
}
```

## Working with BED Records

`BedRecord` preserves all columns from the input:

```rust
use grit_genomics::bed::BedReader;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = BedReader::from_path("genes.bed")?;

    while let Some(record) = reader.read_record()? {
        // Core fields
        println!("{}:{}-{}", record.chrom, record.start, record.end);

        // Optional fields (if present)
        if let Some(name) = &record.name {
            println!("  Name: {}", name);
        }
        if let Some(score) = record.score {
            println!("  Score: {}", score);
        }
        if let Some(strand) = &record.strand {
            println!("  Strand: {:?}", strand);
        }

        // Additional fields preserved as raw string
        if let Some(rest) = &record.rest {
            println!("  Extra: {}", rest);
        }
    }

    Ok(())
}
```

## Error Handling

```rust
use grit_genomics::bed::{BedError, BedReader};

fn main() {
    match BedReader::from_path("nonexistent.bed") {
        Ok(reader) => { /* use reader */ }
        Err(BedError::Io(e)) => {
            eprintln!("File error: {}", e);
        }
        Err(BedError::Parse { line, message }) => {
            eprintln!("Parse error at line {}: {}", line, message);
        }
        Err(BedError::InvalidFormat(msg)) => {
            eprintln!("Invalid format: {}", msg);
        }
    }
}
```

## Parallel Processing

For maximum speed on multi-core systems:

```rust
use grit_genomics::prelude::*;
use grit_genomics::bed::read_intervals;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a = read_intervals("a.bed")?;
    let b = read_intervals("b.bed")?;

    let cmd = IntersectCommand::new();

    // Parallel processing using Rayon
    let results = cmd.find_intersections_parallel(a, b);

    println!("Found {} intersections", results.len());

    Ok(())
}
```

## Complete Example: ChIP-seq Peak Analysis

```rust
use grit_genomics::prelude::*;
use grit_genomics::bed::read_intervals;
use grit_genomics::commands::SlopCommand;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Load peaks and genes
    let peaks = read_intervals("peaks.bed")?;
    let genes = read_intervals("genes.bed")?;

    // Extend gene promoters (2kb upstream)
    let slop_cmd = SlopCommand::new();
    let promoters = slop_cmd.slop_intervals(&genes, 2000, 0, "genome.txt")?;

    // Build index for promoters
    let promoter_index = IntervalIndex::from_intervals(promoters);

    // Find peaks in promoters
    let intersect_cmd = IntersectCommand::new();
    let peaks_in_promoters = intersect_cmd.find_intersections(&peaks, &promoter_index);

    println!("Found {} peaks in promoter regions", peaks_in_promoters.len());

    // Merge nearby peaks
    let mut merge_cmd = MergeCommand::new();
    merge_cmd.distance = 500;  // Within 500bp
    let merged_peaks = merge_cmd.merge(peaks_in_promoters);

    println!("Merged into {} peak regions", merged_peaks.len());

    Ok(())
}
```

## Tips

1. **Use streaming for large files**: If your files exceed available RAM, use `StreamingIntersectCommand` etc.

2. **Pre-sort for streaming**: Streaming commands require sorted input. Sort once with `SortCommand` or `grit sort`.

3. **Build indices for repeated queries**: If querying the same dataset multiple times, build an `IntervalIndex` once.

4. **Release the GIL in Python**: The library automatically releases Python's GIL during computation, enabling true parallelism.

5. **Check stats**: Many streaming commands return statistics about processing.
