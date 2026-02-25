# GRIT Input Validation

This document describes GRIT's input validation system and how to control it.

## Overview

GRIT validates input files to prevent silent failures from incorrectly sorted or malformed data. This differs from bedtools, which often silently produces incorrect results with unsorted input.

## Validation Types

### 1. Sort Order Validation

GRIT validates that intervals are sorted by:
1. Chromosome (contiguous blocks)
2. Start position (ascending within chromosome)

```
Valid:
chr1    100    200
chr1    150    250    ✓ (150 >= 100)
chr1    300    400    ✓ (300 >= 150)
chr2    100    200    ✓ (new chromosome)

Invalid:
chr1    100    200
chr1    50     150    ✗ (50 < 100, unsorted)
```

### 2. Chromosome Order Validation

With `-g/--genome`, GRIT validates chromosome order matches the genome file:

```
genome.txt:
chr1    248956422
chr2    242193529
chr3    198295559

Valid order: chr1 → chr2 → chr3
Invalid: chr2 → chr1 (wrong order)
Invalid: chr4 (not in genome file)
```

### 3. Format Validation

GRIT validates BED format:
- At least 3 tab-separated fields
- Start and end are valid integers
- Start ≤ end (unless zero-length intervals allowed)

## Validation Flags

| Flag | Effect |
|------|--------|
| (default) | Validate sort order, fail on error |
| `--assume-sorted` | Skip validation entirely |
| `--allow-unsorted` | Load and re-sort in memory |
| `-g/--genome FILE` | Validate chromosome order |

### `--assume-sorted`

Skips all sort validation. Use when:
- You know files are pre-sorted
- Processing output from another GRIT command
- Maximum startup speed needed

```bash
# Piping from sort - safe to assume sorted
grit sort -i raw.bed | grit merge -i - --assume-sorted
```

**Warning**: Using `--assume-sorted` with unsorted input produces **incorrect results silently**.

### `--allow-unsorted`

Loads entire file into memory and sorts it:

```bash
grit intersect -a unsorted.bed -b reference.bed --allow-unsorted
```

**Trade-offs**:
- Uses O(n) memory (entire file)
- Adds sorting overhead
- Enables processing of unsorted input

### `-g/--genome`

Validates chromosome order against a genome file:

```bash
grit intersect -a A.bed -b B.bed -g hg38.genome
```

**Genome file format**:
```
chr1    248956422
chr2    242193529
chr3    198295559
...
```

## Error Messages

### Unsorted Input

```
Error: File A is not sorted: position 50 at line 5 comes after 100 on chr1

Fix: Run 'grit sort -i file.bed > sorted.bed' first.
Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).
```

### Wrong Chromosome Order

```
Error: Chromosome 'chr1' appears after 'chr2' in file A

This violates the expected genome order from hg38.genome.
Fix: Run 'grit sort -i file.bed -g hg38.genome > sorted.bed' first.
```

### Invalid BED Format

```
Error: Invalid BED format at line 7: expected integer for start position, got 'abc'
```

## Implementation Details

### Sort Validator

```rust
pub struct SortValidator {
    last_chrom: Option<String>,
    last_start: u64,
    line_number: usize,
    genome_order: Option<HashMap<String, usize>>,
}

impl SortValidator {
    pub fn validate(&mut self, record: &BedRecord) -> Result<()> {
        self.line_number += 1;

        // Check position within chromosome
        if let Some(ref last_chrom) = self.last_chrom {
            if record.chrom == *last_chrom {
                if record.start < self.last_start {
                    return Err(GritError::UnsortedInput {
                        position: record.start,
                        after: self.last_start,
                        chrom: record.chrom.clone(),
                        line: self.line_number,
                    });
                }
            }
        }

        // Check genome order if specified
        if let Some(ref genome) = self.genome_order {
            // ... genome order validation
        }

        self.last_chrom = Some(record.chrom.clone());
        self.last_start = record.start;
        Ok(())
    }
}
```

### stdin Handling

When reading from stdin, GRIT buffers input for validation:

```rust
fn validate_stdin(reader: impl BufRead) -> Result<Vec<BedRecord>> {
    let mut records = Vec::new();
    let mut validator = SortValidator::new();

    for line in reader.lines() {
        let record = parse_bed_line(&line?)?;
        validator.validate(&record)?;
        records.push(record);
    }

    Ok(records)
}
```

**Memory impact**: stdin validation uses O(n) memory. Use `--assume-sorted` to avoid buffering:

```bash
cat sorted.bed | grit merge -i - --assume-sorted
```

## Validation by Command

| Command | Default | `--assume-sorted` | `--allow-unsorted` | `-g/--genome` |
|---------|---------|-------------------|--------------------| --------------|
| intersect | Validates | Skips | Sorts in memory | Validates order |
| subtract | Validates | Skips | Sorts in memory | Validates order |
| closest | Validates | Skips | Sorts in memory | Validates order |
| merge | Validates | Skips | `--in-memory` | Validates order |
| window | Validates | Skips | N/A | Validates order |
| coverage | Validates | Skips | N/A | Validates order |
| sort | N/A | N/A | N/A | Orders by genome |
| slop | N/A | N/A | N/A | N/A |
| complement | Validates | Skips | N/A | N/A |

## Best Practices

### Pipeline Processing

```bash
# Sort once, process multiple times
grit sort -i raw.bed > sorted.bed
grit merge -i sorted.bed --assume-sorted > merged.bed
grit intersect -a merged.bed -b features.bed --assume-sorted > result.bed
```

### Mixed Pipelines

```bash
# GRIT sort output is guaranteed sorted
grit sort -i raw.bed | grit merge -i - --assume-sorted | grit intersect -a - -b ref.bed --assume-sorted
```

### Genome-Specific Ordering

```bash
# Ensure consistent chromosome order across files
grit sort -i file1.bed -g hg38.genome > file1_sorted.bed
grit sort -i file2.bed -g hg38.genome > file2_sorted.bed
grit intersect -a file1_sorted.bed -b file2_sorted.bed -g hg38.genome --assume-sorted
```

## Comparison with bedtools

| Scenario | bedtools | GRIT |
|----------|----------|------|
| Unsorted input | Silent wrong results | Error with fix suggestion |
| Wrong chrom order | Depends on command | Error with `-g` |
| Zero-length intervals | Treats as 1bp | Strict (or `--bedtools-compatible`) |
