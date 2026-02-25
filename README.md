# GRIT: Genomic Range Interval Toolkit

A high-performance genomic interval toolkit written in Rust. Drop-in replacement for bedtools with **2.8-8.3x faster** performance.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Crates.io](https://img.shields.io/crates/v/grit-genomics.svg)](https://crates.io/crates/grit-genomics)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://manish59.github.io/grit/)

---

## Table of Contents

- [Why GRIT?](#why-grit)
- [Installation](#installation)
- [Documentation](https://manish59.github.io/grit/)
- [Quick Start](#quick-start)
- [Commands](#commands)
  - [intersect](#intersect) - Find overlapping intervals
  - [subtract](#subtract) - Remove overlapping regions
  - [merge](#merge) - Combine overlapping intervals
  - [sort](#sort) - Sort BED files
  - [closest](#closest) - Find nearest intervals
  - [window](#window) - Find intervals within a window
  - [coverage](#coverage) - Calculate interval coverage
  - [slop](#slop) - Extend intervals
  - [complement](#complement) - Find gaps between intervals
  - [genomecov](#genomecov) - Genome-wide coverage
  - [jaccard](#jaccard) - Similarity coefficient
  - [multiinter](#multiinter) - Multi-file intersection
- [Input Validation](#input-validation)
  - [Sort Order Validation](#sort-order-validation)
  - [Genome Order Validation](#genome-order-validation)
  - [stdin Validation](#stdin-validation)
- [Streaming Mode](#streaming-mode)
- [Performance](#performance)
- [Testing](#testing)
- [Contributing](#contributing)
- [License](#license)

---

## Why GRIT?

| Feature | bedtools | GRIT |
|---------|----------|------|
| Speed | Baseline | **2.8-8.3x faster** |
| Memory (streaming) | N/A | **O(k) constant** |
| Parallelization | Single-threaded | Multi-core |
| Large file support | Limited by RAM | Process 50GB+ on 4GB RAM |

GRIT is designed for:
- **High-throughput genomics** - Process millions of intervals efficiently
- **Memory-constrained environments** - Streaming mode uses minimal RAM
- **Drop-in replacement** - Same CLI syntax as bedtools
- **Reproducibility** - Deterministic output regardless of thread count

---

## Installation

### From crates.io (Recommended)

```bash
cargo install grit-genomics
```

### From Homebrew (macOS/Linux)

```bash
brew tap manish59/grit
brew install grit
```

### From Source

```bash
git clone https://github.com/manish59/grit
cd grit
cargo build --release
cargo install --path .
```

### Verify Installation

```bash
grit --version
grit --help
```

### Documentation

Full command documentation with examples: **[https://manish59.github.io/grit/](https://manish59.github.io/grit/)**

---

## Quick Start

```bash
# Find overlapping intervals between two BED files
grit intersect -a regions.bed -b features.bed > overlaps.bed

# Merge overlapping intervals
grit merge -i intervals.bed > merged.bed

# Sort a BED file
grit sort -i unsorted.bed > sorted.bed

# Use streaming mode for large files (minimal memory)
grit intersect -a large_a.bed -b large_b.bed --streaming > result.bed
```

---

## Global Options

All commands support these options:

| Option | Description |
|--------|-------------|
| `-t, --threads <N>` | Number of threads (default: all CPUs) |
| `--bedtools-compatible` | Normalize zero-length intervals to 1bp for bedtools parity |
| `-h, --help` | Show help for any command |
| `-V, --version` | Show version |

```bash
# Run with 4 threads
grit -t 4 intersect -a file1.bed -b file2.bed

# Enable bedtools-compatible mode for zero-length intervals
grit --bedtools-compatible intersect -a snps.bed -b features.bed

# Get help for a specific command
grit intersect --help
```

---

## Commands

---

### intersect

Find overlapping intervals between two BED files.

#### When to Use

- Identify genomic regions that overlap between datasets (e.g., peaks vs. promoters)
- Filter intervals based on overlap with a reference set
- Find regions with NO overlap (exclusion analysis)
- Count how many times each region is covered

#### Why Use GRIT

- **2.8x faster** than bedtools intersect
- **3.7x faster** with streaming mode
- **O(k) memory** in streaming mode (k = max concurrent overlaps)

#### How to Use

```
grit intersect -a <FILE_A> -b <FILE_B> [OPTIONS]
```

**Required:**
- `-a, --file-a <FILE>` - Query intervals (file A)
- `-b, --file-b <FILE>` - Reference intervals (file B)

**Output Modes:**

| Option | Output |
|--------|--------|
| *(default)* | Overlap region only |
| `--wa` | Original A entry |
| `--wb` | Overlap region + B entry |
| `--wa --wb` | Both A and B entries |
| `-c, --count` | A entry + overlap count |
| `-u, --unique` | A entry once if ANY overlap |
| `-v, --no-overlap` | A entries with NO overlap |

**Filtering:**

| Option | Description |
|--------|-------------|
| `-f, --fraction <FLOAT>` | Minimum overlap as fraction of A (0.0-1.0) |
| `-r, --reciprocal` | Require reciprocal fraction overlap |

**Performance & Validation:**

| Option | Description |
|--------|-------------|
| `--streaming` | O(k) memory mode (requires sorted input) |
| `--assume-sorted` | Skip sort validation (faster for pre-sorted files) |
| `--allow-unsorted` | Allow unsorted input (loads and re-sorts in memory) |
| `-g, --genome <FILE>` | Validate chromosome order against genome file |
| `--stats` | Print statistics to stderr |

#### Examples

```bash
# Basic: find overlap regions
grit intersect -a peaks.bed -b promoters.bed > overlaps.bed

# Get original entries from both files
grit intersect -a a.bed -b b.bed --wa --wb > both.bed

# Find peaks NOT in blacklist regions
grit intersect -a peaks.bed -b blacklist.bed -v > filtered_peaks.bed

# Require 50% overlap of query interval
grit intersect -a a.bed -b b.bed -f 0.5 > overlap_50pct.bed

# Require 50% reciprocal overlap (both directions)
grit intersect -a a.bed -b b.bed -f 0.5 -r > reciprocal.bed

# Count overlaps per interval
grit intersect -a genes.bed -b variants.bed -c > gene_variant_counts.bed

# Report each query interval once (if it has any overlap)
grit intersect -a a.bed -b b.bed -u > has_overlap.bed

# Large files with minimal memory
grit intersect -a huge_a.bed -b huge_b.bed --streaming > result.bed
```

---

### subtract

Remove portions of A that overlap with B.

#### When to Use

- Remove blacklist regions from your intervals
- Exclude known features from analysis regions
- Clean up interval sets by removing specific regions

#### Why Use GRIT

- **6.6x faster** than bedtools subtract
- Streaming mode for large files
- Precise interval arithmetic

#### How to Use

```
grit subtract -a <FILE_A> -b <FILE_B> [OPTIONS]
```

**Required:**
- `-a, --file-a <FILE>` - Intervals to modify
- `-b, --file-b <FILE>` - Intervals to remove

**Options:**

| Option | Description |
|--------|-------------|
| `-A, --remove-entire` | Remove entire A interval if ANY overlap |
| `-f, --fraction <FLOAT>` | Minimum overlap fraction required |
| `-r, --reciprocal` | Require reciprocal fraction |
| `--streaming` | O(k) memory mode (requires sorted input) |
| `--assume-sorted` | Skip sort validation (faster for pre-sorted files) |
| `--allow-unsorted` | Allow unsorted input (loads and re-sorts in memory) |
| `-g, --genome <FILE>` | Validate chromosome order against genome file |
| `--stats` | Print statistics to stderr |

#### Examples

```bash
# Remove blacklist regions (keeps non-overlapping portions)
grit subtract -a peaks.bed -b blacklist.bed > clean_peaks.bed

# Remove entire interval if ANY overlap with blacklist
grit subtract -a peaks.bed -b blacklist.bed -A > strict_clean.bed

# Only subtract if >50% overlap
grit subtract -a a.bed -b b.bed -f 0.5 > result.bed

# Large files with streaming
grit subtract -a large_a.bed -b large_b.bed --streaming > result.bed
```

---

### merge

Combine overlapping and adjacent intervals into single intervals.

#### When to Use

- Collapse redundant overlapping intervals
- Create non-overlapping interval sets
- Simplify interval data before downstream analysis
- Combine intervals within a certain distance

#### Why Use GRIT

- **7.3x faster** than bedtools merge
- **2 MB memory** regardless of file size
- Streaming by default (no `--streaming` flag needed)

#### How to Use

```
grit merge -i <INPUT> [OPTIONS]
```

**Required:**
- `-i, --input <FILE>` - Input BED file (use `-` for stdin)

**Options:**

| Option | Description |
|--------|-------------|
| `-d, --distance <INT>` | Merge intervals within this distance (default: 0) |
| `-s, --strand` | Only merge intervals on same strand |
| `-c, --count` | Report count of merged intervals |
| `--in-memory` | Load all records (for unsorted input) |
| `--assume-sorted` | Skip sort validation (faster for pre-sorted files) |
| `-g, --genome <FILE>` | Validate chromosome order against genome file |
| `--stats` | Print statistics to stderr |

#### Examples

```bash
# Basic merge (overlapping and adjacent)
grit merge -i intervals.bed > merged.bed

# Merge intervals within 100bp of each other
grit merge -i intervals.bed -d 100 > merged_100bp.bed

# Strand-specific merging
grit merge -i stranded.bed -s > merged_stranded.bed

# Count how many intervals were merged
grit merge -i intervals.bed -c > merged_counts.bed

# Read from stdin (piping)
cat intervals.bed | grit merge -i - > merged.bed

# Handle unsorted input
grit merge -i unsorted.bed --in-memory > merged.bed
```

---

### sort

Sort BED files by chromosome and position.

#### When to Use

- Prepare files for streaming operations
- Ensure consistent ordering for reproducibility
- Sort by interval size for analysis
- Use custom chromosome ordering (genome file)

#### Why Use GRIT

- **O(n) radix sort** vs O(n log n) comparison sort
- Memory-mapped I/O for large files
- Stable sort preserves input order for ties

#### How to Use

```
grit sort -i <INPUT> [OPTIONS]
```

**Required:**
- `-i, --input <FILE>` - Input BED file (use `-` for stdin)

**Options:**

| Option | Description |
|--------|-------------|
| `-g, --genome <FILE>` | Custom chromosome order from genome file |
| `--sizeA` | Sort by interval size (ascending) |
| `--sizeD` | Sort by interval size (descending) |
| `-r, --reverse` | Reverse final sort order |
| `--chrThenSizeA` | Sort by chromosome name only |
| `--stats` | Print statistics to stderr |

#### Examples

```bash
# Default sort (chromosome lexicographic, then start position)
grit sort -i unsorted.bed > sorted.bed

# Custom chromosome order from genome file
grit sort -i input.bed -g genome.txt > sorted.bed

# Sort by interval size (smallest first)
grit sort -i input.bed --sizeA > by_size.bed

# Sort by interval size (largest first)
grit sort -i input.bed --sizeD > by_size_desc.bed

# Reverse sort order
grit sort -i input.bed -r > reversed.bed

# Read from stdin
cat input.bed | grit sort -i - > sorted.bed
```

**Genome File Format:**
```
chr1    248956422
chr2    242193529
chr3    198295559
```

---

### closest

Find the nearest interval in B for each interval in A.

#### When to Use

- Find nearest gene for each variant
- Identify closest regulatory element to each peak
- Distance-to-feature analysis
- Nearest neighbor genomic analysis

#### Why Use GRIT

- Efficient O(n log m) binary search algorithm
- Flexible tie-breaking options
- Direction-aware searching (upstream/downstream)

#### How to Use

```
grit closest -a <FILE_A> -b <FILE_B> [OPTIONS]
```

**Required:**
- `-a, --file-a <FILE>` - Query intervals
- `-b, --file-b <FILE>` - Reference intervals to search

**Options:**

| Option | Description |
|--------|-------------|
| `-d, --distance` | Report distance in output |
| `-t, --tie <MODE>` | Handle ties: `all`, `first`, `last` |
| `--io` | Ignore overlapping intervals |
| `--iu` | Ignore upstream intervals |
| `--id` | Ignore downstream intervals |
| `-D, --max-distance <INT>` | Maximum search distance |
| `--streaming` | O(k) memory mode (requires sorted input) |
| `--assume-sorted` | Skip sort validation (faster for pre-sorted files) |
| `--allow-unsorted` | Allow unsorted input (loads and re-sorts in memory) |
| `-g, --genome <FILE>` | Validate chromosome order against genome file |

#### Examples

```bash
# Find closest gene for each variant
grit closest -a variants.bed -b genes.bed > nearest_genes.bed

# Include distance in output
grit closest -a a.bed -b b.bed -d > closest_with_distance.bed

# Only report first tie
grit closest -a a.bed -b b.bed -t first > closest_first.bed

# Find nearest non-overlapping interval
grit closest -a a.bed -b b.bed --io > nearest_nonoverlap.bed

# Only look downstream
grit closest -a a.bed -b b.bed --iu > downstream_only.bed

# Only look upstream
grit closest -a a.bed -b b.bed --id > upstream_only.bed

# Limit search to 10kb
grit closest -a a.bed -b b.bed -D 10000 > closest_10kb.bed
```

---

### window

Find intervals in B within a window around intervals in A.

#### When to Use

- Find features within a distance of query regions
- Identify nearby regulatory elements
- Proximity-based feature association
- Asymmetric distance searches (different upstream/downstream)

#### Why Use GRIT

- Flexible symmetric and asymmetric windows
- Count or report modes
- Efficient interval tree queries

#### How to Use

```
grit window -a <FILE_A> -b <FILE_B> [OPTIONS]
```

**Required:**
- `-a, --file-a <FILE>` - Query intervals
- `-b, --file-b <FILE>` - Reference intervals

**Options:**

| Option | Description |
|--------|-------------|
| `-w, --window <INT>` | Window size both sides (default: 1000) |
| `-l, --left <INT>` | Left/upstream window size |
| `-r, --right <INT>` | Right/downstream window size |
| `-c, --count` | Report count of matches |
| `-v, --no-overlap` | Report A intervals with NO matches |
| `--assume-sorted` | Skip sort validation (faster for pre-sorted files) |
| `-g, --genome <FILE>` | Validate chromosome order against genome file |

#### Examples

```bash
# Find features within 1kb of query regions
grit window -a genes.bed -b enhancers.bed -w 1000 > nearby.bed

# Asymmetric window: 5kb upstream, 1kb downstream
grit window -a tss.bed -b enhancers.bed -l 5000 -r 1000 > nearby.bed

# Count features in window
grit window -a genes.bed -b variants.bed -w 5000 -c > counts.bed

# Find regions with no features nearby
grit window -a genes.bed -b enhancers.bed -v > isolated.bed
```

---

### coverage

Calculate coverage depth of B intervals over A intervals.

#### When to Use

- Count reads overlapping genomic regions
- Calculate what fraction of each region is covered
- Generate coverage statistics for intervals
- Quality control of sequencing data

#### Why Use GRIT

- **8.3x faster** than bedtools coverage
- Memory-efficient O(B) algorithm
- Multiple output formats (counts, histogram, per-base)

#### How to Use

```
grit coverage -a <FILE_A> -b <FILE_B> [OPTIONS]
```

**Required:**
- `-a, --file-a <FILE>` - Target regions
- `-b, --file-b <FILE>` - Features to count (reads, etc.)

**Options:**

| Option | Description |
|--------|-------------|
| `--hist` | Report histogram of coverage depths |
| `-d, --per-base` | Report depth at each position |
| `--mean` | Report mean depth per region |
| `--assume-sorted` | Skip sort validation (faster for pre-sorted files) |
| `-g, --genome <FILE>` | Validate chromosome order against genome file |

**Output Format (default):**
```
chrom  start  end  name  score  strand  count  bases_covered  length  fraction
```

#### Examples

```bash
# Basic coverage (count, covered bases, length, fraction)
grit coverage -a regions.bed -b reads.bed > coverage.bed

# Mean depth per region
grit coverage -a regions.bed -b reads.bed --mean > mean_depth.bed

# Per-base depth
grit coverage -a regions.bed -b reads.bed -d > per_base.bed

# Histogram of coverage depths
grit coverage -a regions.bed -b reads.bed --hist > histogram.txt

# Streaming mode for large files
grit coverage -a regions.bed -b reads.bed --streaming > coverage.bed
```

---

### slop

Extend intervals by a specified number of bases.

#### When to Use

- Expand peaks to include flanking regions
- Create promoter regions from TSS coordinates
- Add padding around features
- Strand-aware extension (upstream/downstream)

#### Why Use GRIT

- Respects chromosome boundaries
- Strand-aware extension
- Percentage-based extension option

#### How to Use

```
grit slop -i <INPUT> -g <GENOME> [OPTIONS]
```

**Required:**
- `-i, --input <FILE>` - Input BED file
- `-g, --genome <FILE>` - Chromosome sizes file

**Options:**

| Option | Description |
|--------|-------------|
| `-b, --both <INT>` | Extend both sides by N bases |
| `-l, --left <INT>` | Extend left/upstream |
| `-r, --right <INT>` | Extend right/downstream |
| `-s, --strand` | Use strand for upstream/downstream |
| `--pct` | Values are fractions of interval size |

#### Examples

```bash
# Extend 100bp on both sides
grit slop -i peaks.bed -g genome.txt -b 100 > extended.bed

# Create 500bp upstream + 100bp downstream regions
grit slop -i tss.bed -g genome.txt -l 500 -r 100 > promoters.bed

# Strand-aware extension (upstream/downstream relative to strand)
grit slop -i genes.bed -g genome.txt -l 1000 -r 0 -s > upstream_1kb.bed

# Extend by 10% of interval size on each side
grit slop -i peaks.bed -g genome.txt -b 0.1 --pct > extended_10pct.bed
```

**Genome File Format:**
```
chr1    248956422
chr2    242193529
chr3    198295559
```

---

### complement

Find genomic regions NOT covered by input intervals.

#### When to Use

- Find gaps between features
- Identify intergenic regions
- Create inverse of an interval set
- Find uncovered portions of chromosomes

#### Why Use GRIT

- O(n) single-pass streaming algorithm
- Memory efficient
- Simple, focused operation

#### How to Use

```
grit complement -i <INPUT> -g <GENOME>
```

**Required:**
- `-i, --input <FILE>` - Input BED file
- `-g, --genome <FILE>` - Chromosome sizes file

#### Examples

```bash
# Find gaps between intervals
grit complement -i genes.bed -g genome.txt > intergenic.bed

# Find uncovered regions
grit complement -i covered.bed -g genome.txt > gaps.bed
```

---

### genomecov

Compute genome-wide coverage statistics.

#### When to Use

- Generate coverage tracks for visualization
- Compute depth distribution across genome
- Create BedGraph files for genome browsers
- Normalize coverage (scaling)

#### Why Use GRIT

- Multiple output formats (histogram, BedGraph)
- Coverage scaling for normalization
- Efficient whole-genome processing

#### How to Use

```
grit genomecov -i <INPUT> -g <GENOME> [OPTIONS]
```

**Required:**
- `-i, --input <FILE>` - Input BED file
- `-g, --genome <FILE>` - Chromosome sizes file

**Options:**

| Option | Description |
|--------|-------------|
| `-d, --per-base` | Report depth at each position (1-based) |
| `--bg` | BedGraph format (non-zero regions only) |
| `--bga` | BedGraph format (including zero coverage) |
| `--scale <FLOAT>` | Scale depth by factor (default: 1.0) |

#### Examples

```bash
# Default histogram output
grit genomecov -i reads.bed -g genome.txt > histogram.txt

# BedGraph for visualization (non-zero only)
grit genomecov -i reads.bed -g genome.txt --bg > coverage.bedgraph

# BedGraph including zero coverage regions
grit genomecov -i reads.bed -g genome.txt --bga > coverage_all.bedgraph

# Per-base depth (large output)
grit genomecov -i reads.bed -g genome.txt -d > per_base.txt

# Scale coverage (e.g., RPM normalization)
grit genomecov -i reads.bed -g genome.txt --bg --scale 0.5 > scaled.bedgraph
```

**Histogram Output Format:**
```
chrom  depth  bases_at_depth  chrom_size  fraction
```

---

### jaccard

Calculate Jaccard similarity coefficient between two interval sets.

#### When to Use

- Compare similarity of two interval sets
- Measure overlap between experiments
- Quality control: compare replicates
- Quantify agreement between methods

#### Why Use GRIT

- O(n + m) efficient sweep-line algorithm
- Single-pass computation
- Standard Jaccard metric

#### How to Use

```
grit jaccard -a <FILE_A> -b <FILE_B>
```

**Required:**
- `-a, --file-a <FILE>` - First BED file
- `-b, --file-b <FILE>` - Second BED file

**Output Format:**
```
intersection    union    jaccard    n_intersections
15000           45000    0.333333   150
```

#### Examples

```bash
# Compare two peak sets
grit jaccard -a peaks_rep1.bed -b peaks_rep2.bed

# Compare methods
grit jaccard -a method1.bed -b method2.bed
```

---

### multiinter

Identify intervals and which files contain them across multiple BED files.

#### When to Use

- Find common intervals across multiple samples
- Identify sample-specific intervals
- Multi-way intersection analysis
- Consensus peak calling

#### Why Use GRIT

- Handles arbitrary number of files
- Reports which files contain each interval
- Cluster mode for strict consensus

#### How to Use

```
grit multiinter -i <FILE1> <FILE2> [FILE3...] [OPTIONS]
```

**Required:**
- `-i, --input <FILES>` - Two or more input BED files

**Options:**

| Option | Description |
|--------|-------------|
| `--cluster` | Only output intervals in ALL files |

#### Examples

```bash
# Find intervals across 3 files (reports which files contain each)
grit multiinter -i rep1.bed rep2.bed rep3.bed > multi.bed

# Find intervals present in ALL files (consensus)
grit multiinter -i rep1.bed rep2.bed rep3.bed --cluster > consensus.bed
```

---

## Input Validation

GRIT validates input files to prevent silent failures from incorrectly sorted data. This section explains the validation behavior and how to control it.

### Sort Order Validation

By default, GRIT validates that input files are sorted before processing. Most commands require sorted input (by chromosome, then by start position).

**If files are unsorted, you'll see a helpful error:**

```
Error: File A is not sorted: position 100 at line 5 comes after 200 on chr1

Fix: Run 'grit sort -i a.bed > sorted_a.bed' first.
Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).
```

**How to sort files:**

```bash
# Sort with GRIT (recommended)
grit sort -i unsorted.bed > sorted.bed

# Or use standard Unix sort
sort -k1,1 -k2,2n unsorted.bed > sorted.bed
```

### Validation Flags

| Flag | Description | Memory Impact |
|------|-------------|---------------|
| `--assume-sorted` | Skip validation entirely | No change |
| `--allow-unsorted` | Load and re-sort in memory | O(n) |
| `-g, --genome <FILE>` | Validate genome chromosome order | No change |

#### `--assume-sorted`

Skip validation when you know files are pre-sorted:

```bash
# Skip validation for faster startup
grit intersect --streaming --assume-sorted -a sorted_a.bed -b sorted_b.bed

# Useful in pipelines where files are guaranteed sorted
grit merge -i - --assume-sorted < sorted.bed
```

**Warning:** Using `--assume-sorted` with unsorted files produces incorrect results silently.

#### `--allow-unsorted`

For non-streaming commands (`intersect`, `subtract`, `closest`), explicitly allow unsorted input:

```bash
# Load and re-sort in memory (uses O(n) memory)
grit intersect --allow-unsorted -a unsorted_a.bed -b unsorted_b.bed

# Without this flag, unsorted input fails with a clear error
```

This flag is not available for streaming commands, which require pre-sorted input.

### Genome Order Validation

Use `-g, --genome` to validate that chromosomes appear in a specific order (e.g., hg38, mm10):

```bash
# Validate chromosome order against genome file
grit intersect --streaming -a a.bed -b b.bed -g hg38.genome

# Merge with genome order validation
grit merge -i input.bed -g hg38.genome

# Sort files to match genome order
grit sort -i input.bed -g hg38.genome > sorted.bed
```

**Genome file format** (tab-separated: chromosome name and size):

```
chr1    248956422
chr2    242193529
chr3    198295559
chrX    156040895
chrY    57227415
```

**When `-g` is provided:**
- Chromosomes must appear in the genome file order
- Chromosomes not in the genome file cause an error
- Error messages suggest how to fix: `grit sort -i file.bed -g genome.txt`

**Without `-g`:** Any contiguous chromosome order is accepted (lexicographic, natural, etc.)

### stdin Validation

When reading from stdin, GRIT buffers the input to validate sort order:

```bash
# stdin is validated by default (buffers entire input)
cat sorted.bed | grit merge -i - > merged.bed

# Skip stdin validation with --assume-sorted (no buffering)
cat sorted.bed | grit merge -i - --assume-sorted > merged.bed
```

**Note:** stdin validation uses O(n) memory to buffer input. For large piped inputs where data is guaranteed sorted, use `--assume-sorted` to skip buffering.

### Validation Summary by Command

| Command | Requires Sorted | `--allow-unsorted` | `-g, --genome` |
|---------|-----------------|--------------------| ---------------|
| `intersect` | Yes (streaming) / Validates (default) | Yes | Yes |
| `subtract` | Yes (streaming) / Validates (default) | Yes | Yes |
| `closest` | Yes (streaming) / Validates (default) | Yes | Yes |
| `merge` | Yes | No (use `--in-memory`) | Yes |
| `window` | Yes | No | Yes |
| `coverage` | Yes | No | Yes |
| `sort` | No | N/A | Yes (for ordering) |
| `slop` | No | N/A | No |
| `complement` | Yes | No | No |

---

## Streaming Mode

For very large files, streaming mode processes data with constant O(k) memory, where k is the maximum number of overlapping intervals at any position (typically < 100).

### When to Use Streaming

- Files larger than available RAM
- Processing 50GB+ files on laptops
- Memory-constrained environments
- When files are already sorted

### Streaming Commands

Commands that support `--streaming` mode:

```bash
# Intersect
grit intersect -a a.bed -b b.bed --streaming > result.bed
grit intersect -a a.bed -b b.bed --streaming --assume-sorted > result.bed

# Subtract
grit subtract -a a.bed -b b.bed --streaming > result.bed
grit subtract -a a.bed -b b.bed --streaming --assume-sorted > result.bed

# Closest
grit closest -a a.bed -b b.bed --streaming > result.bed
grit closest -a a.bed -b b.bed --streaming --assume-sorted > result.bed

# Window (always uses streaming internally)
grit window -a a.bed -b b.bed > result.bed
grit window -a a.bed -b b.bed --assume-sorted > result.bed

# Coverage (always uses streaming internally)
grit coverage -a a.bed -b b.bed > result.bed
grit coverage -a a.bed -b b.bed --assume-sorted > result.bed

# Merge (streaming by default)
grit merge -i sorted.bed > result.bed
grit merge -i sorted.bed --assume-sorted > result.bed
```

### Memory Comparison

| Mode | Memory Usage | Best For |
|------|--------------|----------|
| Default (parallel) | O(n + m) | Maximum speed |
| Streaming | O(k) ≈ 2 MB | Large files, low RAM |

---

## Zero-Length Interval Semantics

GRIT uses strict half-open interval semantics by default, which differs from bedtools in handling zero-length intervals.

### What Are Zero-Length Intervals?

Zero-length intervals have `start == end`, such as:
```
chr1    100    100
```

These represent point positions (e.g., SNP locations from VCF-to-BED conversion) rather than regions.

### Default Behavior (Strict Mode)

In strict half-open semantics, a zero-length interval `[100, 100)` contains no bases:
- It does **not** overlap with itself
- It does **not** overlap with adjacent intervals like `[100, 101)`

This follows the mathematical definition of half-open intervals.

### Bedtools Behavior

Bedtools treats zero-length intervals as if they were 1bp intervals:
- `[100, 100)` overlaps with `[100, 101)`
- Self-intersection of zero-length intervals produces output

### Enabling Bedtools Compatibility

Use `--bedtools-compatible` to match bedtools behavior:

```bash
# Default: strict semantics (zero-length intervals don't overlap)
grit intersect -a snps.bed -b features.bed

# Bedtools-compatible: zero-length intervals normalized to 1bp
grit --bedtools-compatible intersect -a snps.bed -b features.bed
```

When enabled, zero-length intervals are normalized to 1bp during parsing:
```
chr1    100    100  →  chr1    100    101
```

### When to Use Each Mode

| Mode | Use Case |
|------|----------|
| **Strict (default)** | Mathematical correctness, new projects |
| **Bedtools-compatible** | Reproducing bedtools results, dbSNP data |

### Performance Impact

The `--bedtools-compatible` flag has **negligible performance impact** (<1%). Normalization occurs once during parsing, not in inner loops.

---

## Performance

### Benchmarks

Tested on 500K intervals per file (AMD Ryzen 9, 32GB RAM):

| Command | bedtools | GRIT | Speedup |
|---------|----------|------|---------|
| intersect | 0.67s | 0.24s | **2.8x** |
| intersect --streaming | - | 0.18s | **3.7x** |
| subtract | 1.46s | 0.22s | **6.6x** |
| merge | 0.29s | 0.04s | **7.3x** |
| coverage | 2.08s | 0.25s | **8.3x** |
| closest | 0.53s | 0.50s | **1.1x** |

### Performance Tips

1. **Use streaming for large files** - Constant memory, often faster
2. **Pre-sort your files** - Use `--assume-sorted` to skip validation
3. **Adjust thread count** - Default uses all CPUs, tune with `-t`
4. **Use merge first** - Reduce interval count before expensive operations

---

## Testing

### Quick Start

```bash
# Run all tests
cargo test

# Build release for testing
cargo build --release
```

### Unit Tests

```bash
# All unit tests (290+ tests)
cargo test

# Specific module tests
cargo test streaming       # Streaming infrastructure tests
cargo test intersect       # Intersect command tests
cargo test merge           # Merge command tests
cargo test sort            # Sort command tests
cargo test coverage        # Coverage command tests
cargo test closest         # Closest command tests

# With verbose output
cargo test -- --nocapture
```

### Integration Tests

```bash
# Fast sort integration tests (requires bedtools installed)
cargo test --release --test fast_sort_integration
```

### Test Sorted Input Validation

GRIT validates that input files are sorted for streaming operations:

```bash
# Create test files
cat > sorted.bed << 'EOF'
chr1	100	200
chr1	300	400
chr2	100	200
EOF

cat > unsorted.bed << 'EOF'
chr2	100	200
chr1	300	400
EOF

# This should succeed (sorted input)
grit intersect --streaming -a sorted.bed -b sorted.bed

# This should fail with error (unsorted input)
grit intersect --streaming -a unsorted.bed -b sorted.bed
# Error: File A is not sorted...

# Skip validation with --assume-sorted (faster for pre-sorted files)
grit intersect --streaming --assume-sorted -a sorted.bed -b sorted.bed
```

### Test Commands Individually

```bash
# Intersect
grit intersect -a a.bed -b b.bed > result.bed
grit intersect --streaming -a a.bed -b b.bed > result.bed
grit intersect --streaming --assume-sorted -a a.bed -b b.bed > result.bed

# Subtract
grit subtract --streaming -a a.bed -b b.bed > result.bed
grit subtract --streaming --assume-sorted -a a.bed -b b.bed > result.bed

# Merge
grit merge -i sorted.bed > merged.bed
grit merge --assume-sorted -i sorted.bed > merged.bed

# Closest
grit closest --streaming -a a.bed -b b.bed > result.bed
grit closest --streaming --assume-sorted -a a.bed -b b.bed > result.bed

# Window
grit window -a a.bed -b b.bed -w 1000 > result.bed
grit window --assume-sorted -a a.bed -b b.bed -w 1000 > result.bed

# Coverage
grit coverage -a a.bed -b b.bed > result.bed
grit coverage --assume-sorted -a a.bed -b b.bed > result.bed
```

### Verify Against bedtools

```bash
# Compare intersect output
diff <(bedtools intersect -a a.bed -b b.bed | sort) \
     <(grit intersect -a a.bed -b b.bed | sort)

# Compare sort output
diff <(bedtools sort -i input.bed) <(grit sort -i input.bed)

# Compare merge output
diff <(bedtools merge -i sorted.bed) <(grit merge -i sorted.bed)

# Compare subtract output
diff <(bedtools subtract -a a.bed -b b.bed | sort) \
     <(grit subtract --streaming -a a.bed -b b.bed | sort)

# SHA256 parity check (for large files)
sha256sum <(bedtools intersect -a a.bed -b b.bed | sort) \
          <(grit intersect -a a.bed -b b.bed | sort)
```

### Run Benchmarks

```bash
# Run all benchmarks
cargo bench

# Specific benchmarks
cargo bench intersect
cargo bench sort
cargo bench merge

# Run benchmark script (if available)
./benchmarks/bench.sh run 1M 500K coverage subtract closest merge intersect
./benchmarks/grit-only.sh 10M 5M
```

### Performance Testing

```bash
# Generate large test files
for i in $(seq 1 1000000); do
  echo -e "chr$((RANDOM % 22 + 1))\t$((RANDOM * 100))\t$((RANDOM * 100 + 1000))"
done > large_a.bed

# Sort the test file
grit sort -i large_a.bed > large_a_sorted.bed

# Time streaming vs parallel mode
time grit intersect -a large_a_sorted.bed -b large_a_sorted.bed --streaming > /dev/null
time grit intersect -a large_a_sorted.bed -b large_a_sorted.bed > /dev/null

# Memory usage (on Linux)
/usr/bin/time -v grit intersect --streaming -a large_a_sorted.bed -b large_a_sorted.bed > /dev/null
```

### Test Coverage

```bash
# Install cargo-tarpaulin for coverage
cargo install cargo-tarpaulin

# Run coverage report
cargo tarpaulin --out Html
```

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
- [Rayon](https://github.com/rayon-rs/rayon) for parallel processing
