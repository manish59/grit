# GRIT Tutorial

GRIT (Genomic Range Interval Toolkit) is a high-performance tool for genomic interval operations. This tutorial covers installation, basic usage, and advanced features.

## 1. Installation

### From Source

```bash
git clone https://github.com/manish59/grit.git
cd grit
cargo build --release
```

The binary is located at `target/release/grit`.

### Verify Installation

```bash
./target/release/grit --version
./target/release/grit --help
```

## 2. Basic Usage

### Creating Test Data

Generate synthetic BED files for testing:

```bash
grit generate --output ./test_data --sizes 10K --mode balanced --force
```

This creates:
- `test_data/balanced/10K/A.bed` - 10,000 intervals
- `test_data/balanced/10K/B.bed` - 10,000 intervals

### Sorting

BED files must be sorted for streaming operations:

```bash
grit sort -i input.bed > sorted.bed
```

With genome-specified chromosome order:

```bash
grit sort -i input.bed -g genome.txt > sorted.bed
```

### Finding Overlaps

Basic intersection:

```bash
grit intersect -a regions.bed -b features.bed
```

Report both entries:

```bash
grit intersect -a regions.bed -b features.bed --wa --wb
```

Count overlaps per interval:

```bash
grit intersect -a regions.bed -b features.bed -c
```

### Merging Intervals

```bash
grit sort -i input.bed | grit merge --assume-sorted
```

Merge with distance tolerance:

```bash
grit merge -i sorted.bed -d 100 --assume-sorted
```

## 3. Streaming vs Non-Streaming

### Non-Streaming Mode (Default)

Loads data into memory, uses parallel processing:

```bash
grit intersect -a A.bed -b B.bed
```

**Best for:**
- Unsorted input
- Small to medium files
- Multi-threaded systems

### Streaming Mode

Processes data sequentially with constant memory:

```bash
grit intersect -a A.bed -b B.bed --streaming --assume-sorted
```

**Best for:**
- Large sorted files
- Memory-constrained systems
- Pipeline workflows

### Memory Comparison

| File Size | Non-Streaming | Streaming |
|-----------|---------------|-----------|
| 1M intervals | ~400 MB | ~1 MB |
| 10M intervals | ~4 GB | ~1 MB |
| 100M intervals | ~40 GB | ~1 MB |

## 4. Sorted vs Assume-Sorted

### With Validation

GRIT validates sort order by default:

```bash
grit merge -i input.bed
```

If input is unsorted, GRIT reports an error with guidance.

### Skip Validation

For pre-sorted files, skip validation for faster startup:

```bash
grit merge -i input.bed --assume-sorted
```

### Pipeline Usage

When piping from sort, always use `--assume-sorted`:

```bash
grit sort -i input.bed | grit merge --assume-sorted
```

## 5. Bedtools Compatibility Mode

### Zero-Length Intervals

GRIT uses strict half-open interval semantics:
- Interval `[100, 100)` has length 0
- Zero-length intervals don't overlap with themselves

Bedtools treats zero-length intervals as 1bp:
- Interval `[100, 100)` becomes `[100, 101)`

### Enabling Compatibility

```bash
grit --bedtools-compatible intersect -a snps.bed -b regions.bed
```

### When to Use

Use `--bedtools-compatible` when:
- Working with SNP files (single nucleotide positions)
- Comparing results with bedtools
- Processing data from tools that use bedtools conventions

## 6. Synthetic Data Generation

### Quick Benchmark

```bash
grit generate --sizes 100K --mode balanced
```

### Full Benchmark Suite

```bash
grit generate --sizes 1M,5M,10M --mode all --seed 42
```

### Generation Modes

| Mode | A:B Ratio | Distribution |
|------|-----------|--------------|
| `balanced` | 1:1 | Uniform |
| `skewed-a-gt-b` | 10:1 | Uniform |
| `skewed-b-gt-a` | 1:10 | Uniform |
| `identical` | 1:1 | Same intervals |
| `clustered` | 1:1 | Hotspot regions |

### Custom Sizes

```bash
grit generate --a 1M --b 10M --mode balanced
```

### Clustered Data

Simulate ChIP-seq-like peaks:

```bash
grit generate --mode clustered --hotspot-frac 0.05 --hotspot-weight 0.80
```

## 7. Real Dataset Benchmarking

### Prepare Data

Download real genomic data:

```bash
# dbSNP variants
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/BED/bed_chr_1.bed.gz

# ENCODE peaks
wget https://www.encodeproject.org/.../peaks.bed.gz

# Sort for streaming
zcat peaks.bed.gz | grit sort | gzip > peaks.sorted.bed.gz
```

### Benchmark Commands

```bash
# Time intersect operation
time grit intersect -a regions.bed -b peaks.bed > /dev/null

# Time streaming mode
time grit intersect -a regions.bed -b peaks.bed --streaming --assume-sorted > /dev/null
```

### Compare with Bedtools

```bash
# GRIT
time grit intersect -a A.bed -b B.bed > grit_out.bed

# bedtools
time bedtools intersect -a A.bed -b B.bed > bedtools_out.bed

# Verify identical output
diff <(sort grit_out.bed) <(sort bedtools_out.bed)
```

## 8. Performance Philosophy

### Design Principles

1. **Zero-allocation hot paths**: Critical loops avoid heap allocation
2. **SIMD-accelerated parsing**: Uses memchr for fast field detection
3. **Memory-mapped I/O**: Large files processed without full loading
4. **Radix sort**: O(n) sorting instead of O(n log n)
5. **Streaming algorithms**: O(k) memory where k = max overlaps

### Performance Tips

**Use streaming for large files:**
```bash
grit intersect -a large_a.bed -b large_b.bed --streaming --assume-sorted
```

**Pre-sort once, use many times:**
```bash
grit sort -i raw.bed > sorted.bed
grit intersect -a sorted.bed -b features.bed --assume-sorted
grit coverage -a sorted.bed -b reads.bed --assume-sorted
```

**Skip validation for trusted input:**
```bash
grit merge -i sorted.bed --assume-sorted
```

**Use parallel mode for small files:**
```bash
grit intersect -a small.bed -b small.bed  # Parallel by default
```

### Expected Performance

| Operation | 10M intervals | Memory |
|-----------|---------------|--------|
| sort | 2-5 seconds | O(n) |
| merge (streaming) | 1-3 seconds | O(1) |
| intersect (streaming) | 2-6 seconds | O(k) |
| coverage (streaming) | 2-6 seconds | O(k) |

Actual performance depends on:
- Hardware (CPU, memory speed, SSD vs HDD)
- Data characteristics (overlap density, chromosome distribution)
- Output size

## Quick Reference

### Common Workflows

**Sort and merge:**
```bash
grit sort -i input.bed | grit merge --assume-sorted
```

**Intersect with output:**
```bash
grit intersect -a regions.bed -b features.bed --wa --wb > overlaps.bed
```

**Find non-overlapping:**
```bash
grit intersect -a regions.bed -b features.bed -v > unique.bed
```

**Calculate coverage:**
```bash
grit coverage -a genes.bed -b reads.bed --assume-sorted
```

**Compare datasets:**
```bash
grit jaccard -a sample1.bed -b sample2.bed
```

### See Also

- [Command Reference](COMMANDS.md)
- [Example Documentation](EXAMPLES/)
