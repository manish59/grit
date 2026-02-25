---
layout: default
title: Benchmarks - GRIT
---

# Benchmarks

GRIT is benchmarked against bedtools to ensure correctness and measure performance improvements. All benchmarks verify SHA256 hash parity with bedtools output.

## Summary Results

### 10M x 5M Intervals (Uniform Distribution)

| Command | bedtools | GRIT | Speedup | BT Memory | GRIT Memory | Memory Reduction |
|---------|----------|------|---------|-----------|-------------|------------------|
| window | 32.18s | 2.10s | **15.3x** | 1.5 GB | 11 MB | 137x less |
| merge | 3.68s | 0.34s | **10.8x** | 2.6 MB | 2.8 MB | ~same |
| coverage | 16.53s | 1.84s | **9.0x** | 1.4 GB | 11 MB | 134x less |
| subtract | 9.49s | 1.47s | **6.5x** | 208 MB | 11 MB | 19x less |
| sort | 4.52s | 0.84s | **5.4x** | N/A | 1.1 GB | N/A |
| closest | 9.70s | 1.95s | **5.0x** | 670 MB | 11 MB | 59x less |
| intersect | 6.77s | 1.54s | **4.4x** | 208 MB | 11 MB | 19x less |
| jaccard | 4.98s | 1.59s | **3.1x** | 3.4 GB | 2.8 MB | 1230x less |
| complement | 2.72s | 0.85s | **3.2x** | 2.7 MB | 2.7 MB | ~same |
| slop | 5.56s | 1.77s | **3.1x** | 2.2 MB | 2.6 MB | ~same |

### 10M x 5M Intervals (Clustered Distribution)

Clustered data simulates real-world hotspots (e.g., promoter regions, repeat elements):

| Command | bedtools | GRIT | Speedup | BT Memory | GRIT Memory |
|---------|----------|------|---------|-----------|-------------|
| window | 28.80s | 1.97s | **14.6x** | 1.4 GB | 12 MB |
| subtract | 14.72s | 1.22s | **12.1x** | 1.3 GB | 11 MB |
| coverage | 14.59s | 1.50s | **9.7x** | 1.4 GB | 11 MB |
| merge | 2.17s | 0.31s | **7.0x** | 55 MB | 2.8 MB |
| sort | 6.09s | 0.87s | **7.0x** | 2.7 GB | 1.2 GB |
| closest | 9.51s | 1.80s | **5.3x** | 583 MB | 12 MB |
| intersect | 6.27s | 1.44s | **4.4x** | 207 MB | 11 MB |
| complement | 2.44s | 0.85s | **2.9x** | 54 MB | 2.7 MB |
| slop | 5.49s | 1.78s | **3.1x** | 2.2 MB | 2.6 MB |
| jaccard | 4.51s | 1.95s | **2.3x** | 3.4 GB | 3.4 MB |

## Key Observations

### Memory Efficiency

GRIT's streaming algorithms use **O(k) memory** where k is the maximum number of overlapping intervals at any position (typically < 100). This enables:

- Processing 50GB+ files on machines with 4GB RAM
- Constant memory regardless of file size
- No memory spikes during processing

### Correctness Verification

All benchmarks verify correctness by comparing SHA256 hashes of sorted output:

```
PASS: bedtools hash == GRIT hash
FAIL: outputs differ (investigated and fixed)
```

For commands with non-deterministic output order (e.g., `window`), outputs are sorted before comparison.

## Methodology

### Test Data Generation

Synthetic BED files are generated using `grit generate`:

```bash
# Uniform distribution (random intervals across genome)
./benchmarks/bench.sh data 10M 5M uniform

# Clustered distribution (simulates hotspots)
./benchmarks/bench.sh data 10M 5M clustered
```

**Data characteristics:**
- Genome: hg38 (24 chromosomes)
- Interval sizes: 100-10,000 bp (uniform random)
- Positions: Uniform or clustered distribution
- Files are sorted in genome order (required by bedtools)

### Benchmark Execution

Each benchmark:

1. Clears filesystem cache between runs
2. Runs bedtools with `-sorted` flag where applicable
3. Runs GRIT with `--assume-sorted --streaming` flags
4. Captures wall-clock time and peak RSS memory
5. Verifies output correctness via SHA256 hash

### Hardware

Benchmarks were run on:
- **CPU**: Apple M-series / AMD Ryzen 9
- **RAM**: 16-32 GB
- **Storage**: NVMe SSD
- **OS**: macOS / Ubuntu Linux

### Software Versions

- **GRIT**: 0.1.0
- **bedtools**: 2.31.0
- **Rust**: 1.75+

## Reproducing Benchmarks

### Prerequisites

```bash
# Install bedtools
brew install bedtools  # macOS
# or
conda install -c bioconda bedtools  # conda

# Build GRIT
cargo build --release
```

### Running Benchmarks

```bash
# Check environment
./benchmarks/bench.sh env

# Generate test data (cached for reuse)
./benchmarks/bench.sh data 10M 5M

# Run all benchmarks
./benchmarks/bench.sh run 10M 5M all

# Run specific commands
./benchmarks/bench.sh run 10M 5M coverage intersect merge

# Run with clustered data
./benchmarks/bench.sh run 10M 5M clustered all
```

### Using Pre-computed Truth Data

For faster benchmarking, pre-compute bedtools results:

```bash
# Generate truth data (run once, reuse many times)
./benchmarks/bench.sh truth 10M 5M all

# Subsequent runs use cached bedtools output
./benchmarks/bench.sh run 10M 5M all
# Commands marked with * use truth data
```

### CSV Output

Benchmark results are saved as CSV for analysis:

```
benchmarks/results/bench_10M_5M_YYYYMMDD_HHMMSS.csv
```

## Real-World Datasets

GRIT can be benchmarked against real genomic data:

```bash
# List available datasets
./benchmarks/bench.sh real-list

# Download and prepare dbSNP data
./benchmarks/bench.sh real-download dbsnp --yes
./benchmarks/bench.sh real-prepare dbsnp

# Generate bedtools baseline
./benchmarks/bench.sh real-truth dbsnp all

# Run GRIT benchmark
./benchmarks/bench.sh real-run dbsnp all
```

Available datasets:
- **dbsnp**: dbSNP variant positions
- **encode_peaks**: ENCODE ChIP-seq peaks
- **gencode**: GENCODE gene annotations
- **sv**: Structural variant calls

## Performance Tips

For maximum performance:

1. **Pre-sort files**: Use `grit sort` once, then `--assume-sorted`
2. **Use streaming mode**: `--streaming` for large files
3. **Adjust threads**: `-t N` to control parallelism
4. **Pipeline operations**: Chain commands with pipes

```bash
# Optimal pipeline for large files
grit sort -i raw.bed | grit merge -i - --assume-sorted | grit intersect ...
```
