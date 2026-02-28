# Multi-Tool Benchmark Guide

This guide explains how to run comprehensive benchmarks comparing GRIT against other genomic interval tools.

## Tools Compared

| Tool | Language | Type | Notes |
|------|----------|------|-------|
| **GRIT** | Rust | CLI | This project - streaming, memory-efficient |
| **bedtools** | C++ | CLI | The standard, widely used |
| **bedops** | C++ | CLI | Known for memory efficiency |
| **granges** | Rust | CLI | Another Rust implementation |
| **pyranges** | Python | Library | Popular Python API |
| **polars-bio** | Python/Rust | Library | Built on Polars, gaining traction |

## Quick Start

```bash
cd benchmarks

# Check what's installed
./multi_tool_bench.sh status

# Install missing tools
./multi_tool_bench.sh install bedops
./multi_tool_bench.sh install granges
./multi_tool_bench.sh install python  # installs pyranges + polars-bio

# Run a quick benchmark (100K intervals)
./multi_tool_bench.sh quick

# Run full benchmark suite
./multi_tool_bench.sh run 10M
```

## Installation

### Prerequisites

- **GRIT**: Build from source with `cargo build --release`
- **bedtools**: Required for data generation and baseline
  ```bash
  # macOS
  brew install bedtools

  # Ubuntu/Debian
  apt-get install bedtools
  ```

### Installing Comparison Tools

```bash
cd benchmarks

# bedops (C++)
./multi_tool_bench.sh install bedops

# granges (Rust) - requires cargo
./multi_tool_bench.sh install granges

# Python tools (creates venv in benchmarks/.venv)
./multi_tool_bench.sh install pyranges
./multi_tool_bench.sh install polars-bio

# Or install all Python tools at once
./multi_tool_bench.sh install python
```

### Manual Installation

If automatic installation fails:

```bash
# bedops
# macOS
brew install bedops
# Linux - download from https://github.com/bedops/bedops/releases

# granges
cargo install granges

# Python tools
python3 -m venv benchmarks/.venv
source benchmarks/.venv/bin/activate
pip install pyranges polars-bio pandas polars
```

## Running Benchmarks

### Basic Usage

```bash
# Benchmark with 1 million intervals
./multi_tool_bench.sh run 1M

# Benchmark with specific sizes (A × B)
./multi_tool_bench.sh run 10M 5M

# Benchmark specific operations
./multi_tool_bench.sh run 1M 500K intersect merge closest

# Quick sanity check
./multi_tool_bench.sh quick
```

### Full Benchmark Suite

```bash
# Run benchmarks at multiple scales
./multi_tool_bench.sh full

# This runs: 100K, 1M, 10M
```

### Available Operations

| Operation | Description | Two files? |
|-----------|-------------|:----------:|
| `intersect` | Find overlapping intervals | ✓ |
| `merge` | Combine overlapping intervals | |
| `subtract` | Remove overlapping regions | ✓ |
| `closest` | Find nearest intervals | ✓ |
| `coverage` | Calculate interval coverage | ✓ |
| `sort` | Sort intervals | |

### Size Notation

- `100K` = 100,000 intervals
- `1M` = 1,000,000 intervals
- `10M` = 10,000,000 intervals
- `100M` = 100,000,000 intervals

## Results

Results are saved to `benchmarks/results/multi_tool/`:

```
results/multi_tool/
├── 1M_500K_20240115_143022/
│   ├── results.csv           # Raw benchmark data
│   ├── results.md            # Markdown table
│   ├── grit_intersect.out    # Output files
│   ├── grit_intersect.time   # Timing data
│   └── ...
```

### Generating Visualizations

```bash
# Generate charts from results
python scripts/plot_multi_tool.py results/multi_tool/*/results.csv -o charts/

# Generate markdown table
python scripts/plot_multi_tool.py results.csv --markdown
```

### Example Output

```
=== Multi-Tool Benchmark: 1M × 500K ===

Available tools: grit bedtools bedops pyranges polars_bio

Operation: intersect
  grit:        0.45s  (12.3 MB)
  bedtools:    1.82s  (156.4 MB)
  bedops:      0.92s  (45.2 MB)
  pyranges:    3.21s  (890.5 MB)
  polars_bio:  0.78s  (234.1 MB)

Operation: merge
  grit:        0.12s  (8.1 MB)
  bedtools:    0.89s  (234.5 MB)
  bedops:      0.34s  (23.4 MB)
  pyranges:    1.45s  (456.7 MB)
```

## Understanding Results

### Time

Wall-clock execution time in seconds. Lower is better.

### Memory

Peak resident set size (RSS) in megabytes. Lower is better.

GRIT's streaming architecture typically shows:
- **O(k) memory** where k = max overlapping intervals
- Significantly lower memory than in-memory tools
- Memory advantage increases with file size

### Speedup

`speedup = bedtools_time / tool_time`

- `> 1.0`: Faster than bedtools
- `= 1.0`: Same as bedtools
- `< 1.0`: Slower than bedtools

## Tool-Specific Notes

### bedops

- Requires pre-sorted input (use `.sorted.bed` files)
- Uses different output format for some operations
- `--intersect` produces ranges, not overlapping pairs

### granges

- Limited operation support (mainly filtering/overlap)
- Still under active development
- Good for intersection operations

### pyranges

- Pure Python with Cython optimizations
- Higher memory usage (loads all data into memory)
- Comprehensive API with many operations

### polars-bio

- Built on Polars (Rust-based DataFrame library)
- Good performance for medium-sized datasets
- Growing feature set

## Troubleshooting

### Tool not found

```bash
# Check installation status
./multi_tool_bench.sh status

# Re-install
./multi_tool_bench.sh install <tool>
```

### Python import errors

```bash
# Recreate virtual environment
rm -rf benchmarks/.venv
./multi_tool_bench.sh install python
```

### Memory issues with large datasets

```bash
# Use streaming mode with pre-sorted files
# GRIT's --assume-sorted flag enables O(k) memory
grit intersect -a large.sorted.bed -b other.sorted.bed --assume-sorted
```

### Verifying correctness

The benchmark compares output line counts. For exact verification:

```bash
# Use the original bench.sh for SHA256 verification
./bench.sh verify 1M 500K
```

## Adding New Tools

To add a new comparison tool:

1. Add detection function in `scripts/install_tools.sh`
2. Add runner function in `multi_tool_bench.sh`
3. Update operation mappings
4. Add tool color in `scripts/plot_multi_tool.py`

## CI Integration

For automated benchmarks in CI:

```yaml
- name: Run multi-tool benchmarks
  run: |
    cd benchmarks
    ./multi_tool_bench.sh quick

- name: Upload results
  uses: actions/upload-artifact@v4
  with:
    name: benchmark-results
    path: benchmarks/results/
```
