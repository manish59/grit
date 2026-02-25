# GRIT Benchmark Suite

Comprehensive benchmarking tools for comparing GRIT performance against bedtools.

## Quick Start

```bash
# Check environment (grit built, bedtools installed)
./bench.sh env

# Generate test data
./bench.sh data 10M 5M

# Run all benchmarks
./bench.sh run 10M 5M all
```

## Commands

### Environment Setup

```bash
./bench.sh env
```

Verifies:
- GRIT binary is built (`target/release/grit`)
- bedtools is installed
- Rust toolchain is available

### Data Generation

```bash
./bench.sh data <A_size> <B_size> [distribution]
```

**Size notation:** `1K`, `10K`, `100K`, `1M`, `10M`, `100M`

**Distributions:**
- `uniform` (default): Random intervals across genome
- `clustered`: Simulates real-world hotspots

**Examples:**
```bash
./bench.sh data 10M 5M           # 10M A intervals, 5M B intervals
./bench.sh data 10M 5M clustered # Clustered distribution
```

Data is cached in `benchmarks/data/<size>/<distribution>/` and reused.

### Running Benchmarks

```bash
./bench.sh run <A_size> <B_size> [distribution] <commands...>
```

**Examples:**
```bash
./bench.sh run 10M 5M all                    # All commands
./bench.sh run 10M 5M coverage intersect     # Specific commands
./bench.sh run 10M 5M clustered all          # Clustered data
```

**Available commands:**
- `coverage`, `intersect`, `merge`, `subtract`
- `closest`, `window`, `jaccard`
- `complement`, `slop`, `sort`
- `all` (runs all commands)

### Truth Data (Pre-computed Baselines)

For faster repeated benchmarking, pre-compute bedtools results:

```bash
# Generate truth data (slow, but only once)
./bench.sh truth 10M 5M all

# Subsequent runs use cached results (fast)
./bench.sh run 10M 5M all
# Commands marked with * use truth data

# List available truth data
./bench.sh truth-list

# Verify truth data integrity
./bench.sh truth-verify 10M 5M
```

### Real-World Datasets

Benchmark against actual genomic data:

```bash
# List available datasets
./bench.sh real-list

# Download dataset
./bench.sh real-download dbsnp --yes

# Prepare for benchmarking
./bench.sh real-prepare dbsnp

# Generate bedtools baseline
./bench.sh real-truth dbsnp all

# Run benchmark
./bench.sh real-run dbsnp all
```

**Available datasets:**
- `dbsnp`: dbSNP variant positions
- `encode_peaks`: ENCODE ChIP-seq peaks
- `gencode`: GENCODE gene annotations
- `sv`: Structural variant calls

### Cleanup

```bash
./bench.sh clean   # Remove all generated data and results
```

## Output

### Console Output

```
GRIT Benchmark Results
═══════════════════════════════════════════════════════════════════════════════
Data: 10M A intervals, 5M B intervals
Distribution: uniform

Command      bedtools       GRIT    Speedup       BT_RAM     GRIT_RAM   Status
───────────────────────────────────────────────────────────────────────────────
coverage*      16.53s      1.84s       8.98x       1.4 GB       10.7 MB   PASS
intersect*      6.77s      1.54s       4.40x     207.5 MB       11.2 MB   PASS
merge*          3.68s      0.34s      10.82x       2.6 MB        2.8 MB   PASS
═══════════════════════════════════════════════════════════════════════════════
```

Commands marked with `*` use pre-computed truth data.

### CSV Output

Results are saved to `benchmarks/results/`:

```
bench_10M_5M_20260225_130345.csv
```

**CSV columns:**
- `command`: The benchmark command
- `distribution`: uniform or clustered
- `tool`: bedtools or grit
- `mode`: default or streaming
- `time_s`: Execution time in seconds
- `rss_mb`: Peak memory in MB
- `speedup_vs_bedtools`: Speed improvement factor
- `correctness`: PASS or FAIL
- `source`: (truth) if using pre-computed data

## Directory Structure

```
benchmarks/
├── bench.sh              # Main benchmark script
├── data/                 # Generated test data (gitignored)
│   └── 10M_5M/
│       ├── uniform/
│       │   ├── A.sorted.bed
│       │   ├── B.sorted.bed
│       │   └── genome.txt
│       └── clustered/
├── truth/                # Pre-computed bedtools results
│   └── 10M_5M/
│       └── uniform/
│           ├── coverage.out
│           ├── coverage.time
│           ├── coverage.sha256
│           └── ...
├── results/              # Benchmark run results
│   └── bench_10M_5M_*.csv
├── real/                 # Real-world dataset configs
│   ├── dbsnp/
│   ├── encode_peaks/
│   └── ...
└── scripts/              # Helper scripts
    ├── common.sh
    ├── bench_all.sh
    └── verify_all.sh
```

## Correctness Verification

All benchmarks verify that GRIT output matches bedtools:

1. **Line count**: Same number of output lines
2. **SHA256 hash**: Identical content (after normalization)

**Normalization rules:**
- `coverage`: Compare first 6 columns only (bedtools adds extra fields)
- `window`: Sort output before comparison (order may vary)
- Others: Direct byte-for-byte comparison

## Adding New Benchmarks

To add a new command to the benchmark suite:

1. Add command to `ALL_COMMANDS` in `bench.sh`
2. Add bedtools/grit command mapping in `run_single_benchmark()`
3. Add hash normalization if needed (coverage, window style)
4. Generate truth data: `./bench.sh truth 10M 5M <new_command>`

## Troubleshooting

### "GRIT binary not found"
```bash
cargo build --release
```

### "bedtools not found"
```bash
brew install bedtools      # macOS
conda install -c bioconda bedtools  # conda
```

### "Data not found"
```bash
./bench.sh data 10M 5M     # Generate data first
```

### Benchmark shows FAIL
Check `benchmarks/results/<run>/` for `.out` files to diff outputs:
```bash
diff results/<run>/coverage_bedtools.out results/<run>/coverage_grit.out
```

## CI Integration

The benchmark suite can be run in CI:

```yaml
# .github/workflows/bench.yml
- name: Run benchmarks
  run: |
    ./benchmarks/bench.sh env
    ./benchmarks/bench.sh data 100K 50K
    ./benchmarks/bench.sh run 100K 50K all
```

Use smaller sizes (100K) for CI to keep runtime reasonable.
