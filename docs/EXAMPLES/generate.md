# generate

## Description

Generate synthetic BED datasets for benchmarking. Creates reproducible test data with various distribution patterns.

## Command

```bash
grit generate --output ./benchmark_data --sizes 1M,5M,10M --seed 42
```

## Output

Creates directory structure with A.bed and B.bed files:

```
benchmark_data/
├── balanced/
│   ├── 1M/
│   │   ├── A.bed
│   │   └── B.bed
│   ├── 5M/
│   └── 10M/
├── clustered/
│   └── ...
├── identical/
│   └── ...
├── skewed-a-gt-b/
│   └── ...
└── skewed-b-gt-a/
    └── ...
```

## Options

| Flag | Description |
|------|-------------|
| `-o, --output` | Output directory (default: `./grit_bench_data`) |
| `--sizes` | Sizes to generate (default: `1M,5M,10M,25M,50M`) |
| `--seed` | Random seed for reproducibility (default: 42) |
| `--mode` | Generation mode (default: `all`) |
| `--sorted` | Sorting behavior: `yes`, `no`, `auto` (default: `auto`) |
| `--no-sort` | Alias for `--sorted no` |
| `--a` | Custom A file size |
| `--b` | Custom B file size |
| `--hotspot-frac` | Genome fraction for hotspots (default: 0.05) |
| `--hotspot-weight` | Interval fraction in hotspots (default: 0.80) |
| `--len-min` | Minimum interval length (default: 50) |
| `--len-max` | Maximum interval length (default: 1000) |
| `--force` | Overwrite existing files |

## Generation Modes

| Mode | Description |
|------|-------------|
| `balanced` | Equal-sized A and B files with uniform distribution |
| `skewed-a-gt-b` | A file 10x larger than B |
| `skewed-b-gt-a` | B file 10x larger than A |
| `identical` | A and B are identical |
| `clustered` | Intervals concentrated in hotspot regions |
| `all` | Generate all modes |

## Size Specifications

| Format | Example | Value |
|--------|---------|-------|
| Number | `1000` | 1,000 intervals |
| K suffix | `10K` | 10,000 intervals |
| M suffix | `5M` | 5,000,000 intervals |

## Examples

### Quick Benchmark Data

```bash
grit generate --sizes 100K --mode balanced --force
```

### Large-Scale Testing

```bash
grit generate --sizes 10M,50M,100M --mode all --seed 123
```

### Custom A/B Sizes

```bash
grit generate --a 1M --b 10M --mode balanced
```

### Clustered Data

```bash
grit generate --mode clustered --hotspot-frac 0.1 --hotspot-weight 0.9
```

### Unsorted Output

```bash
grit generate --sizes 1M --no-sort
```

## Notes

- All generated files use human genome chromosome sizes (hg38)
- Seed ensures reproducibility across runs
- Auto sort mode sorts when size <= 10M
- Generated data suitable for benchmarking intersect, coverage, etc.
