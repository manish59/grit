---
layout: default
title: generate
parent: Commands
nav_order: 13
---

# grit generate

Generate synthetic BED datasets for benchmarking and testing.

## Usage

```bash
grit generate [OPTIONS]
```

## Options

| Option | Description |
|--------|-------------|
| `-o, --output <DIR>` | Output directory (default: `./grit_bench_data`) |
| `--sizes <SIZES>` | Comma-separated sizes to generate (default: `1M,5M,10M,25M,50M`) |
| `--seed <INT>` | Random seed for reproducibility (default: 42) |
| `--mode <MODE>` | Generation mode: `balanced`, `clustered`, `identical`, `skewed-a-gt-b`, `skewed-b-gt-a`, `all` |
| `--sorted <yes\|no\|auto>` | Sorting behavior (default: `auto`) |
| `--no-sort` | Alias for `--sorted no` |
| `--a <SIZE>` | Custom A file size |
| `--b <SIZE>` | Custom B file size |
| `--hotspot-frac <FLOAT>` | Genome fraction for hotspots (default: 0.05) |
| `--hotspot-weight <FLOAT>` | Interval fraction in hotspots (default: 0.80) |
| `--len-min <INT>` | Minimum interval length (default: 50) |
| `--len-max <INT>` | Maximum interval length (default: 1000) |
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

### Quick benchmark data

```bash
# Generate 100K intervals for quick testing
grit generate --sizes 100K --mode balanced --force
```

### Large-scale testing

```bash
# Generate multiple sizes for comprehensive benchmarking
grit generate --sizes 10M,50M,100M --mode all --seed 123
```

### Custom A/B sizes

```bash
# Generate asymmetric datasets
grit generate --a 1M --b 10M --mode balanced
```

### Clustered data

```bash
# Simulate real-world hotspots (e.g., promoter regions)
grit generate --mode clustered --hotspot-frac 0.1 --hotspot-weight 0.9
```

### Unsorted output

```bash
# Generate unsorted data for testing sort validation
grit generate --sizes 1M --no-sort
```

## Output Structure

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

## Notes

- All generated files use human genome chromosome sizes (hg38)
- Seed ensures reproducibility across runs
- Auto sort mode sorts when size <= 10M
- Generated data is suitable for benchmarking intersect, coverage, merge, etc.

[← Back to Commands](../index.html)
