# genomecov

## Description

Compute genome-wide coverage statistics. Reports coverage depth across the genome in various output formats.

## Example Input

```bash
cat example_a.bed
```
```
chr1	100	200	gene1	100	+
chr1	150	250	gene2	200	-
chr1	400	500	gene3	300	+
chr2	100	300	gene4	400	+
chr2	500	700	gene5	500	-
```

```bash
cat genome.txt
```
```
chr1	1000
chr2	1000
chr3	1000
```

## Command

```bash
grit genomecov -i example_a.bed -g genome.txt
```

## Output

```
chr1	0	750	1000	0.75
chr1	1	200	1000	0.2
chr1	2	50	1000	0.05
chr2	0	600	1000	0.6
chr2	1	400	1000	0.4
chr3	0	1000	1000	1
genome	0	2350	3000	0.783333
genome	1	600	3000	0.2
genome	2	50	3000	0.0166667
```

Output columns: chromosome, depth, bases at depth, chromosome size, fraction at depth.

## Options

| Flag | Description |
|------|-------------|
| `-i, --input` | Input BED file |
| `-g, --genome` | Genome file (chrom sizes) |
| `-d, --per-base` | Report depth at each position (1-based) |
| `--bg` | Report BedGraph format (non-zero only) |
| `--bga` | Report BedGraph format (including zero coverage) |
| `--scale` | Scale depth by factor (default: 1.0) |

## BedGraph Output

Non-zero coverage regions only:

```bash
grit genomecov -i example_a.bed -g genome.txt --bg
```

## BedGraph with Zero Coverage

Include zero-coverage regions:

```bash
grit genomecov -i example_a.bed -g genome.txt --bga
```

## Per-Base Depth

```bash
grit genomecov -i example_a.bed -g genome.txt -d
```

## Scaled Coverage

Apply a scaling factor to depth values:

```bash
grit genomecov -i example_a.bed -g genome.txt --scale 0.5
```

## Notes

- Default output is histogram format
- BedGraph formats are useful for genome browsers
- Per-base mode generates large output for whole genomes
