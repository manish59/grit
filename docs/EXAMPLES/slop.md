# slop

## Description

Extend intervals by a given number of bases. Respects chromosome boundaries defined in a genome file.

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
grit slop -i example_a.bed -g genome.txt -b 50
```

## Output

```
chr1	50	250	gene1	100	+
chr1	100	300	gene2	200	-
chr1	350	550	gene3	300	+
chr2	50	350	gene4	400	+
chr2	450	750	gene5	500	-
```

## Options

| Flag | Description |
|------|-------------|
| `-i, --input` | Input BED file |
| `-g, --genome` | Genome file (chrom sizes) |
| `-b, --both` | Extend both sides by this many bases |
| `-l, --left` | Extend left/upstream |
| `-r, --right` | Extend right/downstream |
| `-s, --strand` | Use strand info (left=upstream, right=downstream) |
| `--pct` | Interpret values as fraction of interval size |

## Asymmetric Extension

```bash
grit slop -i example_a.bed -g genome.txt -l 25 -r 75
```

## Strand-Aware Extension

Extend upstream and downstream relative to strand:

```bash
grit slop -i example_a.bed -g genome.txt -l 100 -r 50 -s
```

## Percentage-Based Extension

Extend by 50% of interval length on each side:

```bash
grit slop -i example_a.bed -g genome.txt -b 0.5 --pct
```

## Notes

- Coordinates are clipped to chromosome boundaries
- Requires a genome file with chromosome sizes
- Strand-aware mode uses strand column (column 6) for direction
