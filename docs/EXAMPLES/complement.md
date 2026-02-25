# complement

## Description

Return intervals NOT covered by the input BED file. Reports gaps between intervals and uncovered chromosome regions.

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
grit complement -i example_a.bed -g genome.txt --assume-sorted
```

## Output

```
chr1	0	100
chr1	250	400
chr1	500	1000
chr2	0	100
chr2	300	500
chr2	700	1000
chr3	0	1000
```

## Options

| Flag | Description |
|------|-------------|
| `-i, --input` | Input BED file |
| `-g, --genome` | Genome file (chrom sizes) |
| `--assume-sorted` | Assume input is sorted (enables O(1) memory streaming) |

## Notes

- Input must be sorted by chromosome and position
- Reports regions from genome file that have no coverage
- Chromosomes in genome file but not in input are fully reported
- Uses O(1) memory when `--assume-sorted` is specified
