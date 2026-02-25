# sort

## Description

Sort a BED file by chromosome and position. Uses radix sort with memory-mapped I/O for optimal performance.

## Example Input

```bash
cat unsorted.bed
```
```
chr2	100	200	feat1
chr1	300	400	feat2
chr1	100	200	feat3
chr3	50	150	feat4
chr1	250	350	feat5
```

## Command

```bash
grit sort -i unsorted.bed
```

## Output

```
chr1	100	200	feat3
chr1	250	350	feat5
chr1	300	400	feat2
chr2	100	200	feat1
chr3	50	150	feat4
```

## Options

| Flag | Description |
|------|-------------|
| `-i, --input` | Input BED file (use `-` for stdin) |
| `-g, --genome` | Genome file for chromosome ordering |
| `--sizeA` | Sort by interval size (ascending) |
| `--sizeD` | Sort by interval size (descending) |
| `-r, --reverse` | Reverse the sort order |
| `--chrThenSizeA` | Sort by chromosome name only |
| `--stats` | Print sorting statistics to stderr |

## Sorting by Genome Order

```bash
grit sort -i unsorted.bed -g genome.txt
```

## Reverse Sort

```bash
grit sort -i unsorted.bed --reverse
```

## Notes

- GRIT sort uses radix sort by default, providing 2-5x speedup over comparison-based sorts
- Memory-mapped I/O enables efficient processing of large files
- Output is deterministic for identical inputs
