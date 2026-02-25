# merge

## Description

Merge overlapping or adjacent intervals into single intervals. Uses streaming algorithm with O(1) memory.

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

## Command

```bash
grit merge -i example_a.bed --assume-sorted
```

## Output

```
chr1	100	250
chr1	400	500
chr2	100	300
chr2	500	700
```

## Options

| Flag | Description |
|------|-------------|
| `-i, --input` | Input BED file (use `-` for stdin) |
| `-d, --distance` | Maximum distance between intervals to merge (default: 0) |
| `-s, --strand` | Require strand to match for merging |
| `-c, --count` | Report count of merged intervals |
| `--in-memory` | Use in-memory mode (handles unsorted input) |
| `--assume-sorted` | Skip sorted validation |
| `--stats` | Print streaming statistics to stderr |

## Merge with Distance

Merge intervals within 100bp of each other:

```bash
grit merge -i example_a.bed -d 100 --assume-sorted
```

## Count Merged Intervals

```bash
grit merge -i example_a.bed -c --assume-sorted
```

Output:
```
chr1	100	250	2
chr1	400	500	1
chr2	100	300	1
chr2	500	700	1
```

## Strand-Specific Merge

```bash
grit merge -i example_a.bed -s --assume-sorted
```

## Notes

- Input must be sorted by chromosome and position
- Use `--assume-sorted` for pre-sorted input to skip validation
- Uses O(1) memory in streaming mode
