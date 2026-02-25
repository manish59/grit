# intersect

## Description

Find overlapping intervals between two BED files. Supports both parallel in-memory mode and streaming mode.

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
cat example_b.bed
```
```
chr1	120	180	feat1	50	+
chr1	220	280	feat2	60	-
chr1	450	480	feat3	70	+
chr2	150	250	feat4	80	+
chr2	600	650	feat5	90	-
```

## Command

```bash
grit intersect -a example_a.bed -b example_b.bed
```

## Output

```
chr1	120	180	gene1	100	+
chr1	150	180	gene2	200	-
chr1	220	250	gene2	200	-
chr1	450	480	gene3	300	+
chr2	150	250	gene4	400	+
chr2	600	650	gene5	500	-
```

## Options

| Flag | Description |
|------|-------------|
| `-a, --file-a` | Input BED file A |
| `-b, --file-b` | Input BED file B |
| `--wa` | Write original A entry |
| `--wb` | Write original B entry |
| `-u, --unique` | Only report unique A intervals |
| `-v, --no-overlap` | Only report A intervals with NO overlap |
| `-f, --fraction` | Minimum overlap fraction for A |
| `-r, --reciprocal` | Require reciprocal fraction overlap |
| `-c, --count` | Report the number of overlaps |
| `--streaming` | Use streaming mode (constant memory) |
| `--assume-sorted` | Skip sorted validation |
| `--stats` | Print streaming statistics |

## Write Both Entries

```bash
grit intersect -a example_a.bed -b example_b.bed --wa --wb
```

Output:
```
chr1	100	200	gene1	100	+	chr1	120	180	feat1	50	+
chr1	150	250	gene2	200	-	chr1	120	180	feat1	50	+
chr1	150	250	gene2	200	-	chr1	220	280	feat2	60	-
chr1	400	500	gene3	300	+	chr1	450	480	feat3	70	+
chr2	100	300	gene4	400	+	chr2	150	250	feat4	80	+
chr2	500	700	gene5	500	-	chr2	600	650	feat5	90	-
```

## Count Overlaps

```bash
grit intersect -a example_a.bed -b example_b.bed -c
```

Output:
```
chr1	100	200	gene1	100	+	1
chr1	150	250	gene2	200	-	2
chr1	400	500	gene3	300	+	1
chr2	100	300	gene4	400	+	1
chr2	500	700	gene5	500	-	1
```

## Streaming Mode

For large files with sorted input:

```bash
grit intersect -a example_a.bed -b example_b.bed --streaming --assume-sorted
```

## Notes

- Streaming mode requires sorted input and uses O(k) memory where k = max overlapping intervals
- Parallel mode loads files into memory and uses interval trees
- Use `--assume-sorted` to skip validation for pre-sorted files
