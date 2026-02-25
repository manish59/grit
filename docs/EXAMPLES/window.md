# window

## Description

Find intervals in B that are within a window of A. Extends each A interval by the window size before checking for overlaps.

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
grit window -a example_a.bed -b example_b.bed -w 100 --assume-sorted
```

## Output

```
chr1	100	200	gene1	100	+	chr1	120	180	feat1	50	+
chr1	100	200	gene1	100	+	chr1	220	280	feat2	60	-
chr1	150	250	gene2	200	-	chr1	120	180	feat1	50	+
chr1	150	250	gene2	200	-	chr1	220	280	feat2	60	-
chr1	400	500	gene3	300	+	chr1	450	480	feat3	70	+
chr2	100	300	gene4	400	+	chr2	150	250	feat4	80	+
chr2	500	700	gene5	500	-	chr2	600	650	feat5	90	-
```

## Options

| Flag | Description |
|------|-------------|
| `-a, --file-a` | Input BED file A |
| `-b, --file-b` | Input BED file B |
| `-w, --window` | Window size (both sides, default: 1000) |
| `-l, --left` | Left window size |
| `-r, --right` | Right window size |
| `-c, --count` | Report number of overlaps |
| `-v, --no-overlap` | Only report A intervals with no matches |
| `--assume-sorted` | Skip sorted validation |

## Asymmetric Window

Different window sizes for left and right:

```bash
grit window -a example_a.bed -b example_b.bed -l 50 -r 200 --assume-sorted
```

## Count Matches

```bash
grit window -a example_a.bed -b example_b.bed -w 100 -c --assume-sorted
```

## Notes

- Requires sorted input
- Window extends in both directions unless `-l` and `-r` specified
- Uses streaming implementation for memory efficiency
