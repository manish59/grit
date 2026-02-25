# coverage

## Description

Calculate coverage of A intervals by B intervals. Reports overlap count, bases covered, interval length, and fraction covered.

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
grit coverage -a example_a.bed -b example_b.bed --assume-sorted
```

## Output

```
chr1	100	200	gene1	100	+	1	60	100	0.6000000
chr1	150	250	gene2	200	-	2	60	100	0.6000000
chr1	400	500	gene3	300	+	1	30	100	0.3000000
chr2	100	300	gene4	400	+	1	100	200	0.5000000
chr2	500	700	gene5	500	-	1	50	200	0.2500000
```

Output columns: A interval fields, overlap count, bases covered, interval length, fraction covered.

## Options

| Flag | Description |
|------|-------------|
| `-a, --file-a` | Input BED file A (regions) |
| `-b, --file-b` | Input BED file B (reads/features) |
| `--hist` | Report a histogram of coverage |
| `-d, --per-base` | Report depth at each position |
| `--mean` | Report mean depth |
| `--assume-sorted` | Skip sorted validation |

## Histogram Mode

```bash
grit coverage -a example_a.bed -b example_b.bed --hist --assume-sorted
```

## Mean Depth

```bash
grit coverage -a example_a.bed -b example_b.bed --mean --assume-sorted
```

## Notes

- Uses streaming mode by default for memory efficiency
- Requires sorted input
- Memory usage: O(k) where k = max overlapping intervals
