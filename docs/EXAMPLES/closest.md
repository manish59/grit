# closest

## Description

Find the closest interval in B for each interval in A. Reports both overlapping and nearest non-overlapping intervals.

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
grit closest -a example_a.bed -b example_b.bed
```

## Output

```
chr1	100	200	gene1	100	+	chr1	120	180	feat1	50	+
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
| `-d, --distance` | Report distance in output |
| `-t, --tie` | Tie handling: `all`, `first`, `last` |
| `--io` | Ignore overlapping intervals |
| `--iu` | Ignore upstream intervals |
| `--id` | Ignore downstream intervals |
| `-D, --max-distance` | Maximum distance to report |
| `--streaming` | Use streaming mode (O(k) memory) |
| `--assume-sorted` | Skip sorted validation |

## Report Distance

```bash
grit closest -a example_a.bed -b example_b.bed -d
```

## Ignore Overlaps

Find closest non-overlapping interval:

```bash
grit closest -a example_a.bed -b example_b.bed --io
```

## Streaming Mode

```bash
grit closest -a example_a.bed -b example_b.bed --streaming --assume-sorted
```

## Notes

- By default, reports all ties
- Use `-t first` or `-t last` to report only first/last tie
- Streaming mode requires sorted input
