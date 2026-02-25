# subtract

## Description

Remove intervals in A that overlap with B. Can either trim overlapping portions or remove entire intervals.

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
grit subtract -a example_a.bed -b example_b.bed
```

## Output

```
chr1	100	120	gene1	100	+
chr1	180	200	gene1	100	+
chr1	180	220	gene2	200	-
chr1	400	450	gene3	300	+
chr1	480	500	gene3	300	+
chr2	100	150	gene4	400	+
chr2	250	300	gene4	400	+
chr2	500	600	gene5	500	-
chr2	650	700	gene5	500	-
```

## Options

| Flag | Description |
|------|-------------|
| `-a, --file-a` | Input BED file A |
| `-b, --file-b` | Input BED file B |
| `-A, --remove-entire` | Remove entire A feature if any overlap |
| `-f, --fraction` | Minimum overlap fraction required |
| `-r, --reciprocal` | Require reciprocal fraction overlap |
| `--streaming` | Use streaming mode (O(k) memory) |
| `--assume-sorted` | Skip sorted validation |
| `--stats` | Print streaming statistics |

## Remove Entire Features

Remove complete intervals that have any overlap:

```bash
grit subtract -a example_a.bed -b example_b.bed -A
```

## Streaming Mode

```bash
grit subtract -a example_a.bed -b example_b.bed --streaming --assume-sorted
```

## Notes

- Default behavior trims overlapping portions
- Use `-A` to remove entire intervals with any overlap
- Streaming mode uses O(k) memory where k = max overlapping intervals
