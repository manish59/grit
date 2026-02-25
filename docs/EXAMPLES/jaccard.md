# jaccard

## Description

Calculate Jaccard similarity between two BED files. Measures the overlap between two interval sets.

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
grit jaccard -a example_a.bed -b example_b.bed
```

## Output

```
intersection	union	jaccard	n_intersections
270	680	0.397059	5
```

## Options

| Flag | Description |
|------|-------------|
| `-a, --file-a` | Input BED file A |
| `-b, --file-b` | Input BED file B |

## Output Columns

| Column | Description |
|--------|-------------|
| intersection | Total bases in intersection |
| union | Total bases in union |
| jaccard | Jaccard index (intersection / union) |
| n_intersections | Number of intersecting interval pairs |

## Notes

- Jaccard index ranges from 0 (no overlap) to 1 (identical)
- Useful for comparing interval set similarity
- Both files are merged internally before comparison
