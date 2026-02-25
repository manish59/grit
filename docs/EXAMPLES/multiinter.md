# multiinter

## Description

Identify common intervals across multiple BED files. Reports regions and which files contain them.

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
grit multiinter -i example_a.bed example_b.bed
```

## Output

```
chr1	100	120	1	1	1	0
chr1	120	150	2	1,2	1	1
chr1	150	180	2	1,2	1	1
chr1	180	200	1	1	1	0
chr1	200	220	1	1	1	0
chr1	220	250	2	1,2	1	1
chr1	250	280	1	2	0	1
chr1	400	450	1	1	1	0
chr1	450	480	2	1,2	1	1
chr1	480	500	1	1	1	0
chr2	100	150	1	1	1	0
chr2	150	250	2	1,2	1	1
chr2	250	300	1	1	1	0
chr2	500	600	1	1	1	0
chr2	600	650	2	1,2	1	1
chr2	650	700	1	1	1	0
```

Output columns: chrom, start, end, file count, file list, presence flags for each file.

## Options

| Flag | Description |
|------|-------------|
| `-i, --input` | Input BED files (multiple allowed) |
| `--cluster` | Only output intervals found in all files |

## Cluster Mode

Only report intervals present in ALL files:

```bash
grit multiinter -i example_a.bed example_b.bed --cluster
```

## Multiple Files

```bash
grit multiinter -i file1.bed file2.bed file3.bed file4.bed
```

## Notes

- Useful for finding consensus regions across samples
- File numbering starts at 1
- Cluster mode filters to intervals present in every input file
