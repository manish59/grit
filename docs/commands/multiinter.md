---
layout: default
title: multiinter
parent: Commands
nav_order: 12
---

# grit multiinter

Identify common intervals across multiple BED files.

## Usage

```bash
grit multiinter [OPTIONS] -i <FILE1> <FILE2> [FILE3...]
```

## Options

| Option | Description |
|--------|-------------|
| `-i, --input <FILES>` | Input BED files (2 or more) |
| `--cluster` | Only output intervals found in all files |
| `--streaming` | Use streaming mode (O(k) memory) |
| `--assume-sorted` | Skip sorted validation |

## Examples

### Basic multi-intersection

```bash
# Find common regions across 3 files
grit multiinter -i file1.bed file2.bed file3.bed > common.bed

# Streaming mode for large files
grit multiinter -i *.bed --streaming --assume-sorted > common.bed
```

### Intervals in ALL files

```bash
# Only report intervals present in every file
grit multiinter -i rep1.bed rep2.bed rep3.bed --cluster > consensus.bed
```

### Compare replicates

```bash
# Find consensus peaks across replicates
grit multiinter -i rep1_peaks.bed rep2_peaks.bed rep3_peaks.bed --cluster > consensus_peaks.bed
```

### Compare conditions

```bash
# Find regions present in all samples
grit multiinter -i sample1.bed sample2.bed sample3.bed sample4.bed > multi.bed
```

## Output

**Default output:**
```
chr1    100    150    2    1,2
chr1    150    200    3    1,2,3
chr1    200    250    1    3
```

| Column | Description |
|--------|-------------|
| 1-3 | Chromosome, start, end |
| 4 | Number of files with this interval |
| 5 | List of file indices (1-based) |

**With --cluster:**
```
chr1    150    200    3    1,2,3
```
Only intervals present in ALL files are reported.

## Visual Example

```
File 1:  |--------|
File 2:      |--------|
File 3:          |--------|

Output:
         |--|              (in file 1 only)
             |--|          (in files 1,2)
                 |--|      (in files 1,2,3)  ← --cluster reports only this
                     |--|  (in files 2,3)
                         |-| (in file 3 only)
```

## Use Cases

### Consensus peaks
```bash
# Require peak in at least 2 of 3 replicates
grit multiinter -i rep1.bed rep2.bed rep3.bed | awk '$4 >= 2' > consensus.bed
```

### Universal binding sites
```bash
# Sites bound in ALL cell types
grit multiinter -i cellA.bed cellB.bed cellC.bed --cluster > universal.bed
```

### Sample overlap analysis
```bash
# Analyze overlap patterns
grit multiinter -i *.bed > overlap_matrix.bed
```

## Performance

```bash
# Fastest with streaming mode
grit multiinter -i *.bed --streaming --assume-sorted > result.bed
```

[← Back to Commands](../index.html)
