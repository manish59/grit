---
layout: default
title: Troubleshooting
nav_order: 10
---

# Troubleshooting Guide

This guide covers common errors and how to resolve them.

## Quick Reference

| Error | Cause | Fix |
|-------|-------|-----|
| "Input is not sorted" | Unsorted BED file | Run `grit sort -i file.bed` |
| "Chromosome not found in genome" | Missing chromosome in genome file | Update genome file or check chromosome names |
| "Expected at least 3 fields" | Malformed BED line | Check file format (tab-separated) |
| "start > end" | Invalid interval | Fix coordinates or check for parsing issues |
| High memory usage | Large file in parallel mode | Use `--streaming` mode |

---

## Sorted Input Errors

### Error: "Input is not sorted"

**Full message**:
```
Error: Input is not sorted: record at line 42 out of order

Fix: Run 'grit sort -i input.bed' first.
Or use '--assume-sorted' if you know the input is sorted.
```

**Cause**: Streaming operations require input sorted by chromosome, then by start position.

**Solutions**:

1. **Sort your file**:
   ```bash
   grit sort -i input.bed > sorted.bed
   grit intersect -a sorted.bed -b other.bed --streaming
   ```

2. **Use parallel mode** (handles unsorted input):
   ```bash
   # Remove --streaming flag
   grit intersect -a unsorted.bed -b other.bed > output.bed
   ```

3. **Skip validation** (only if you're sure input is sorted):
   ```bash
   grit intersect -a input.bed -b other.bed --streaming --assume-sorted
   ```

### Error: "Chromosome revisited"

**Full message**:
```
Error: Chromosome 'chr1' revisited at line 500 (previously seen, then other chromosomes appeared)
```

**Cause**: In streaming mode, all intervals for a chromosome must be contiguous. This error means:
```
chr1    100    200    ✓
chr1    300    400    ✓
chr2    100    200    ✓
chr1    500    600    ✗ ERROR: chr1 after chr2
```

**Solutions**:

1. **Sort your file**:
   ```bash
   grit sort -i input.bed > sorted.bed
   ```

2. **Use parallel mode**:
   ```bash
   grit intersect -a input.bed -b other.bed > output.bed
   ```

---

## Genome File Errors

### Error: "Chromosome not found in genome file"

**Full message**:
```
Error: Chromosome 'chrUn_gl000220' at line 42 not found in genome file
```

**Cause**: The BED file contains a chromosome not listed in the genome file.

**Solutions**:

1. **Add missing chromosome to genome file**:
   ```bash
   echo -e "chrUn_gl000220\t168386" >> genome.txt
   ```

2. **Filter out unknown chromosomes**:
   ```bash
   # Get chromosomes from genome file
   cut -f1 genome.txt > valid_chroms.txt
   grep -f valid_chroms.txt input.bed > filtered.bed
   ```

3. **Use a complete genome file**:
   - Human hg38: Download from UCSC
   - Download chromosome sizes:
     ```bash
     mysql --host=genome-mysql.soe.ucsc.edu --user=genome -N -A \
       -e "SELECT chrom,size FROM hg38.chromInfo" > hg38.genome
     ```

### Error: "Genome file required"

**Full message**:
```
Error: genome file is required for this operation
```

**Cause**: Commands like `slop`, `complement`, and `genomecov` need chromosome sizes.

**Solution**: Provide a genome file (tab-separated chromosome and size):

```bash
# Create genome file
echo -e "chr1\t248956422\nchr2\t242193529\nchr3\t198295559" > genome.txt

# Use it
grit slop -i input.bed -g genome.txt -b 1000 > extended.bed
```

---

## BED Format Errors

### Error: "Expected at least 3 fields"

**Full message**:
```
Error: Parse error at line 5: Expected at least 3 fields, got 2
```

**Cause**: BED format requires at least 3 tab-separated fields (chrom, start, end).

**Common causes**:
- Space-separated instead of tab-separated
- Missing columns
- Header lines without `#` prefix

**Solutions**:

1. **Convert spaces to tabs**:
   ```bash
   sed 's/ \+/\t/g' input.txt > input.bed
   ```

2. **Check file format**:
   ```bash
   # Show tabs as visible characters
   cat -A input.bed | head
   # Tabs appear as ^I
   ```

3. **Add header comment**:
   ```bash
   # If your file has a header, prefix with #
   sed '1s/^/#/' input.bed > fixed.bed
   ```

### Error: "start > end"

**Full message**:
```
Error: Parse error at line 10: start (500) > end (100)
```

**Cause**: BED format requires start ≤ end.

**Solutions**:

1. **Swap coordinates**:
   ```bash
   awk -F'\t' 'BEGIN{OFS="\t"} {if($2>$3){t=$2;$2=$3;$3=t} print}' input.bed > fixed.bed
   ```

2. **Check for parsing issues** (wrong column order):
   ```bash
   head input.bed
   # Verify: chrom<TAB>start<TAB>end
   ```

### Error: "Invalid coordinate"

**Full message**:
```
Error: Parse error at line 3: invalid coordinate 'abc'
```

**Cause**: Non-numeric value in start or end column.

**Solutions**:

1. **Check for header lines**:
   ```bash
   # Remove header or add # prefix
   tail -n +2 input.bed > no_header.bed
   # Or:
   sed '1s/^/#/' input.bed > fixed.bed
   ```

2. **Find problematic lines**:
   ```bash
   awk -F'\t' '$2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/' input.bed
   ```

---

## Memory Issues

### High Memory Usage

**Symptom**: Process uses excessive RAM, system becomes slow or swaps.

**Cause**: Large files loaded entirely into memory in parallel mode.

**Solutions**:

1. **Use streaming mode**:
   ```bash
   grit intersect -a large.bed -b large.bed --streaming --assume-sorted
   ```

2. **Process chromosomes separately**:
   ```bash
   for chr in chr1 chr2 chr3; do
       grep "^$chr\t" large.bed > ${chr}.bed
       grit merge -i ${chr}.bed > ${chr}_merged.bed
   done
   cat *_merged.bed > all_merged.bed
   ```

3. **Limit threads** (reduces memory):
   ```bash
   RAYON_NUM_THREADS=2 grit intersect -a a.bed -b b.bed
   ```

### Out of Memory (OOM)

**Symptom**: Process killed by system, "out of memory" error.

**Solutions**:

1. **Switch to streaming mode** (see above)

2. **Increase swap space** (temporary fix):
   ```bash
   # Linux: Add swap file
   sudo fallocate -l 8G /swapfile
   sudo chmod 600 /swapfile
   sudo mkswap /swapfile
   sudo swapon /swapfile
   ```

3. **Use a machine with more RAM**

---

## Performance Issues

### Slow Processing

**Possible causes and solutions**:

1. **Unsorted input being re-sorted**:
   ```bash
   # Pre-sort once, reuse
   grit sort -i input.bed > sorted.bed
   ```

2. **Wrong mode for data size**:
   - Large files: Use `--streaming`
   - Small files: Use parallel mode (default)

3. **High overlap density**:
   ```bash
   # Use filtering to reduce output
   grit intersect -a a.bed -b b.bed -u  # unique only
   grit intersect -a a.bed -b b.bed -f 0.5  # minimum overlap
   ```

4. **I/O bottleneck**:
   ```bash
   # Use SSD if available
   # Avoid network filesystems for large files
   ```

### No Output

**Possible causes**:

1. **No overlaps exist**: Verify intervals actually overlap
   ```bash
   # Check chromosome names match
   cut -f1 a.bed | sort -u
   cut -f1 b.bed | sort -u
   ```

2. **Coordinates don't overlap**: Check your data
   ```bash
   # Sample from each file
   head a.bed b.bed
   ```

3. **Filtering too strict**:
   ```bash
   # Remove fraction requirement to test
   grit intersect -a a.bed -b b.bed  # without -f flag
   ```

---

## Python-Specific Issues

### Import Error: "No module named 'pygrit'"

**Solution**: Install the package:
```bash
pip install grit-genomics
```

Note: Package name is `grit-genomics`, import name is `pygrit`.

### Incorrect Results in Python

**Cause**: Input files not sorted (Python API uses streaming internally).

**Solution**: Always sort input files:
```python
import pygrit

# Sort first
pygrit.sort("unsorted.bed", output="sorted.bed")

# Then process
result = pygrit.intersect("sorted.bed", "other_sorted.bed")
```

### TypeError: "expected str, got bytes"

**Cause**: Passing bytes instead of file path string.

**Solution**:
```python
# Wrong
pygrit.intersect(b"a.bed", b"b.bed")

# Correct
pygrit.intersect("a.bed", "b.bed")
```

---

## File Issues

### Error: "No such file or directory"

**Solutions**:

1. **Check file exists**:
   ```bash
   ls -la input.bed
   ```

2. **Use absolute path**:
   ```bash
   grit intersect -a /full/path/to/a.bed -b /full/path/to/b.bed
   ```

3. **Check permissions**:
   ```bash
   chmod +r input.bed
   ```

### Error: "Permission denied"

**Solutions**:

1. **Check read permissions**:
   ```bash
   ls -la input.bed
   chmod +r input.bed
   ```

2. **Check output directory is writable**:
   ```bash
   ls -la output_directory/
   chmod +w output_directory/
   ```

---

## Getting Help

If your issue isn't covered here:

1. **Check the documentation**:
   - [STREAMING_MODEL.md](STREAMING_MODEL.md) - Memory and algorithm details
   - [PERFORMANCE.md](PERFORMANCE.md) - Performance tuning
   - [COMMANDS.md](COMMANDS.md) - Command reference

2. **Search existing issues**:
   - [GitHub Issues](https://github.com/manish59/grit/issues)

3. **Open a new issue** with:
   - GRIT version (`grit --version`)
   - Command that failed
   - Full error message
   - Sample input (if possible)
   - Operating system

---

## Common Fixes Summary

```bash
# Fix: Unsorted input
grit sort -i input.bed > sorted.bed

# Fix: Memory issues
grit intersect -a a.bed -b b.bed --streaming --assume-sorted

# Fix: Chromosome naming mismatch
# Check if files use "chr1" vs "1"
cut -f1 a.bed | head
cut -f1 b.bed | head

# Fix: Space-separated file
sed 's/ \+/\t/g' input.txt > input.bed

# Fix: Header in file
sed '1s/^/#/' input.bed > fixed.bed
# Or remove it:
tail -n +2 input.bed > no_header.bed
```
