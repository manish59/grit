# Common Workflows

Real-world examples of using pygrit for genomic analysis.

## ChIP-seq Peak Analysis

### Find Peaks Overlapping Promoters

```python
import pygrit

# Find peaks that overlap promoter regions
overlapping_peaks = pygrit.intersect(
    "peaks.bed",
    "promoters.bed",
    unique=True,  # Each peak once
)

print(f"Found {len(overlapping_peaks)} peaks in promoters")
```

### Require Significant Overlap

```python
# Peaks must have at least 50% overlap with promoter
significant_peaks = pygrit.intersect(
    "peaks.bed",
    "promoters.bed",
    fraction=0.5,
    unique=True,
)
```

### Find Peaks NOT in Known Regions

```python
# Find novel peaks (not overlapping known enhancers)
novel_peaks = pygrit.intersect(
    "peaks.bed",
    "known_enhancers.bed",
    no_overlap=True,
)

print(f"Found {len(novel_peaks)} novel peaks")
```

---

## RNA-seq Analysis

### Count Reads Per Gene

```python
# Coverage of gene regions by aligned reads
pygrit.coverage(
    "genes.bed",
    "aligned_reads.bed",
    mean=True,
    output="gene_coverage.bed"
)
```

### Find Reads in Exons vs Introns

```python
# Reads overlapping exons
exon_reads = pygrit.intersect(
    "reads.bed",
    "exons.bed",
    count=True,
    output="exon_counts.bed"
)

# Reads in introns (gene body minus exons)
intron_regions = pygrit.subtract(
    "gene_bodies.bed",
    "exons.bed",
    output="introns.bed"
)

pygrit.intersect(
    "reads.bed",
    "introns.bed",
    count=True,
    output="intron_counts.bed"
)
```

---

## Variant Analysis

### Filter Variants by Region

```python
# Keep variants in coding regions
coding_variants = pygrit.intersect(
    "variants.bed",
    "coding_regions.bed",
)

# Remove variants in repetitive regions
filtered_variants = pygrit.subtract(
    "variants.bed",
    "repeats.bed",
)
```

### Find Nearest Gene

```python
# For each variant, find the closest gene
pygrit.closest(
    "variants.bed",
    "genes.bed",
    output="variant_nearest_gene.bed"
)
```

### Variants in Promoter Windows

```python
# Variants within 2kb upstream of genes
pygrit.window(
    "gene_starts.bed",
    "variants.bed",
    left=2000,
    right=0,
    output="promoter_variants.bed"
)
```

---

## Genome Annotation

### Merge Overlapping Annotations

```python
# Merge overlapping gene annotations
merged_genes = pygrit.merge("all_gene_annotations.bed")

# Merge nearby features (within 1kb)
merged_features = pygrit.merge(
    "features.bed",
    distance=1000,
)
```

### Create Non-Overlapping Bins

```python
import numpy as np
import pygrit

# Create 10kb bins across a chromosome
chrom_length = 250_000_000  # 250 Mb
bin_size = 10_000

starts = np.arange(0, chrom_length, bin_size)
ends = np.minimum(starts + bin_size, chrom_length)

bins = pygrit.from_numpy("chr1", np.column_stack([starts, ends]))
```

### Subtract Known Regions

```python
# Find intergenic regions
intergenic = pygrit.subtract(
    "chromosome.bed",  # Whole chromosome
    "genes.bed",
    output="intergenic.bed"
)
```

---

## Pipeline Example

### Complete Analysis Pipeline

```python
import pygrit
from pathlib import Path

def analyze_chip_peaks(
    peaks_file: str,
    genes_file: str,
    enhancers_file: str,
    output_dir: str,
):
    """Comprehensive ChIP-seq peak analysis."""

    out = Path(output_dir)
    out.mkdir(exist_ok=True)

    # 1. Merge nearby peaks
    pygrit.merge(
        peaks_file,
        distance=500,
        output=str(out / "merged_peaks.bed")
    )

    # 2. Find promoter peaks
    pygrit.intersect(
        str(out / "merged_peaks.bed"),
        genes_file,
        fraction=0.5,
        unique=True,
        output=str(out / "promoter_peaks.bed")
    )

    # 3. Find enhancer peaks
    pygrit.intersect(
        str(out / "merged_peaks.bed"),
        enhancers_file,
        unique=True,
        output=str(out / "enhancer_peaks.bed")
    )

    # 4. Find novel peaks (not promoter or enhancer)
    pygrit.intersect(
        str(out / "merged_peaks.bed"),
        genes_file,
        no_overlap=True,
        output=str(out / "temp_not_promoter.bed")
    )

    pygrit.intersect(
        str(out / "temp_not_promoter.bed"),
        enhancers_file,
        no_overlap=True,
        output=str(out / "novel_peaks.bed")
    )

    # 5. Count results
    promoter_peaks = pygrit.read_bed(str(out / "promoter_peaks.bed"))
    enhancer_peaks = pygrit.read_bed(str(out / "enhancer_peaks.bed"))
    novel_peaks = pygrit.read_bed(str(out / "novel_peaks.bed"))

    return {
        "promoter": len(promoter_peaks),
        "enhancer": len(enhancer_peaks),
        "novel": len(novel_peaks),
    }

# Run analysis
results = analyze_chip_peaks(
    "chip_peaks.bed",
    "gene_promoters.bed",
    "enhancer_regions.bed",
    "analysis_output",
)

print(f"Promoter peaks: {results['promoter']}")
print(f"Enhancer peaks: {results['enhancer']}")
print(f"Novel peaks: {results['novel']}")
```

---

## Parallel Processing

### Process Multiple Chromosomes

```python
from concurrent.futures import ThreadPoolExecutor
import pygrit

chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

def process_chromosome(chrom):
    """Process a single chromosome."""
    pygrit.intersect(
        f"peaks_{chrom}.bed",
        f"genes_{chrom}.bed",
        output=f"results_{chrom}.bed"
    )
    return chrom

# Process in parallel (pygrit releases GIL)
with ThreadPoolExecutor(max_workers=8) as executor:
    completed = list(executor.map(process_chromosome, chromosomes))

print(f"Processed {len(completed)} chromosomes")
```

### Batch Processing

```python
from pathlib import Path

def batch_merge(input_dir: str, output_dir: str, distance: int = 100):
    """Merge all BED files in a directory."""

    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    for bed_file in input_path.glob("*.bed"):
        output_file = output_path / f"merged_{bed_file.name}"
        pygrit.merge(
            str(bed_file),
            distance=distance,
            output=str(output_file)
        )
        print(f"Merged: {bed_file.name}")

batch_merge("raw_peaks", "merged_peaks", distance=500)
```
