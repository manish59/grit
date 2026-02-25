#!/bin/bash
# Generate test datasets for bedtools-rs validation

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATASETS="$PROJECT_ROOT/benches/datasets"

mkdir -p "$DATASETS"

generate_bed() {
    local count=$1
    local output=$2
    local min_size=${3:-50}
    local max_size=${4:-5000}

    echo "Generating $count intervals -> $output"

    awk -v count="$count" -v min_size="$min_size" -v max_size="$max_size" '
    BEGIN {
        srand()
        # Chromosome names
        split("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY", chroms)
        # hg19 chromosome sizes (must use actual sizes for complement/genomecov)
        sizes["chr1"]=249250621; sizes["chr2"]=243199373; sizes["chr3"]=198022430
        sizes["chr4"]=191154276; sizes["chr5"]=180915260; sizes["chr6"]=171115067
        sizes["chr7"]=159138663; sizes["chr8"]=146364022; sizes["chr9"]=141213431
        sizes["chr10"]=135534747; sizes["chr11"]=135006516; sizes["chr12"]=133851895
        sizes["chr13"]=115169878; sizes["chr14"]=107349540; sizes["chr15"]=102531392
        sizes["chr16"]=90354753; sizes["chr17"]=81195210; sizes["chr18"]=78077248
        sizes["chr19"]=59128983; sizes["chr20"]=63025520; sizes["chr21"]=48129895
        sizes["chr22"]=51304566; sizes["chrX"]=155270560; sizes["chrY"]=59373566
        strands[1]="+"; strands[2]="-"

        for (i = 1; i <= count; i++) {
            chrom_idx = int(rand() * 24) + 1
            chrom = chroms[chrom_idx]
            chrom_size = sizes[chrom]
            interval_size = int(rand() * (max_size - min_size)) + min_size
            start = int(rand() * (chrom_size - interval_size))
            end = start + interval_size
            strand_idx = int(rand() * 2) + 1
            strand = strands[strand_idx]
            score = int(rand() * 1000)
            printf "%s\t%d\t%d\tinterval_%d\t%d\t%s\n", chrom, start, end, i, score, strand
        }
    }' > "$output"

    sort -k1,1V -k2,2n "$output" -o "$output"
    echo "  Done: $(wc -l < "$output") intervals"
}

generate_high_overlap() {
    local count=$1
    local output=$2

    echo "Generating $count high-overlap intervals -> $output"

    awk -v count="$count" '
    BEGIN {
        srand()
        split("chr1 chr2 chr3 chr4 chr5", chroms)
        # Use chr1-5 sizes (all > 180M so hotspots up to 170M are safe)
        sizes["chr1"]=249250621; sizes["chr2"]=243199373; sizes["chr3"]=198022430
        sizes["chr4"]=191154276; sizes["chr5"]=180915260

        for (i = 1; i <= count; i++) {
            chrom_idx = int(rand() * 5) + 1
            chrom = chroms[chrom_idx]
            chrom_size = sizes[chrom]
            # Hotspots at 10M, 20M, ... up to 170M (safe for all chr1-5)
            hotspot = int(rand() * 17) * 10000000
            offset = int(rand() * 100000)
            start = hotspot + offset
            size = int(rand() * 1000) + 100
            end = start + size
            # Clamp to chromosome bounds
            if (end > chrom_size) {
                end = chrom_size
                start = end - size
                if (start < 0) start = 0
            }
            strand = (rand() > 0.5) ? "+" : "-"
            score = int(rand() * 1000)
            printf "%s\t%d\t%d\tho_%d\t%d\t%s\n", chrom, start, end, i, score, strand
        }
    }' > "$output"

    sort -k1,1V -k2,2n "$output" -o "$output"
    echo "  Done: $(wc -l < "$output") intervals"
}

echo "=== Generating benchmark datasets ==="
echo ""

# Small datasets
generate_bed 1000 "$DATASETS/1k_a.bed"
generate_bed 1000 "$DATASETS/1k_b.bed"

# Medium datasets
generate_bed 10000 "$DATASETS/10k_a.bed"
generate_bed 10000 "$DATASETS/10k_b.bed"

# Large datasets
generate_bed 100000 "$DATASETS/100k_a.bed"
generate_bed 100000 "$DATASETS/100k_b.bed"

# 1M datasets (the real test)
echo ""
echo "Generating 1M interval datasets..."
generate_bed 1000000 "$DATASETS/1m_a.bed"
generate_bed 1000000 "$DATASETS/1m_b.bed"

# High overlap datasets for stress testing
echo ""
echo "Generating high-overlap datasets..."
generate_high_overlap 100000 "$DATASETS/100k_high_overlap_a.bed"
generate_high_overlap 100000 "$DATASETS/100k_high_overlap_b.bed"

echo ""
echo "=== Dataset generation complete ==="
echo ""
ls -lh "$DATASETS"/*.bed
