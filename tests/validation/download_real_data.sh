#!/bin/bash
# Download real genomic datasets for benchmarking
# Sources: UCSC, direct BED files

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATA_DIR="$PROJECT_ROOT/benches/real_data"

mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

echo "=============================================="
echo "  Downloading Real Genomic Datasets"
echo "=============================================="
echo ""

download_ucsc() {
    local table="$1"
    local output="$2"
    local description="$3"

    if [ -f "$output" ]; then
        echo "✓ $description already exists"
        return
    fi

    echo "Downloading: $description"
    local url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/${table}.txt.gz"
    
    if curl -sL --fail "$url" | gunzip -c > "${output}.tmp" 2>/dev/null; then
        mv "${output}.tmp" "$output"
        lines=$(wc -l < "$output" | tr -d ' ')
        echo "  Downloaded: $lines lines"
    else
        echo "  Failed to download $description"
        rm -f "${output}.tmp"
        return 1
    fi
}

echo "=== UCSC Genome Annotations (hg38) ==="
echo ""

# CpG Islands - reliable UCSC source
download_ucsc "cpgIslandExt" "cpg_islands_raw.txt" "CpG Islands"

if [ -f "cpg_islands_raw.txt" ] && [ ! -f "cpg_islands.bed" ]; then
    echo "Converting CpG Islands to BED format..."
    awk -F'\t' 'BEGIN{OFS="\t"} $2 ~ /^chr[0-9XYM]+$/ {print $2, $3, $4, "cpg_"NR, $6, "+"}' cpg_islands_raw.txt 2>/dev/null | \
        sort -k1,1V -k2,2n > cpg_islands.bed
    lines=$(wc -l < "cpg_islands.bed" | tr -d ' ')
    echo "  Converted: $lines CpG islands"
fi

# Simple Repeats - large dataset
download_ucsc "simpleRepeat" "simple_repeats_raw.txt" "Simple Repeats"

if [ -f "simple_repeats_raw.txt" ] && [ ! -f "simple_repeats.bed" ]; then
    echo "Converting Simple Repeats to BED format..."
    awk -F'\t' 'BEGIN{OFS="\t"} $2 ~ /^chr[0-9XYM]+$/ {print $2, $3, $4, $5, 0, "+"}' simple_repeats_raw.txt 2>/dev/null | \
        sort -k1,1V -k2,2n > simple_repeats.bed
    lines=$(wc -l < "simple_repeats.bed" | tr -d ' ')
    echo "  Converted: $lines repeat intervals"
fi

# DNase clusters
download_ucsc "wgEncodeRegDnaseClustered" "dnase_clusters_raw.txt" "ENCODE DNase Clusters"

if [ -f "dnase_clusters_raw.txt" ] && [ ! -f "dnase_clusters.bed" ]; then
    echo "Converting DNase Clusters to BED format..."
    awk -F'\t' 'BEGIN{OFS="\t"} $2 ~ /^chr[0-9XYM]+$/ {print $2, $3, $4, "dnase_"NR, $5, "+"}' dnase_clusters_raw.txt 2>/dev/null | \
        sort -k1,1V -k2,2n > dnase_clusters.bed
    lines=$(wc -l < "dnase_clusters.bed" | tr -d ' ')
    echo "  Converted: $lines DNase clusters"
fi

# Transcription Factor Binding Sites (hg38 uses encRegTfbsClustered)
download_ucsc "encRegTfbsClustered" "tfbs_raw.txt" "ENCODE TF Binding Sites"

if [ -f "tfbs_raw.txt" ] && [ ! -f "tfbs_clusters.bed" ]; then
    echo "Converting TFBS to BED format..."
    awk -F'\t' 'BEGIN{OFS="\t"} $2 ~ /^chr[0-9XYM]+$/ {print $2, $3, $4, $5, $6, "+"}' tfbs_raw.txt 2>/dev/null | \
        sort -k1,1V -k2,2n > tfbs_clusters.bed
    lines=$(wc -l < "tfbs_clusters.bed" | tr -d ' ')
    echo "  Converted: $lines TFBS clusters"
fi

# RefSeq genes
download_ucsc "ncbiRefSeq" "refseq_raw.txt" "RefSeq Genes"

if [ -f "refseq_raw.txt" ] && [ ! -f "refseq_genes.bed" ]; then
    echo "Converting RefSeq to BED format..."
    awk -F'\t' 'BEGIN{OFS="\t"} $3 ~ /^chr[0-9XYM]+$/ {print $3, $5, $6, $2, 0, $4}' refseq_raw.txt 2>/dev/null | \
        sort -k1,1V -k2,2n > refseq_genes.bed
    lines=$(wc -l < "refseq_genes.bed" | tr -d ' ')
    echo "  Converted: $lines gene records"
fi

# Download hg38 chromosome sizes (genome file)
echo ""
echo "=== Downloading Genome File ==="
echo ""

if [ ! -f "hg38.genome" ]; then
    echo "Downloading hg38 chromosome sizes..."
    curl -sL "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes" | \
        grep "^chr" | grep -v "_" | sort -k1,1V > hg38.genome
    lines=$(wc -l < "hg38.genome" | tr -d ' ')
    echo "  Downloaded: $lines chromosomes"
else
    echo "✓ hg38.genome already exists"
fi

# Create derived datasets
echo ""
echo "=== Creating Derived Datasets ==="
echo ""

# Create promoter regions from genes
if [ -f "refseq_genes.bed" ] && [ ! -f "promoters.bed" ]; then
    echo "Creating promoter regions (TSS ± 1kb)..."
    awk -F'\t' 'BEGIN{OFS="\t"} {
        if ($6 == "+") { start = ($2 > 1000) ? $2 - 1000 : 0; end = $2 + 1000 }
        else { start = ($3 > 1000) ? $3 - 1000 : 0; end = $3 + 1000 }
        print $1, start, end, $4"_promoter", 0, $6
    }' refseq_genes.bed | sort -k1,1V -k2,2n > promoters.bed
    lines=$(wc -l < "promoters.bed" | tr -d ' ')
    echo "  Created: $lines promoter regions"
fi

echo ""
echo "=== Dataset Summary ==="
echo ""

printf "%-30s %15s %10s\n" "Dataset" "Intervals" "Size"
printf "%-30s %15s %10s\n" "-------" "---------" "----"

for f in *.bed; do
    if [ -f "$f" ]; then
        lines=$(wc -l < "$f" | tr -d ' ')
        size=$(du -h "$f" | cut -f1)
        printf "%-30s %15s %10s\n" "$f" "$lines" "$size"
    fi
done

# Cleanup
echo ""
rm -f *_raw.txt 2>/dev/null || true
echo "Cleaned up temporary files"

echo ""
echo "=============================================="
echo "  Download Complete"
echo "=============================================="
