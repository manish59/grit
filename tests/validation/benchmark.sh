#!/bin/bash
# Performance benchmark: bedtools-rs vs bedtools
# Records runtime, memory, and throughput
#
# Usage:
#   ./benchmark.sh              # Run all benchmarks
#   ./benchmark.sh complement   # Run only complement benchmarks
#   ./benchmark.sh intersect merge  # Run intersect and merge benchmarks
#
# Available commands: intersect, merge, sort, subtract, closest, coverage,
#                     slop, complement, genomecov, jaccard, multiinter, window

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATASETS="$PROJECT_ROOT/benches/datasets"
RESULTS="$PROJECT_ROOT/tests/validation/results"
BEDTOOLS_RS="$PROJECT_ROOT/target/release/rbedtools"

mkdir -p "$RESULTS"

# Parse command line arguments for selective benchmark running
RUN_ALL=true
COMMANDS_TO_RUN=""
if [ $# -gt 0 ]; then
    RUN_ALL=false
    for cmd in "$@"; do
        COMMANDS_TO_RUN="$COMMANDS_TO_RUN $cmd "
    done
fi

should_run() {
    local cmd="$1"
    if [ "$RUN_ALL" = true ]; then
        return 0
    fi
    case "$COMMANDS_TO_RUN" in
        *" $cmd "*) return 0 ;;
        *) return 1 ;;
    esac
}

# Check for hyperfine
if ! command -v hyperfine &> /dev/null; then
    echo "hyperfine not found. Installing via brew..."
    brew install hyperfine
fi

# Build release binary with optimizations
echo "Building bedtools-rs with full optimizations..."
RUSTFLAGS="-C target-cpu=native" cargo build --release --manifest-path "$PROJECT_ROOT/Cargo.toml"

# Check if datasets exist
if [ ! -f "$DATASETS/1m_a.bed" ]; then
    echo "Datasets not found. Generating..."
    bash "$SCRIPT_DIR/generate_data.sh"
fi

echo ""
echo "=============================================="
echo "  BEDTOOLS-RS PERFORMANCE BENCHMARKS"
echo "=============================================="
echo ""
echo "System info:"
echo "  CPU: $(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo 'Unknown')"
echo "  Cores: $(sysctl -n hw.ncpu 2>/dev/null || nproc)"
echo "  RAM: $(sysctl -n hw.memsize 2>/dev/null | awk '{print $1/1024/1024/1024 " GB"}')"
echo "  bedtools: $(bedtools --version)"
echo "  bedtools-rs: $($BEDTOOLS_RS --version 2>/dev/null || echo 'dev')"
echo ""

BENCHMARK_FILE="$RESULTS/benchmark_$(date +%Y%m%d_%H%M%S).md"

cat > "$BENCHMARK_FILE" << 'EOF'
# Bedtools-rs Performance Benchmark Results

## System Information
EOF

echo "- **Date**: $(date)" >> "$BENCHMARK_FILE"
echo "- **CPU**: $(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo 'Unknown')" >> "$BENCHMARK_FILE"
echo "- **Cores**: $(sysctl -n hw.ncpu 2>/dev/null || nproc)" >> "$BENCHMARK_FILE"
echo "- **RAM**: $(sysctl -n hw.memsize 2>/dev/null | awk '{print $1/1024/1024/1024 " GB"}')" >> "$BENCHMARK_FILE"
echo "- **OS**: $(sw_vers -productName 2>/dev/null) $(sw_vers -productVersion 2>/dev/null)" >> "$BENCHMARK_FILE"
echo "- **bedtools version**: $(bedtools --version)" >> "$BENCHMARK_FILE"
echo "" >> "$BENCHMARK_FILE"

run_benchmark() {
    local name="$1"
    local bedtools_cmd="$2"
    local rs_cmd="$3"
    local warmup=${4:-3}
    local runs=${5:-10}

    echo "Benchmarking: $name"
    echo ""

    echo "## $name" >> "$BENCHMARK_FILE"
    echo "" >> "$BENCHMARK_FILE"

    # Run hyperfine comparison
    hyperfine \
        --warmup "$warmup" \
        --runs "$runs" \
        --export-markdown "$RESULTS/temp_bench.md" \
        --command-name "bedtools" "$bedtools_cmd" \
        --command-name "bedtools-rs" "$rs_cmd"

    cat "$RESULTS/temp_bench.md" >> "$BENCHMARK_FILE"
    echo "" >> "$BENCHMARK_FILE"
}

if should_run "intersect"; then
echo "=============================================="
echo "  INTERSECT BENCHMARKS"
echo "=============================================="
echo ""

# 10K intersect
run_benchmark "Intersect 10K x 10K" \
    "bedtools intersect -sorted -a $DATASETS/10k_a.bed -b $DATASETS/10k_b.bed > /dev/null" \
    "$BEDTOOLS_RS intersect -a $DATASETS/10k_a.bed -b $DATASETS/10k_b.bed > /dev/null"

# 100K intersect
run_benchmark "Intersect 100K x 100K" \
    "bedtools intersect -sorted -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed > /dev/null" \
    "$BEDTOOLS_RS intersect -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed > /dev/null"

# 1M intersect
run_benchmark "Intersect 1M x 1M" \
    "bedtools intersect -sorted -a $DATASETS/1m_a.bed -b $DATASETS/1m_b.bed > /dev/null" \
    "$BEDTOOLS_RS intersect -a $DATASETS/1m_a.bed -b $DATASETS/1m_b.bed > /dev/null" \
    2 5

# With flags
run_benchmark "Intersect 100K -wa -wb" \
    "bedtools intersect -sorted -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed -wa -wb > /dev/null" \
    "$BEDTOOLS_RS intersect -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed --wa --wb > /dev/null"

# High overlap stress test
run_benchmark "Intersect 100K High Overlap" \
    "bedtools intersect -sorted -a $DATASETS/100k_high_overlap_a.bed -b $DATASETS/100k_high_overlap_b.bed > /dev/null" \
    "$BEDTOOLS_RS intersect -a $DATASETS/100k_high_overlap_a.bed -b $DATASETS/100k_high_overlap_b.bed > /dev/null"
fi

if should_run "merge"; then
echo "=============================================="
echo "  MERGE BENCHMARKS"
echo "=============================================="
echo ""

run_benchmark "Merge 100K" \
    "bedtools merge -i $DATASETS/100k_a.bed > /dev/null" \
    "$BEDTOOLS_RS merge -i $DATASETS/100k_a.bed > /dev/null"

run_benchmark "Merge 1M" \
    "bedtools merge -i $DATASETS/1m_a.bed > /dev/null" \
    "$BEDTOOLS_RS merge -i $DATASETS/1m_a.bed > /dev/null" \
    2 5

run_benchmark "Merge 100K High Overlap" \
    "bedtools merge -i $DATASETS/100k_high_overlap_a.bed > /dev/null" \
    "$BEDTOOLS_RS merge -i $DATASETS/100k_high_overlap_a.bed > /dev/null"
fi

if should_run "sort"; then
echo "=============================================="
echo "  SORT BENCHMARKS"
echo "=============================================="
echo ""

# Create unsorted version of 100k
if [ ! -f "$DATASETS/100k_unsorted.bed" ]; then
    shuf "$DATASETS/100k_a.bed" > "$DATASETS/100k_unsorted.bed"
fi

if [ ! -f "$DATASETS/1m_unsorted.bed" ]; then
    shuf "$DATASETS/1m_a.bed" > "$DATASETS/1m_unsorted.bed"
fi

run_benchmark "Sort 100K" \
    "bedtools sort -i $DATASETS/100k_unsorted.bed > /dev/null" \
    "$BEDTOOLS_RS sort -i $DATASETS/100k_unsorted.bed > /dev/null"

run_benchmark "Sort 1M" \
    "bedtools sort -i $DATASETS/1m_unsorted.bed > /dev/null" \
    "$BEDTOOLS_RS sort -i $DATASETS/1m_unsorted.bed > /dev/null" \
    2 5
fi

if should_run "subtract"; then
echo "=============================================="
echo "  SUBTRACT BENCHMARKS"
echo "=============================================="
echo ""

run_benchmark "Subtract 100K x 100K" \
    "bedtools subtract -sorted -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed > /dev/null" \
    "$BEDTOOLS_RS subtract -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed > /dev/null"

run_benchmark "Subtract 1M x 1M" \
    "bedtools subtract -sorted -a $DATASETS/1m_a.bed -b $DATASETS/1m_b.bed > /dev/null" \
    "$BEDTOOLS_RS subtract -a $DATASETS/1m_a.bed -b $DATASETS/1m_b.bed > /dev/null" \
    2 5
fi

if should_run "closest"; then
echo "=============================================="
echo "  CLOSEST BENCHMARKS"
echo "=============================================="
echo ""

run_benchmark "Closest 10K x 10K" \
    "bedtools closest -a $DATASETS/10k_a.bed -b $DATASETS/10k_b.bed > /dev/null" \
    "$BEDTOOLS_RS closest -a $DATASETS/10k_a.bed -b $DATASETS/10k_b.bed > /dev/null"

run_benchmark "Closest 100K x 100K" \
    "bedtools closest -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed > /dev/null" \
    "$BEDTOOLS_RS closest -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed > /dev/null"
fi

if should_run "coverage"; then
echo "=============================================="
echo "  COVERAGE BENCHMARKS"
echo "=============================================="
echo ""

run_benchmark "Coverage 100K x 100K" \
    "bedtools coverage -sorted -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed > /dev/null" \
    "$BEDTOOLS_RS coverage -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed > /dev/null"

run_benchmark "Coverage 1M x 1M" \
    "bedtools coverage -sorted -a $DATASETS/1m_a.bed -b $DATASETS/1m_b.bed > /dev/null" \
    "$BEDTOOLS_RS coverage -a $DATASETS/1m_a.bed -b $DATASETS/1m_b.bed > /dev/null" \
    2 5

run_benchmark "Coverage 100K -hist" \
    "bedtools coverage -sorted -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed -hist > /dev/null" \
    "$BEDTOOLS_RS coverage -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed --hist > /dev/null"

run_benchmark "Coverage 100K -mean" \
    "bedtools coverage -sorted -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed -mean > /dev/null" \
    "$BEDTOOLS_RS coverage -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed --mean > /dev/null"
fi

# Create genome file if it doesn't exist (hg19 chromosome sizes)
GENOME_FILE="$DATASETS/test.genome"
if [ ! -f "$GENOME_FILE" ] || [ $(wc -l < "$GENOME_FILE") -lt 22 ]; then
    echo "Creating test genome file..."
    cat > "$GENOME_FILE" << 'GENOME_EOF'
chr1	249250621
chr2	243199373
chr3	198022430
chr4	191154276
chr5	180915260
chr6	171115067
chr7	159138663
chr8	146364022
chr9	141213431
chr10	135534747
chr11	135006516
chr12	133851895
chr13	115169878
chr14	107349540
chr15	102531392
chr16	90354753
chr17	81195210
chr18	78077248
chr19	59128983
chr20	63025520
chr21	48129895
chr22	51304566
chrX	155270560
chrY	59373566
chrM	16571
GENOME_EOF
fi

if should_run "slop"; then
echo "=============================================="
echo "  SLOP BENCHMARKS"
echo "=============================================="
echo ""

run_benchmark "Slop 100K (extend 100bp)" \
    "bedtools slop -i $DATASETS/100k_a.bed -g $GENOME_FILE -b 100 > /dev/null" \
    "$BEDTOOLS_RS slop -i $DATASETS/100k_a.bed -g $GENOME_FILE -b 100 > /dev/null"

run_benchmark "Slop 1M (extend 100bp)" \
    "bedtools slop -i $DATASETS/1m_a.bed -g $GENOME_FILE -b 100 > /dev/null" \
    "$BEDTOOLS_RS slop -i $DATASETS/1m_a.bed -g $GENOME_FILE -b 100 > /dev/null" \
    2 5

run_benchmark "Slop 100K (extend l/r)" \
    "bedtools slop -i $DATASETS/100k_a.bed -g $GENOME_FILE -l 50 -r 150 > /dev/null" \
    "$BEDTOOLS_RS slop -i $DATASETS/100k_a.bed -g $GENOME_FILE -l 50 -r 150 > /dev/null"
fi

if should_run "complement"; then
echo "=============================================="
echo "  COMPLEMENT BENCHMARKS"
echo "=============================================="
echo ""

run_benchmark "Complement 100K" \
    "bedtools complement -i $DATASETS/100k_a.bed -g $GENOME_FILE > /dev/null" \
    "$BEDTOOLS_RS complement -i $DATASETS/100k_a.bed -g $GENOME_FILE > /dev/null"

run_benchmark "Complement 1M" \
    "bedtools complement -i $DATASETS/1m_a.bed -g $GENOME_FILE > /dev/null" \
    "$BEDTOOLS_RS complement -i $DATASETS/1m_a.bed -g $GENOME_FILE > /dev/null" \
    2 5
fi

if should_run "genomecov"; then
echo "=============================================="
echo "  GENOMECOV BENCHMARKS"
echo "=============================================="
echo ""

run_benchmark "Genomecov 100K (histogram)" \
    "bedtools genomecov -i $DATASETS/100k_a.bed -g $GENOME_FILE > /dev/null" \
    "$BEDTOOLS_RS genomecov -i $DATASETS/100k_a.bed -g $GENOME_FILE > /dev/null"

run_benchmark "Genomecov 100K (bedgraph)" \
    "bedtools genomecov -i $DATASETS/100k_a.bed -g $GENOME_FILE -bg > /dev/null" \
    "$BEDTOOLS_RS genomecov -i $DATASETS/100k_a.bed -g $GENOME_FILE -bg > /dev/null"

run_benchmark "Genomecov 1M (histogram)" \
    "bedtools genomecov -i $DATASETS/1m_a.bed -g $GENOME_FILE > /dev/null" \
    "$BEDTOOLS_RS genomecov -i $DATASETS/1m_a.bed -g $GENOME_FILE > /dev/null" \
    2 5

run_benchmark "Genomecov 1M (bedgraph)" \
    "bedtools genomecov -i $DATASETS/1m_a.bed -g $GENOME_FILE -bg > /dev/null" \
    "$BEDTOOLS_RS genomecov -i $DATASETS/1m_a.bed -g $GENOME_FILE -bg > /dev/null" \
    2 5
fi

if should_run "jaccard"; then
echo "=============================================="
echo "  JACCARD BENCHMARKS"
echo "=============================================="
echo ""

run_benchmark "Jaccard 10K x 10K" \
    "bedtools jaccard -a $DATASETS/10k_a.bed -b $DATASETS/10k_b.bed > /dev/null" \
    "$BEDTOOLS_RS jaccard -a $DATASETS/10k_a.bed -b $DATASETS/10k_b.bed > /dev/null"

run_benchmark "Jaccard 100K x 100K" \
    "bedtools jaccard -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed > /dev/null" \
    "$BEDTOOLS_RS jaccard -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed > /dev/null"

run_benchmark "Jaccard 1M x 1M" \
    "bedtools jaccard -a $DATASETS/1m_a.bed -b $DATASETS/1m_b.bed > /dev/null" \
    "$BEDTOOLS_RS jaccard -a $DATASETS/1m_a.bed -b $DATASETS/1m_b.bed > /dev/null" \
    2 5
fi

if should_run "multiinter"; then
echo "=============================================="
echo "  MULTIINTER BENCHMARKS"
echo "=============================================="
echo ""

run_benchmark "Multiinter 3 files (10K each)" \
    "bedtools multiinter -i $DATASETS/10k_a.bed $DATASETS/10k_b.bed $DATASETS/10k_a.bed > /dev/null" \
    "$BEDTOOLS_RS multiinter -i $DATASETS/10k_a.bed $DATASETS/10k_b.bed $DATASETS/10k_a.bed > /dev/null"

run_benchmark "Multiinter 3 files (100K each)" \
    "bedtools multiinter -i $DATASETS/100k_a.bed $DATASETS/100k_b.bed $DATASETS/100k_a.bed > /dev/null" \
    "$BEDTOOLS_RS multiinter -i $DATASETS/100k_a.bed $DATASETS/100k_b.bed $DATASETS/100k_a.bed > /dev/null"

run_benchmark "Multiinter 3 files (1M each)" \
    "bedtools multiinter -i $DATASETS/1m_a.bed $DATASETS/1m_b.bed $DATASETS/1m_a.bed > /dev/null" \
    "$BEDTOOLS_RS multiinter -i $DATASETS/1m_a.bed $DATASETS/1m_b.bed $DATASETS/1m_a.bed > /dev/null" \
    2 5
fi

if should_run "window"; then
echo "=============================================="
echo "  WINDOW BENCHMARKS"
echo "=============================================="
echo ""

run_benchmark "Window 10K x 10K (w=1000)" \
    "bedtools window -a $DATASETS/10k_a.bed -b $DATASETS/10k_b.bed -w 1000 > /dev/null" \
    "$BEDTOOLS_RS window -a $DATASETS/10k_a.bed -b $DATASETS/10k_b.bed -w 1000 > /dev/null"

run_benchmark "Window 100K x 100K (w=1000)" \
    "bedtools window -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed -w 1000 > /dev/null" \
    "$BEDTOOLS_RS window -a $DATASETS/100k_a.bed -b $DATASETS/100k_b.bed -w 1000 > /dev/null"
fi

# ============================================
# REAL GENOMIC DATA BENCHMARKS
# ============================================

REAL_DATA="$PROJECT_ROOT/benches/real_data"

# Check if real data exists, download if not
if [ ! -f "$REAL_DATA/encode_h3k4me3_k562.bed" ]; then
    echo ""
    echo "=============================================="
    echo "  DOWNLOADING REAL GENOMIC DATASETS"
    echo "=============================================="
    echo ""
    bash "$SCRIPT_DIR/download_real_data.sh"
fi

# Only run real data benchmarks if files exist
echo "# Real Genomic Data Benchmarks" >> "$BENCHMARK_FILE"
echo "" >> "$BENCHMARK_FILE"

if [ -f "$REAL_DATA/cpg_islands.bed" ] && [ -f "$REAL_DATA/dnase_clusters.bed" ]; then
    if should_run "real" || should_run "intersect"; then
    echo ""
    echo "=============================================="
    echo "  REAL DATA: CpG Islands vs DNase (32K x 2.1M)"
    echo "=============================================="
    echo ""

    # CpG islands vs DNase clusters
    run_benchmark "Real: CpG Islands ∩ DNase Clusters" \
        "bedtools intersect -sorted -a $REAL_DATA/cpg_islands.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null" \
        "$BEDTOOLS_RS intersect -a $REAL_DATA/cpg_islands.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null"
    fi

    if should_run "real" || should_run "coverage"; then
    # Coverage of CpG islands by DNase
    run_benchmark "Real: DNase coverage of CpG Islands" \
        "bedtools coverage -sorted -a $REAL_DATA/cpg_islands.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null" \
        "$BEDTOOLS_RS coverage -a $REAL_DATA/cpg_islands.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null"
    fi

    if should_run "real" || should_run "subtract"; then
    # Subtract DNase from CpG islands
    run_benchmark "Real: CpG Islands - DNase" \
        "bedtools subtract -sorted -a $REAL_DATA/cpg_islands.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null" \
        "$BEDTOOLS_RS subtract -a $REAL_DATA/cpg_islands.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null"
    fi

    if should_run "real" || should_run "closest"; then
    # Closest DNase to CpG islands
    run_benchmark "Real: Closest DNase to CpG Islands" \
        "bedtools closest -a $REAL_DATA/cpg_islands.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null" \
        "$BEDTOOLS_RS closest -a $REAL_DATA/cpg_islands.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null"
    fi
fi

if [ -f "$REAL_DATA/simple_repeats.bed" ]; then
    if should_run "real" || should_run "merge"; then
    echo ""
    echo "=============================================="
    echo "  REAL DATA: Simple Repeats (1M+ intervals)"
    echo "=============================================="
    echo ""

    # Merge overlapping repeats
    run_benchmark "Real: Merge Simple Repeats (1M)" \
        "bedtools merge -i $REAL_DATA/simple_repeats.bed > /dev/null" \
        "$BEDTOOLS_RS merge -i $REAL_DATA/simple_repeats.bed > /dev/null" \
        2 5
    fi

    if should_run "real" || should_run "sort"; then
    # Sort repeats (shuffled)
    if [ ! -f "$REAL_DATA/simple_repeats_shuffled.bed" ]; then
        shuf "$REAL_DATA/simple_repeats.bed" > "$REAL_DATA/simple_repeats_shuffled.bed"
    fi
    run_benchmark "Real: Sort Simple Repeats (1M)" \
        "bedtools sort -i $REAL_DATA/simple_repeats_shuffled.bed > /dev/null" \
        "$BEDTOOLS_RS sort -i $REAL_DATA/simple_repeats_shuffled.bed > /dev/null" \
        2 5
    fi
fi

if [ -f "$REAL_DATA/simple_repeats.bed" ] && [ -f "$REAL_DATA/dnase_clusters.bed" ]; then
    if should_run "real" || should_run "intersect"; then
    echo ""
    echo "=============================================="
    echo "  REAL DATA: Repeats vs DNase (1M x 2.1M)"
    echo "=============================================="
    echo ""

    # Large-scale intersect stress test
    run_benchmark "Real: Simple Repeats ∩ DNase (1Mx2.1M)" \
        "bedtools intersect -sorted -a $REAL_DATA/simple_repeats.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null" \
        "$BEDTOOLS_RS intersect -a $REAL_DATA/simple_repeats.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null" \
        1 3
    fi

    if should_run "real" || should_run "subtract"; then
    # Large-scale subtract
    run_benchmark "Real: Simple Repeats - DNase (1Mx2.1M)" \
        "bedtools subtract -sorted -a $REAL_DATA/simple_repeats.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null" \
        "$BEDTOOLS_RS subtract -a $REAL_DATA/simple_repeats.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null" \
        1 3
    fi
fi

if [ -f "$REAL_DATA/dnase_clusters.bed" ]; then
    if should_run "real" || should_run "merge"; then
    echo ""
    echo "=============================================="
    echo "  REAL DATA: DNase Clusters (2.1M intervals)"
    echo "=============================================="
    echo ""

    # Merge DNase clusters
    run_benchmark "Real: Merge DNase Clusters (2.1M)" \
        "bedtools merge -i $REAL_DATA/dnase_clusters.bed > /dev/null" \
        "$BEDTOOLS_RS merge -i $REAL_DATA/dnase_clusters.bed > /dev/null" \
        2 5
    fi

    if should_run "real" || should_run "intersect"; then
    # Self-intersect DNase (find overlapping peaks)
    run_benchmark "Real: DNase self-intersect (2.1M)" \
        "bedtools intersect -sorted -a $REAL_DATA/dnase_clusters.bed -b $REAL_DATA/dnase_clusters.bed -u > /dev/null" \
        "$BEDTOOLS_RS intersect -a $REAL_DATA/dnase_clusters.bed -b $REAL_DATA/dnase_clusters.bed -u > /dev/null" \
        1 3
    fi
fi

# Real data benchmarks for new commands
REAL_GENOME="$REAL_DATA/hg38.genome"

if [ -f "$REAL_GENOME" ] && [ -f "$REAL_DATA/cpg_islands.bed" ]; then
    if should_run "real" || should_run "slop"; then
    echo ""
    echo "=============================================="
    echo "  REAL DATA: New Commands Benchmarks"
    echo "=============================================="
    echo ""

    # Slop on CpG islands
    run_benchmark "Real: Slop CpG Islands ±500bp" \
        "bedtools slop -i $REAL_DATA/cpg_islands.bed -g $REAL_GENOME -b 500 > /dev/null" \
        "$BEDTOOLS_RS slop -i $REAL_DATA/cpg_islands.bed -g $REAL_GENOME -b 500 > /dev/null"
    fi

    if should_run "real" || should_run "complement"; then
    # Complement of CpG islands
    run_benchmark "Real: Complement CpG Islands" \
        "bedtools complement -i $REAL_DATA/cpg_islands.bed -g $REAL_GENOME > /dev/null" \
        "$BEDTOOLS_RS complement -i $REAL_DATA/cpg_islands.bed -g $REAL_GENOME > /dev/null"
    fi

    # Genomecov on DNase clusters
    if [ -f "$REAL_DATA/dnase_clusters.bed" ]; then
        if should_run "real" || should_run "genomecov"; then
        run_benchmark "Real: Genomecov DNase (histogram)" \
            "bedtools genomecov -i $REAL_DATA/dnase_clusters.bed -g $REAL_GENOME > /dev/null" \
            "$BEDTOOLS_RS genomecov -i $REAL_DATA/dnase_clusters.bed -g $REAL_GENOME > /dev/null" \
            2 5

        run_benchmark "Real: Genomecov DNase (bedgraph)" \
            "bedtools genomecov -i $REAL_DATA/dnase_clusters.bed -g $REAL_GENOME -bg > /dev/null" \
            "$BEDTOOLS_RS genomecov -i $REAL_DATA/dnase_clusters.bed -g $REAL_GENOME -bg > /dev/null" \
            2 5
        fi
    fi

    # Jaccard on real data
    if [ -f "$REAL_DATA/dnase_clusters.bed" ]; then
        if should_run "real" || should_run "jaccard"; then
        run_benchmark "Real: Jaccard CpG vs DNase" \
            "bedtools jaccard -a $REAL_DATA/cpg_islands.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null" \
            "$BEDTOOLS_RS jaccard -a $REAL_DATA/cpg_islands.bed -b $REAL_DATA/dnase_clusters.bed > /dev/null"
        fi
    fi

    # Multiinter on real data
    if [ -f "$REAL_DATA/dnase_clusters.bed" ] && [ -f "$REAL_DATA/promoters.bed" ]; then
        if should_run "real" || should_run "multiinter"; then
        run_benchmark "Real: Multiinter CpG+DNase+Promoters" \
            "bedtools multiinter -i $REAL_DATA/cpg_islands.bed $REAL_DATA/dnase_clusters.bed $REAL_DATA/promoters.bed > /dev/null" \
            "$BEDTOOLS_RS multiinter -i $REAL_DATA/cpg_islands.bed $REAL_DATA/dnase_clusters.bed $REAL_DATA/promoters.bed > /dev/null" \
            2 5
        fi
    fi
fi

echo ""
echo "=============================================="
echo "  BENCHMARK COMPLETE"
echo "=============================================="
echo ""
echo "Results saved to: $BENCHMARK_FILE"
echo ""

# Calculate and display summary
echo "## Summary" >> "$BENCHMARK_FILE"
echo "" >> "$BENCHMARK_FILE"
echo "Benchmark completed at $(date)" >> "$BENCHMARK_FILE"

cat "$BENCHMARK_FILE"
