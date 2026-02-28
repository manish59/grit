#!/usr/bin/env bash
# Rigorous benchmark: GRIT vs bedops
# Runs each operation N times, reports median

set -euo pipefail

RUNS=5
DATA_DIR="benchmarks/data/100M_10M"
GRIT_BIN="./target/release/grit"

echo "=============================================="
echo "GRIT vs bedops Benchmark"
echo "=============================================="
echo ""
echo "Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Runs per test: $RUNS"
echo ""

echo "## System"
echo ""
echo "- OS: $(uname -s) $(uname -r)"
echo "- CPU: $(sysctl -n machdep.cpu.brand_string 2>/dev/null || lscpu | grep 'Model name' | cut -d: -f2 | xargs)"
echo "- Cores: $(sysctl -n hw.ncpu 2>/dev/null || nproc)"
echo "- Memory: $(( $(sysctl -n hw.memsize 2>/dev/null || grep MemTotal /proc/meminfo | awk '{print $2 * 1024}') / 1024 / 1024 / 1024 )) GB"
echo ""

echo "## Versions"
echo ""
echo "- grit: $($GRIT_BIN --version 2>&1 | head -1)"
echo "- bedops: $(bedops --version 2>&1 | head -1)"
echo ""

echo "## Data"
echo ""
echo "- File A: $(wc -l < "$DATA_DIR/a.sorted.bed" | xargs) intervals"
echo "- File B: $(wc -l < "$DATA_DIR/b.sorted.bed" | xargs) intervals"
echo ""

echo "## Methodology"
echo ""
echo "- Single-threaded execution"
echo "- Same sorted input files for all tools"
echo "- Output redirected to /dev/null"
echo "- Filesystem cache warmed before timing"
echo ""

median() {
    local arr=("$@")
    local n=${#arr[@]}
    IFS=$'\n' sorted=($(sort -n <<<"${arr[*]}")); unset IFS
    local mid=$((n / 2))
    if (( n % 2 == 1 )); then
        echo "${sorted[$mid]}"
    else
        echo "scale=2; (${sorted[$mid-1]} + ${sorted[$mid]}) / 2" | bc
    fi
}

run_bench() {
    local name="$1"
    shift
    local cmd=("$@")
    
    local times=()
    local mems=()
    
    for i in $(seq 1 $RUNS); do
        [[ $i -eq 1 ]] && "${cmd[@]}" > /dev/null 2>&1 || true  # warmup
        
        result=$({ /usr/bin/time -l "${cmd[@]}" > /dev/null; } 2>&1)
        t=$(echo "$result" | grep "real" | awk '{print $1}')
        m=$(echo "$result" | grep "maximum resident" | awk '{print $1}')
        times+=("$t")
        mems+=("$m")
    done
    
    local med_time=$(median "${times[@]}")
    local med_mem=$(median "${mems[@]}")
    local med_mem_mb=$(echo "scale=1; $med_mem / 1024 / 1024" | bc)
    
    IFS=$'\n' sorted_times=($(sort -n <<<"${times[*]}")); unset IFS
    local min_time="${sorted_times[0]}"
    local max_time="${sorted_times[$((RUNS-1))]}"
    
    printf "| %-8s | %7.2fs | %5.2fs - %5.2fs | %6.1f MB |\n" \
        "$name" "$med_time" "$min_time" "$max_time" "$med_mem_mb"
}

echo "## intersect"
echo ""
echo "| Tool     |  Median |     Range     | Memory |"
echo "|----------|---------|---------------|--------|"
run_bench "grit" $GRIT_BIN intersect -a "$DATA_DIR/a.sorted.bed" -b "$DATA_DIR/b.sorted.bed" --streaming --assume-sorted
run_bench "bedops" bedops --intersect "$DATA_DIR/a.sorted.bed" "$DATA_DIR/b.sorted.bed"

echo ""
echo "## subtract"
echo ""
echo "| Tool     |  Median |     Range     | Memory |"
echo "|----------|---------|---------------|--------|"
run_bench "grit" $GRIT_BIN subtract -a "$DATA_DIR/a.sorted.bed" -b "$DATA_DIR/b.sorted.bed" --streaming --assume-sorted
run_bench "bedops" bedops --difference "$DATA_DIR/a.sorted.bed" "$DATA_DIR/b.sorted.bed"

echo ""
echo "## closest"
echo ""
echo "| Tool     |  Median |     Range     | Memory |"
echo "|----------|---------|---------------|--------|"
run_bench "grit" $GRIT_BIN closest -a "$DATA_DIR/a.sorted.bed" -b "$DATA_DIR/b.sorted.bed" --streaming --assume-sorted
run_bench "bedops" closest-features --closest "$DATA_DIR/a.sorted.bed" "$DATA_DIR/b.sorted.bed"
