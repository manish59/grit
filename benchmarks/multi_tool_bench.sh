#!/usr/bin/env bash
#############################################################################
# GRIT Multi-Tool Benchmark Suite
#
# Comprehensive benchmarking against:
#   - bedtools (C++, the standard)
#   - bedops (C++, memory-efficient)
#   - granges (Rust)
#   - pyranges (Python)
#   - polars-bio (Python/Rust)
#
# Usage:
#   ./multi_tool_bench.sh status           Check installed tools
#   ./multi_tool_bench.sh run <size>       Run benchmarks (1M, 10M, etc.)
#   ./multi_tool_bench.sh compare <size>   Generate comparison tables
#
#############################################################################

set -euo pipefail

BENCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$BENCH_DIR/data"
RESULTS_DIR="$BENCH_DIR/results/multi_tool"
VENV_DIR="$BENCH_DIR/.venv"

source "$BENCH_DIR/scripts/install_tools.sh" 2>/dev/null || true

# Restore our directory reference after sourcing
SCRIPT_DIR="$BENCH_DIR"

# Colors
if [[ -t 1 ]]; then
  GREEN='\033[0;32m'
  RED='\033[0;31m'
  YELLOW='\033[0;33m'
  CYAN='\033[0;36m'
  BOLD='\033[1m'
  NC='\033[0m'
else
  GREEN='' RED='' YELLOW='' CYAN='' BOLD='' NC=''
fi

# Supported operations (bash 3 compatible - no associative arrays)
ALL_OPERATIONS="intersect merge subtract closest coverage sort"

# Operations that need two files
TWO_FILE_OPS="intersect subtract closest coverage"

# Default benchmark sizes
DEFAULT_SIZES="100K 1M 10M"

#############################################################################
# Utility Functions
#############################################################################

parse_size() {
  local size="$1"
  local num="${size%[KkMm]}"
  local suffix="${size: -1}"

  case "$suffix" in
    K|k) echo $((num * 1000)) ;;
    M|m) echo $((num * 1000000)) ;;
    *)   echo "$size" ;;
  esac
}

format_time() {
  local seconds="$1"
  if (( $(echo "$seconds < 0.01" | bc -l) )); then
    printf "%.4fs" "$seconds"
  elif (( $(echo "$seconds < 1" | bc -l) )); then
    printf "%.3fs" "$seconds"
  elif (( $(echo "$seconds < 60" | bc -l) )); then
    printf "%.2fs" "$seconds"
  else
    local mins=$(echo "$seconds / 60" | bc)
    local secs=$(echo "$seconds - $mins * 60" | bc)
    printf "%dm%.1fs" "$mins" "$secs"
  fi
}

format_memory() {
  local mb="$1"
  if (( $(echo "$mb < 1" | bc -l) )); then
    printf "%.2f MB" "$mb"
  elif (( $(echo "$mb < 1024" | bc -l) )); then
    printf "%.1f MB" "$mb"
  else
    printf "%.2f GB" "$(echo "$mb / 1024" | bc -l)"
  fi
}

detect_time_cmd() {
  if [[ "$(uname)" == "Darwin" ]]; then
    echo "/usr/bin/time -l"
  else
    echo "/usr/bin/time -v"
  fi
}

extract_metrics() {
  local time_file="$1"

  if [[ "$(uname)" == "Darwin" ]]; then
    TIME_SEC=$(grep "real" "$time_file" 2>/dev/null | awk '{print $1}' || echo "0")
    local rss_bytes=$(grep "maximum resident set size" "$time_file" 2>/dev/null | awk '{print $1}' || echo "0")
    TIME_RSS_MB=$(echo "scale=2; $rss_bytes / 1048576" | bc)
  else
    TIME_SEC=$(grep "Elapsed" "$time_file" 2>/dev/null | sed 's/.*: //' | \
      awk -F: '{if(NF==3)print $1*3600+$2*60+$3; else if(NF==2)print $1*60+$2; else print $1}' || echo "0")
    local rss_kb=$(grep "Maximum resident set size" "$time_file" 2>/dev/null | awk '{print $NF}' || echo "0")
    TIME_RSS_MB=$(echo "scale=2; $rss_kb / 1024" | bc)
  fi
}

#############################################################################
# Data Generation
#############################################################################

generate_data() {
  local a_size="$1"
  local b_size="${2:-$(( $(parse_size "$a_size") / 2 ))}"

  local a_count=$(parse_size "$a_size")
  local b_count=$(parse_size "$b_size")

  local data_subdir="$DATA_DIR/${a_size}_${b_size}"
  mkdir -p "$data_subdir"

  local genome_file="$data_subdir/genome.txt"
  local file_a="$data_subdir/a.bed"
  local file_b="$data_subdir/b.bed"
  local file_a_sorted="$data_subdir/a.sorted.bed"
  local file_b_sorted="$data_subdir/b.sorted.bed"

  # Create genome file if needed
  if [[ ! -f "$genome_file" ]]; then
    cat > "$genome_file" <<'EOF'
chr1	248956422
chr2	242193529
chr3	198295559
chr4	190214555
chr5	181538259
chr6	170805979
chr7	159345973
chr8	145138636
chr9	138394717
chr10	133797422
chr11	135086622
chr12	133275309
chr13	114364328
chr14	107043718
chr15	101991189
chr16	90338345
chr17	83257441
chr18	80373285
chr19	58617616
chr20	64444167
chr21	46709983
chr22	50818468
chrX	156040895
chrY	57227415
EOF
  fi

  # Generate A file (messages to stderr so path is only stdout)
  if [[ ! -f "$file_a" ]]; then
    echo -e "${CYAN}Generating $a_size intervals for A...${NC}" >&2
    bedtools random -l 150 -n "$a_count" -seed 42 -g "$genome_file" > "$file_a"
  fi

  # Generate B file
  if [[ ! -f "$file_b" ]]; then
    echo -e "${CYAN}Generating $b_size intervals for B...${NC}" >&2
    bedtools random -l 500 -n "$b_count" -seed 123 -g "$genome_file" > "$file_b"
  fi

  # Sort files (use lowercase to match existing convention)
  if [[ ! -f "$file_a_sorted" ]] && [[ ! -f "$data_subdir/A.sorted.bed" ]]; then
    echo -e "${CYAN}Sorting A...${NC}" >&2
    LC_ALL=C sort -k1,1 -k2,2n "$file_a" > "$file_a_sorted"
  fi

  if [[ ! -f "$file_b_sorted" ]] && [[ ! -f "$data_subdir/B.sorted.bed" ]]; then
    echo -e "${CYAN}Sorting B...${NC}" >&2
    LC_ALL=C sort -k1,1 -k2,2n "$file_b" > "$file_b_sorted"
  fi

  echo "$data_subdir"
}

#############################################################################
# Tool Runners
#############################################################################

run_grit() {
  local op="$1"
  local file_a="$2"
  local file_b="$3"
  local output="$4"
  local time_file="$5"

  local grit_bin
  grit_bin=$(detect_grit 2>/dev/null) || return 1

  local TIME_CMD=$(detect_time_cmd)

  # Use streaming/assume-sorted flags where supported for O(k) memory
  case "$op" in
    intersect)
      $TIME_CMD "$grit_bin" intersect -a "$file_a" -b "$file_b" --streaming --assume-sorted > "$output" 2> "$time_file"
      ;;
    merge)
      # merge doesn't have --streaming, uses streaming by default with --assume-sorted
      $TIME_CMD "$grit_bin" merge -i "$file_a" --assume-sorted > "$output" 2> "$time_file"
      ;;
    subtract)
      $TIME_CMD "$grit_bin" subtract -a "$file_a" -b "$file_b" --streaming --assume-sorted > "$output" 2> "$time_file"
      ;;
    closest)
      $TIME_CMD "$grit_bin" closest -a "$file_a" -b "$file_b" --streaming --assume-sorted > "$output" 2> "$time_file"
      ;;
    coverage)
      $TIME_CMD "$grit_bin" coverage -a "$file_a" -b "$file_b" --streaming --assume-sorted > "$output" 2> "$time_file"
      ;;
    sort)
      $TIME_CMD "$grit_bin" sort -i "$file_a" > "$output" 2> "$time_file"
      ;;
    window)
      $TIME_CMD "$grit_bin" window -a "$file_a" -b "$file_b" --streaming --assume-sorted > "$output" 2> "$time_file"
      ;;
    *)
      echo "Unknown operation: $op" >&2
      return 1
      ;;
  esac
}

run_bedtools() {
  local op="$1"
  local file_a="$2"
  local file_b="$3"
  local output="$4"
  local time_file="$5"

  command -v bedtools &>/dev/null || return 1

  local TIME_CMD=$(detect_time_cmd)

  case "$op" in
    intersect)
      $TIME_CMD bedtools intersect -a "$file_a" -b "$file_b" -sorted > "$output" 2> "$time_file"
      ;;
    merge)
      $TIME_CMD bedtools merge -i "$file_a" > "$output" 2> "$time_file"
      ;;
    subtract)
      $TIME_CMD bedtools subtract -a "$file_a" -b "$file_b" -sorted > "$output" 2> "$time_file"
      ;;
    closest)
      $TIME_CMD bedtools closest -a "$file_a" -b "$file_b" > "$output" 2> "$time_file"
      ;;
    coverage)
      $TIME_CMD bedtools coverage -a "$file_a" -b "$file_b" -sorted > "$output" 2> "$time_file"
      ;;
    sort)
      $TIME_CMD bedtools sort -i "$file_a" > "$output" 2> "$time_file"
      ;;
    *)
      return 1
      ;;
  esac
}

run_bedops() {
  local op="$1"
  local file_a="$2"
  local file_b="$3"
  local output="$4"
  local time_file="$5"

  command -v bedops &>/dev/null || return 1

  local TIME_CMD=$(detect_time_cmd)

  case "$op" in
    intersect)
      $TIME_CMD bedops --intersect "$file_a" "$file_b" > "$output" 2> "$time_file"
      ;;
    merge)
      $TIME_CMD bedops --merge "$file_a" > "$output" 2> "$time_file"
      ;;
    subtract)
      $TIME_CMD bedops --difference "$file_a" "$file_b" > "$output" 2> "$time_file"
      ;;
    closest)
      $TIME_CMD closest-features --closest "$file_a" "$file_b" > "$output" 2> "$time_file"
      ;;
    sort)
      $TIME_CMD sort-bed "$file_a" > "$output" 2> "$time_file"
      ;;
    *)
      return 1
      ;;
  esac
}

run_granges() {
  local op="$1"
  local file_a="$2"
  local file_b="$3"
  local output="$4"
  local time_file="$5"

  command -v granges &>/dev/null || return 1

  local TIME_CMD=$(detect_time_cmd)

  # Derive genome file path from input file directory
  local data_dir=$(dirname "$file_a")
  local genome_file="$data_dir/genome.txt"
  [[ -f "$genome_file" ]] || return 1

  case "$op" in
    intersect)
      $TIME_CMD granges filter --genome "$genome_file" --left "$file_a" --right "$file_b" -o "$output" 2> "$time_file"
      ;;
    *)
      # granges has limited operations
      return 1
      ;;
  esac
}

run_pyranges() {
  local op="$1"
  local file_a="$2"
  local file_b="$3"
  local output="$4"
  local time_file="$5"

  [[ -f "$VENV_DIR/bin/python" ]] || return 1
  "$VENV_DIR/bin/python" -c "import pyranges" &>/dev/null || return 1

  local TIME_CMD=$(detect_time_cmd)
  local py_script="$SCRIPT_DIR/scripts/python_benchmark.py"

  # Use system time for consistent timing across all tools
  if [[ " $TWO_FILE_OPS " == *" $op "* ]]; then
    $TIME_CMD "$VENV_DIR/bin/python" "$py_script" pyranges "$op" "$file_a" "$file_b" -o "$output" --csv > /dev/null 2> "$time_file"
  else
    $TIME_CMD "$VENV_DIR/bin/python" "$py_script" pyranges "$op" "$file_a" -o "$output" --csv > /dev/null 2> "$time_file"
  fi
}

run_polars_bio() {
  local op="$1"
  local file_a="$2"
  local file_b="$3"
  local output="$4"
  local time_file="$5"

  [[ -f "$VENV_DIR/bin/python" ]] || return 1
  "$VENV_DIR/bin/python" -c "import polars_bio" &>/dev/null || return 1

  local TIME_CMD=$(detect_time_cmd)
  local py_script="$SCRIPT_DIR/scripts/python_benchmark.py"

  # Use system time for consistent timing across all tools
  if [[ " $TWO_FILE_OPS " == *" $op "* ]]; then
    $TIME_CMD "$VENV_DIR/bin/python" "$py_script" polars-bio "$op" "$file_a" "$file_b" -o "$output" --csv > /dev/null 2> "$time_file"
  else
    $TIME_CMD "$VENV_DIR/bin/python" "$py_script" polars-bio "$op" "$file_a" -o "$output" --csv > /dev/null 2> "$time_file"
  fi
}

#############################################################################
# Benchmark Runner
#############################################################################

run_single_benchmark() {
  local tool="$1"
  local op="$2"
  local file_a="$3"
  local file_b="$4"
  local results_dir="$5"

  local output="$results_dir/${tool}_${op}.out"
  local time_file="$results_dir/${tool}_${op}.time"

  local run_func="run_$tool"

  # Run benchmark
  if $run_func "$op" "$file_a" "$file_b" "$output" "$time_file" 2>/dev/null; then
    extract_metrics "$time_file"
    echo "$TIME_SEC,$TIME_RSS_MB"
  else
    echo "NA,NA"
  fi
}

run_benchmarks() {
  local a_size="$1"
  local b_size="${2:-}"
  local operations="${3:-intersect merge subtract closest}"

  # Calculate b_size if not provided
  if [[ -z "$b_size" ]]; then
    local a_count=$(parse_size "$a_size")
    b_size=$(echo "$a_count / 2" | bc)
    if (( b_size >= 1000000 )); then
      b_size="$((b_size / 1000000))M"
    elif (( b_size >= 1000 )); then
      b_size="$((b_size / 1000))K"
    fi
  fi

  echo -e "${BOLD}${CYAN}=== Multi-Tool Benchmark: ${a_size} × ${b_size} ===${NC}"
  echo ""

  # Generate data
  local data_dir
  data_dir=$(generate_data "$a_size" "$b_size")

  # Support both naming conventions (lowercase and uppercase)
  local file_a file_b
  if [[ -f "$data_dir/a.sorted.bed" ]]; then
    file_a="$data_dir/a.sorted.bed"
    file_b="$data_dir/b.sorted.bed"
  elif [[ -f "$data_dir/A.sorted.bed" ]]; then
    file_a="$data_dir/A.sorted.bed"
    file_b="$data_dir/B.sorted.bed"
  else
    echo -e "${RED}Error: Sorted BED files not found in $data_dir${NC}"
    return 1
  fi

  # Create results directory
  local timestamp=$(date +%Y%m%d_%H%M%S)
  local run_results_dir="$RESULTS_DIR/${a_size}_${b_size}_${timestamp}"
  mkdir -p "$run_results_dir"

  # Detect available tools
  local available_tools=""
  detect_grit &>/dev/null && available_tools="$available_tools grit"
  detect_bedtools &>/dev/null && available_tools="$available_tools bedtools"
  detect_bedops &>/dev/null && available_tools="$available_tools bedops"
  detect_granges &>/dev/null && available_tools="$available_tools granges"
  detect_pyranges && available_tools="$available_tools pyranges"
  detect_polars_bio && available_tools="$available_tools polars_bio"

  echo -e "${CYAN}Available tools:${NC} $available_tools"
  echo ""

  # CSV header
  local csv_file="$run_results_dir/results.csv"
  echo "operation,tool,time_s,memory_mb,speedup_vs_bedtools" > "$csv_file"

  # Run benchmarks for each operation
  for op in $operations; do
    echo -e "${BOLD}Operation: $op${NC}"

    # Temporary file to store results for this operation
    local op_results_file="$run_results_dir/${op}_results.tmp"
    > "$op_results_file"

    # Run benchmarks and collect results
    for tool in $available_tools; do
      local result
      result=$(run_single_benchmark "$tool" "$op" "$file_a" "$file_b" "$run_results_dir")

      local time_s=$(echo "$result" | cut -d, -f1)
      local mem_mb=$(echo "$result" | cut -d, -f2)

      # Store result
      echo "$tool,$time_s,$mem_mb" >> "$op_results_file"

      if [[ "$time_s" != "NA" ]]; then
        printf "  %-12s %s  (%.1f MB)\n" "$tool:" "$(format_time "$time_s")" "$mem_mb"
      else
        printf "  %-12s %s\n" "$tool:" "N/A"
      fi
    done

    # Get bedtools time for speedup calculation
    local bedtools_time="NA"
    if grep -q "^bedtools," "$op_results_file"; then
      bedtools_time=$(grep "^bedtools," "$op_results_file" | cut -d, -f2)
    fi

    # Calculate speedups and write to CSV
    while IFS=, read -r tool time_s mem_mb; do
      local speedup="NA"

      if [[ "$time_s" != "NA" && "$bedtools_time" != "NA" ]]; then
        # Use awk for safer division (handles zero/small numbers)
        speedup=$(awk -v bt="$bedtools_time" -v tt="$time_s" 'BEGIN {
          if (tt > 0) printf "%.2f", bt / tt; else print "NA"
        }' 2>/dev/null || echo "NA")
      fi

      echo "$op,$tool,$time_s,$mem_mb,$speedup" >> "$csv_file"
    done < "$op_results_file"

    rm -f "$op_results_file"
    echo ""
  done

  echo -e "${GREEN}Results saved to: $csv_file${NC}"
  return 0
}

#############################################################################
# Comparison Table Generation
#############################################################################

generate_comparison_table() {
  local csv_file="$1"

  echo ""
  echo -e "${BOLD}=== Performance Comparison ===${NC}"
  echo ""

  # Header
  printf "%-12s" "Operation"
  for tool in grit bedtools bedops granges pyranges polars_bio; do
    printf " | %-12s" "$tool"
  done
  echo ""
  printf "%s" "------------"
  for i in {1..6}; do
    printf "%s" "-+--------------"
  done
  echo ""

  # Parse CSV and print table
  local current_op=""
  while IFS=, read -r op tool time_s mem_mb speedup; do
    [[ "$op" == "operation" ]] && continue  # Skip header

    if [[ "$op" != "$current_op" ]]; then
      if [[ -n "$current_op" ]]; then
        echo ""
      fi
      current_op="$op"
      printf "%-12s" "$op"
    fi

    if [[ "$time_s" != "NA" ]]; then
      printf " | %7s" "$(format_time "$time_s")"
    else
      printf " | %12s" "-"
    fi
  done < "$csv_file"
  echo ""
}

generate_markdown_table() {
  local csv_file="$1"
  local output_file="$2"

  cat > "$output_file" <<'EOF'
# Multi-Tool Benchmark Results

## Performance Comparison

EOF

  echo "| Operation | Tool | Time | Memory | Speedup vs bedtools |" >> "$output_file"
  echo "|-----------|------|------|--------|---------------------|" >> "$output_file"

  while IFS=, read -r op tool time_s mem_mb speedup; do
    [[ "$op" == "operation" ]] && continue
    if [[ "$time_s" != "NA" ]]; then
      echo "| $op | $tool | ${time_s}s | ${mem_mb} MB | ${speedup}x |" >> "$output_file"
    else
      echo "| $op | $tool | - | - | - |" >> "$output_file"
    fi
  done < "$csv_file"

  echo "" >> "$output_file"
  echo "_Generated on $(date)_" >> "$output_file"

  echo -e "${GREEN}Markdown table saved to: $output_file${NC}"
}

#############################################################################
# Main
#############################################################################

usage() {
  cat <<EOF
GRIT Multi-Tool Benchmark Suite

Compare GRIT against bedtools, bedops, granges, pyranges, and polars-bio.

Usage:
  $0 status                         Show installed tools
  $0 install <tool>                 Install a benchmark tool
  $0 data <a_size> [b_size]         Generate benchmark data
  $0 run <a_size> [b_size] [ops]    Run benchmarks
  $0 quick                          Quick benchmark (100K intervals)
  $0 full                           Full benchmark suite (100K, 1M, 10M)

Tools: bedops, granges, pyranges, polars-bio

Examples:
  $0 status                         # Check what's installed
  $0 install bedops                 # Install bedops
  $0 run 1M                         # Benchmark with 1M intervals
  $0 run 10M 5M intersect merge     # Specific operations
  $0 quick                          # Quick sanity check

Size notation: 100K, 1M, 10M, 100M
Operations: intersect, merge, subtract, closest, coverage, sort
EOF
}

main() {
  local cmd="${1:-help}"
  shift || true

  case "$cmd" in
    status)
      print_status
      ;;
    install)
      local tool="${1:-}"
      case "$tool" in
        bedops) install_bedops ;;
        granges) install_granges ;;
        pyranges) install_pyranges ;;
        polars-bio|polars_bio) install_polars_bio ;;
        python|all-python) install_all_python ;;
        "") echo "Usage: $0 install <tool>"; exit 1 ;;
        *) echo "Unknown tool: $tool"; exit 1 ;;
      esac
      ;;
    data)
      local a_size="${1:-1M}"
      local b_size="${2:-}"
      generate_data "$a_size" "$b_size"
      ;;
    run)
      local a_size="${1:-1M}"
      local b_size="${2:-}"
      shift 2 2>/dev/null || shift 1 2>/dev/null || true
      local ops="${*:-intersect merge subtract closest}"
      run_benchmarks "$a_size" "$b_size" "$ops"
      ;;
    quick)
      run_benchmarks "100K" "50K" "intersect merge"
      ;;
    full)
      for size in 100K 1M 10M; do
        run_benchmarks "$size"
        echo ""
      done
      ;;
    compare)
      local results_dir="${1:-$RESULTS_DIR}"
      local latest=$(ls -t "$results_dir"/*/results.csv 2>/dev/null | head -1)
      if [[ -n "$latest" ]]; then
        generate_comparison_table "$latest"
        generate_markdown_table "$latest" "${latest%.csv}.md"
      else
        echo "No benchmark results found in $results_dir"
        exit 1
      fi
      ;;
    help|--help|-h)
      usage
      ;;
    *)
      usage
      exit 1
      ;;
  esac
}

main "$@"
