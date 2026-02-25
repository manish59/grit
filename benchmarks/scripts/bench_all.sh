#!/usr/bin/env bash
set -euo pipefail

#############################################################################
# bench_all.sh - Performance benchmarking harness for GRIT vs bedtools
#
# Usage:
#   ./bench_all.sh <A_size> <B_size> [distribution] [commands...] [--trials N]
#   ./bench_all.sh 10M 5M                        # All commands, uniform, 3 trials
#   ./bench_all.sh 10M 5M clustered              # All commands, clustered
#   ./bench_all.sh 10M 5M all all --trials 5    # All distributions, 5 trials
#
# Output:
#   - CSV: benchmarks/results/bench_<size>_<timestamp>.csv
#   - Console: formatted table with speedup calculations
#############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$SCRIPT_DIR/.."
DATA_DIR="$BENCH_DIR/data"
RESULTS_DIR="$BENCH_DIR/results"
TRUTH_DIR="$BENCH_DIR/truth"

# Source common functions
source "$SCRIPT_DIR/common.sh" 2>/dev/null || true

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

# All supported commands
ALL_COMMANDS="coverage intersect merge subtract closest window jaccard complement slop sort"

# All distributions
ALL_DISTRIBUTIONS="uniform clustered"

# Default trials
TRIALS=3

#############################################################################
# Utilities
#############################################################################

# Detect GRIT binary
detect_grit() {
  if [[ -x "$BENCH_DIR/../target/release/grit" ]]; then
    echo "$BENCH_DIR/../target/release/grit"
  elif command -v grit &>/dev/null; then
    command -v grit
  else
    echo -e "${RED}ERROR: GRIT binary not found.${NC}" >&2
    exit 1
  fi
}

# Platform-aware time command
detect_time_cmd() {
  if [[ "$(uname)" == "Darwin" ]]; then
    echo "/usr/bin/time -l"
  else
    echo "/usr/bin/time -v"
  fi
}

# Extract metrics from time output
# Sets: TIME_SEC, TIME_RSS_MB
extract_metrics() {
  local time_file="$1"

  if [[ "$(uname)" == "Darwin" ]]; then
    TIME_SEC=$(grep "real" "$time_file" 2>/dev/null | awk '{print $1}' || echo "0")
    local rss_bytes=$(grep "maximum resident set size" "$time_file" 2>/dev/null | awk '{print $1}' || echo "0")
    TIME_RSS_MB=$(awk "BEGIN {printf \"%.1f\", $rss_bytes/1048576}")
  else
    TIME_SEC=$(grep "Elapsed" "$time_file" 2>/dev/null | sed 's/.*: //' | \
      awk -F: '{if(NF==3)print $1*3600+$2*60+$3; else if(NF==2)print $1*60+$2; else print $1}' || echo "0")
    local rss_kb=$(grep "Maximum resident set size" "$time_file" 2>/dev/null | awk '{print $NF}' || echo "0")
    TIME_RSS_MB=$(awk "BEGIN {printf \"%.1f\", $rss_kb/1024}")
  fi
}

# Platform-aware SHA256
sha256_hash() {
  if command -v sha256sum &>/dev/null; then
    sha256sum "$@" | awk '{print $1}'
  else
    shasum -a 256 "$@" | awk '{print $1}'
  fi
}

# SHA256 hash from stdin
sha256_hash_stdin() {
  if command -v sha256sum &>/dev/null; then
    sha256sum | awk '{print $1}'
  else
    shasum -a 256 | awk '{print $1}'
  fi
}

# Calculate median of array
calc_median() {
  printf '%s\n' "$@" | sort -n | awk '{a[NR]=$1} END {
    if (NR % 2 == 1) print a[(NR+1)/2]
    else print (a[NR/2] + a[NR/2+1]) / 2
  }'
}

# Clear cache
clear_cache() {
  sync 2>/dev/null || true
  sleep 0.1
}

#############################################################################
# Benchmarking
#############################################################################

run_trials() {
  local cmd_str="$1"
  local trials="$2"
  local time_cmd="$3"

  local times=()
  local rss_values=()

  for i in $(seq 1 "$trials"); do
    local tmp_out=$(mktemp)
    local tmp_time=$(mktemp)

    clear_cache

    # Run command
    eval "$time_cmd $cmd_str" > "$tmp_out" 2> "$tmp_time" || true

    extract_metrics "$tmp_time"
    times+=("$TIME_SEC")
    rss_values+=("$TIME_RSS_MB")

    rm -f "$tmp_out" "$tmp_time"
  done

  MEDIAN_TIME=$(calc_median "${times[@]}")
  MEDIAN_RSS=$(calc_median "${rss_values[@]}")
}

bench_command() {
  local cmd="$1"
  local a_file="$2"
  local b_file="$3"
  local genome="$4"
  local grit="$5"
  local time_cmd="$6"
  local trials="$7"
  local csv_file="$8"
  local distribution="$9"
  local size_spec="${10}"

  # Check for truth data
  local truth_base="$TRUTH_DIR/${size_spec}/${distribution}"
  local use_truth=false
  if [[ -f "${truth_base}/${cmd}.out" && -f "${truth_base}/${cmd}.time" && -f "${truth_base}/${cmd}.sha256" ]]; then
    use_truth=true
  fi

  # Build commands
  local bt_cmd gr_cmd
  case "$cmd" in
    coverage)
      bt_cmd="bedtools coverage -a $a_file -b $b_file -sorted"
      gr_cmd="$grit coverage -a $a_file -b $b_file --assume-sorted"
      ;;
    intersect)
      bt_cmd="bedtools intersect -a $a_file -b $b_file -sorted"
      gr_cmd="$grit intersect -a $a_file -b $b_file --streaming --assume-sorted"
      ;;
    merge)
      bt_cmd="bedtools merge -i $a_file"
      gr_cmd="$grit merge -i $a_file --assume-sorted"
      ;;
    subtract)
      bt_cmd="bedtools subtract -a $a_file -b $b_file -sorted"
      gr_cmd="$grit subtract -a $a_file -b $b_file --streaming --assume-sorted"
      ;;
    closest)
      bt_cmd="bedtools closest -a $a_file -b $b_file"
      gr_cmd="$grit closest -a $a_file -b $b_file --streaming --assume-sorted"
      ;;
    window)
      bt_cmd="bedtools window -a $a_file -b $b_file -w 1000"
      gr_cmd="$grit window -a $a_file -b $b_file -w 1000 --assume-sorted"
      ;;
    jaccard)
      bt_cmd="bedtools jaccard -a $a_file -b $b_file"
      gr_cmd="$grit jaccard -a $a_file -b $b_file"
      ;;
    complement)
      bt_cmd="bedtools complement -i $a_file -g $genome"
      gr_cmd="$grit complement -i $a_file -g $genome --assume-sorted"
      ;;
    slop)
      bt_cmd="bedtools slop -i $a_file -g $genome -b 100"
      gr_cmd="$grit slop -i $a_file -g $genome -b 100"
      ;;
    sort)
      bt_cmd="LC_ALL=C sort -k1,1 -k2,2n -k3,3n $a_file"
      gr_cmd="$grit sort -i $a_file"
      ;;
    *)
      echo "Unknown command: $cmd"
      return 1
      ;;
  esac

  local bt_time bt_rss bt_lines bt_hash
  local cmd_display="$cmd"

  if $use_truth; then
    # Use pre-computed truth data
    bt_time=$(cat "${truth_base}/${cmd}.time")
    local bt_rss_bytes=$(cat "${truth_base}/${cmd}.rss" 2>/dev/null || echo "0")
    bt_rss=$(awk "BEGIN {printf \"%.1f\", $bt_rss_bytes/1048576}")
    bt_lines=$(cat "${truth_base}/${cmd}.lines" 2>/dev/null || wc -l < "${truth_base}/${cmd}.out" | tr -d ' ')
    bt_hash=$(cat "${truth_base}/${cmd}.sha256")
    cmd_display="${cmd}*"
  else
    # Benchmark bedtools
    run_trials "$bt_cmd" "$trials" "$time_cmd"
    bt_time="$MEDIAN_TIME"
    bt_rss="$MEDIAN_RSS"

    # Get bedtools output for correctness check
    local bt_out_tmp=$(mktemp)
    eval "$bt_cmd" > "$bt_out_tmp" 2>/dev/null || true
    bt_lines=$(wc -l < "$bt_out_tmp" | tr -d ' ')

    # Calculate hash
    if [[ "$cmd" == "coverage" ]]; then
      bt_hash=$(cut -f1-6 "$bt_out_tmp" | sha256_hash_stdin 2>/dev/null || shasum -a 256 | awk '{print $1}')
    elif [[ "$cmd" == "window" ]]; then
      bt_hash=$(sort "$bt_out_tmp" | sha256_hash_stdin 2>/dev/null || shasum -a 256 | awk '{print $1}')
    else
      bt_hash=$(sha256_hash "$bt_out_tmp")
    fi
    rm -f "$bt_out_tmp"
  fi

  # Benchmark grit
  run_trials "$gr_cmd" "$trials" "$time_cmd"
  local gr_time="$MEDIAN_TIME"
  local gr_rss="$MEDIAN_RSS"

  # Calculate speedup
  local speedup
  if [[ "$gr_time" != "0" && "$gr_time" != "0.00" ]]; then
    speedup=$(awk "BEGIN {printf \"%.2f\", $bt_time/$gr_time}")
  else
    speedup="inf"
  fi

  # Correctness check using hash
  local gr_out=$(mktemp)
  eval "$gr_cmd" > "$gr_out" 2>/dev/null || true
  local gr_lines=$(wc -l < "$gr_out" | tr -d ' ')

  local gr_hash
  if [[ "$cmd" == "coverage" ]]; then
    gr_hash=$(cut -f1-6 "$gr_out" | sha256_hash_stdin 2>/dev/null || shasum -a 256 | awk '{print $1}')
  elif [[ "$cmd" == "window" ]]; then
    gr_hash=$(sort "$gr_out" | sha256_hash_stdin 2>/dev/null || shasum -a 256 | awk '{print $1}')
  else
    gr_hash=$(sha256_hash "$gr_out")
  fi
  rm -f "$gr_out"

  local correctness="PASS"
  if [[ "$bt_lines" != "$gr_lines" || "$bt_hash" != "$gr_hash" ]]; then
    correctness="FAIL"
  fi

  # Write CSV rows
  local source_indicator=""
  if $use_truth; then
    source_indicator="(truth)"
  fi
  echo "$cmd,$distribution,bedtools,default,$bt_time,$bt_rss,1.00,$correctness,$source_indicator" >> "$csv_file"
  echo "$cmd,$distribution,grit,streaming,$gr_time,$gr_rss,$speedup,$correctness," >> "$csv_file"

  # Color output
  local speedup_colored
  if [[ "$speedup" == "inf" ]]; then
    speedup_colored="${GREEN}inf${NC}"
  elif (( $(echo "$speedup >= 5" | bc -l) )); then
    speedup_colored="${GREEN}${speedup}x${NC}"
  elif (( $(echo "$speedup >= 2" | bc -l) )); then
    speedup_colored="${CYAN}${speedup}x${NC}"
  elif (( $(echo "$speedup >= 1" | bc -l) )); then
    speedup_colored="${YELLOW}${speedup}x${NC}"
  else
    speedup_colored="${RED}${speedup}x${NC}"
  fi

  local correct_colored
  if [[ "$correctness" == "PASS" ]]; then
    correct_colored="${GREEN}PASS${NC}"
  else
    correct_colored="${RED}FAIL${NC}"
  fi

  printf "  %-12s %10.2fs %10.2fs   %b   %8s MB   %b\n" \
    "$cmd_display" "$bt_time" "$gr_time" "$speedup_colored" "$gr_rss" "$correct_colored"
}

bench_distribution() {
  local a_size="$1"
  local b_size="$2"
  local distribution="$3"
  local commands="$4"
  local grit="$5"
  local time_cmd="$6"
  local trials="$7"
  local csv_file="$8"

  local size_spec="${a_size}_${b_size}"
  local data_path="$DATA_DIR/${size_spec}/${distribution}"
  local a_file="$data_path/A.sorted.bed"
  local b_file="$data_path/B.sorted.bed"
  local genome="$data_path/genome.txt"

  # Check data exists
  if [[ ! -f "$a_file" || ! -f "$b_file" ]]; then
    echo -e "${YELLOW}Data not found for $distribution. Run: ./bench.sh data $a_size $b_size $distribution${NC}"
    return 1
  fi

  local a_lines=$(wc -l < "$a_file" | tr -d ' ')
  local b_lines=$(wc -l < "$b_file" | tr -d ' ')

  # Check for truth data
  local truth_base="$TRUTH_DIR/${size_spec}/${distribution}"
  local truth_status="${YELLOW}not available${NC}"
  if [[ -d "$truth_base" ]]; then
    truth_status="${GREEN}available${NC} (* = using cached bedtools)"
  fi

  echo -e "${BOLD}Distribution: $distribution${NC}"
  echo "  Data: A=$a_lines intervals, B=$b_lines intervals"
  echo -e "  Truth data: $truth_status"
  echo ""
  printf "  %-12s %11s %11s %11s %12s   %s\n" \
    "Command" "bedtools" "GRIT" "Speedup" "GRIT_RSS" "Correct"
  echo "  ─────────────────────────────────────────────────────────────────────────"

  for cmd in $commands; do
    bench_command "$cmd" "$a_file" "$b_file" "$genome" "$grit" "$time_cmd" "$trials" "$csv_file" "$distribution" "$size_spec"
  done

  echo ""
}

#############################################################################
# Main
#############################################################################

main() {
  # Parse arguments
  local a_size=""
  local b_size=""
  local distribution="uniform"
  local commands=""
  local trials=$TRIALS

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --trials)
        trials="$2"
        shift 2
        ;;
      *)
        if [[ -z "$a_size" ]]; then
          a_size="$1"
        elif [[ -z "$b_size" ]]; then
          b_size="$1"
        elif [[ "$1" == "uniform" || "$1" == "clustered" || "$1" == "all" ]]; then
          distribution="$1"
        else
          commands="$commands $1"
        fi
        shift
        ;;
    esac
  done

  if [[ -z "$a_size" || -z "$b_size" ]]; then
    echo "Usage: $0 <A_size> <B_size> [distribution] [commands...] [--trials N]"
    echo "Example: $0 10M 5M"
    echo "         $0 10M 5M clustered --trials 5"
    echo "         $0 10M 5M all all"
    exit 1
  fi

  # Handle 'all' keywords
  if [[ "$distribution" == "all" ]]; then
    distribution="$ALL_DISTRIBUTIONS"
  fi
  commands="${commands:-$ALL_COMMANDS}"
  if [[ "$commands" == *"all"* ]]; then
    commands="$ALL_COMMANDS"
  fi
  commands=$(echo "$commands" | xargs)  # trim

  local grit=$(detect_grit)
  local time_cmd=$(detect_time_cmd)

  mkdir -p "$RESULTS_DIR"
  local timestamp=$(date +%Y%m%d_%H%M%S)
  local csv_file="$RESULTS_DIR/bench_${a_size}_${b_size}_${timestamp}.csv"

  echo -e "${BOLD}=== GRIT Performance Benchmark ===${NC}"
  echo "Size: ${a_size} A × ${b_size} B"
  echo "Distributions: $distribution"
  echo "Commands: $commands"
  echo "Trials: $trials"
  echo "Output: $csv_file"
  echo ""

  # Write CSV header
  echo "command,distribution,tool,mode,time_s,rss_mb,speedup_vs_bedtools,correctness,source" > "$csv_file"

  for dist in $distribution; do
    bench_distribution "$a_size" "$b_size" "$dist" "$commands" "$grit" "$time_cmd" "$trials" "$csv_file"
  done

  echo "═══════════════════════════════════════════════════════════════════════════"
  echo "Results saved: $csv_file"

  # Summary
  echo ""
  echo -e "${BOLD}Summary:${NC}"
  local total=$(grep -c "grit" "$csv_file" || echo 0)
  local fast=$(awk -F, '$3 == "grit" && $7 >= 2 {count++} END {print count+0}' "$csv_file")
  local target=$(awk -F, '$3 == "grit" && $7 >= 5 {count++} END {print count+0}' "$csv_file")
  echo "  Total benchmarks: $total"
  echo -e "  ≥2x speedup: ${GREEN}$fast${NC}/$total"
  echo -e "  ≥5x speedup: ${GREEN}$target${NC}/$total"
}

main "$@"
