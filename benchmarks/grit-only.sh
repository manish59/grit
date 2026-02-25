#!/usr/bin/env bash
set -euo pipefail

#############################################################################
# GRIT-Only Benchmark
#
# Benchmarks GRIT performance without bedtools comparison.
# Uses existing sorted data for fast iteration.
#
# Usage:
#   ./grit-only.sh <A_size> <B_size> [commands...]
#
# Examples:
#   ./grit-only.sh 10M 5M                    # Run all 10 commands
#   ./grit-only.sh 100M 50M coverage merge   # Run specific commands
#   ./grit-only.sh 50M 25M all               # Explicit all
#
# Available commands:
#   coverage, intersect, merge, subtract, closest,
#   window, jaccard, complement, slop, sort
#############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data"

# Colors
if [[ -t 1 ]]; then
  GREEN='\033[0;32m'
  CYAN='\033[0;36m'
  BOLD='\033[1m'
  NC='\033[0m'
else
  GREEN='' CYAN='' BOLD='' NC=''
fi

# All 10 commands
ALL_COMMANDS="coverage intersect merge subtract closest window jaccard complement slop sort"

usage() {
  cat <<EOF
GRIT-Only Benchmark

Usage: $0 <A_size> <B_size> [commands...]

Examples:
  $0 10M 5M                    # Run all 10 commands
  $0 100M 50M coverage merge   # Run specific commands

Available commands:
  coverage, intersect, merge, subtract, closest,
  window, jaccard, complement, slop, sort

Available data sizes:
EOF
  for dir in "$DATA_DIR"/*/; do
    if [[ -d "$dir" ]]; then
      basename "$dir"
    fi
  done
  exit 1
}

# Format bytes to human readable
format_bytes() {
  local bytes="$1"
  if (( bytes >= 1073741824 )); then
    awk "BEGIN {printf \"%.1f GB\", $bytes/1073741824}"
  elif (( bytes >= 1048576 )); then
    awk "BEGIN {printf \"%.1f MB\", $bytes/1048576}"
  elif (( bytes >= 1024 )); then
    awk "BEGIN {printf \"%.1f KB\", $bytes/1024}"
  else
    echo "${bytes} B"
  fi
}

# Detect GRIT binary
detect_grit() {
  if [[ -x "$SCRIPT_DIR/../target/release/grit" ]]; then
    echo "$SCRIPT_DIR/../target/release/grit"
  elif command -v grit &>/dev/null; then
    command -v grit
  else
    echo "ERROR: GRIT binary not found. Build with: cargo build --release" >&2
    exit 1
  fi
}

# Platform-aware time command
get_time_cmd() {
  if [[ "$(uname)" == "Darwin" ]]; then
    echo "/usr/bin/time -l"
  else
    echo "/usr/bin/time -v"
  fi
}

# Extract metrics from time output
extract_metrics() {
  local time_file="$1"

  if [[ "$(uname)" == "Darwin" ]]; then
    TIME_SEC=$(grep "real" "$time_file" 2>/dev/null | awk '{print $1}' || echo "0")
    TIME_RSS=$(grep "maximum resident set size" "$time_file" 2>/dev/null | awk '{print $1}' || echo "0")
  else
    TIME_SEC=$(grep "Elapsed" "$time_file" 2>/dev/null | sed 's/.*: //' | \
      awk -F: '{if(NF==3)print $1*3600+$2*60+$3; else if(NF==2)print $1*60+$2; else print $1}' || echo "0")
    local rss_kb=$(grep "Maximum resident set size" "$time_file" 2>/dev/null | awk '{print $NF}' || echo "0")
    TIME_RSS=$((rss_kb * 1024))
  fi
}

# Run single GRIT benchmark
run_grit_benchmark() {
  local cmd="$1"
  local a_file="$2"
  local b_file="$3"
  local genome="$4"
  local grit="$5"
  local time_cmd="$6"

  local tmp_out=$(mktemp)
  local tmp_time=$(mktemp)

  # Build GRIT command
  local grit_cmd
  case "$cmd" in
    coverage)
      grit_cmd="$grit coverage -a $a_file -b $b_file --assume-sorted"
      ;;
    intersect)
      grit_cmd="$grit intersect -a $a_file -b $b_file --assume-sorted --streaming"
      ;;
    merge)
      grit_cmd="$grit merge -i $a_file --assume-sorted"
      ;;
    subtract)
      grit_cmd="$grit subtract -a $a_file -b $b_file --assume-sorted --streaming"
      ;;
    closest)
      grit_cmd="$grit closest -a $a_file -b $b_file --assume-sorted --streaming"
      ;;
    window)
      grit_cmd="$grit window -a $a_file -b $b_file -w 1000 --assume-sorted"
      ;;
    jaccard)
      grit_cmd="$grit jaccard -a $a_file -b $b_file"
      ;;
    complement)
      grit_cmd="$grit complement -i $a_file -g $genome --assume-sorted"
      ;;
    slop)
      grit_cmd="$grit slop -i $a_file -g $genome -b 100"
      ;;
    sort)
      grit_cmd="$grit sort -i $a_file"
      ;;
    *)
      echo "Unknown command: $cmd" >&2
      return 1
      ;;
  esac

  # Warm up filesystem cache
  cat "$a_file" > /dev/null 2>&1
  # Only warm up B file for commands that use it
  case "$cmd" in
    coverage|intersect|subtract|closest|window|jaccard)
      cat "$b_file" > /dev/null 2>&1
      ;;
  esac

  # Run benchmark
  eval "$time_cmd $grit_cmd" > "$tmp_out" 2> "$tmp_time"

  extract_metrics "$tmp_time"
  local time_sec="$TIME_SEC"
  local mem_rss="$TIME_RSS"
  local mem_fmt=$(format_bytes "$mem_rss")
  local out_lines=$(wc -l < "$tmp_out" | tr -d ' ')

  # Print result
  printf "%-12s %10.3fs %12s %12s\n" "$cmd" "$time_sec" "$mem_fmt" "$out_lines"

  # Cleanup
  rm -f "$tmp_out" "$tmp_time"
}

#############################################################################
# Main
#############################################################################

if [[ $# -lt 2 ]]; then
  usage
fi

A_SIZE="$1"
B_SIZE="$2"
shift 2

# Parse commands
if [[ $# -eq 0 ]] || [[ "$1" == "all" ]]; then
  COMMANDS="$ALL_COMMANDS"
else
  COMMANDS="$*"
fi

# Validate data exists (use uniform distribution by default)
DATA_PATH="$DATA_DIR/${A_SIZE}_${B_SIZE}/uniform"
A_FILE="$DATA_PATH/A.sorted.bed"
B_FILE="$DATA_PATH/B.sorted.bed"
GENOME="$DATA_PATH/genome.txt"

if [[ ! -f "$A_FILE" || ! -f "$B_FILE" ]]; then
  echo "ERROR: Data not found at $DATA_PATH"
  echo ""
  echo "Available data sizes:"
  for dir in "$DATA_DIR"/*/; do
    [[ -d "$dir" ]] && echo "  $(basename "$dir")"
  done
  echo ""
  echo "Generate data with: ./bench.sh data $A_SIZE $B_SIZE"
  exit 1
fi

# Setup
GRIT=$(detect_grit)
TIME_CMD=$(get_time_cmd)

A_COUNT=$(wc -l < "$A_FILE" | tr -d ' ')
B_COUNT=$(wc -l < "$B_FILE" | tr -d ' ')

# Print header
echo ""
echo -e "${BOLD}GRIT-Only Benchmark${NC}"
echo "═══════════════════════════════════════════════════════"
echo -e "Data: ${CYAN}$A_COUNT${NC} A intervals, ${CYAN}$B_COUNT${NC} B intervals"
echo "Commands: $COMMANDS"
echo ""
printf "%-12s %11s %12s %12s\n" "Command" "Time" "Memory" "Output"
echo "─────────────────────────────────────────────────────────"

# Run benchmarks
for cmd in $COMMANDS; do
  run_grit_benchmark "$cmd" "$A_FILE" "$B_FILE" "$GENOME" "$GRIT" "$TIME_CMD"
done

echo "═══════════════════════════════════════════════════════"
echo -e "${GREEN}Done${NC}"
