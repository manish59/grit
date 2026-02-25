#!/usr/bin/env bash
set -euo pipefail

#############################################################################
# verify_all.sh - Correctness verification harness for GRIT vs bedtools
#
# Usage:
#   ./verify_all.sh <A_size> <B_size> [distribution] [commands...]
#   ./verify_all.sh 10M 5M                    # All commands, uniform
#   ./verify_all.sh 10M 5M clustered          # All commands, clustered
#   ./verify_all.sh 10M 5M uniform intersect  # Single command
#   ./verify_all.sh 10M 5M all all            # All distributions, all commands
#
# Output:
#   - Pass/fail status per command/distribution/mode
#   - Diffs saved to benchmarks/diffs/<command>/<distribution>/
#############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$SCRIPT_DIR/.."
DATA_DIR="$BENCH_DIR/data"
DIFF_DIR="$BENCH_DIR/diffs"
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

# Platform-aware SHA256
sha256_hash() {
  if command -v sha256sum &>/dev/null; then
    sha256sum "$@" | awk '{print $1}'
  else
    shasum -a 256 "$@" | awk '{print $1}'
  fi
}

# Canonicalize output for comparison
# Sort, normalize whitespace
canonicalize() {
  local input="$1"
  local output="$2"
  LC_ALL=C sort "$input" | sed 's/[[:space:]][[:space:]]*/	/g' | grep -v '^$' > "$output" || true
}

# Canonicalize coverage output (first 6 columns only due to float precision)
canonicalize_coverage() {
  local input="$1"
  local output="$2"
  cut -f1-6 "$input" | LC_ALL=C sort | sed 's/[[:space:]][[:space:]]*/	/g' | grep -v '^$' > "$output" || true
}

#############################################################################
# Verification
#############################################################################

verify_command() {
  local cmd="$1"
  local a_file="$2"
  local b_file="$3"
  local genome="$4"
  local diff_path="$5"
  local grit="$6"
  local size_spec="$7"
  local distribution="$8"

  mkdir -p "$diff_path"

  # Check for truth data
  local truth_base="$TRUTH_DIR/${size_spec}/${distribution}"
  local use_truth=false
  if [[ -f "${truth_base}/${cmd}.out" && -f "${truth_base}/${cmd}.sha256" ]]; then
    use_truth=true
  fi

  local bt_out=$(mktemp)
  local gr_out=$(mktemp)
  local gr_stream_out=$(mktemp)
  local bt_canon=$(mktemp)
  local gr_canon=$(mktemp)
  local gr_stream_canon=$(mktemp)

  # Build commands based on type
  local bt_cmd gr_cmd gr_stream_cmd has_streaming=true
  case "$cmd" in
    coverage)
      bt_cmd="bedtools coverage -a $a_file -b $b_file -sorted"
      gr_cmd="$grit coverage -a $a_file -b $b_file --assume-sorted"
      gr_stream_cmd="$gr_cmd"  # Coverage is always streaming
      has_streaming=false
      ;;
    intersect)
      bt_cmd="bedtools intersect -a $a_file -b $b_file -sorted"
      gr_cmd="$grit intersect -a $a_file -b $b_file"
      gr_stream_cmd="$grit intersect -a $a_file -b $b_file --streaming --assume-sorted"
      ;;
    merge)
      bt_cmd="bedtools merge -i $a_file"
      gr_cmd="$grit merge -i $a_file"
      gr_stream_cmd="$grit merge -i $a_file --assume-sorted"
      has_streaming=false  # Default is already streaming
      ;;
    subtract)
      bt_cmd="bedtools subtract -a $a_file -b $b_file -sorted"
      gr_cmd="$grit subtract -a $a_file -b $b_file"
      gr_stream_cmd="$grit subtract -a $a_file -b $b_file --streaming --assume-sorted"
      ;;
    closest)
      bt_cmd="bedtools closest -a $a_file -b $b_file"
      gr_cmd="$grit closest -a $a_file -b $b_file"
      gr_stream_cmd="$grit closest -a $a_file -b $b_file --streaming --assume-sorted"
      ;;
    window)
      bt_cmd="bedtools window -a $a_file -b $b_file -w 1000"
      gr_cmd="$grit window -a $a_file -b $b_file -w 1000 --assume-sorted"
      gr_stream_cmd="$gr_cmd"  # Window is always streaming
      has_streaming=false
      ;;
    jaccard)
      bt_cmd="bedtools jaccard -a $a_file -b $b_file"
      gr_cmd="$grit jaccard -a $a_file -b $b_file"
      gr_stream_cmd="$gr_cmd"  # Jaccard doesn't have streaming mode
      has_streaming=false
      ;;
    complement)
      bt_cmd="bedtools complement -i $a_file -g $genome"
      gr_cmd="$grit complement -i $a_file -g $genome"
      gr_stream_cmd="$grit complement -i $a_file -g $genome --assume-sorted"
      has_streaming=false  # Default handles sorting
      ;;
    slop)
      bt_cmd="bedtools slop -i $a_file -g $genome -b 100"
      gr_cmd="$grit slop -i $a_file -g $genome -b 100"
      gr_stream_cmd="$gr_cmd"  # Slop doesn't have streaming mode
      has_streaming=false
      ;;
    sort)
      bt_cmd="LC_ALL=C sort -k1,1 -k2,2n -k3,3n $a_file"
      gr_cmd="$grit sort -i $a_file"
      gr_stream_cmd="$gr_cmd"  # Sort doesn't have streaming mode
      has_streaming=false
      ;;
    *)
      echo "Unknown command: $cmd"
      return 1
      ;;
  esac

  local cmd_display="$cmd"
  if $use_truth; then
    # Use pre-computed truth data
    cp "${truth_base}/${cmd}.out" "$bt_out"
    cmd_display="${cmd}*"
  else
    # Run bedtools
    eval "$bt_cmd" > "$bt_out" 2>/dev/null || true
  fi

  # Run grit default
  eval "$gr_cmd" > "$gr_out" 2>/dev/null || true

  # Run grit streaming (if different)
  if $has_streaming; then
    eval "$gr_stream_cmd" > "$gr_stream_out" 2>/dev/null || true
  else
    cp "$gr_out" "$gr_stream_out"
  fi

  # Canonicalize all outputs
  if [[ "$cmd" == "coverage" ]]; then
    canonicalize_coverage "$bt_out" "$bt_canon"
    canonicalize_coverage "$gr_out" "$gr_canon"
    canonicalize_coverage "$gr_stream_out" "$gr_stream_canon"
  else
    canonicalize "$bt_out" "$bt_canon"
    canonicalize "$gr_out" "$gr_canon"
    canonicalize "$gr_stream_out" "$gr_stream_canon"
  fi

  # Compare hashes
  local bt_hash=$(sha256_hash "$bt_canon")
  local gr_hash=$(sha256_hash "$gr_canon")
  local gr_stream_hash=$(sha256_hash "$gr_stream_canon")

  local status_default="PASS"
  local status_stream="PASS"

  if [[ "$bt_hash" != "$gr_hash" ]]; then
    status_default="FAIL"
    diff -u "$bt_canon" "$gr_canon" > "$diff_path/default.diff" 2>/dev/null || true
    # Save line counts for debugging
    echo "bedtools: $(wc -l < "$bt_canon") lines" > "$diff_path/default.info"
    echo "grit:     $(wc -l < "$gr_canon") lines" >> "$diff_path/default.info"
  fi

  if [[ "$bt_hash" != "$gr_stream_hash" ]]; then
    status_stream="FAIL"
    diff -u "$bt_canon" "$gr_stream_canon" > "$diff_path/streaming.diff" 2>/dev/null || true
    echo "bedtools: $(wc -l < "$bt_canon") lines" > "$diff_path/streaming.info"
    echo "grit:     $(wc -l < "$gr_stream_canon") lines" >> "$diff_path/streaming.info"
  fi

  # Color output
  local default_colored stream_colored
  if [[ "$status_default" == "PASS" ]]; then
    default_colored="${GREEN}PASS${NC}"
  else
    default_colored="${RED}FAIL${NC}"
  fi
  if [[ "$status_stream" == "PASS" ]]; then
    stream_colored="${GREEN}PASS${NC}"
  else
    stream_colored="${RED}FAIL${NC}"
  fi

  # Print result
  if $has_streaming; then
    printf "  %-12s default=%b  streaming=%b\n" "$cmd_display" "$default_colored" "$stream_colored"
  else
    printf "  %-12s default=%b\n" "$cmd_display" "$default_colored"
  fi

  # Cleanup
  rm -f "$bt_out" "$gr_out" "$gr_stream_out" "$bt_canon" "$gr_canon" "$gr_stream_canon"

  # Return status
  [[ "$status_default" == "PASS" && "$status_stream" == "PASS" ]]
}

verify_distribution() {
  local a_size="$1"
  local b_size="$2"
  local distribution="$3"
  local commands="$4"
  local grit="$5"

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

  # Check for truth data
  local truth_base="$TRUTH_DIR/${size_spec}/${distribution}"
  local truth_status="${YELLOW}not available${NC}"
  if [[ -d "$truth_base" ]]; then
    truth_status="${GREEN}available${NC} (* = using cached bedtools)"
  fi

  echo -e "${BOLD}Distribution: $distribution${NC}"
  echo "  Data: $data_path"
  echo -e "  Truth data: $truth_status"
  echo ""

  local all_pass=true
  for cmd in $commands; do
    local diff_path="$DIFF_DIR/${size_spec}/$cmd/$distribution"
    if ! verify_command "$cmd" "$a_file" "$b_file" "$genome" "$diff_path" "$grit" "$size_spec" "$distribution"; then
      all_pass=false
    fi
  done

  echo ""
  $all_pass
}

#############################################################################
# Main
#############################################################################

main() {
  if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <A_size> <B_size> [distribution] [commands...]"
    echo "Example: $0 10M 5M"
    echo "         $0 10M 5M clustered"
    echo "         $0 10M 5M uniform intersect merge"
    echo "         $0 10M 5M all all"
    exit 1
  fi

  local a_size="$1"
  local b_size="$2"
  local distribution="${3:-uniform}"
  shift 2
  shift || true  # shift distribution if provided
  local commands="${*:-$ALL_COMMANDS}"

  # Handle 'all' keywords
  if [[ "$distribution" == "all" ]]; then
    distribution="$ALL_DISTRIBUTIONS"
  fi
  if [[ "$commands" == "all" || -z "$commands" ]]; then
    commands="$ALL_COMMANDS"
  fi

  local grit=$(detect_grit)

  echo -e "${BOLD}=== GRIT Correctness Verification ===${NC}"
  echo "Size: ${a_size} A × ${b_size} B"
  echo "Distributions: $distribution"
  echo "Commands: $commands"
  echo ""

  local all_pass=true
  for dist in $distribution; do
    if ! verify_distribution "$a_size" "$b_size" "$dist" "$commands" "$grit"; then
      all_pass=false
    fi
  done

  echo "═══════════════════════════════════════════════════════════════════"
  if $all_pass; then
    echo -e "${GREEN}All verifications PASSED${NC}"
  else
    echo -e "${RED}Some verifications FAILED${NC}"
    echo "Check diffs at: $DIFF_DIR/${a_size}_${b_size}/"
    exit 1
  fi
}

main "$@"
