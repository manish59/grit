#!/usr/bin/env bash
#############################################################################
# Focused Benchmark: GRIT vs bedtools vs bedops
#
# Clear comparison with explicit flags for each tool
#############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
RESULTS_DIR="$SCRIPT_DIR/results/focused"

# Colors
if [[ -t 1 ]]; then
  GREEN='\033[0;32m'
  RED='\033[0;31m'
  YELLOW='\033[1;33m'
  CYAN='\033[0;36m'
  BOLD='\033[1m'
  NC='\033[0m'
else
  GREEN='' RED='' YELLOW='' CYAN='' BOLD='' NC=''
fi

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
  if (( $(echo "$seconds < 1" | bc -l) )); then
    printf "%.3fs" "$seconds"
  elif (( $(echo "$seconds < 60" | bc -l) )); then
    printf "%.2fs" "$seconds"
  else
    local mins=$(echo "$seconds / 60" | bc)
    local secs=$(echo "$seconds - $mins * 60" | bc)
    printf "%dm%.1fs" "$mins" "$secs"
  fi
}

format_memory_mb() {
  local mb="$1"
  if (( $(echo "$mb < 1024" | bc -l) )); then
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
  local b_size="$2"

  local a_count=$(parse_size "$a_size")
  local b_count=$(parse_size "$b_size")

  local data_subdir="$DATA_DIR/${a_size}_${b_size}"
  mkdir -p "$data_subdir"

  local genome_file="$data_subdir/genome.txt"
  local file_a="$data_subdir/a.bed"
  local file_b="$data_subdir/b.bed"
  local file_a_sorted="$data_subdir/a.sorted.bed"
  local file_b_sorted="$data_subdir/b.sorted.bed"

  # Create genome file
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

  # Generate and sort A
  if [[ ! -f "$file_a_sorted" ]]; then
    if [[ ! -f "$file_a" ]]; then
      echo -e "${CYAN}Generating $a_size intervals for A...${NC}" >&2
      bedtools random -l 150 -n "$a_count" -seed 42 -g "$genome_file" > "$file_a"
    fi
    echo -e "${CYAN}Sorting A ($a_size)...${NC}" >&2
    LC_ALL=C sort -k1,1 -k2,2n -S 4G --parallel=4 "$file_a" > "$file_a_sorted"
  fi

  # Generate and sort B
  if [[ ! -f "$file_b_sorted" ]]; then
    if [[ ! -f "$file_b" ]]; then
      echo -e "${CYAN}Generating $b_size intervals for B...${NC}" >&2
      bedtools random -l 500 -n "$b_count" -seed 123 -g "$genome_file" > "$file_b"
    fi
    echo -e "${CYAN}Sorting B ($b_size)...${NC}" >&2
    LC_ALL=C sort -k1,1 -k2,2n -S 4G --parallel=4 "$file_b" > "$file_b_sorted"
  fi

  echo "$data_subdir"
}

#############################################################################
# Benchmark Functions
#############################################################################

run_benchmark() {
  local tool="$1"
  local op="$2"
  local file_a="$3"
  local file_b="$4"
  local output="$5"
  local time_file="$6"

  local TIME_CMD=$(detect_time_cmd)
  local grit_bin="$SCRIPT_DIR/../target/release/grit"

  case "$tool" in
    grit)
      case "$op" in
        intersect)
          # GRIT: --streaming for O(k) memory, --assume-sorted to skip validation
          $TIME_CMD "$grit_bin" intersect -a "$file_a" -b "$file_b" \
            --streaming --assume-sorted > "$output" 2> "$time_file"
          ;;
        merge)
          # GRIT merge: streaming by default, --assume-sorted for speed
          $TIME_CMD "$grit_bin" merge -i "$file_a" \
            --assume-sorted > "$output" 2> "$time_file"
          ;;
        subtract)
          # GRIT: --streaming for O(k) memory
          $TIME_CMD "$grit_bin" subtract -a "$file_a" -b "$file_b" \
            --streaming --assume-sorted > "$output" 2> "$time_file"
          ;;
        closest)
          # GRIT: --streaming for O(k) memory
          $TIME_CMD "$grit_bin" closest -a "$file_a" -b "$file_b" \
            --streaming --assume-sorted > "$output" 2> "$time_file"
          ;;
        coverage)
          # GRIT coverage
          $TIME_CMD "$grit_bin" coverage -a "$file_a" -b "$file_b" \
            --assume-sorted > "$output" 2> "$time_file"
          ;;
      esac
      ;;
    bedtools)
      case "$op" in
        intersect)
          # bedtools: -sorted for streaming on pre-sorted input
          $TIME_CMD bedtools intersect -a "$file_a" -b "$file_b" \
            -sorted > "$output" 2> "$time_file"
          ;;
        merge)
          # bedtools merge: input must be sorted (no flag needed)
          $TIME_CMD bedtools merge -i "$file_a" > "$output" 2> "$time_file"
          ;;
        subtract)
          # bedtools: -sorted for streaming
          $TIME_CMD bedtools subtract -a "$file_a" -b "$file_b" \
            -sorted > "$output" 2> "$time_file"
          ;;
        closest)
          # bedtools closest: no -sorted flag available
          $TIME_CMD bedtools closest -a "$file_a" -b "$file_b" > "$output" 2> "$time_file"
          ;;
        coverage)
          # bedtools: -sorted for streaming
          $TIME_CMD bedtools coverage -a "$file_a" -b "$file_b" \
            -sorted > "$output" 2> "$time_file"
          ;;
      esac
      ;;
    bedops)
      case "$op" in
        intersect)
          # bedops: requires pre-sorted input (no flag, it's mandatory)
          $TIME_CMD bedops --intersect "$file_a" "$file_b" > "$output" 2> "$time_file"
          ;;
        merge)
          # bedops merge
          $TIME_CMD bedops --merge "$file_a" > "$output" 2> "$time_file"
          ;;
        subtract)
          # bedops: --difference is subtract
          $TIME_CMD bedops --difference "$file_a" "$file_b" > "$output" 2> "$time_file"
          ;;
        closest)
          # bedops: closest-features
          $TIME_CMD closest-features --closest "$file_a" "$file_b" > "$output" 2> "$time_file"
          ;;
        *)
          return 1
          ;;
      esac
      ;;
  esac
}

#############################################################################
# Main Benchmark Runner
#############################################################################

run_full_benchmark() {
  local a_size="$1"
  local b_size="$2"
  local operations="${3:-intersect merge subtract closest}"

  echo ""
  echo -e "${BOLD}${CYAN}╔════════════════════════════════════════════════════════════════╗${NC}"
  echo -e "${BOLD}${CYAN}║  GRIT vs bedtools vs bedops Benchmark: ${a_size} × ${b_size}${NC}"
  echo -e "${BOLD}${CYAN}╚════════════════════════════════════════════════════════════════╝${NC}"
  echo ""

  # Generate data
  local data_dir
  data_dir=$(generate_data "$a_size" "$b_size")

  local file_a="$data_dir/a.sorted.bed"
  local file_b="$data_dir/b.sorted.bed"

  # Verify files exist
  if [[ ! -f "$file_a" ]] || [[ ! -f "$file_b" ]]; then
    echo -e "${RED}Error: Data files not found${NC}"
    return 1
  fi

  # Show data info
  local a_lines=$(wc -l < "$file_a" | tr -d ' ')
  local b_lines=$(wc -l < "$file_b" | tr -d ' ')
  local a_size_bytes=$(stat -f%z "$file_a" 2>/dev/null || stat -c%s "$file_a")
  local b_size_bytes=$(stat -f%z "$file_b" 2>/dev/null || stat -c%s "$file_b")

  echo -e "${YELLOW}Dataset:${NC}"
  echo "  A: $a_lines intervals ($(format_memory_mb $(echo "$a_size_bytes / 1048576" | bc -l)))"
  echo "  B: $b_lines intervals ($(format_memory_mb $(echo "$b_size_bytes / 1048576" | bc -l)))"
  echo ""

  echo -e "${YELLOW}Commands used:${NC}"
  echo "  GRIT:     --streaming --assume-sorted (O(k) memory)"
  echo "  bedtools: -sorted (streaming on sorted input)"
  echo "  bedops:   requires pre-sorted input (default)"
  echo ""

  # Create results directory
  local timestamp=$(date +%Y%m%d_%H%M%S)
  local results_subdir="$RESULTS_DIR/${a_size}_${b_size}_${timestamp}"
  mkdir -p "$results_subdir"

  # CSV output
  local csv_file="$results_subdir/results.csv"
  echo "operation,tool,time_s,memory_mb,output_lines" > "$csv_file"

  # Run benchmarks
  local tools="grit bedtools bedops"

  for op in $operations; do
    echo -e "${BOLD}┌─────────────────────────────────────────────────────────────────┐${NC}"
    echo -e "${BOLD}│ Operation: ${op}${NC}"
    echo -e "${BOLD}└─────────────────────────────────────────────────────────────────┘${NC}"

    for tool in $tools; do
      local output="$results_subdir/${tool}_${op}.bed"
      local time_file="$results_subdir/${tool}_${op}.time"

      printf "  %-10s " "$tool:"

      if run_benchmark "$tool" "$op" "$file_a" "$file_b" "$output" "$time_file" 2>/dev/null; then
        extract_metrics "$time_file"
        local out_lines=$(wc -l < "$output" 2>/dev/null | tr -d ' ' || echo "0")

        printf "${GREEN}%10s${NC}  " "$(format_time "$TIME_SEC")"
        printf "Memory: ${CYAN}%s${NC}  " "$(format_memory_mb "$TIME_RSS_MB")"
        printf "Output: %s lines\n" "$out_lines"

        echo "$op,$tool,$TIME_SEC,$TIME_RSS_MB,$out_lines" >> "$csv_file"
      else
        printf "${RED}FAILED${NC}\n"
        echo "$op,$tool,NA,NA,NA" >> "$csv_file"
      fi
    done
    echo ""
  done

  # Summary table
  echo -e "${BOLD}${CYAN}═══════════════════════════════════════════════════════════════════${NC}"
  echo -e "${BOLD}Summary: Speedup vs bedtools${NC}"
  echo -e "${CYAN}═══════════════════════════════════════════════════════════════════${NC}"

  printf "%-12s" "Operation"
  printf "%12s" "GRIT"
  printf "%12s" "bedops"
  printf "%15s" "GRIT Memory"
  printf "%17s\n" "bedtools Memory"
  echo "─────────────────────────────────────────────────────────────────────"

  for op in $operations; do
    local grit_time=$(grep "^$op,grit," "$csv_file" | cut -d, -f3)
    local bedtools_time=$(grep "^$op,bedtools," "$csv_file" | cut -d, -f3)
    local bedops_time=$(grep "^$op,bedops," "$csv_file" | cut -d, -f3)
    local grit_mem=$(grep "^$op,grit," "$csv_file" | cut -d, -f4)
    local bedtools_mem=$(grep "^$op,bedtools," "$csv_file" | cut -d, -f4)

    printf "%-12s" "$op"

    if [[ "$grit_time" != "NA" && "$bedtools_time" != "NA" ]]; then
      local grit_speedup=$(awk "BEGIN {printf \"%.1f\", $bedtools_time / $grit_time}")
      printf "${GREEN}%10sx${NC}  " "$grit_speedup"
    else
      printf "%12s" "-"
    fi

    if [[ "$bedops_time" != "NA" && "$bedtools_time" != "NA" ]]; then
      local bedops_speedup=$(awk "BEGIN {printf \"%.1f\", $bedtools_time / $bedops_time}")
      printf "%10sx  " "$bedops_speedup"
    else
      printf "%12s" "-"
    fi

    if [[ "$grit_mem" != "NA" ]]; then
      printf "%12s MB" "$grit_mem"
    else
      printf "%15s" "-"
    fi

    if [[ "$bedtools_mem" != "NA" ]]; then
      printf "%14s MB" "$bedtools_mem"
    else
      printf "%17s" "-"
    fi

    echo ""
  done

  echo ""
  echo -e "${GREEN}Results saved to: $csv_file${NC}"
}

#############################################################################
# Main
#############################################################################

usage() {
  cat <<EOF
Focused Benchmark: GRIT vs bedtools vs bedops

Usage:
  $0 <a_size> <b_size> [operations]

Examples:
  $0 10M 5M                    # All operations
  $0 50M 5M intersect merge    # Specific operations
  $0 100M 10M                  # Large scale test

Operations: intersect, merge, subtract, closest
EOF
}

if [[ $# -lt 2 ]]; then
  usage
  exit 1
fi

a_size="$1"
b_size="$2"
shift 2
operations="${*:-intersect merge subtract closest}"

run_full_benchmark "$a_size" "$b_size" "$operations"
