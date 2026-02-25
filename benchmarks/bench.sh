#!/usr/bin/env bash
set -euo pipefail

#############################################################################
# GRIT Benchmark Suite
#
# Unified benchmarking tool for comparing GRIT vs bedtools performance.
#
# Usage:
#   ./bench.sh data <A_size> <B_size>     Generate test data
#   ./bench.sh run <A_size> <B_size> <commands...>  Run benchmarks
#   ./bench.sh list                        List available commands
#   ./bench.sh clean                       Remove all generated data
#
# Examples:
#   ./bench.sh data 10M 5M                 # Generate 10M A, 5M B intervals
#   ./bench.sh run 10M 5M coverage         # Benchmark coverage only
#   ./bench.sh run 10M 5M all              # Benchmark all commands
#   ./bench.sh run 100M 50M coverage intersect merge
#
# Size notation:
#   1K = 1,000  |  1M = 1,000,000  |  100M = 100,000,000
#
# Data is cached in benchmarks/data/<size>/ and reused across runs.
#############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
RESULTS_DIR="$SCRIPT_DIR/results"
TRUTH_DIR="$SCRIPT_DIR/truth"

# Source common functions
source "$SCRIPT_DIR/scripts/common.sh" 2>/dev/null || true

# Colors (disabled if not tty)
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

# All dataset distributions
ALL_DISTRIBUTIONS="uniform clustered"

#############################################################################
# Utility Functions
#############################################################################

usage() {
  cat <<EOF
GRIT Benchmark Suite

Usage:
  $0 env                                    Verify environment (grit, bedtools, rustc)
  $0 data <A_size> <B_size> [distribution]  Generate test data (uniform|clustered)
  $0 run <A_size> <B_size> <cmd...>         Run benchmarks (uses truth data if available)
  $0 verify <A_size> <B_size> [cmd...]      Verify correctness vs bedtools
  $0 bench <A_size> <B_size> [cmd...]       Run performance benchmarks with CSV output
  $0 list                                   List available commands
  $0 clean                                  Remove all generated data

Truth Data (pre-computed bedtools results):
  $0 truth <A_size> <B_size> [cmd...]       Generate truth data from bedtools
  $0 truth-list [A_size_B_size]             List available truth data
  $0 truth-verify <A_size> <B_size>         Verify truth data integrity

Real-World Datasets:
  $0 real-list                              List available real-world datasets
  $0 real-download <dataset> [--yes]        Download dataset from source
  $0 real-prepare <dataset>                 Prepare dataset for benchmarking
  $0 real-truth <dataset> <cmd|all>         Generate bedtools truth baseline
  $0 real-run <dataset> <cmd|all>           Run GRIT benchmark vs bedtools

Examples:
  $0 env                              # Check environment
  $0 data 10M 5M                      # Generate 10M A, 5M B uniform data
  $0 data 10M 5M clustered            # Generate clustered (adversarial) data
  $0 truth 10M 5M all                 # Pre-compute bedtools outputs
  $0 truth 10M 5M coverage --force    # Regenerate specific command
  $0 run 10M 5M coverage intersect    # Quick benchmark (uses truth if available)
  $0 verify 10M 5M all                # Verify all commands
  $0 bench 10M 5M all                 # Full benchmark with CSV

Real-World Examples:
  $0 real-list                        # Show available datasets
  $0 real-download dbsnp --yes        # Download dbSNP data
  $0 real-prepare dbsnp               # Convert and sort for benchmarking
  $0 real-truth dbsnp all             # Generate bedtools baselines
  $0 real-run dbsnp all               # Benchmark GRIT vs bedtools

Size notation: 1K, 10K, 100K, 1M, 10M, 100M
Distributions: uniform (default), clustered
Datasets: dbsnp, encode_peaks, gencode, sv
EOF
  exit 1
}

# Parse size notation: 10M -> 10000000, 500K -> 500000
parse_size() {
  local size="$1"
  local num="${size%[KkMm]}"
  local suffix="${size: -1}"

  case "$suffix" in
    K|k) echo $((num * 1000)) ;;
    M|m) echo $((num * 1000000)) ;;
    *)   echo "$size" ;;  # Already a number
  esac
}

# Format number with suffix: 10000000 -> 10M
format_size() {
  local n="$1"
  if (( n >= 1000000 )); then
    echo "$((n / 1000000))M"
  elif (( n >= 1000 )); then
    echo "$((n / 1000))K"
  else
    echo "$n"
  fi
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
    echo -e "${RED}ERROR: GRIT binary not found.${NC}" >&2
    echo "Build with: cargo build --release" >&2
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

# Platform-aware SHA256 (macOS uses shasum, Ubuntu uses sha256sum)
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

# Extract time and memory from time output
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

# Create genome file (hg38)
create_genome() {
  local genome_file="$1"
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
}

# Clear filesystem cache by reading files
clear_cache() {
  local a_file="$1"
  local b_file="$2"

  # Read files to clear any OS-level caching effects
  cat "$a_file" > /dev/null 2>&1 || true
  cat "$b_file" > /dev/null 2>&1 || true

  # Sync to flush buffers
  sync 2>/dev/null || true

  # Small delay to let caches settle
  sleep 0.1
}

#############################################################################
# Data Generation
#############################################################################

cmd_data() {
  if [[ $# -lt 2 ]]; then
    echo "Usage: $0 data <A_size> <B_size> [distribution]"
    echo "Example: $0 data 10M 5M"
    echo "         $0 data 10M 5M clustered"
    exit 1
  fi

  local a_size="$1"
  local b_size="$2"
  local distribution="${3:-uniform}"

  # Validate distribution
  if [[ "$distribution" != "uniform" && "$distribution" != "clustered" ]]; then
    echo -e "${RED}ERROR: Invalid distribution '$distribution'. Use 'uniform' or 'clustered'.${NC}"
    exit 1
  fi

  local data_path="$DATA_DIR/${a_size}_${b_size}/${distribution}"

  # Check if data already exists
  if [[ -f "$data_path/A.sorted.bed" && -f "$data_path/B.sorted.bed" ]]; then
    echo -e "${GREEN}Data already exists:${NC} $data_path"
    echo "  A: $(wc -l < "$data_path/A.sorted.bed" | tr -d ' ') intervals"
    echo "  B: $(wc -l < "$data_path/B.sorted.bed" | tr -d ' ') intervals"
    return 0
  fi

  mkdir -p "$data_path"
  local grit=$(detect_grit)
  local tmp_dir="$data_path/tmp"
  local genome="$data_path/genome.txt"

  echo -e "${BOLD}Generating benchmark data with GRIT${NC}"
  echo "  A intervals: $a_size"
  echo "  B intervals: $b_size"
  echo "  Distribution: $distribution"
  echo "  Output: $data_path"
  echo ""

  # Create genome file first (needed for re-sorting to genome order)
  echo -n "[1/4] Creating genome file... "
  create_genome "$genome"
  echo "done"

  # Generate using grit generate (creates lexicographically sorted files)
  echo -n "[2/4] Generating intervals with GRIT... "
  local start_time=$(date +%s)

  # Build generate command based on distribution
  local gen_cmd="$grit generate --output $tmp_dir --a $a_size --b $b_size --sorted yes --seed 42 --force"
  if [[ "$distribution" == "clustered" ]]; then
    gen_cmd="$gen_cmd --mode clustered"
  fi

  eval "$gen_cmd" >/dev/null 2>&1
  local elapsed=$(($(date +%s) - start_time))
  echo "done (${elapsed}s)"

  # grit generate outputs to: <output>/custom_A<size>_B<size>/A.bed
  local gen_dir="$tmp_dir/custom_A${a_size}_B${b_size}"

  # Re-sort files to genome order (required for bedtools complement/closest)
  echo -n "[3/4] Re-sorting to genome order... "
  start_time=$(date +%s)
  bedtools sort -i "$gen_dir/A.bed" -g "$genome" > "$data_path/A.sorted.bed"
  bedtools sort -i "$gen_dir/B.bed" -g "$genome" > "$data_path/B.sorted.bed"
  elapsed=$(($(date +%s) - start_time))
  echo "done (${elapsed}s)"

  # Cleanup temp files
  echo -n "[4/4] Cleaning up... "
  rm -rf "$tmp_dir"
  echo "done"

  echo ""
  echo -e "${GREEN}Data ready:${NC} $data_path"
  echo "  A.sorted.bed: $(wc -l < "$data_path/A.sorted.bed" | tr -d ' ') intervals"
  echo "  B.sorted.bed: $(wc -l < "$data_path/B.sorted.bed" | tr -d ' ') intervals"
}

#############################################################################
# Benchmark Execution
#############################################################################

run_single_benchmark() {
  local cmd="$1"
  local a_file="$2"
  local b_file="$3"
  local genome="$4"
  local results_path="$5"
  local grit="$6"
  local time_cmd="$7"
  local size_spec="$8"
  local distribution="$9"

  local bt_out="$results_path/${cmd}_bedtools.out"
  local gr_out="$results_path/${cmd}_grit.out"
  local bt_time="$results_path/${cmd}_bedtools.time"
  local gr_time="$results_path/${cmd}_grit.time"

  local bt_sec bt_rss bt_lines bt_hash
  local use_truth=false

  # Check for truth data
  local truth_base="$TRUTH_DIR/${size_spec}/${distribution}"
  if [[ -f "${truth_base}/${cmd}.out" && -f "${truth_base}/${cmd}.time" && -f "${truth_base}/${cmd}.sha256" ]]; then
    use_truth=true
  fi

  # Build commands based on type
  local bt_cmd gr_cmd
  case "$cmd" in
    coverage)
      bt_cmd="bedtools coverage -a $a_file -b $b_file -sorted"
      gr_cmd="$grit coverage -a $a_file -b $b_file --assume-sorted"
      ;;
    intersect)
      bt_cmd="bedtools intersect -a $a_file -b $b_file -sorted"
      gr_cmd="$grit intersect -a $a_file -b $b_file --assume-sorted --streaming"
      ;;
    merge)
      bt_cmd="bedtools merge -i $a_file"
      gr_cmd="$grit merge -i $a_file --assume-sorted"
      ;;
    subtract)
      bt_cmd="bedtools subtract -a $a_file -b $b_file -sorted"
      gr_cmd="$grit subtract -a $a_file -b $b_file --assume-sorted --streaming"
      ;;
    closest)
      bt_cmd="bedtools closest -a $a_file -b $b_file"
      gr_cmd="$grit closest -a $a_file -b $b_file --assume-sorted --streaming"
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

  if $use_truth; then
    # Use pre-computed truth data
    bt_out="${truth_base}/${cmd}.out"
    bt_sec=$(cat "${truth_base}/${cmd}.time")
    bt_rss=$(cat "${truth_base}/${cmd}.rss" 2>/dev/null || echo "0")
    bt_lines=$(cat "${truth_base}/${cmd}.lines" 2>/dev/null || wc -l < "$bt_out" | tr -d ' ')
    bt_hash=$(cat "${truth_base}/${cmd}.sha256")
  else
    # Clear cache before bedtools
    clear_cache "$a_file" "$b_file"

    # Run bedtools
    eval "$time_cmd $bt_cmd" > "$bt_out" 2> "$bt_time"
    extract_metrics "$bt_time"
    bt_sec="$TIME_SEC"
    bt_rss="$TIME_RSS"
    bt_lines=$(wc -l < "$bt_out" | tr -d ' ')

    # Calculate hash
    if [[ "$cmd" == "coverage" ]]; then
      bt_hash=$(cut -f1-6 "$bt_out" | sha256_hash_stdin)
    elif [[ "$cmd" == "window" ]]; then
      bt_hash=$(sort "$bt_out" | sha256_hash_stdin)
    else
      bt_hash=$(sha256_hash "$bt_out")
    fi
  fi

  # Clear cache before grit
  clear_cache "$a_file" "$b_file"

  # Run grit
  eval "$time_cmd $gr_cmd" > "$gr_out" 2> "$gr_time"
  extract_metrics "$gr_time"
  local gr_sec="$TIME_SEC"
  local gr_rss="$TIME_RSS"

  # Check correctness
  local status="FAIL"
  local gr_lines=$(wc -l < "$gr_out" | tr -d ' ')

  if [[ "$bt_lines" == "$gr_lines" ]]; then
    # Calculate grit hash using same method
    local gr_hash
    if [[ "$cmd" == "coverage" ]]; then
      gr_hash=$(cut -f1-6 "$gr_out" | sha256_hash_stdin)
    elif [[ "$cmd" == "window" ]]; then
      gr_hash=$(sort "$gr_out" | sha256_hash_stdin)
    else
      gr_hash=$(sha256_hash "$gr_out")
    fi

    if [[ "$bt_hash" == "$gr_hash" ]]; then
      status="PASS"
    fi
  fi

  # Calculate speedup
  local speedup=$(awk "BEGIN {printf \"%.1f\", $bt_sec/$gr_sec}")

  # Format output
  local bt_rss_fmt=$(format_bytes "$bt_rss")
  local gr_rss_fmt=$(format_bytes "$gr_rss")

  # Color status and add truth indicator
  local cmd_display="$cmd"
  if $use_truth; then
    cmd_display="${cmd}*"
  fi

  local status_colored
  if [[ "$status" == "PASS" ]]; then
    status_colored="${GREEN}PASS${NC}"
  else
    status_colored="${RED}FAIL${NC}"
  fi

  # Print row
  printf "%-12s %10.2fs %10.2fs %10sx %12s %12s   %b\n" \
    "$cmd_display" "$bt_sec" "$gr_sec" "$speedup" "$bt_rss_fmt" "$gr_rss_fmt" "$status_colored"

  # Return status for summary
  [[ "$status" == "PASS" ]]
}

cmd_run() {
  if [[ $# -lt 3 ]]; then
    echo "Usage: $0 run <A_size> <B_size> [distribution] <commands...>"
    echo "Example: $0 run 10M 5M coverage intersect merge"
    echo "         $0 run 10M 5M all"
    echo "         $0 run 10M 5M clustered all"
    exit 1
  fi

  local a_size="$1"
  local b_size="$2"
  shift 2

  # Parse distribution and commands
  local distribution="uniform"
  local commands=""

  while [[ $# -gt 0 ]]; do
    case "$1" in
      uniform|clustered)
        distribution="$1"
        shift
        ;;
      all)
        commands="$ALL_COMMANDS"
        shift
        ;;
      *)
        commands="$commands $1"
        shift
        ;;
    esac
  done

  commands=$(echo "$commands" | xargs)  # trim
  if [[ -z "$commands" ]]; then
    commands="$ALL_COMMANDS"
  fi

  local size_spec="${a_size}_${b_size}"
  local data_path="$DATA_DIR/${size_spec}/${distribution}"
  local a_file="$data_path/A.sorted.bed"
  local b_file="$data_path/B.sorted.bed"
  local genome="$data_path/genome.txt"

  # Check data exists
  if [[ ! -f "$a_file" || ! -f "$b_file" ]]; then
    echo -e "${YELLOW}Data not found. Generating...${NC}"
    cmd_data "$a_size" "$b_size" "$distribution"
    echo ""
  fi

  # Setup
  local grit=$(detect_grit)
  local time_cmd=$(detect_time_cmd)
  local results_path="$RESULTS_DIR/${size_spec}_${distribution}_$(date +%Y%m%d_%H%M%S)"
  mkdir -p "$results_path"

  local a_n=$(wc -l < "$a_file" | tr -d ' ')
  local b_n=$(wc -l < "$b_file" | tr -d ' ')

  # Check if truth data is available
  local truth_base="$TRUTH_DIR/${size_spec}/${distribution}"
  local using_truth=false
  if [[ -d "$truth_base" ]]; then
    using_truth=true
  fi

  # Print header
  echo ""
  echo -e "${BOLD}GRIT Benchmark Results${NC}"
  echo "═══════════════════════════════════════════════════════════════════════════════"
  echo "Data: $(format_size $a_n) A intervals, $(format_size $b_n) B intervals"
  echo "Distribution: $distribution"
  echo "Commands: $commands"
  if $using_truth; then
    echo -e "Truth data: ${GREEN}available${NC} (commands marked * use cached bedtools output)"
  else
    echo -e "Truth data: ${YELLOW}not found${NC} (generate with: ./bench.sh truth $a_size $b_size $distribution all)"
  fi
  echo ""
  printf "%-12s %11s %11s %11s %13s %13s   %s\n" \
    "Command" "bedtools" "GRIT" "Speedup" "BT_RAM" "GRIT_RAM" "Status"
  echo "───────────────────────────────────────────────────────────────────────────────"

  # Run benchmarks
  local all_passed=true
  for cmd in $commands; do
    if ! run_single_benchmark "$cmd" "$a_file" "$b_file" "$genome" "$results_path" "$grit" "$time_cmd" "$size_spec" "$distribution"; then
      all_passed=false
    fi
  done

  echo "═══════════════════════════════════════════════════════════════════════════════"
  echo "Results saved: $results_path"

  if $all_passed; then
    echo -e "Status: ${GREEN}All benchmarks passed${NC}"
  else
    echo -e "Status: ${RED}Some benchmarks failed${NC}"
  fi
}

cmd_list() {
  echo "Available benchmark commands:"
  echo ""
  for cmd in $ALL_COMMANDS; do
    echo "  $cmd"
  done
  echo ""
  echo "Use 'all' to run all commands: $0 run 10M 5M all"
}

cmd_clean() {
  echo "Removing generated data and results..."
  rm -rf "$DATA_DIR" "$RESULTS_DIR"
  echo "Done."
}

#############################################################################
# Truth Data Management
#############################################################################

# Generate truth data for a single command
generate_truth_single() {
  local cmd="$1"
  local a_file="$2"
  local b_file="$3"
  local genome="$4"
  local truth_path="$5"
  local time_cmd="$6"
  local force="${7:-false}"

  # Check if already exists and not forcing
  if [[ -f "${truth_path}.out" && -f "${truth_path}.time" && "$force" != "true" ]]; then
    echo -e "  ${YELLOW}$cmd${NC}: already exists (use --force to regenerate)"
    return 0
  fi

  # Build bedtools command
  local bt_cmd
  case "$cmd" in
    coverage)
      bt_cmd="bedtools coverage -a $a_file -b $b_file -sorted"
      ;;
    intersect)
      bt_cmd="bedtools intersect -a $a_file -b $b_file -sorted"
      ;;
    merge)
      bt_cmd="bedtools merge -i $a_file"
      ;;
    subtract)
      bt_cmd="bedtools subtract -a $a_file -b $b_file -sorted"
      ;;
    closest)
      bt_cmd="bedtools closest -a $a_file -b $b_file"
      ;;
    window)
      bt_cmd="bedtools window -a $a_file -b $b_file -w 1000"
      ;;
    jaccard)
      bt_cmd="bedtools jaccard -a $a_file -b $b_file"
      ;;
    complement)
      bt_cmd="bedtools complement -i $a_file -g $genome"
      ;;
    slop)
      bt_cmd="bedtools slop -i $a_file -g $genome -b 100"
      ;;
    sort)
      bt_cmd="LC_ALL=C sort -k1,1 -k2,2n -k3,3n $a_file"
      ;;
    *)
      echo -e "  ${RED}$cmd${NC}: unknown command"
      return 1
      ;;
  esac

  echo -n "  $cmd: running bedtools... "

  # Clear cache
  clear_cache "$a_file" "$b_file"

  # Run and capture timing
  local tmp_time=$(mktemp)
  eval "$time_cmd $bt_cmd" > "${truth_path}.out" 2> "$tmp_time"

  # Extract metrics
  extract_metrics "$tmp_time"
  echo "$TIME_SEC" > "${truth_path}.time"
  echo "$TIME_RSS" > "${truth_path}.rss"
  rm -f "$tmp_time"

  # Calculate line count
  local lines=$(wc -l < "${truth_path}.out" | tr -d ' ')
  echo "$lines" > "${truth_path}.lines"

  # Calculate SHA256
  # For coverage, hash first 6 columns; for window, hash sorted output
  if [[ "$cmd" == "coverage" ]]; then
    cut -f1-6 "${truth_path}.out" | sha256_hash_stdin > "${truth_path}.sha256"
  elif [[ "$cmd" == "window" ]]; then
    sort "${truth_path}.out" | sha256_hash_stdin > "${truth_path}.sha256"
  else
    sha256_hash "${truth_path}.out" > "${truth_path}.sha256"
  fi

  local rss_mb=$(awk "BEGIN {printf \"%.1f\", $TIME_RSS/1048576}")
  echo "done (${TIME_SEC}s, ${rss_mb}MB, ${lines} lines)"
}

cmd_truth() {
  if [[ $# -lt 2 ]]; then
    echo "Usage: $0 truth <A_size> <B_size> [distribution] [commands...] [--force]"
    echo "Example: $0 truth 10M 5M all"
    echo "         $0 truth 10M 5M clustered coverage intersect"
    echo "         $0 truth 10M 5M all --force"
    exit 1
  fi

  local a_size="$1"
  local b_size="$2"
  shift 2

  local distribution="uniform"
  local commands=""
  local force="false"

  # Parse remaining arguments
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --force)
        force="true"
        shift
        ;;
      uniform|clustered)
        distribution="$1"
        shift
        ;;
      all)
        commands="$ALL_COMMANDS"
        shift
        ;;
      *)
        commands="$commands $1"
        shift
        ;;
    esac
  done

  # Default to all commands if none specified
  commands="${commands:-$ALL_COMMANDS}"
  commands=$(echo "$commands" | xargs)  # trim

  local size_spec="${a_size}_${b_size}"
  local data_path="$DATA_DIR/${size_spec}/${distribution}"
  local a_file="$data_path/A.sorted.bed"
  local b_file="$data_path/B.sorted.bed"
  local genome="$data_path/genome.txt"

  # Check data exists
  if [[ ! -f "$a_file" || ! -f "$b_file" ]]; then
    echo -e "${YELLOW}Data not found. Generating...${NC}"
    cmd_data "$a_size" "$b_size" "$distribution"
    echo ""
  fi

  local truth_base="$TRUTH_DIR/${size_spec}/${distribution}"
  mkdir -p "$truth_base"

  local time_cmd=$(detect_time_cmd)

  echo -e "${BOLD}Generating Truth Data${NC}"
  echo "═══════════════════════════════════════════════════════════════════════════════"
  echo "Size: ${a_size} A × ${b_size} B"
  echo "Distribution: $distribution"
  echo "Commands: $commands"
  echo "Output: $truth_base"
  echo ""

  for cmd in $commands; do
    generate_truth_single "$cmd" "$a_file" "$b_file" "$genome" "$truth_base/$cmd" "$time_cmd" "$force"
  done

  echo ""
  echo -e "${GREEN}Truth data generated successfully.${NC}"
  echo "Use './bench.sh run $a_size $b_size <cmd>' to benchmark against stored truth."
}

cmd_truth_list() {
  local filter="${1:-}"

  echo -e "${BOLD}Available Truth Data${NC}"
  echo "═══════════════════════════════════════════════════════════════════════════════"

  if [[ ! -d "$TRUTH_DIR" ]]; then
    echo "No truth data found. Generate with: ./bench.sh truth <A_size> <B_size> all"
    return 0
  fi

  local count=0
  for size_dir in "$TRUTH_DIR"/*/; do
    [[ -d "$size_dir" ]] || continue
    local size_spec=$(basename "$size_dir")

    # Apply filter if specified
    if [[ -n "$filter" && "$size_spec" != "$filter" ]]; then
      continue
    fi

    for dist_dir in "$size_dir"*/; do
      [[ -d "$dist_dir" ]] || continue
      local distribution=$(basename "$dist_dir")

      echo -e "\n${CYAN}${size_spec}/${distribution}${NC}:"
      printf "  %-12s %10s %10s %12s\n" "Command" "Time(s)" "RSS(MB)" "Lines"
      echo "  ─────────────────────────────────────────────────"

      for out_file in "$dist_dir"*.out; do
        [[ -f "$out_file" ]] || continue
        local cmd=$(basename "$out_file" .out)
        local time_file="${out_file%.out}.time"
        local rss_file="${out_file%.out}.rss"
        local lines_file="${out_file%.out}.lines"

        local time_s="?"
        local rss_mb="?"
        local lines="?"

        [[ -f "$time_file" ]] && time_s=$(cat "$time_file")
        [[ -f "$rss_file" ]] && rss_mb=$(awk "BEGIN {printf \"%.1f\", $(cat "$rss_file")/1048576}")
        [[ -f "$lines_file" ]] && lines=$(cat "$lines_file")

        printf "  %-12s %10s %10s %12s\n" "$cmd" "$time_s" "$rss_mb" "$lines"
        ((count++))
      done
    done
  done

  if [[ $count -eq 0 ]]; then
    echo "No truth data found. Generate with: ./bench.sh truth <A_size> <B_size> all"
  else
    echo ""
    echo "Total: $count truth entries"
  fi
}

cmd_truth_verify() {
  if [[ $# -lt 2 ]]; then
    echo "Usage: $0 truth-verify <A_size> <B_size> [distribution]"
    echo "Example: $0 truth-verify 10M 5M"
    echo "         $0 truth-verify 10M 5M clustered"
    exit 1
  fi

  local a_size="$1"
  local b_size="$2"
  local distribution="${3:-uniform}"
  local size_spec="${a_size}_${b_size}"
  local truth_base="$TRUTH_DIR/${size_spec}/${distribution}"

  echo -e "${BOLD}Verifying Truth Data Integrity${NC}"
  echo "═══════════════════════════════════════════════════════════════════════════════"
  echo "Path: $truth_base"
  echo ""

  if [[ ! -d "$truth_base" ]]; then
    echo -e "${RED}No truth data found for ${size_spec}/${distribution}${NC}"
    echo "Generate with: ./bench.sh truth $a_size $b_size $distribution all"
    exit 1
  fi

  local all_valid=true
  local count=0
  local valid=0

  for out_file in "$truth_base"/*.out; do
    [[ -f "$out_file" ]] || continue
    local cmd=$(basename "$out_file" .out)
    local sha_file="${out_file%.out}.sha256"

    ((count++))
    echo -n "  $cmd: "

    if [[ ! -f "$sha_file" ]]; then
      echo -e "${RED}MISSING SHA256${NC}"
      all_valid=false
      continue
    fi

    local stored_hash=$(cat "$sha_file")

    # Calculate actual hash (same logic as generation)
    local actual_hash
    if [[ "$cmd" == "coverage" ]]; then
      actual_hash=$(cut -f1-6 "$out_file" | sha256_hash_stdin)
    elif [[ "$cmd" == "window" ]]; then
      actual_hash=$(sort "$out_file" | sha256_hash_stdin)
    else
      actual_hash=$(sha256_hash "$out_file")
    fi

    if [[ "$stored_hash" == "$actual_hash" ]]; then
      echo -e "${GREEN}VALID${NC}"
      ((valid++))
    else
      echo -e "${RED}CORRUPTED${NC}"
      echo "    Expected: $stored_hash"
      echo "    Actual:   $actual_hash"
      all_valid=false
    fi
  done

  echo ""
  echo "═══════════════════════════════════════════════════════════════════════════════"
  echo "Verified: $valid/$count"

  if $all_valid; then
    echo -e "${GREEN}All truth data is valid.${NC}"
  else
    echo -e "${RED}Some truth data is corrupted. Regenerate with --force.${NC}"
    exit 1
  fi
}

#############################################################################
# Environment Verification
#############################################################################

cmd_env() {
  echo -e "${BOLD}=== GRIT Audit Environment ===${NC}"
  echo ""

  # GRIT version
  local grit
  if [[ -x "$SCRIPT_DIR/../target/release/grit" ]]; then
    grit="$SCRIPT_DIR/../target/release/grit"
    echo -e "GRIT:     $("$grit" --version 2>/dev/null || echo 'version unknown')"
  else
    echo -e "GRIT:     ${RED}NOT BUILT${NC}"
    echo "          Run: cargo build --release"
  fi

  # bedtools version
  if command -v bedtools &>/dev/null; then
    echo -e "bedtools: $(bedtools --version 2>/dev/null)"
  else
    echo -e "bedtools: ${RED}NOT FOUND${NC}"
  fi

  # Rust version
  if command -v rustc &>/dev/null; then
    echo -e "rustc:    $(rustc --version)"
  else
    echo -e "rustc:    ${RED}NOT FOUND${NC}"
  fi

  echo ""
  echo "Platform: $(uname -s) $(uname -r)"
  echo "Time cmd: $(detect_time_cmd)"
  echo ""

  # Verify requirements
  local ok=true

  if ! command -v bedtools &>/dev/null; then
    echo -e "${RED}ERROR: bedtools not found on PATH${NC}"
    ok=false
  fi

  if [[ ! -x "$SCRIPT_DIR/../target/release/grit" ]]; then
    echo -e "${YELLOW}Building GRIT release...${NC}"
    if (cd "$SCRIPT_DIR/.." && cargo build --release); then
      echo -e "${GREEN}Build successful${NC}"
    else
      echo -e "${RED}Build failed${NC}"
      ok=false
    fi
  fi

  if $ok; then
    echo -e "${GREEN}Environment OK${NC}"
  else
    echo -e "${RED}Environment check failed${NC}"
    exit 1
  fi
}

#############################################################################
# Real-World Dataset Management
#############################################################################

# Available real-world datasets and commands
REAL_DATASETS="dbsnp encode_peaks gencode sv"
REAL_COMMANDS="coverage intersect merge subtract closest"
REAL_DIR="$SCRIPT_DIR/real"

# Validate dataset name
validate_real_dataset() {
  local dataset="$1"
  for d in $REAL_DATASETS; do
    if [[ "$d" == "$dataset" ]]; then
      return 0
    fi
  done
  echo -e "${RED}ERROR: Unknown dataset '$dataset'${NC}" >&2
  echo "Available datasets: $REAL_DATASETS" >&2
  return 1
}

# Get source URL from source_url.txt
get_source_url() {
  local dataset="$1"
  local url_file="$REAL_DIR/$dataset/source_url.txt"
  if [[ -f "$url_file" ]]; then
    grep -v '^#' "$url_file" | grep -v '^$' | grep -E '^https?://' | head -1
  fi
}

# Read JSON field (simple)
json_get() {
  local file="$1"
  local key="$2"
  grep "\"$key\"" "$file" 2>/dev/null | sed 's/.*: *"\{0,1\}\([^",}]*\)"\{0,1\}.*/\1/' | head -1
}

# Update JSON field
json_set() {
  local file="$1"
  local key="$2"
  local value="$3"
  if [[ "$value" =~ ^[0-9]+$ ]] || [[ "$value" == "null" ]]; then
    sed -i.bak "s/\"$key\": *[^,}]*/\"$key\": $value/" "$file"
  else
    sed -i.bak "s|\"$key\": *[^,}]*|\"$key\": \"$value\"|" "$file"
  fi
  rm -f "${file}.bak"
}

# Check if dataset is downloaded
is_downloaded() {
  local dataset="$1"
  local raw_dir="$REAL_DIR/$dataset/raw"
  [[ -d "$raw_dir" && -n "$(ls -A "$raw_dir" 2>/dev/null)" ]]
}

# Check if dataset is prepared
is_prepared() {
  local dataset="$1"
  [[ -f "$REAL_DIR/$dataset/processed/data.sorted.bed" ]]
}

# File size in bytes (cross-platform)
file_size() {
  local file="$1"
  if [[ "$(uname)" == "Darwin" ]]; then
    stat -f%z "$file"
  else
    stat -c%s "$file"
  fi
}

#############################################################################
# real-download: Download real-world dataset
#############################################################################

cmd_real_download() {
  local dataset=""
  local yes_flag=false

  # Parse arguments
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --yes|-y)
        yes_flag=true
        shift
        ;;
      *)
        dataset="$1"
        shift
        ;;
    esac
  done

  if [[ -z "$dataset" ]]; then
    echo "Usage: $0 real-download <dataset> [--yes]"
    echo ""
    echo "Available datasets:"
    for d in $REAL_DATASETS; do
      local status="not downloaded"
      is_downloaded "$d" && status="downloaded"
      echo "  $d ($status)"
    done
    exit 1
  fi

  validate_real_dataset "$dataset" || exit 1

  local ds_dir="$REAL_DIR/$dataset"
  local raw_dir="$ds_dir/raw"
  local metadata="$ds_dir/metadata.json"
  local url=$(get_source_url "$dataset")

  # Check if URL is configured
  if [[ -z "$url" ]]; then
    echo -e "${RED}ERROR: No download URL configured for '$dataset'${NC}"
    echo ""
    echo "Configure the URL in: $ds_dir/source_url.txt"
    echo "Add a line starting with http:// or https://"
    exit 1
  fi

  # Check if already downloaded
  if is_downloaded "$dataset"; then
    echo -e "${YELLOW}Dataset '$dataset' already downloaded.${NC}"
    echo "Contents of $raw_dir:"
    ls -lh "$raw_dir"
    echo ""
    read -p "Re-download? [y/N] " confirm
    if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
      echo "Skipping download."
      return 0
    fi
    rm -rf "$raw_dir"/*
  fi

  # Confirmation prompt (unless --yes)
  if ! $yes_flag; then
    echo -e "${BOLD}Download Real-World Dataset${NC}"
    echo "═══════════════════════════════════════════════════════════════════════════════"
    echo "Dataset: $dataset"
    echo "URL: $url"
    echo "Target: $raw_dir/"
    echo ""
    echo -e "${YELLOW}WARNING: This will download data from an external source.${NC}"
    read -p "Continue? [y/N] " confirm
    if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
      echo "Download cancelled."
      exit 0
    fi
  fi

  mkdir -p "$raw_dir"

  # Extract filename from URL
  local filename=$(basename "$url")

  echo ""
  echo -n "Downloading $filename... "

  # Download using curl (prefer) or wget
  if command -v curl &>/dev/null; then
    if ! curl -L -f -o "$raw_dir/$filename" "$url" 2>/dev/null; then
      echo -e "${RED}FAILED${NC}"
      echo "Download failed. Check the URL and network connection."
      exit 1
    fi
  elif command -v wget &>/dev/null; then
    if ! wget -q -O "$raw_dir/$filename" "$url" 2>/dev/null; then
      echo -e "${RED}FAILED${NC}"
      echo "Download failed. Check the URL and network connection."
      exit 1
    fi
  else
    echo -e "${RED}FAILED${NC}"
    echo "ERROR: Neither curl nor wget found. Install one to download files."
    exit 1
  fi

  echo -e "${GREEN}done${NC}"

  # Verify checksum if provided
  local expected_sha256=$(json_get "$metadata" "expected_sha256")
  if [[ -n "$expected_sha256" && "$expected_sha256" != "null" ]]; then
    echo -n "Verifying checksum... "
    local actual_sha256=$(sha256_hash "$raw_dir/$filename")
    if [[ "$actual_sha256" == "$expected_sha256" ]]; then
      echo -e "${GREEN}PASS${NC}"
    else
      echo -e "${RED}FAIL${NC}"
      echo "Expected: $expected_sha256"
      echo "Actual:   $actual_sha256"
      exit 1
    fi
  fi

  # Update metadata
  json_set "$metadata" "filename" "$filename"

  local size=$(file_size "$raw_dir/$filename")
  echo ""
  echo -e "${GREEN}Download complete.${NC}"
  echo "File: $raw_dir/$filename ($(format_bytes $size))"
  echo ""
  echo "Next step: ./bench.sh real-prepare $dataset"
}

#############################################################################
# real-prepare: Prepare real-world dataset for benchmarking
#############################################################################

cmd_real_prepare() {
  if [[ $# -lt 1 ]]; then
    echo "Usage: $0 real-prepare <dataset>"
    echo ""
    echo "Available datasets:"
    for d in $REAL_DATASETS; do
      local status="not downloaded"
      is_downloaded "$d" && status="downloaded"
      is_prepared "$d" && status="prepared"
      echo "  $d ($status)"
    done
    exit 1
  fi

  local dataset="$1"
  validate_real_dataset "$dataset" || exit 1

  local ds_dir="$REAL_DIR/$dataset"
  local raw_dir="$ds_dir/raw"
  local proc_dir="$ds_dir/processed"
  local metadata="$ds_dir/metadata.json"

  # Check if downloaded
  if ! is_downloaded "$dataset"; then
    echo -e "${RED}ERROR: Dataset '$dataset' not downloaded.${NC}"
    echo "Run: ./bench.sh real-download $dataset"
    exit 1
  fi

  # Check if already prepared
  if is_prepared "$dataset"; then
    echo -e "${YELLOW}Dataset '$dataset' already prepared.${NC}"
    local bed_file="$proc_dir/data.sorted.bed"
    local n_intervals=$(wc -l < "$bed_file" | tr -d ' ')
    local size=$(file_size "$bed_file")
    echo "  Intervals: $n_intervals"
    echo "  Size: $(format_bytes $size)"
    echo ""
    read -p "Re-prepare? [y/N] " confirm
    if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
      return 0
    fi
  fi

  mkdir -p "$proc_dir"

  echo -e "${BOLD}Preparing Real-World Dataset${NC}"
  echo "═══════════════════════════════════════════════════════════════════════════════"
  echo "Dataset: $dataset"
  echo ""

  # Find raw file
  local raw_file=$(ls "$raw_dir"/* 2>/dev/null | head -1)
  if [[ -z "$raw_file" ]]; then
    echo -e "${RED}ERROR: No raw file found in $raw_dir${NC}"
    exit 1
  fi

  local tmp_dir=$(mktemp -d)
  local tmp_decompressed="$tmp_dir/decompressed"
  local tmp_bed="$tmp_dir/converted.bed"
  local tmp_sorted="$tmp_dir/sorted.bed"
  local output="$proc_dir/data.sorted.bed"
  local genome="$proc_dir/genome.txt"

  # Step 1: Decompress if needed
  echo -n "[1/4] Decompressing... "
  case "$raw_file" in
    *.gz)
      gunzip -c "$raw_file" > "$tmp_decompressed"
      echo "done (gzip)"
      ;;
    *.bz2)
      bunzip2 -c "$raw_file" > "$tmp_decompressed"
      echo "done (bzip2)"
      ;;
    *.xz)
      xz -dc "$raw_file" > "$tmp_decompressed"
      echo "done (xz)"
      ;;
    *.zip)
      unzip -p "$raw_file" > "$tmp_decompressed"
      echo "done (zip)"
      ;;
    *)
      cp "$raw_file" "$tmp_decompressed"
      echo "done (no compression)"
      ;;
  esac

  # Step 2: Detect format and convert to BED
  echo -n "[2/4] Converting to BED... "

  # Auto-detect format from content
  local format="bed"
  local first_line=$(head -1 "$tmp_decompressed")

  if [[ "$first_line" =~ ^##fileformat=VCF ]]; then
    format="vcf"
  elif [[ "$first_line" =~ ^##gff-version ]] || [[ "$first_line" =~ ^# ]] && grep -q '\tgene\t' "$tmp_decompressed"; then
    format="gtf"
  elif [[ "$raw_file" == *.gtf* ]] || [[ "$raw_file" == *.gff* ]]; then
    format="gtf"
  elif [[ "$raw_file" == *.vcf* ]]; then
    format="vcf"
  fi

  case "$format" in
    vcf)
      awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && NF>=2 {print $1, $2-1, $2}' "$tmp_decompressed" > "$tmp_bed"
      ;;
    gtf)
      awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && NF>=5 {print $1, $4-1, $5}' "$tmp_decompressed" > "$tmp_bed"
      ;;
    *)
      awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && !/^track/ && !/^browser/ && NF>=3 {print $1,$2,$3}' "$tmp_decompressed" > "$tmp_bed"
      ;;
  esac

  local n_raw=$(wc -l < "$tmp_bed" | tr -d ' ')
  echo "done ($format format, $n_raw intervals)"

  # Step 3: Create genome file and sort
  echo -n "[3/4] Sorting (genome order)... "
  create_genome "$genome"

  # Sort using bedtools if available, otherwise GNU sort
  if command -v bedtools &>/dev/null; then
    bedtools sort -i "$tmp_bed" -g "$genome" > "$tmp_sorted"
  else
    # GNU sort fallback
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n "$tmp_bed" > "$tmp_sorted"
  fi

  # Remove any lines with chromosomes not in genome file
  local valid_chroms=$(cut -f1 "$genome" | tr '\n' '|' | sed 's/|$//')
  grep -E "^($valid_chroms)\t" "$tmp_sorted" > "$output" || cp "$tmp_sorted" "$output"

  echo "done"

  # Step 4: Calculate statistics and update metadata
  echo -n "[4/4] Updating metadata... "

  local n_intervals=$(wc -l < "$output" | tr -d ' ')
  local file_size_val=$(file_size "$output")
  local sha256=$(sha256_hash "$output")
  local timestamp=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

  # Update metadata.json with processed file info
  # Using sed to update nested JSON (limited but works for this structure)
  sed -i.bak \
    -e "s|\"file\": *[^,}]*|\"file\": \"data.sorted.bed\"|" \
    -e "s|\"intervals\": *[^,}]*|\"intervals\": $n_intervals|" \
    -e "s|\"file_size_bytes\": *[^,}]*|\"file_size_bytes\": $file_size_val|" \
    -e "s|\"prepared_at\": *[^,}]*|\"prepared_at\": \"$timestamp\"|" \
    "$metadata"
  rm -f "${metadata}.bak"

  # Update SHA256 (need to handle quoted string)
  sed -i.bak "s|\"sha256\": *[^,}]*|\"sha256\": \"$sha256\"|" "$metadata"
  rm -f "${metadata}.bak"

  echo "done"

  # Cleanup
  rm -rf "$tmp_dir"

  echo ""
  echo "═══════════════════════════════════════════════════════════════════════════════"
  echo -e "${GREEN}Dataset prepared successfully.${NC}"
  echo ""
  echo "  Output:    $output"
  echo "  Intervals: $n_intervals"
  echo "  Size:      $(format_bytes $file_size_val)"
  echo "  SHA256:    $sha256"
  echo ""
  echo "Next step: ./bench.sh real-truth $dataset all"
}

#############################################################################
# real-truth: Generate bedtools truth baseline for real-world dataset
#############################################################################

cmd_real_truth() {
  local dataset=""
  local commands=""
  local force=false

  # Parse arguments
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --force)
        force=true
        shift
        ;;
      all)
        commands="$REAL_COMMANDS"
        shift
        ;;
      *)
        if [[ -z "$dataset" ]]; then
          dataset="$1"
        else
          commands="$commands $1"
        fi
        shift
        ;;
    esac
  done

  if [[ -z "$dataset" ]]; then
    echo "Usage: $0 real-truth <dataset> <command|all> [--force]"
    echo ""
    echo "Available datasets:"
    for d in $REAL_DATASETS; do
      local status="not prepared"
      is_prepared "$d" && status="prepared"
      echo "  $d ($status)"
    done
    echo ""
    echo "Commands: $REAL_COMMANDS"
    exit 1
  fi

  validate_real_dataset "$dataset" || exit 1

  # Default to all commands if none specified
  commands="${commands:-$REAL_COMMANDS}"
  commands=$(echo "$commands" | xargs)  # trim

  # Check if prepared
  if ! is_prepared "$dataset"; then
    echo -e "${RED}ERROR: Dataset '$dataset' not prepared.${NC}"
    echo "Run: ./bench.sh real-prepare $dataset"
    exit 1
  fi

  local ds_dir="$REAL_DIR/$dataset"
  local proc_dir="$ds_dir/processed"
  local truth_dir="$ds_dir/truth"
  local bed_file="$proc_dir/data.sorted.bed"
  local genome="$proc_dir/genome.txt"
  local metadata="$ds_dir/metadata.json"

  mkdir -p "$truth_dir"

  local time_cmd=$(detect_time_cmd)
  local n_intervals=$(wc -l < "$bed_file" | tr -d ' ')

  echo -e "${BOLD}Generating Truth Baseline${NC}"
  echo "═══════════════════════════════════════════════════════════════════════════════"
  echo "Dataset: $dataset ($n_intervals intervals)"
  echo "Commands: $commands"
  echo "Output: $truth_dir"
  echo ""

  for cmd in $commands; do
    local truth_base="$truth_dir/$cmd"

    # Check if already exists
    if [[ -f "${truth_base}.out" && -f "${truth_base}.sha256" && "$force" != "true" ]]; then
      echo -e "  ${YELLOW}$cmd${NC}: already exists (use --force to regenerate)"
      continue
    fi

    # Build bedtools command
    # For real-world datasets, we use the same file for -a and -b where applicable
    local bt_cmd
    case "$cmd" in
      coverage)
        bt_cmd="bedtools coverage -a $bed_file -b $bed_file -sorted"
        ;;
      intersect)
        bt_cmd="bedtools intersect -a $bed_file -b $bed_file -sorted"
        ;;
      merge)
        bt_cmd="bedtools merge -i $bed_file"
        ;;
      subtract)
        # Subtract from self produces empty output, so skip
        echo -e "  ${YELLOW}$cmd${NC}: skipped (self-subtract not meaningful)"
        continue
        ;;
      closest)
        bt_cmd="bedtools closest -a $bed_file -b $bed_file"
        ;;
      *)
        echo -e "  ${RED}$cmd${NC}: unknown command"
        continue
        ;;
    esac

    echo -n "  $cmd: running bedtools... "

    # Clear cache
    clear_cache "$bed_file" "$bed_file"

    # Run with timing
    local tmp_time=$(mktemp)
    eval "$time_cmd $bt_cmd" > "${truth_base}.out" 2> "$tmp_time"

    # Extract metrics
    extract_metrics "$tmp_time"
    echo "$TIME_SEC" > "${truth_base}.time"
    echo "$TIME_RSS" > "${truth_base}.rss"
    rm -f "$tmp_time"

    # Line count
    local lines=$(wc -l < "${truth_base}.out" | tr -d ' ')
    echo "$lines" > "${truth_base}.lines"

    # SHA256 (with canonicalization for coverage)
    if [[ "$cmd" == "coverage" ]]; then
      cut -f1-6 "${truth_base}.out" | sha256_hash_stdin > "${truth_base}.sha256"
    else
      sha256_hash "${truth_base}.out" > "${truth_base}.sha256"
    fi

    local rss_mb=$(awk "BEGIN {printf \"%.1f\", $TIME_RSS/1048576}")
    echo "done (${TIME_SEC}s, ${rss_mb}MB, $lines lines)"
  done

  echo ""
  echo -e "${GREEN}Truth baseline generated.${NC}"
  echo "Run benchmarks with: ./bench.sh real-run $dataset all"
}

#############################################################################
# real-run: Benchmark GRIT against bedtools on real-world dataset
#############################################################################

cmd_real_run() {
  local dataset=""
  local commands=""

  # Parse arguments
  while [[ $# -gt 0 ]]; do
    case "$1" in
      all)
        commands="$REAL_COMMANDS"
        shift
        ;;
      *)
        if [[ -z "$dataset" ]]; then
          dataset="$1"
        else
          commands="$commands $1"
        fi
        shift
        ;;
    esac
  done

  if [[ -z "$dataset" ]]; then
    echo "Usage: $0 real-run <dataset> <command|all>"
    echo ""
    echo "Available datasets:"
    for d in $REAL_DATASETS; do
      local status="not prepared"
      is_prepared "$d" && status="prepared"
      echo "  $d ($status)"
    done
    echo ""
    echo "Commands: $REAL_COMMANDS"
    exit 1
  fi

  validate_real_dataset "$dataset" || exit 1

  # Default to all commands
  commands="${commands:-$REAL_COMMANDS}"
  commands=$(echo "$commands" | xargs)

  # Check if prepared
  if ! is_prepared "$dataset"; then
    echo -e "${RED}ERROR: Dataset '$dataset' not prepared.${NC}"
    echo "Run: ./bench.sh real-prepare $dataset"
    exit 1
  fi

  local ds_dir="$REAL_DIR/$dataset"
  local proc_dir="$ds_dir/processed"
  local truth_dir="$ds_dir/truth"
  local bed_file="$proc_dir/data.sorted.bed"
  local genome="$proc_dir/genome.txt"

  local grit=$(detect_grit)
  local time_cmd=$(detect_time_cmd)
  local n_intervals=$(wc -l < "$bed_file" | tr -d ' ')
  local file_size_val=$(file_size "$bed_file")

  # Create results directory
  local results_dir="$RESULTS_DIR/real_${dataset}_$(date +%Y%m%d_%H%M%S)"
  mkdir -p "$results_dir"

  echo ""
  echo -e "${BOLD}GRIT Real-World Benchmark${NC}"
  echo "═══════════════════════════════════════════════════════════════════════════════"
  echo "Dataset: $dataset"
  echo "Intervals: $n_intervals"
  echo "File size: $(format_bytes $file_size_val)"
  echo "Commands: $commands"
  echo ""

  # Check for truth data
  local has_any_truth=false
  for cmd in $commands; do
    if [[ -f "$truth_dir/${cmd}.out" && -f "$truth_dir/${cmd}.sha256" ]]; then
      has_any_truth=true
      break
    fi
  done

  if $has_any_truth; then
    echo -e "Truth data: ${GREEN}available${NC} (commands marked * use cached bedtools output)"
  else
    echo -e "Truth data: ${YELLOW}not found${NC} (generate with: ./bench.sh real-truth $dataset all)"
  fi
  echo ""

  printf "%-12s %11s %11s %11s %13s %13s   %s\n" \
    "Command" "bedtools" "GRIT" "Speedup" "BT_RAM" "GRIT_RAM" "Status"
  echo "───────────────────────────────────────────────────────────────────────────────"

  local all_passed=true

  for cmd in $commands; do
    # Skip subtract (not meaningful for self-comparison)
    if [[ "$cmd" == "subtract" ]]; then
      continue
    fi

    local truth_base="$truth_dir/$cmd"
    local use_truth=false
    local bt_sec bt_rss bt_hash bt_lines

    # Check for truth data
    if [[ -f "${truth_base}.out" && -f "${truth_base}.sha256" ]]; then
      use_truth=true
      bt_sec=$(cat "${truth_base}.time")
      bt_rss=$(cat "${truth_base}.rss" 2>/dev/null || echo "0")
      bt_lines=$(cat "${truth_base}.lines" 2>/dev/null || wc -l < "${truth_base}.out" | tr -d ' ')
      bt_hash=$(cat "${truth_base}.sha256")
    fi

    # Build commands
    local bt_cmd gr_cmd
    case "$cmd" in
      coverage)
        bt_cmd="bedtools coverage -a $bed_file -b $bed_file -sorted"
        gr_cmd="$grit coverage -a $bed_file -b $bed_file --assume-sorted"
        ;;
      intersect)
        bt_cmd="bedtools intersect -a $bed_file -b $bed_file -sorted"
        gr_cmd="$grit intersect -a $bed_file -b $bed_file --assume-sorted --streaming"
        ;;
      merge)
        bt_cmd="bedtools merge -i $bed_file"
        gr_cmd="$grit merge -i $bed_file --assume-sorted"
        ;;
      closest)
        bt_cmd="bedtools closest -a $bed_file -b $bed_file"
        gr_cmd="$grit closest -a $bed_file -b $bed_file --assume-sorted --streaming"
        ;;
      *)
        continue
        ;;
    esac

    # Run bedtools if no truth data
    if ! $use_truth; then
      clear_cache "$bed_file" "$bed_file"

      local bt_out="$results_dir/${cmd}_bedtools.out"
      local bt_time_file="$results_dir/${cmd}_bedtools.time"

      eval "$time_cmd $bt_cmd" > "$bt_out" 2> "$bt_time_file"
      extract_metrics "$bt_time_file"
      bt_sec="$TIME_SEC"
      bt_rss="$TIME_RSS"
      bt_lines=$(wc -l < "$bt_out" | tr -d ' ')

      if [[ "$cmd" == "coverage" ]]; then
        bt_hash=$(cut -f1-6 "$bt_out" | sha256_hash_stdin)
      else
        bt_hash=$(sha256_hash "$bt_out")
      fi
    fi

    # Run GRIT
    clear_cache "$bed_file" "$bed_file"

    local gr_out="$results_dir/${cmd}_grit.out"
    local gr_time_file="$results_dir/${cmd}_grit.time"

    eval "$time_cmd $gr_cmd" > "$gr_out" 2> "$gr_time_file"
    extract_metrics "$gr_time_file"
    local gr_sec="$TIME_SEC"
    local gr_rss="$TIME_RSS"
    local gr_lines=$(wc -l < "$gr_out" | tr -d ' ')

    # Calculate GRIT hash
    local gr_hash
    if [[ "$cmd" == "coverage" ]]; then
      gr_hash=$(cut -f1-6 "$gr_out" | sha256_hash_stdin)
    else
      gr_hash=$(sha256_hash "$gr_out")
    fi

    # Check correctness
    local status="FAIL"
    if [[ "$bt_lines" == "$gr_lines" && "$bt_hash" == "$gr_hash" ]]; then
      status="PASS"
    fi

    # Calculate speedup
    local speedup="N/A"
    if [[ "$gr_sec" != "0" ]]; then
      speedup=$(awk "BEGIN {printf \"%.1f\", $bt_sec/$gr_sec}")
    fi

    # Calculate memory reduction
    local mem_reduction="N/A"
    if [[ "$gr_rss" != "0" && "$bt_rss" != "0" ]]; then
      mem_reduction=$(awk "BEGIN {printf \"%.1f\", $bt_rss/$gr_rss}")
    fi

    # Format output
    local bt_rss_fmt=$(format_bytes "$bt_rss")
    local gr_rss_fmt=$(format_bytes "$gr_rss")

    # Command display (with truth indicator)
    local cmd_display="$cmd"
    if $use_truth; then
      cmd_display="${cmd}*"
    fi

    # Colored status
    local status_colored
    if [[ "$status" == "PASS" ]]; then
      status_colored="${GREEN}PASS${NC}"
    else
      status_colored="${RED}FAIL${NC}"
      all_passed=false
    fi

    # Print row
    printf "%-12s %10.2fs %10.2fs %10sx %12s %12s   %b\n" \
      "$cmd_display" "$bt_sec" "$gr_sec" "$speedup" "$bt_rss_fmt" "$gr_rss_fmt" "$status_colored"

    # Calculate and report throughput
    if [[ "$status" == "PASS" && "$gr_sec" != "0" ]]; then
      local throughput=$(awk "BEGIN {printf \"%.0f\", $n_intervals/$gr_sec}")
      # Store extended metrics for later reporting
    fi
  done

  echo "═══════════════════════════════════════════════════════════════════════════════"
  echo "Results saved: $results_dir"

  if $all_passed; then
    echo -e "Status: ${GREEN}All benchmarks passed${NC}"
  else
    echo -e "Status: ${RED}Some benchmarks failed${NC}"
  fi

  # Print throughput summary
  echo ""
  echo "Throughput Summary:"
  for cmd in $commands; do
    [[ "$cmd" == "subtract" ]] && continue
    local gr_out="$results_dir/${cmd}_grit.out"
    local gr_time_file="$results_dir/${cmd}_grit.time"
    if [[ -f "$gr_time_file" ]]; then
      extract_metrics "$gr_time_file"
      if [[ "$TIME_SEC" != "0" ]]; then
        local throughput=$(awk "BEGIN {printf \"%.2f\", $n_intervals/$TIME_SEC}")
        echo "  $cmd: ${throughput} intervals/sec"
      fi
    fi
  done
}

#############################################################################
# real-list: List available real-world datasets
#############################################################################

cmd_real_list() {
  echo -e "${BOLD}Available Real-World Datasets${NC}"
  echo "═══════════════════════════════════════════════════════════════════════════════"
  echo ""

  for dataset in $REAL_DATASETS; do
    local ds_dir="$REAL_DIR/$dataset"
    local metadata="$ds_dir/metadata.json"

    # Read metadata
    local name=$(json_get "$metadata" "name")
    local desc=$(json_get "$metadata" "description")

    # Determine status
    local status="${RED}not configured${NC}"
    local details=""

    if [[ -f "$ds_dir/source_url.txt" ]]; then
      local url=$(get_source_url "$dataset")
      if [[ -n "$url" ]]; then
        status="${YELLOW}URL configured${NC}"
      fi
    fi

    if is_downloaded "$dataset"; then
      status="${YELLOW}downloaded${NC}"
    fi

    if is_prepared "$dataset"; then
      status="${GREEN}prepared${NC}"
      local bed_file="$ds_dir/processed/data.sorted.bed"
      if [[ -f "$bed_file" ]]; then
        local intervals=$(wc -l < "$bed_file" | tr -d ' ')
        local size=$(file_size "$bed_file")
        details="  ($intervals intervals, $(format_bytes $size))"
      fi
    fi

    echo -e "${CYAN}$dataset${NC}: $name"
    echo "  Description: $desc"
    echo -e "  Status: $status$details"

    # Check for truth data
    local truth_dir="$ds_dir/truth"
    if [[ -d "$truth_dir" ]]; then
      local truth_cmds=$(ls "$truth_dir"/*.out 2>/dev/null | xargs -n1 basename 2>/dev/null | sed 's/.out$//' | tr '\n' ' ')
      if [[ -n "$truth_cmds" ]]; then
        echo "  Truth data: $truth_cmds"
      fi
    fi
    echo ""
  done
}

#############################################################################
# Main
#############################################################################

case "${1:-}" in
  env)
    cmd_env
    ;;
  data)
    shift
    cmd_data "$@"
    ;;
  run)
    shift
    cmd_run "$@"
    ;;
  verify)
    shift
    "$SCRIPT_DIR/scripts/verify_all.sh" "$@"
    ;;
  bench)
    shift
    "$SCRIPT_DIR/scripts/bench_all.sh" "$@"
    ;;
  truth)
    shift
    cmd_truth "$@"
    ;;
  truth-list)
    shift
    cmd_truth_list "$@"
    ;;
  truth-verify)
    shift
    cmd_truth_verify "$@"
    ;;
  real-list)
    cmd_real_list
    ;;
  real-download)
    shift
    cmd_real_download "$@"
    ;;
  real-prepare)
    shift
    cmd_real_prepare "$@"
    ;;
  real-truth)
    shift
    cmd_real_truth "$@"
    ;;
  real-run)
    shift
    cmd_real_run "$@"
    ;;
  list)
    cmd_list
    ;;
  clean)
    cmd_clean
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    usage
    ;;
esac
