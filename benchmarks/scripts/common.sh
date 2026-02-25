#!/usr/bin/env bash
#############################################
# GRIT Benchmark Common Functions
# Source this file in individual benchmark scripts
#############################################

# Default values
DEFAULT_A_INTERVALS=1000000
DEFAULT_B_INTERVALS=500000
DEFAULT_SINGLE_INTERVALS=1000000
DEFAULT_LEN_A=150
DEFAULT_LEN_B=500
DEFAULT_SEED=42

# Benchmark directory (set early for all functions)
BENCH_DIR="${BENCH_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)}"

# GRIT binary detection
detect_grit_binary() {
  if [[ -n "${GRIT_BIN:-}" && -x "$GRIT_BIN" ]]; then
    echo "$GRIT_BIN"
  elif [[ -x "$HOME/Projects/bedtools-rs/target/release/grit" ]]; then
    echo "$HOME/Projects/bedtools-rs/target/release/grit"
  elif command -v grit >/dev/null 2>&1; then
    command -v grit
  else
    echo "ERROR: GRIT binary not found." >&2
    echo "Build with: cargo build --release" >&2
    exit 1
  fi
}

# Platform-aware time command detection
detect_time_cmd() {
  if [[ "$(uname)" == "Darwin" ]]; then
    # macOS: use /usr/bin/time -l for RSS in bytes
    echo "/usr/bin/time -l"
  else
    # Linux: use /usr/bin/time -v for RSS in KB
    echo "/usr/bin/time -v"
  fi
}

# Extract timing metrics from time output
# Usage: extract_time_metrics <time_file>
# Sets: TIME_REAL, TIME_RSS_BYTES
extract_time_metrics() {
  local time_file="$1"

  if [[ "$(uname)" == "Darwin" ]]; then
    # macOS format: "X.XX real" and "NNNN maximum resident set size"
    TIME_REAL=$(grep " real" "$time_file" | awk '{print $1}')
    TIME_RSS_BYTES=$(grep "maximum resident set size" "$time_file" | awk '{print $1}')
  else
    # Linux format: "Elapsed (wall clock) time: X:XX.XX" and "Maximum resident set size (kbytes): NNNN"
    TIME_REAL=$(grep "Elapsed" "$time_file" | sed 's/.*: //' | awk -F: '{if (NF==3) print $1*3600+$2*60+$3; else if (NF==2) print $1*60+$2; else print $1}')
    local rss_kb=$(grep "Maximum resident set size" "$time_file" | awk '{print $NF}')
    TIME_RSS_BYTES=$((rss_kb * 1024))
  fi
}

# Create genome file for hg38
create_genome_file() {
  local genome_file="${1:-genome.txt}"

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
  echo "$genome_file"
}

# Generate BED file using bedtools random
# Usage: generate_bed_random <output_file> <num_intervals> <interval_length> <genome_file> [seed]
generate_bed_random() {
  local output_file="$1"
  local num_intervals="$2"
  local interval_length="$3"
  local genome_file="$4"
  local seed="${5:-$DEFAULT_SEED}"

  bedtools random -l "$interval_length" -n "$num_intervals" -seed "$seed" -g "$genome_file" > "$output_file"
}

# Generate BED file using AWK (no bedtools dependency)
# Usage: generate_bed_awk <output_file> <num_intervals> [min_length] [max_length] [seed]
generate_bed_awk() {
  local output_file="$1"
  local num_intervals="$2"
  local min_length="${3:-50}"
  local max_length="${4:-500}"
  local seed="${5:-$DEFAULT_SEED}"

  awk -v N="$num_intervals" -v MIN="$min_length" -v MAX="$max_length" -v SEED="$seed" 'BEGIN {
    srand(SEED);
    for (i = 0; i < N; i++) {
      chr = "chr" int(1 + rand()*22);
      start = int(rand() * 100000000);
      len = int(MIN + rand()*(MAX-MIN));
      end = start + len;
      print chr "\t" start "\t" end;
    }
  }' > "$output_file"
}

# GNU sort BED file
# Usage: sort_bed <input_file> <output_file>
sort_bed() {
  local input_file="$1"
  local output_file="$2"

  LC_ALL=C sort -k1,1 -k2,2n "$input_file" > "$output_file"
}

# Compare SHA256 hashes
# Usage: compare_sha256 <file1> <file2>
# Returns: 0 if match, 1 if mismatch
compare_sha256() {
  local file1="$1"
  local file2="$2"

  local hash1=$(shasum -a 256 "$file1" | awk '{print $1}')
  local hash2=$(shasum -a 256 "$file2" | awk '{print $1}')

  if [[ "$hash1" == "$hash2" ]]; then
    echo "PASS"
    return 0
  else
    echo "MISMATCH"
    return 1
  fi
}

# Compare line counts
# Usage: compare_line_count <file1> <file2>
compare_line_count() {
  local file1="$1"
  local file2="$2"

  local count1=$(wc -l < "$file1")
  local count2=$(wc -l < "$file2")

  if [[ "$count1" == "$count2" ]]; then
    echo "PASS ($count1 lines)"
    return 0
  else
    echo "MISMATCH (bedtools: $count1, grit: $count2)"
    return 1
  fi
}

# Calculate throughput
# Usage: calc_throughput <intervals> <time_seconds>
calc_throughput() {
  local intervals="$1"
  local time_seconds="$2"

  awk -v n="$intervals" -v t="$time_seconds" 'BEGIN { printf "%.2f", n/t }'
}

# Calculate bytes per interval
# Usage: calc_bytes_per_interval <bytes> <intervals>
calc_bytes_per_interval() {
  local bytes="$1"
  local intervals="$2"

  echo $((bytes / intervals))
}

# Print benchmark header
# Usage: print_header <benchmark_name> <description>
print_header() {
  local name="$1"
  local desc="$2"

  echo "=========================================="
  echo "$name"
  echo "$desc"
  echo "=========================================="
}

# Print results table
# Usage: print_results <bt_time> <bt_throughput> <bt_rss> <gr_time> <gr_throughput> <gr_rss> <correctness>
print_results() {
  local bt_time="$1"
  local bt_thr="$2"
  local bt_rss="$3"
  local gr_time="$4"
  local gr_thr="$5"
  local gr_rss="$6"
  local correct="$7"

  echo ""
  echo "============= RESULTS ============="
  printf "%-15s %-10s %-15s %-15s\n" "Tool" "Time(s)" "Throughput/s" "Max RSS(bytes)"
  printf "%-15s %-10s %-15s %-15s\n" "bedtools" "$bt_time" "$bt_thr" "$bt_rss"
  printf "%-15s %-10s %-15s %-15s\n" "grit" "$gr_time" "$gr_thr" "$gr_rss"
  echo ""
  echo "Correctness: $correct"
  echo "===================================="
}

# Print extended results table (with bytes per interval)
# Usage: print_results_extended <bt_time> <bt_throughput> <bt_rss> <bt_bpi> <gr_time> <gr_throughput> <gr_rss> <gr_bpi> <correctness>
print_results_extended() {
  local bt_time="$1"
  local bt_thr="$2"
  local bt_rss="$3"
  local bt_bpi="$4"
  local gr_time="$5"
  local gr_thr="$6"
  local gr_rss="$7"
  local gr_bpi="$8"
  local correct="$9"

  echo ""
  echo "============= RESULTS ============="
  printf "%-12s %-10s %-15s %-15s %-15s\n" "Tool" "Time(s)" "Throughput/s" "Max RSS(bytes)" "Bytes/interval"
  printf "%-12s %-10s %-15s %-15s %-15s\n" "bedtools" "$bt_time" "$bt_thr" "$bt_rss" "$bt_bpi"
  printf "%-12s %-10s %-15s %-15s %-15s\n" "grit" "$gr_time" "$gr_thr" "$gr_rss" "$gr_bpi"
  echo ""
  echo "Correctness: $correct"
  echo "===================================="
}

# Create results directory with timestamp
# Usage: create_results_dir <prefix>
create_results_dir() {
  local prefix="$1"
  local timestamp=$(date +%Y%m%d_%H%M%S)
  local dir="${prefix}_${timestamp}"
  mkdir -p "$dir"
  echo "$dir"
}

# Run timed command and capture output
# Usage: run_timed <output_file> <time_file> <command...>
run_timed() {
  local output_file="$1"
  local time_file="$2"
  shift 2

  local time_cmd=$(detect_time_cmd)

  $time_cmd "$@" > "$output_file" 2> "$time_file"
}

# Run timed command with hash output (no file written)
# Usage: run_timed_hash <hash_file> <time_file> <command...>
run_timed_hash() {
  local hash_file="$1"
  local time_file="$2"
  shift 2

  local time_cmd=$(detect_time_cmd)

  $time_cmd sh -c "$* | shasum -a 256" > "$hash_file" 2> "$time_file"
}

#############################################
# Truth Data Management Functions
#
# Truth data is pre-computed bedtools output stored
# for faster benchmarking iterations.
#############################################

# Get the truth directory for a given size and distribution
# Usage: truth_dir <size_spec> <distribution>
# Example: truth_dir "10M_5M" "uniform"
truth_dir() {
  local size_spec="$1"
  local distribution="$2"
  local bench_dir="${BENCH_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
  echo "$bench_dir/truth/${size_spec}/${distribution}"
}

# Get the path to a truth data file (without extension)
# Usage: truth_path <size_spec> <distribution> <command>
truth_path() {
  local size_spec="$1"
  local distribution="$2"
  local command="$3"
  echo "$(truth_dir "$size_spec" "$distribution")/${command}"
}

# Check if truth data exists for a command
# Usage: has_truth <size_spec> <distribution> <command>
# Returns: 0 if exists, 1 if not
has_truth() {
  local size_spec="$1"
  local distribution="$2"
  local command="$3"
  local path=$(truth_path "$size_spec" "$distribution" "$command")
  [[ -f "${path}.out" && -f "${path}.time" && -f "${path}.sha256" ]]
}

# Get truth timing (seconds)
# Usage: get_truth_time <size_spec> <distribution> <command>
get_truth_time() {
  local size_spec="$1"
  local distribution="$2"
  local command="$3"
  local path=$(truth_path "$size_spec" "$distribution" "$command")
  if [[ -f "${path}.time" ]]; then
    cat "${path}.time"
  else
    echo "0"
  fi
}

# Get truth RSS (MB)
# Usage: get_truth_rss <size_spec> <distribution> <command>
get_truth_rss() {
  local size_spec="$1"
  local distribution="$2"
  local command="$3"
  local path=$(truth_path "$size_spec" "$distribution" "$command")
  if [[ -f "${path}.rss" ]]; then
    cat "${path}.rss"
  else
    echo "0"
  fi
}

# Get truth output file path
# Usage: get_truth_output <size_spec> <distribution> <command>
get_truth_output() {
  local size_spec="$1"
  local distribution="$2"
  local command="$3"
  echo "$(truth_path "$size_spec" "$distribution" "$command").out"
}

# Get truth SHA256 hash
# Usage: get_truth_sha256 <size_spec> <distribution> <command>
get_truth_sha256() {
  local size_spec="$1"
  local distribution="$2"
  local command="$3"
  local path=$(truth_path "$size_spec" "$distribution" "$command")
  if [[ -f "${path}.sha256" ]]; then
    cat "${path}.sha256"
  else
    echo ""
  fi
}

# Get truth line count
# Usage: get_truth_lines <size_spec> <distribution> <command>
get_truth_lines() {
  local size_spec="$1"
  local distribution="$2"
  local command="$3"
  local path=$(truth_path "$size_spec" "$distribution" "$command")
  if [[ -f "${path}.lines" ]]; then
    cat "${path}.lines"
  else
    echo "0"
  fi
}

# Platform-aware SHA256 for files
sha256_file() {
  if command -v sha256sum &>/dev/null; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

# Platform-aware SHA256 for stdin
sha256_stdin() {
  if command -v sha256sum &>/dev/null; then
    sha256sum | awk '{print $1}'
  else
    shasum -a 256 | awk '{print $1}'
  fi
}

# Verify truth data integrity
# Usage: verify_truth <size_spec> <distribution> <command>
# Returns: 0 if valid, 1 if invalid
verify_truth() {
  local size_spec="$1"
  local distribution="$2"
  local command="$3"
  local path=$(truth_path "$size_spec" "$distribution" "$command")

  if ! has_truth "$size_spec" "$distribution" "$command"; then
    return 1
  fi

  local stored_hash=$(cat "${path}.sha256")
  local actual_hash=$(sha256_file "${path}.out")

  [[ "$stored_hash" == "$actual_hash" ]]
}

# List all available truth data
# Usage: list_truth [size_spec]
list_truth() {
  local size_spec="${1:-}"
  local bench_dir="${BENCH_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
  local truth_base="$bench_dir/truth"

  if [[ ! -d "$truth_base" ]]; then
    echo "No truth data found."
    return 0
  fi

  if [[ -n "$size_spec" ]]; then
    # List for specific size
    if [[ -d "$truth_base/$size_spec" ]]; then
      find "$truth_base/$size_spec" -name "*.out" -type f | while read -r f; do
        local rel_path="${f#$truth_base/}"
        local dir=$(dirname "$rel_path")
        local cmd=$(basename "$f" .out)
        echo "$dir/$cmd"
      done | sort
    fi
  else
    # List all
    find "$truth_base" -name "*.out" -type f 2>/dev/null | while read -r f; do
      local rel_path="${f#$truth_base/}"
      local dir=$(dirname "$rel_path")
      local cmd=$(basename "$f" .out)
      echo "$dir/$cmd"
    done | sort
  fi
}

#############################################
# Real-World Dataset Management Functions
#
# Functions for downloading, preparing, and
# benchmarking real-world genomic datasets.
#############################################

# Available real-world datasets
REAL_DATASETS="dbsnp encode_peaks gencode sv"

# Commands supported for real-world benchmarks
REAL_COMMANDS="coverage intersect merge subtract closest"

# Get the real dataset directory
# Usage: real_dataset_dir <dataset>
real_dataset_dir() {
  local dataset="$1"
  local bench_dir="${BENCH_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
  echo "$bench_dir/real/${dataset}"
}

# Validate dataset name
# Usage: validate_real_dataset <dataset>
# Returns: 0 if valid, 1 if invalid
validate_real_dataset() {
  local dataset="$1"
  local valid=false
  for d in $REAL_DATASETS; do
    if [[ "$d" == "$dataset" ]]; then
      valid=true
      break
    fi
  done
  $valid
}

# Get URL from source_url.txt (first non-comment, non-empty line starting with http)
# Usage: get_real_source_url <dataset>
get_real_source_url() {
  local dataset="$1"
  local url_file="$(real_dataset_dir "$dataset")/source_url.txt"
  if [[ -f "$url_file" ]]; then
    grep -v '^#' "$url_file" | grep -v '^$' | grep -E '^https?://' | head -1
  fi
}

# Read JSON field value (simple parser for flat JSON)
# Usage: json_get <file> <key>
json_get() {
  local file="$1"
  local key="$2"
  grep "\"$key\"" "$file" 2>/dev/null | sed 's/.*: *"\{0,1\}\([^",}]*\)"\{0,1\}.*/\1/' | head -1
}

# Read nested JSON field value
# Usage: json_get_nested <file> <parent_key> <child_key>
json_get_nested() {
  local file="$1"
  local parent="$2"
  local child="$3"
  awk -v p="\"$parent\"" -v c="\"$child\"" '
    $0 ~ p { in_parent=1 }
    in_parent && $0 ~ c {
      gsub(/.*:[ ]*"?/, "");
      gsub(/"?[ ]*,?[ ]*$/, "");
      print;
      exit
    }
    in_parent && /\}/ { in_parent=0 }
  ' "$file"
}

# Update JSON field (creates backup)
# Usage: json_set <file> <key> <value>
json_set() {
  local file="$1"
  local key="$2"
  local value="$3"

  if [[ "$value" =~ ^[0-9]+$ ]] || [[ "$value" == "null" ]]; then
    # Numeric or null value - no quotes
    sed -i.bak "s/\"$key\": *[^,}]*/\"$key\": $value/" "$file"
  else
    # String value - with quotes
    sed -i.bak "s|\"$key\": *[^,}]*|\"$key\": \"$value\"|" "$file"
  fi
  rm -f "${file}.bak"
}

# Check if dataset is downloaded (raw file exists)
# Usage: is_real_downloaded <dataset>
is_real_downloaded() {
  local dataset="$1"
  local raw_dir="$(real_dataset_dir "$dataset")/raw"
  [[ -d "$raw_dir" && -n "$(ls -A "$raw_dir" 2>/dev/null)" ]]
}

# Check if dataset is prepared (processed BED exists)
# Usage: is_real_prepared <dataset>
is_real_prepared() {
  local dataset="$1"
  local proc_dir="$(real_dataset_dir "$dataset")/processed"
  [[ -f "$proc_dir/data.sorted.bed" ]]
}

# Get processed BED file path
# Usage: get_real_bed <dataset>
get_real_bed() {
  local dataset="$1"
  echo "$(real_dataset_dir "$dataset")/processed/data.sorted.bed"
}

# Get real dataset metadata file path
# Usage: get_real_metadata <dataset>
get_real_metadata() {
  local dataset="$1"
  echo "$(real_dataset_dir "$dataset")/metadata.json"
}

# Verify SHA256 checksum of a file
# Usage: verify_checksum <file> <expected_hash>
# Returns: 0 if match, 1 if mismatch
verify_checksum() {
  local file="$1"
  local expected="$2"

  if [[ -z "$expected" || "$expected" == "null" ]]; then
    return 0  # No checksum to verify
  fi

  local actual=$(sha256_file "$file")
  [[ "$actual" == "$expected" ]]
}

# Decompress file based on extension
# Usage: decompress_file <input_file> <output_file>
decompress_file() {
  local input="$1"
  local output="$2"

  case "$input" in
    *.gz)
      gunzip -c "$input" > "$output"
      ;;
    *.bz2)
      bunzip2 -c "$input" > "$output"
      ;;
    *.xz)
      xz -dc "$input" > "$output"
      ;;
    *.zip)
      unzip -p "$input" > "$output"
      ;;
    *)
      cp "$input" "$output"
      ;;
  esac
}

# Convert various formats to BED (3+ column)
# Usage: convert_to_bed <input_file> <output_file> <format>
# Formats: bed, vcf, gtf, gff, narrowpeak, broadpeak
convert_to_bed() {
  local input="$1"
  local output="$2"
  local format="${3:-bed}"

  case "$format" in
    bed|narrowpeak|broadpeak)
      # Already BED-like: extract first 3 columns minimum
      awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && !/^track/ && !/^browser/ && NF>=3 {print $1,$2,$3}' "$input" > "$output"
      ;;
    vcf)
      # VCF: chr, pos-1, pos (0-based BED)
      awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && NF>=2 {print $1, $2-1, $2}' "$input" > "$output"
      ;;
    gtf|gff|gff3)
      # GTF/GFF: chr, start-1, end (convert 1-based to 0-based)
      awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && NF>=5 {print $1, $4-1, $5}' "$input" > "$output"
      ;;
    ucsc)
      # UCSC table format: skip header, extract chrom, chromStart, chromEnd
      awk -F'\t' 'BEGIN{OFS="\t"} NR>1 && NF>=4 {print $2,$3,$4}' "$input" > "$output"
      ;;
    *)
      echo "ERROR: Unknown format '$format'" >&2
      return 1
      ;;
  esac
}

# Sort BED file in genome order (cross-platform deterministic)
# Usage: sort_bed_genome <input_file> <output_file> [genome_file]
sort_bed_genome() {
  local input="$1"
  local output="$2"
  local genome="${3:-}"

  if [[ -n "$genome" ]] && command -v bedtools &>/dev/null; then
    # Use bedtools sort with genome file for proper chromosome order
    bedtools sort -i "$input" -g "$genome" > "$output"
  else
    # Fallback: GNU sort with locale-independent ordering
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n "$input" > "$output"
  fi
}

# Get file size in bytes (cross-platform)
# Usage: file_size_bytes <file>
file_size_bytes() {
  local file="$1"
  if [[ "$(uname)" == "Darwin" ]]; then
    stat -f%z "$file"
  else
    stat -c%s "$file"
  fi
}

# Count intervals in BED file
# Usage: count_intervals <file>
count_intervals() {
  local file="$1"
  wc -l < "$file" | tr -d ' '
}

# Get current timestamp in ISO format
# Usage: timestamp_iso
timestamp_iso() {
  date -u +"%Y-%m-%dT%H:%M:%SZ"
}

# Real dataset truth directory
# Usage: real_truth_dir <dataset>
real_truth_dir() {
  local dataset="$1"
  echo "$(real_dataset_dir "$dataset")/truth"
}

# Check if real dataset has truth data for a command
# Usage: has_real_truth <dataset> <command>
has_real_truth() {
  local dataset="$1"
  local command="$2"
  local truth_dir="$(real_truth_dir "$dataset")"
  [[ -f "${truth_dir}/${command}.out" && -f "${truth_dir}/${command}.sha256" ]]
}

# Get real truth file path (without extension)
# Usage: real_truth_path <dataset> <command>
real_truth_path() {
  local dataset="$1"
  local command="$2"
  echo "$(real_truth_dir "$dataset")/${command}"
}

# List available real datasets
# Usage: list_real_datasets
list_real_datasets() {
  local bench_dir="${BENCH_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
  local real_dir="$bench_dir/real"

  if [[ ! -d "$real_dir" ]]; then
    echo "No real datasets configured."
    return 0
  fi

  for dataset in $REAL_DATASETS; do
    local ds_dir="$real_dir/$dataset"
    if [[ -d "$ds_dir" ]]; then
      local status=""
      if is_real_prepared "$dataset"; then
        status="[prepared]"
      elif is_real_downloaded "$dataset"; then
        status="[downloaded]"
      else
        status="[not downloaded]"
      fi
      echo "$dataset $status"
    fi
  done
}
