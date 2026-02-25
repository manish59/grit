#!/usr/bin/env bash
set -uo pipefail
# Note: -e removed to handle SIGPIPE from pipes gracefully

#############################################################################
# validate_dataset.sh - Statistical validation of synthetic BED datasets
#
# Usage:
#   ./validate_dataset.sh <bed_file> <genome_sizes.txt>
#   ./validate_dataset.sh data/10M_5M/uniform/A.sorted.bed data/10M_5M/uniform/genome.txt
#
# Output:
#   Comprehensive validation report with distribution classification
#
# Requirements:
#   - bedtools (for overlap analysis)
#   - awk, sort, uniq (standard unix tools)
#   - Works on macOS and Linux
#############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#############################################################################
# Colors (disabled if not tty)
#############################################################################
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

#############################################################################
# Usage and argument parsing
#############################################################################
usage() {
  cat <<EOF
Usage: $0 <bed_file> <genome_sizes.txt> [options]

Options:
  --bins N        Number of bins for uniformity test (default: 100)
  --hist-bins N   Number of bins for length histogram (default: 10)
  --top-dense N   Number of densest regions to report (default: 10)
  --window W      Window size for density analysis (default: 1000000)
  --quiet         Suppress progress messages

Arguments:
  bed_file           Path to BED file to validate
  genome_sizes.txt   Tab-separated file: chrom<tab>size

Output:
  Validation report with distribution classification

Examples:
  $0 A.sorted.bed genome.txt
  $0 clustered.bed genome.txt --bins 50 --window 500000
EOF
  exit 1
}

if [[ $# -lt 2 ]]; then
  usage
fi

BED_FILE="$1"
GENOME_FILE="$2"
shift 2

# Default parameters
NUM_BINS=100
HIST_BINS=10
TOP_DENSE=10
WINDOW_SIZE=1000000
QUIET=false

# Parse optional arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bins) NUM_BINS="$2"; shift 2 ;;
    --hist-bins) HIST_BINS="$2"; shift 2 ;;
    --top-dense) TOP_DENSE="$2"; shift 2 ;;
    --window) WINDOW_SIZE="$2"; shift 2 ;;
    --quiet) QUIET=true; shift ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# Validate inputs
if [[ ! -f "$BED_FILE" ]]; then
  echo -e "${RED}ERROR: BED file not found: $BED_FILE${NC}" >&2
  exit 1
fi

if [[ ! -f "$GENOME_FILE" ]]; then
  echo -e "${RED}ERROR: Genome file not found: $GENOME_FILE${NC}" >&2
  exit 1
fi

# Create temp directory for intermediate files
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

log() {
  if ! $QUIET; then
    echo -e "$@" >&2
  fi
}

#############################################################################
# SECTION 1: Basic Statistics
#############################################################################
log "${BOLD}[1/6] Computing basic statistics...${NC}"

# Total intervals
TOTAL_INTERVALS=$(wc -l < "$BED_FILE" | tr -d ' ')

# Load genome sizes into associative structure via temp file
awk '{print $1, $2}' "$GENOME_FILE" > "$TEMP_DIR/genome_sizes.txt"
TOTAL_GENOME_SIZE=$(awk '{sum += $2} END {print sum}' "$GENOME_FILE")

#############################################################################
# SECTION 2: Chromosome Distribution Analysis
#############################################################################
log "${BOLD}[2/6] Analyzing chromosome distribution...${NC}"

# Count intervals per chromosome
cut -f1 "$BED_FILE" | sort | uniq -c | awk '{print $2, $1}' > "$TEMP_DIR/chrom_counts.txt"

# Calculate expected vs observed proportions
# Expected: proportional to chromosome size
# Observed: actual counts
awk -v total_intervals="$TOTAL_INTERVALS" -v total_genome="$TOTAL_GENOME_SIZE" '
  NR==FNR {
    genome_size[$1] = $2
    next
  }
  {
    chrom = $1
    observed = $2
    expected_frac = genome_size[chrom] / total_genome
    observed_frac = observed / total_intervals
    expected_count = expected_frac * total_intervals

    # Chi-square contribution
    if (expected_count > 0) {
      chi_sq = ((observed - expected_count)^2) / expected_count
    } else {
      chi_sq = 0
    }

    # Percent difference
    if (expected_frac > 0) {
      pct_diff = ((observed_frac - expected_frac) / expected_frac) * 100
    } else {
      pct_diff = 0
    }

    printf "%s\t%d\t%.4f\t%.4f\t%.2f\t%.4f\n", chrom, observed, expected_frac, observed_frac, pct_diff, chi_sq
  }
' "$TEMP_DIR/genome_sizes.txt" "$TEMP_DIR/chrom_counts.txt" > "$TEMP_DIR/chrom_analysis.txt"

# Total chi-square for chromosome distribution
CHROM_CHI_SQ=$(awk '{sum += $6} END {printf "%.2f", sum}' "$TEMP_DIR/chrom_analysis.txt")

# Degrees of freedom = num_chromosomes - 1
NUM_CHROMS=$(wc -l < "$TEMP_DIR/chrom_counts.txt" | tr -d ' ')
CHROM_DF=$((NUM_CHROMS - 1))

# Critical chi-square at p=0.05 (approximate lookup table)
# For df > 30, use normal approximation
get_chi_sq_critical() {
  local df=$1
  case $df in
    1) echo "3.84" ;;
    2) echo "5.99" ;;
    3) echo "7.81" ;;
    4) echo "9.49" ;;
    5) echo "11.07" ;;
    10) echo "18.31" ;;
    15) echo "25.00" ;;
    20) echo "31.41" ;;
    22) echo "33.92" ;;
    *)
      # Approximation for larger df
      awk -v df="$df" 'BEGIN {printf "%.2f", df + 2*sqrt(2*df)}'
      ;;
  esac
}

CHROM_CRITICAL=$(get_chi_sq_critical $CHROM_DF)

# Calculate proportionality score (0-100)
# 100 = perfect, 0 = completely wrong
CHROM_PROP_SCORE=$(awk -v chi="$CHROM_CHI_SQ" -v crit="$CHROM_CRITICAL" '
  BEGIN {
    ratio = chi / crit
    if (ratio < 1) {
      score = 100
    } else {
      score = 100 / (1 + log(ratio))
    }
    printf "%.1f", (score > 0 ? score : 0)
  }
')

#############################################################################
# SECTION 3: Start Position Uniformity Test
#############################################################################
log "${BOLD}[3/6] Testing start position uniformity...${NC}"

# For each chromosome, bin positions and compute variance
awk -v bins="$NUM_BINS" '
  NR==FNR {
    chrom_size[$1] = $2
    next
  }
  {
    chrom = $1
    start = $2
    if (chrom_size[chrom] > 0) {
      bin_size = chrom_size[chrom] / bins
      bin_idx = int(start / bin_size)
      if (bin_idx >= bins) bin_idx = bins - 1
      key = chrom ":" bin_idx
      counts[key]++
      chrom_counts[chrom]++
    }
  }
  END {
    total_chi_sq = 0
    total_df = 0

    for (chrom in chrom_counts) {
      n = chrom_counts[chrom]
      expected = n / bins
      chi_sq = 0

      for (i = 0; i < bins; i++) {
        key = chrom ":" i
        obs = (key in counts) ? counts[key] : 0
        if (expected > 0) {
          chi_sq += ((obs - expected)^2) / expected
        }
      }

      total_chi_sq += chi_sq
      total_df += (bins - 1)

      # Variance ratio (coefficient of variation)
      sum_sq = 0
      for (i = 0; i < bins; i++) {
        key = chrom ":" i
        obs = (key in counts) ? counts[key] : 0
        sum_sq += (obs - expected)^2
      }
      variance = sum_sq / bins
      cv = (expected > 0) ? sqrt(variance) / expected : 0

      printf "%s\t%d\t%.2f\t%.4f\n", chrom, n, chi_sq, cv
    }

    # Output totals on special line
    printf "TOTAL\t%d\t%.2f\t%.4f\n", total_chi_sq, total_df, 0
  }
' "$TEMP_DIR/genome_sizes.txt" "$BED_FILE" > "$TEMP_DIR/uniformity_analysis.txt"

# Extract total chi-square and df
UNIFORMITY_CHI_SQ=$(grep "^TOTAL" "$TEMP_DIR/uniformity_analysis.txt" | awk '{print $2}')
UNIFORMITY_DF=$(grep "^TOTAL" "$TEMP_DIR/uniformity_analysis.txt" | awk '{print $3}')

# Calculate uniformity score
UNIFORMITY_CRITICAL=$(awk -v df="$UNIFORMITY_DF" 'BEGIN {printf "%.2f", df + 2*sqrt(2*df)}')
UNIFORMITY_SCORE=$(awk -v chi="$UNIFORMITY_CHI_SQ" -v crit="$UNIFORMITY_CRITICAL" '
  BEGIN {
    ratio = chi / crit
    if (ratio < 1) {
      score = 100
    } else {
      score = 100 / (1 + log(ratio))
    }
    printf "%.1f", (score > 0 ? score : 0)
  }
')

# Average coefficient of variation (excluding TOTAL line)
AVG_CV=$(grep -v "^TOTAL" "$TEMP_DIR/uniformity_analysis.txt" | \
  awk '{sum += $4; n++} END {printf "%.4f", (n > 0 ? sum/n : 0)}')

#############################################################################
# SECTION 4: Interval Length Distribution
#############################################################################
log "${BOLD}[4/6] Analyzing interval length distribution...${NC}"

# Calculate length stats
awk '{len = $3 - $2; print len}' "$BED_FILE" | sort -n > "$TEMP_DIR/lengths.txt"

LENGTH_STATS=$(awk '
  {
    lengths[NR] = $1
    sum += $1
    n++
  }
  END {
    # Min, max
    min = lengths[1]
    max = lengths[n]

    # Mean
    mean = sum / n

    # Median
    if (n % 2 == 1) {
      median = lengths[int(n/2) + 1]
    } else {
      median = (lengths[n/2] + lengths[n/2 + 1]) / 2
    }

    # Variance and stddev
    sum_sq = 0
    for (i = 1; i <= n; i++) {
      sum_sq += (lengths[i] - mean)^2
    }
    variance = sum_sq / n
    stddev = sqrt(variance)

    # Skewness (simplified)
    sum_cube = 0
    for (i = 1; i <= n; i++) {
      sum_cube += ((lengths[i] - mean) / stddev)^3
    }
    skewness = sum_cube / n

    printf "%d\t%d\t%.2f\t%.2f\t%.2f\t%.4f", min, max, mean, median, stddev, skewness
  }
' "$TEMP_DIR/lengths.txt")

LEN_MIN=$(echo "$LENGTH_STATS" | cut -f1)
LEN_MAX=$(echo "$LENGTH_STATS" | cut -f2)
LEN_MEAN=$(echo "$LENGTH_STATS" | cut -f3)
LEN_MEDIAN=$(echo "$LENGTH_STATS" | cut -f4)
LEN_STDDEV=$(echo "$LENGTH_STATS" | cut -f5)
LEN_SKEWNESS=$(echo "$LENGTH_STATS" | cut -f6)

# Create histogram
LEN_RANGE=$((LEN_MAX - LEN_MIN))
if [[ $LEN_RANGE -eq 0 ]]; then
  LEN_RANGE=1
fi

awk -v min="$LEN_MIN" -v range="$LEN_RANGE" -v bins="$HIST_BINS" '
  {
    len = $1
    bin = int((len - min) * bins / range)
    if (bin >= bins) bin = bins - 1
    if (bin < 0) bin = 0
    counts[bin]++
    n++
  }
  END {
    bin_size = range / bins
    for (i = 0; i < bins; i++) {
      lower = min + i * bin_size
      upper = lower + bin_size
      count = (i in counts) ? counts[i] : 0
      pct = (n > 0) ? count * 100 / n : 0
      printf "[%d-%d)\t%d\t%.2f%%\n", lower, upper, count, pct
    }
  }
' "$TEMP_DIR/lengths.txt" > "$TEMP_DIR/length_histogram.txt"

# Determine if length distribution is uniform or skewed
# For uniform: skewness should be near 0, CV should be ~0.29 for U(a,b)
LEN_CV=$(awk -v mean="$LEN_MEAN" -v stddev="$LEN_STDDEV" 'BEGIN {
  cv = (mean > 0) ? stddev / mean : 0
  printf "%.4f", cv
}')

# Expected CV for uniform distribution U(a,b) = (b-a)/(sqrt(12)*((a+b)/2)) ≈ 0.289 * (range/mean)
EXPECTED_UNIFORM_CV=$(awk -v range="$LEN_RANGE" -v mean="$LEN_MEAN" 'BEGIN {
  if (mean > 0) {
    expected = range / (sqrt(12) * mean)
    printf "%.4f", expected
  } else {
    print "0"
  }
}')

LEN_UNIFORMITY=$(awk -v cv="$LEN_CV" -v expected="$EXPECTED_UNIFORM_CV" -v skew="$LEN_SKEWNESS" '
  BEGIN {
    # Check if CV is close to expected uniform CV
    cv_diff = (expected > 0) ? (cv - expected) / expected : 0

    # For uniform distribution: |skewness| < 0.5 and CV close to expected
    if (cv_diff > -0.3 && cv_diff < 0.3 && skew > -0.5 && skew < 0.5) {
      print "UNIFORM"
    } else if (skew > 0.5) {
      print "RIGHT_SKEWED"
    } else if (skew < -0.5) {
      print "LEFT_SKEWED"
    } else {
      print "NON_UNIFORM"
    }
  }
')

#############################################################################
# SECTION 5: Overlap Density Analysis
#############################################################################
log "${BOLD}[5/6] Computing overlap density...${NC}"

# Check if bedtools is available
if ! command -v bedtools &>/dev/null; then
  log "${YELLOW}WARNING: bedtools not found, skipping overlap analysis${NC}"
  AVG_OVERLAPS="N/A"
  MAX_DEPTH="N/A"
  OVERLAP_STATS="bedtools_not_available"
else
  # Self-intersection to count overlaps per interval
  # Use bedtools intersect -c for count
  bedtools intersect -a "$BED_FILE" -b "$BED_FILE" -c -sorted 2>/dev/null > "$TEMP_DIR/overlap_counts.txt" || \
  bedtools intersect -a "$BED_FILE" -b "$BED_FILE" -c > "$TEMP_DIR/overlap_counts.txt"

  # Each interval overlaps itself, so subtract 1
  OVERLAP_STATS=$(awk '
    {
      overlaps = $NF - 1  # Subtract self-overlap
      if (overlaps < 0) overlaps = 0
      sum += overlaps
      n++
      if (overlaps > max) max = overlaps
    }
    END {
      avg = (n > 0) ? sum / n : 0
      printf "%.2f\t%d", avg, max
    }
  ' "$TEMP_DIR/overlap_counts.txt")

  AVG_OVERLAPS=$(echo "$OVERLAP_STATS" | cut -f1)
  MAX_DEPTH=$(echo "$OVERLAP_STATS" | cut -f2)

  # Calculate genome coverage using bedtools genomecov
  bedtools genomecov -i "$BED_FILE" -g "$GENOME_FILE" -max 1 2>/dev/null > "$TEMP_DIR/genomecov.txt" || true

  if [[ -s "$TEMP_DIR/genomecov.txt" ]]; then
    # Get fraction of genome covered
    GENOME_COVERAGE=$(awk '$1 == "genome" && $2 == 1 {print $5}' "$TEMP_DIR/genomecov.txt")
    if [[ -z "$GENOME_COVERAGE" ]]; then
      GENOME_COVERAGE="0"
    fi
  else
    GENOME_COVERAGE="N/A"
  fi
fi

#############################################################################
# SECTION 6: Cluster/Hotspot Detection
#############################################################################
log "${BOLD}[6/6] Detecting cluster hotspots...${NC}"

# Create genome windows and count intervals per window
if command -v bedtools &>/dev/null; then
  # Create windows
  bedtools makewindows -g "$GENOME_FILE" -w "$WINDOW_SIZE" 2>/dev/null > "$TEMP_DIR/windows.bed" || true

  if [[ -s "$TEMP_DIR/windows.bed" ]]; then
    # Count intervals per window
    bedtools intersect -a "$TEMP_DIR/windows.bed" -b "$BED_FILE" -c -sorted 2>/dev/null > "$TEMP_DIR/window_counts.txt" || \
    bedtools intersect -a "$TEMP_DIR/windows.bed" -b "$BED_FILE" -c > "$TEMP_DIR/window_counts.txt"

    # Calculate density statistics
    DENSITY_STATS=$(awk '
      {
        count = $4
        sum += count
        sum_sq += count^2
        n++
        if (count > max) max = count
        counts[NR] = count
      }
      END {
        mean = sum / n
        variance = (sum_sq / n) - (mean^2)
        stddev = sqrt(variance > 0 ? variance : 0)
        cv = (mean > 0) ? stddev / mean : 0

        # Sort counts for percentiles
        for (i = 1; i <= n; i++) {
          for (j = i+1; j <= n; j++) {
            if (counts[i] > counts[j]) {
              tmp = counts[i]
              counts[i] = counts[j]
              counts[j] = tmp
            }
          }
        }

        p50 = counts[int(n*0.5)]
        p90 = counts[int(n*0.9)]
        p99 = counts[int(n*0.99)]

        # Cluster ratio: p99/median indicates how much denser top 1% is
        cluster_ratio = (p50 > 0) ? p99 / p50 : 0

        printf "%.2f\t%.2f\t%.4f\t%d\t%d\t%d\t%.2f", mean, stddev, cv, p50, p90, p99, cluster_ratio
      }
    ' "$TEMP_DIR/window_counts.txt")

    DENSITY_MEAN=$(echo "$DENSITY_STATS" | cut -f1)
    DENSITY_STDDEV=$(echo "$DENSITY_STATS" | cut -f2)
    DENSITY_CV=$(echo "$DENSITY_STATS" | cut -f3)
    DENSITY_P50=$(echo "$DENSITY_STATS" | cut -f4)
    DENSITY_P90=$(echo "$DENSITY_STATS" | cut -f5)
    DENSITY_P99=$(echo "$DENSITY_STATS" | cut -f6)
    CLUSTER_RATIO=$(echo "$DENSITY_STATS" | cut -f7)

    # Get top dense regions
    sort -k4,4nr "$TEMP_DIR/window_counts.txt" | head -n "$TOP_DENSE" > "$TEMP_DIR/top_dense.txt"
  else
    DENSITY_MEAN="N/A"
    DENSITY_CV="N/A"
    CLUSTER_RATIO="N/A"
  fi
else
  DENSITY_MEAN="N/A"
  DENSITY_CV="N/A"
  CLUSTER_RATIO="N/A"
fi

#############################################################################
# SECTION 7: Distribution Classification
#############################################################################

classify_distribution() {
  local chrom_score="$1"
  local uniformity_score="$2"
  local cluster_ratio="$3"
  local avg_cv="$4"

  # Check for clustering FIRST (takes precedence over chromosome bias)
  # Clustered data naturally has chromosome bias due to hotspot placement
  if [[ "$cluster_ratio" != "N/A" ]]; then
    # High cluster ratio (p99/p50 > 5) suggests clustering
    if (( $(echo "$cluster_ratio > 5" | bc -l) )); then
      # Also check uniformity - clustered data should have poor uniformity
      if (( $(echo "$uniformity_score < 50" | bc -l) )); then
        echo "LIKELY_CLUSTERED"
        return
      fi
    fi
  fi

  # Check chromosome bias (only if not clustered)
  if (( $(echo "$chrom_score < 50" | bc -l) )); then
    echo "SKEWED_CHROMOSOME_BIAS"
    return
  fi

  # Check uniformity
  if (( $(echo "$uniformity_score > 80" | bc -l) )); then
    echo "LIKELY_UNIFORM"
    return
  elif (( $(echo "$uniformity_score > 50" | bc -l) )); then
    echo "POSSIBLY_UNIFORM"
    return
  fi

  echo "INCORRECT_DISTRIBUTION"
}

if [[ "$CLUSTER_RATIO" != "N/A" ]]; then
  CLASSIFICATION=$(classify_distribution "$CHROM_PROP_SCORE" "$UNIFORMITY_SCORE" "$CLUSTER_RATIO" "$AVG_CV")
else
  # Without cluster ratio, use simpler classification
  if (( $(echo "$CHROM_PROP_SCORE < 50" | bc -l) )); then
    CLASSIFICATION="SKEWED_CHROMOSOME_BIAS"
  elif (( $(echo "$UNIFORMITY_SCORE > 80" | bc -l) )); then
    CLASSIFICATION="LIKELY_UNIFORM"
  elif (( $(echo "$UNIFORMITY_SCORE > 50" | bc -l) )); then
    CLASSIFICATION="POSSIBLY_UNIFORM"
  else
    CLASSIFICATION="INCORRECT_DISTRIBUTION"
  fi
fi

#############################################################################
# SECTION 8: Generate Report
#############################################################################

echo ""
echo "═══════════════════════════════════════════════════════════════════════════════"
echo "                         DATASET VALIDATION REPORT"
echo "═══════════════════════════════════════════════════════════════════════════════"
echo ""
echo "Input file:       $BED_FILE"
echo "Genome file:      $GENOME_FILE"
echo "Analysis date:    $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "───────────────────────────────────────────────────────────────────────────────"
echo "                           BASIC STATISTICS"
echo "───────────────────────────────────────────────────────────────────────────────"
echo ""
printf "Total intervals:              %'d\n" "$TOTAL_INTERVALS"
printf "Total genome size:            %'d bp\n" "$TOTAL_GENOME_SIZE"
printf "Number of chromosomes:        %d\n" "$NUM_CHROMS"
echo ""
echo "───────────────────────────────────────────────────────────────────────────────"
echo "                       CHROMOSOME DISTRIBUTION"
echo "───────────────────────────────────────────────────────────────────────────────"
echo ""
printf "%-8s %12s %12s %12s %12s\n" "Chrom" "Count" "Expected%" "Observed%" "Diff%"
echo "────────────────────────────────────────────────────────────────────"
head -10 "$TEMP_DIR/chrom_analysis.txt" | while IFS=$'\t' read -r chrom count exp_frac obs_frac diff chi; do
  printf "%-8s %12d %11.2f%% %11.2f%% %+11.2f%%\n" "$chrom" "$count" "$(echo "$exp_frac * 100" | bc -l)" "$(echo "$obs_frac * 100" | bc -l)" "$diff"
done
echo "..."
echo ""
printf "Chi-square statistic:         %.2f (df=%d, critical=%.2f)\n" "$CHROM_CHI_SQ" "$CHROM_DF" "$CHROM_CRITICAL"
printf "Proportionality score:        %.1f / 100\n" "$CHROM_PROP_SCORE"
echo ""
echo "───────────────────────────────────────────────────────────────────────────────"
echo "                       POSITION UNIFORMITY TEST"
echo "───────────────────────────────────────────────────────────────────────────────"
echo ""
printf "Number of bins per chromosome: %d\n" "$NUM_BINS"
printf "Chi-square statistic:         %.2f (df=%.0f)\n" "$UNIFORMITY_CHI_SQ" "$UNIFORMITY_DF"
printf "Average CV across chromosomes: %.4f\n" "$AVG_CV"
printf "Uniformity score:             %.1f / 100\n" "$UNIFORMITY_SCORE"
echo ""
echo "───────────────────────────────────────────────────────────────────────────────"
echo "                     INTERVAL LENGTH DISTRIBUTION"
echo "───────────────────────────────────────────────────────────────────────────────"
echo ""
printf "Minimum length:               %d bp\n" "$LEN_MIN"
printf "Maximum length:               %d bp\n" "$LEN_MAX"
printf "Mean length:                  %.2f bp\n" "$LEN_MEAN"
printf "Median length:                %.2f bp\n" "$LEN_MEDIAN"
printf "Std deviation:                %.2f bp\n" "$LEN_STDDEV"
printf "Skewness:                     %.4f\n" "$LEN_SKEWNESS"
printf "Coefficient of variation:     %.4f (expected uniform: %.4f)\n" "$LEN_CV" "$EXPECTED_UNIFORM_CV"
printf "Length distribution type:     %s\n" "$LEN_UNIFORMITY"
echo ""
echo "Histogram:"
cat "$TEMP_DIR/length_histogram.txt" | while IFS=$'\t' read -r range count pct; do
  bar_len=$(echo "$pct" | awk '{printf "%d", $1/2}')
  bar=$(printf '%*s' "$bar_len" '' | tr ' ' '#')
  printf "  %-20s %8s %6s %s\n" "$range" "$count" "$pct" "$bar"
done
echo ""
echo "───────────────────────────────────────────────────────────────────────────────"
echo "                          OVERLAP DENSITY"
echo "───────────────────────────────────────────────────────────────────────────────"
echo ""
if [[ "$AVG_OVERLAPS" != "N/A" ]]; then
  printf "Average overlaps per interval: %.2f\n" "$AVG_OVERLAPS"
  printf "Maximum overlap depth:         %d\n" "$MAX_DEPTH"
  if [[ "$GENOME_COVERAGE" != "N/A" ]]; then
    printf "Genome coverage:              %.4f%%\n" "$(echo "$GENOME_COVERAGE * 100" | bc -l)"
  fi
else
  echo "Overlap analysis: skipped (bedtools not available)"
fi
echo ""
echo "───────────────────────────────────────────────────────────────────────────────"
echo "                         CLUSTER DETECTION"
echo "───────────────────────────────────────────────────────────────────────────────"
echo ""
if [[ "$DENSITY_MEAN" != "N/A" ]]; then
  printf "Window size:                  %'d bp\n" "$WINDOW_SIZE"
  printf "Mean intervals per window:    %.2f\n" "$DENSITY_MEAN"
  printf "Std deviation:                %.2f\n" "$DENSITY_STDDEV"
  printf "Coefficient of variation:     %.4f\n" "$DENSITY_CV"
  printf "Median (P50):                 %d\n" "$DENSITY_P50"
  printf "90th percentile (P90):        %d\n" "$DENSITY_P90"
  printf "99th percentile (P99):        %d\n" "$DENSITY_P99"
  printf "Cluster density ratio (P99/P50): %.2f\n" "$CLUSTER_RATIO"
  echo ""
  echo "Top $TOP_DENSE densest regions:"
  printf "  %-8s %15s %15s %8s\n" "Chrom" "Start" "End" "Count"
  head -n "$TOP_DENSE" "$TEMP_DIR/top_dense.txt" | while IFS=$'\t' read -r chrom start end count; do
    printf "  %-8s %15d %15d %8d\n" "$chrom" "$start" "$end" "$count"
  done
else
  echo "Cluster detection: skipped (bedtools not available)"
fi
echo ""
echo "═══════════════════════════════════════════════════════════════════════════════"
echo "                              SUMMARY"
echo "═══════════════════════════════════════════════════════════════════════════════"
echo ""
printf "Chromosome proportionality score: %.1f / 100\n" "$CHROM_PROP_SCORE"
printf "Uniformity score:                 %.1f / 100\n" "$UNIFORMITY_SCORE"
if [[ "$CLUSTER_RATIO" != "N/A" ]]; then
  printf "Cluster density ratio:            %.2f\n" "$CLUSTER_RATIO"
fi
printf "Overlap depth stats:              avg=%.2f, max=%s\n" "${AVG_OVERLAPS:-0}" "${MAX_DEPTH:-N/A}"
printf "Length distribution summary:      %s (range: %d-%d, mean: %.1f)\n" "$LEN_UNIFORMITY" "$LEN_MIN" "$LEN_MAX" "$LEN_MEAN"
echo ""
echo "───────────────────────────────────────────────────────────────────────────────"
echo ""

# Color the classification
case "$CLASSIFICATION" in
  LIKELY_UNIFORM)
    echo -e "Final classification:            ${GREEN}${BOLD}$CLASSIFICATION${NC}"
    ;;
  LIKELY_CLUSTERED)
    echo -e "Final classification:            ${CYAN}${BOLD}$CLASSIFICATION${NC}"
    ;;
  POSSIBLY_UNIFORM)
    echo -e "Final classification:            ${YELLOW}${BOLD}$CLASSIFICATION${NC}"
    ;;
  SKEWED_CHROMOSOME_BIAS|INCORRECT_DISTRIBUTION)
    echo -e "Final classification:            ${RED}${BOLD}$CLASSIFICATION${NC}"
    ;;
  *)
    echo "Final classification:            $CLASSIFICATION"
    ;;
esac

echo ""
echo "═══════════════════════════════════════════════════════════════════════════════"

# Exit with appropriate code
case "$CLASSIFICATION" in
  LIKELY_UNIFORM|LIKELY_CLUSTERED)
    exit 0
    ;;
  POSSIBLY_UNIFORM)
    exit 0
    ;;
  *)
    exit 1
    ;;
esac
