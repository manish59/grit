#!/usr/bin/env bash
set -euo pipefail

# Generate complete truth data baseline for GRIT benchmarking
# Usage: ./generate_all_truth.sh [--force] [--skip-large]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

FORCE_FLAG=""
SKIP_LARGE=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --force) FORCE_FLAG="--force"; shift ;;
    --skip-large) SKIP_LARGE=true; shift ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# All target sizes (A_size B_size)
SIZES=(
  "1M 500K"
  "2M 1M"
  "10M 5M"
  "50M 25M"
  "100M 50M"
)

DISTRIBUTIONS="uniform clustered"

echo "=============================================="
echo "GRIT Truth Data Generation"
echo "=============================================="
echo ""

total_start=$(date +%s)

for size_pair in "${SIZES[@]}"; do
  read a_size b_size <<< "$size_pair"
  size_name="${a_size}_${b_size}"

  # Skip large sizes if requested
  if $SKIP_LARGE; then
    case "$a_size" in
      50M|100M)
        echo "Skipping $size_name (--skip-large)"
        continue
        ;;
    esac
  fi

  echo ""
  echo "====== Processing $size_name ======"

  for dist in $DISTRIBUTIONS; do
    echo ""
    echo "--- $size_name / $dist ---"

    # Generate data
    ./bench.sh data "$a_size" "$b_size" "$dist"

    # Generate truth data for all commands
    ./bench.sh truth "$a_size" "$b_size" "$dist" all $FORCE_FLAG
  done
done

total_end=$(date +%s)
total_elapsed=$((total_end - total_start))

echo ""
echo "=============================================="
echo "Truth generation complete!"
echo "Total time: ${total_elapsed}s"
echo ""
echo "Verify with: ./bench.sh truth-list"
echo "Run benchmarks: ./bench.sh run <A_size> <B_size> all"
echo "=============================================="
