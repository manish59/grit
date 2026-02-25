#!/bin/bash
# Bedtools-rs Validation Suite
# Compares bedtools-rs output against bedtools for bit-for-bit correctness

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
FIXTURES="$PROJECT_ROOT/tests/fixtures"
RESULTS="$PROJECT_ROOT/tests/validation/results"
BEDTOOLS_RS="$PROJECT_ROOT/target/release/rbedtools"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

PASSED=0
FAILED=0
SKIPPED=0

mkdir -p "$RESULTS"

# Build release binary
echo "Building bedtools-rs in release mode..."
cargo build --release --manifest-path "$PROJECT_ROOT/Cargo.toml" 2>/dev/null

log_pass() {
    echo -e "${GREEN}[PASS]${NC} $1"
    PASSED=$((PASSED + 1))
}

log_fail() {
    echo -e "${RED}[FAIL]${NC} $1"
    FAILED=$((FAILED + 1))
}

log_skip() {
    echo -e "${YELLOW}[SKIP]${NC} $1"
    SKIPPED=$((SKIPPED + 1))
}

# Compare two outputs, return 0 if identical
compare_output() {
    local name="$1"
    local bedtools_out="$2"
    local rs_out="$3"

    # Sort both outputs to handle order differences
    local bt_sorted=$(echo "$bedtools_out" | sort)
    local rs_sorted=$(echo "$rs_out" | sort)

    if [ "$bt_sorted" = "$rs_sorted" ]; then
        log_pass "$name"
        return 0
    else
        log_fail "$name"
        echo "  Expected (bedtools):"
        echo "$bedtools_out" | head -5 | sed 's/^/    /'
        echo "  Got (bedtools-rs):"
        echo "$rs_out" | head -5 | sed 's/^/    /'

        # Save diff to results
        echo "=== $name ===" >> "$RESULTS/failures.log"
        echo "bedtools output:" >> "$RESULTS/failures.log"
        echo "$bedtools_out" >> "$RESULTS/failures.log"
        echo "bedtools-rs output:" >> "$RESULTS/failures.log"
        echo "$rs_out" >> "$RESULTS/failures.log"
        echo "" >> "$RESULTS/failures.log"
        return 1
    fi
}

# Compare outputs preserving order (for sorted outputs)
compare_ordered() {
    local name="$1"
    local bedtools_out="$2"
    local rs_out="$3"

    if [ "$bedtools_out" = "$rs_out" ]; then
        log_pass "$name"
        return 0
    else
        log_fail "$name"
        echo "  Diff (first 10 lines):"
        diff <(echo "$bedtools_out") <(echo "$rs_out") | head -10 | sed 's/^/    /'
        return 1
    fi
}

echo ""
echo "=============================================="
echo "  PHASE 1: INTERSECT CORRECTNESS"
echo "=============================================="
echo ""

# Basic intersect
bt_out=$(bedtools intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
compare_output "intersect: basic" "$bt_out" "$rs_out"

# -wa flag
bt_out=$(bedtools intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -wa 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" --wa 2>/dev/null || true)
compare_output "intersect: -wa" "$bt_out" "$rs_out"

# -wb flag
bt_out=$(bedtools intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -wb 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" --wb 2>/dev/null || true)
compare_output "intersect: -wb" "$bt_out" "$rs_out"

# -wa -wb flags
bt_out=$(bedtools intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -wa -wb 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" --wa --wb 2>/dev/null || true)
compare_output "intersect: -wa -wb" "$bt_out" "$rs_out"

# -u flag (unique)
bt_out=$(bedtools intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -u 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -u 2>/dev/null || true)
compare_output "intersect: -u (unique)" "$bt_out" "$rs_out"

# -v flag (no overlap)
bt_out=$(bedtools intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -v 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -v 2>/dev/null || true)
compare_output "intersect: -v (no overlap)" "$bt_out" "$rs_out"

# -c flag (count)
bt_out=$(bedtools intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -c 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -c 2>/dev/null || true)
compare_output "intersect: -c (count)" "$bt_out" "$rs_out"

# -f flag (fraction overlap)
for frac in 0.1 0.25 0.5 0.75 0.9 1.0; do
    bt_out=$(bedtools intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -f "$frac" 2>/dev/null || true)
    rs_out=$("$BEDTOOLS_RS" intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -f "$frac" 2>/dev/null || true)
    compare_output "intersect: -f $frac" "$bt_out" "$rs_out"
done

# -f with -r (reciprocal)
for frac in 0.1 0.5 0.9; do
    bt_out=$(bedtools intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -f "$frac" -r 2>/dev/null || true)
    rs_out=$("$BEDTOOLS_RS" intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -f "$frac" -r 2>/dev/null || true)
    compare_output "intersect: -f $frac -r" "$bt_out" "$rs_out"
done

# Combinations
bt_out=$(bedtools intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -wa -u 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" --wa -u 2>/dev/null || true)
compare_output "intersect: -wa -u" "$bt_out" "$rs_out"

bt_out=$(bedtools intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -wa -c 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" --wa -c 2>/dev/null || true)
compare_output "intersect: -wa -c" "$bt_out" "$rs_out"

# Nested intervals
bt_out=$(bedtools intersect -a "$FIXTURES/nested.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" intersect -a "$FIXTURES/nested.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
compare_output "intersect: nested intervals" "$bt_out" "$rs_out"

# Self intersection
bt_out=$(bedtools intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/a.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" intersect -a "$FIXTURES/a.bed" -b "$FIXTURES/a.bed" 2>/dev/null || true)
compare_output "intersect: self" "$bt_out" "$rs_out"

echo ""
echo "=============================================="
echo "  PHASE 2: MERGE CORRECTNESS"
echo "=============================================="
echo ""

# Basic merge
bt_out=$(bedtools merge -i "$FIXTURES/a.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" merge -i "$FIXTURES/a.bed" 2>/dev/null || true)
compare_output "merge: basic" "$bt_out" "$rs_out"

# Merge with distance
for dist in 0 1 10 50 100 1000; do
    bt_out=$(bedtools merge -i "$FIXTURES/a.bed" -d "$dist" 2>/dev/null || true)
    rs_out=$("$BEDTOOLS_RS" merge -i "$FIXTURES/a.bed" -d "$dist" 2>/dev/null || true)
    compare_output "merge: -d $dist" "$bt_out" "$rs_out"
done

# Merge nested intervals
bt_out=$(bedtools merge -i "$FIXTURES/nested.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" merge -i "$FIXTURES/nested.bed" 2>/dev/null || true)
compare_output "merge: nested" "$bt_out" "$rs_out"

# Merge touching intervals
bt_out=$(bedtools merge -i "$FIXTURES/touching.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" merge -i "$FIXTURES/touching.bed" 2>/dev/null || true)
compare_output "merge: touching" "$bt_out" "$rs_out"

# Merge touching with d=1
bt_out=$(bedtools merge -i "$FIXTURES/touching.bed" -d 1 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" merge -i "$FIXTURES/touching.bed" -d 1 2>/dev/null || true)
compare_output "merge: touching -d 1" "$bt_out" "$rs_out"

# Merge unsorted (bedtools requires sorted, we handle unsorted)
sorted_unsorted=$(sort -k1,1 -k2,2n "$FIXTURES/unsorted.bed")
bt_out=$(echo "$sorted_unsorted" | bedtools merge 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" merge -i "$FIXTURES/unsorted.bed" 2>/dev/null || true)
compare_output "merge: unsorted input" "$bt_out" "$rs_out"

# Merge with strand
bt_out=$(bedtools merge -i "$FIXTURES/a.bed" -s 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" merge -i "$FIXTURES/a.bed" -s 2>/dev/null || true)
compare_output "merge: -s (strand)" "$bt_out" "$rs_out"

echo ""
echo "=============================================="
echo "  PHASE 3: SUBTRACT CORRECTNESS"
echo "=============================================="
echo ""

# Basic subtract
bt_out=$(bedtools subtract -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" subtract -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
compare_output "subtract: basic" "$bt_out" "$rs_out"

# Subtract with -A (remove entire)
bt_out=$(bedtools subtract -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -A 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" subtract -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -A 2>/dev/null || true)
compare_output "subtract: -A" "$bt_out" "$rs_out"

# Subtract with fraction
for frac in 0.1 0.5 0.9; do
    bt_out=$(bedtools subtract -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -f "$frac" 2>/dev/null || true)
    rs_out=$("$BEDTOOLS_RS" subtract -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -f "$frac" 2>/dev/null || true)
    compare_output "subtract: -f $frac" "$bt_out" "$rs_out"
done

# Subtract with complex overlaps
bt_out=$(bedtools subtract -a "$FIXTURES/subtract_cases.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" subtract -a "$FIXTURES/subtract_cases.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
compare_output "subtract: complex overlaps" "$bt_out" "$rs_out"

# Self subtraction (should yield nothing or fragments)
bt_out=$(bedtools subtract -a "$FIXTURES/a.bed" -b "$FIXTURES/a.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" subtract -a "$FIXTURES/a.bed" -b "$FIXTURES/a.bed" 2>/dev/null || true)
compare_output "subtract: self" "$bt_out" "$rs_out"

echo ""
echo "=============================================="
echo "  PHASE 3: CLOSEST CORRECTNESS"
echo "=============================================="
echo ""

# Basic closest
bt_out=$(bedtools closest -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" closest -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
compare_output "closest: basic" "$bt_out" "$rs_out"

# Closest with distance
bt_out=$(bedtools closest -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -d 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" closest -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -d 2>/dev/null || true)
compare_output "closest: -d" "$bt_out" "$rs_out"

# Closest tie handling
for tie in all first last; do
    bt_out=$(bedtools closest -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -t "$tie" 2>/dev/null || true)
    rs_out=$("$BEDTOOLS_RS" closest -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -t "$tie" 2>/dev/null || true)
    compare_output "closest: -t $tie" "$bt_out" "$rs_out"
done

# Closest ignore overlapping
bt_out=$(bedtools closest -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -io 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" closest -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" --io 2>/dev/null || true)
compare_output "closest: -io" "$bt_out" "$rs_out"

echo ""
echo "=============================================="
echo "  PHASE 4: COVERAGE CORRECTNESS"
echo "=============================================="
echo ""

# Basic coverage
bt_out=$(bedtools coverage -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" coverage -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
compare_output "coverage: basic" "$bt_out" "$rs_out"

# Coverage -hist
bt_out=$(bedtools coverage -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -hist 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" coverage -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" --hist 2>/dev/null || true)
compare_output "coverage: -hist" "$bt_out" "$rs_out"

# Coverage -mean
bt_out=$(bedtools coverage -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -mean 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" coverage -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" --mean 2>/dev/null || true)
compare_output "coverage: -mean" "$bt_out" "$rs_out"

# Coverage -d (per-base)
bt_out=$(bedtools coverage -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -d 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" coverage -a "$FIXTURES/a.bed" -b "$FIXTURES/b.bed" -d 2>/dev/null || true)
compare_output "coverage: -d (per-base)" "$bt_out" "$rs_out"

# Coverage with nested intervals
bt_out=$(bedtools coverage -a "$FIXTURES/nested.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" coverage -a "$FIXTURES/nested.bed" -b "$FIXTURES/b.bed" 2>/dev/null || true)
compare_output "coverage: nested intervals" "$bt_out" "$rs_out"

# Coverage self
bt_out=$(bedtools coverage -a "$FIXTURES/a.bed" -b "$FIXTURES/a.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" coverage -a "$FIXTURES/a.bed" -b "$FIXTURES/a.bed" 2>/dev/null || true)
compare_output "coverage: self" "$bt_out" "$rs_out"

echo ""
echo "=============================================="
echo "  PHASE 5: SORT CORRECTNESS"
echo "=============================================="
echo ""

# Basic sort
bt_out=$(bedtools sort -i "$FIXTURES/unsorted.bed" 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" sort -i "$FIXTURES/unsorted.bed" 2>/dev/null || true)
compare_ordered "sort: basic" "$bt_out" "$rs_out"

# Sort by size ascending
bt_out=$(bedtools sort -i "$FIXTURES/a.bed" -sizeA 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" sort -i "$FIXTURES/a.bed" --sizeA 2>/dev/null || true)
compare_ordered "sort: -sizeA" "$bt_out" "$rs_out"

# Sort by size descending
bt_out=$(bedtools sort -i "$FIXTURES/a.bed" -sizeD 2>/dev/null || true)
rs_out=$("$BEDTOOLS_RS" sort -i "$FIXTURES/a.bed" --sizeD 2>/dev/null || true)
compare_ordered "sort: -sizeD" "$bt_out" "$rs_out"

echo ""
echo "=============================================="
echo "  SUMMARY"
echo "=============================================="
echo ""
echo -e "Passed: ${GREEN}$PASSED${NC}"
echo -e "Failed: ${RED}$FAILED${NC}"
echo -e "Skipped: ${YELLOW}$SKIPPED${NC}"
echo ""

if [ $FAILED -gt 0 ]; then
    echo "Failures logged to: $RESULTS/failures.log"
    exit 1
else
    echo "All tests passed!"
    exit 0
fi
