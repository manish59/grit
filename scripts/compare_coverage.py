#!/usr/bin/env python3
"""
Compare grit coverage output with bedtools coverage output.

Handles floating-point precision differences (1e-7) that occur due to
different rounding implementations between tools.

Usage:
    ./compare_coverage.py bedtools.out grit.out
    ./compare_coverage.py bedtools.out grit.out --tolerance 1e-6
"""

import sys
import argparse

def parse_line(line):
    """Parse a coverage output line, returning fields with fraction as float."""
    fields = line.rstrip('\n').split('\t')
    if len(fields) >= 10:
        # Standard coverage output: chrom, start, end, name, score, strand, count, bases, len, fraction
        try:
            fraction = float(fields[9])
            return fields[:9], fraction
        except (ValueError, IndexError):
            return fields, None
    return fields, None

def compare_files(file1, file2, tolerance=1e-6):
    """Compare two coverage output files with floating-point tolerance."""
    mismatches = 0
    line_num = 0

    with open(file1) as f1, open(file2) as f2:
        for line1, line2 in zip(f1, f2):
            line_num += 1

            fields1, frac1 = parse_line(line1)
            fields2, frac2 = parse_line(line2)

            # Compare non-fraction fields exactly
            if fields1 != fields2:
                print(f"Line {line_num}: Field mismatch")
                print(f"  File1: {fields1}")
                print(f"  File2: {fields2}")
                mismatches += 1
                continue

            # Compare fractions with tolerance
            if frac1 is not None and frac2 is not None:
                if abs(frac1 - frac2) > tolerance:
                    print(f"Line {line_num}: Fraction difference exceeds tolerance")
                    print(f"  File1: {frac1:.10f}")
                    print(f"  File2: {frac2:.10f}")
                    print(f"  Diff:  {abs(frac1 - frac2):.2e}")
                    mismatches += 1

    # Check for different line counts
    with open(file1) as f1, open(file2) as f2:
        count1 = sum(1 for _ in f1)
        count2 = sum(1 for _ in f2)
        if count1 != count2:
            print(f"Line count mismatch: {count1} vs {count2}")
            mismatches += 1

    return mismatches

def main():
    parser = argparse.ArgumentParser(description='Compare coverage outputs with floating-point tolerance')
    parser.add_argument('file1', help='First file (e.g., bedtools output)')
    parser.add_argument('file2', help='Second file (e.g., grit output)')
    parser.add_argument('--tolerance', '-t', type=float, default=1e-6,
                        help='Floating-point tolerance (default: 1e-6)')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Only print summary')

    args = parser.parse_args()

    mismatches = compare_files(args.file1, args.file2, args.tolerance)

    if mismatches == 0:
        print(f"✓ Files match (tolerance: {args.tolerance:.0e})")
        sys.exit(0)
    else:
        print(f"✗ {mismatches} differences found")
        sys.exit(1)

if __name__ == '__main__':
    main()
