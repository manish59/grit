//! Tests for zero-length interval semantics.
//!
//! These tests verify that GRIT correctly handles zero-length intervals
//! in both strict mode (default) and bedtools-compatible mode.
//!
//! Note: Tests are run serially to avoid global config race conditions.

use grit::bed::{parse_intervals, BedReader, FastBedParser};
use grit::config;
use serial_test::serial;

/// Reset config to default state before each test
fn reset_config() {
    config::set_bedtools_compatible(false);
}

// =============================================================================
// Test 1: Strict Mode Default - Zero-length intervals remain zero-length
// =============================================================================

#[test]
#[serial]
fn test_strict_mode_zero_length_no_overlap() {
    reset_config();

    // In strict mode, zero-length intervals should not overlap with themselves
    let content = "chr1\t100\t100\nchr1\t100\t101\n";
    let intervals = parse_intervals(content).unwrap();

    assert_eq!(intervals.len(), 2);

    // First interval: zero-length (start == end)
    assert_eq!(intervals[0].start, 100);
    assert_eq!(intervals[0].end, 100);
    assert_eq!(intervals[0].len(), 0);
    assert!(intervals[0].is_empty());

    // Second interval: 1bp interval
    assert_eq!(intervals[1].start, 100);
    assert_eq!(intervals[1].end, 101);
    assert_eq!(intervals[1].len(), 1);
    assert!(!intervals[1].is_empty());

    // Zero-length interval should NOT overlap with the 1bp interval
    // because the zero-length interval [100,100) contains no bases
    assert!(!intervals[0].overlaps(&intervals[1]));

    // Zero-length interval should NOT overlap with itself
    assert!(!intervals[0].overlaps(&intervals[0]));
}

#[test]
#[serial]
fn test_strict_mode_preserves_zero_length() {
    reset_config();

    let content = "chr1\t50\t50\nchr2\t100\t100\nchr3\t200\t200\n";
    let intervals = parse_intervals(content).unwrap();

    assert_eq!(intervals.len(), 3);
    for interval in &intervals {
        assert_eq!(interval.len(), 0);
        assert!(interval.is_empty());
    }
}

// =============================================================================
// Test 2: Bedtools-Compatible Mode - Zero-length intervals become 1bp
// =============================================================================

#[test]
#[serial]
fn test_bedtools_compatible_zero_length_overlap() {
    reset_config();
    config::set_bedtools_compatible(true);

    let content = "chr1\t100\t100\nchr1\t100\t101\n";
    let intervals = parse_intervals(content).unwrap();

    assert_eq!(intervals.len(), 2);

    // First interval: was zero-length, now normalized to 1bp
    assert_eq!(intervals[0].start, 100);
    assert_eq!(intervals[0].end, 101); // Normalized!
    assert_eq!(intervals[0].len(), 1);
    assert!(!intervals[0].is_empty());

    // Second interval: unchanged
    assert_eq!(intervals[1].start, 100);
    assert_eq!(intervals[1].end, 101);

    // Both intervals are now [100,101) and should overlap
    assert!(intervals[0].overlaps(&intervals[1]));

    // Each interval should overlap with itself
    assert!(intervals[0].overlaps(&intervals[0]));
    assert!(intervals[1].overlaps(&intervals[1]));

    reset_config();
}

#[test]
#[serial]
fn test_bedtools_compatible_normalizes_all_zero_length() {
    reset_config();
    config::set_bedtools_compatible(true);

    let content = "chr1\t50\t50\nchr2\t100\t100\nchr3\t200\t200\n";
    let intervals = parse_intervals(content).unwrap();

    assert_eq!(intervals.len(), 3);
    for interval in &intervals {
        assert_eq!(interval.len(), 1); // All normalized to 1bp
        assert!(!interval.is_empty());
    }

    // Verify specific values
    assert_eq!(intervals[0].end, 51);
    assert_eq!(intervals[1].end, 101);
    assert_eq!(intervals[2].end, 201);

    reset_config();
}

// =============================================================================
// Test 3: Multiple zero-length intervals across chromosomes
// =============================================================================

#[test]
#[serial]
fn test_strict_mode_multiple_chroms() {
    reset_config();

    let content = "\
chr1\t100\t100
chr1\t200\t200
chr2\t50\t50
chr2\t100\t150
chrX\t0\t0
";
    let intervals = parse_intervals(content).unwrap();

    assert_eq!(intervals.len(), 5);

    // Zero-length intervals stay zero-length
    assert_eq!(intervals[0].len(), 0); // chr1:100-100
    assert_eq!(intervals[1].len(), 0); // chr1:200-200
    assert_eq!(intervals[2].len(), 0); // chr2:50-50
    assert_eq!(intervals[3].len(), 50); // chr2:100-150 (non-zero)
    assert_eq!(intervals[4].len(), 0); // chrX:0-0

    // Non-zero interval should not be affected
    assert_eq!(intervals[3].start, 100);
    assert_eq!(intervals[3].end, 150);
}

#[test]
#[serial]
fn test_bedtools_compatible_multiple_chroms() {
    reset_config();
    config::set_bedtools_compatible(true);

    let content = "\
chr1\t100\t100
chr1\t200\t200
chr2\t50\t50
chr2\t100\t150
chrX\t0\t0
";
    let intervals = parse_intervals(content).unwrap();

    assert_eq!(intervals.len(), 5);

    // Zero-length intervals are normalized
    assert_eq!(intervals[0].len(), 1); // chr1:100-101
    assert_eq!(intervals[1].len(), 1); // chr1:200-201
    assert_eq!(intervals[2].len(), 1); // chr2:50-51
    assert_eq!(intervals[3].len(), 50); // chr2:100-150 (unchanged)
    assert_eq!(intervals[4].len(), 1); // chrX:0-1

    reset_config();
}

// =============================================================================
// Test 4: Non-zero intervals unaffected in both modes
// =============================================================================

#[test]
#[serial]
fn test_nonzero_intervals_unchanged_strict() {
    reset_config();

    let content = "\
chr1\t100\t200
chr1\t0\t1
chr2\t500\t1000
chr3\t1\t2
";
    let intervals = parse_intervals(content).unwrap();

    assert_eq!(intervals.len(), 4);
    assert_eq!(intervals[0].start, 100);
    assert_eq!(intervals[0].end, 200);
    assert_eq!(intervals[1].start, 0);
    assert_eq!(intervals[1].end, 1);
    assert_eq!(intervals[2].start, 500);
    assert_eq!(intervals[2].end, 1000);
    assert_eq!(intervals[3].start, 1);
    assert_eq!(intervals[3].end, 2);
}

#[test]
#[serial]
fn test_nonzero_intervals_unchanged_compatible() {
    reset_config();
    config::set_bedtools_compatible(true);

    let content = "\
chr1\t100\t200
chr1\t0\t1
chr2\t500\t1000
chr3\t1\t2
";
    let intervals = parse_intervals(content).unwrap();

    assert_eq!(intervals.len(), 4);
    // All non-zero intervals remain unchanged
    assert_eq!(intervals[0].start, 100);
    assert_eq!(intervals[0].end, 200);
    assert_eq!(intervals[1].start, 0);
    assert_eq!(intervals[1].end, 1);
    assert_eq!(intervals[2].start, 500);
    assert_eq!(intervals[2].end, 1000);
    assert_eq!(intervals[3].start, 1);
    assert_eq!(intervals[3].end, 2);

    reset_config();
}

// =============================================================================
// Test 5: FastBedParser respects config
// =============================================================================

#[test]
#[serial]
fn test_fast_parser_strict_mode() {
    reset_config();

    let parser = FastBedParser::new();

    let line1 = b"chr1\t100\t100";
    let interval1 = parser.parse_interval(line1).unwrap();
    assert_eq!(interval1.start, 100);
    assert_eq!(interval1.end, 100); // Not normalized

    let line2 = b"chr1\t100\t200";
    let interval2 = parser.parse_interval(line2).unwrap();
    assert_eq!(interval2.start, 100);
    assert_eq!(interval2.end, 200);
}

#[test]
#[serial]
fn test_fast_parser_compatible_mode() {
    reset_config();
    config::set_bedtools_compatible(true);

    let parser = FastBedParser::new();

    let line1 = b"chr1\t100\t100";
    let interval1 = parser.parse_interval(line1).unwrap();
    assert_eq!(interval1.start, 100);
    assert_eq!(interval1.end, 101); // Normalized!

    let line2 = b"chr1\t100\t200";
    let interval2 = parser.parse_interval(line2).unwrap();
    assert_eq!(interval2.start, 100);
    assert_eq!(interval2.end, 200); // Unchanged

    reset_config();
}

// =============================================================================
// Test 6: BedReader respects config
// =============================================================================

#[test]
#[serial]
fn test_bed_reader_strict_mode() {
    reset_config();

    let content = b"chr1\t100\t100\tname\t500\t+\n";
    let reader = BedReader::new(&content[..]);
    let records: Vec<_> = reader.records().collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].interval.start, 100);
    assert_eq!(records[0].interval.end, 100); // Not normalized
}

#[test]
#[serial]
fn test_bed_reader_compatible_mode() {
    reset_config();
    config::set_bedtools_compatible(true);

    let content = b"chr1\t100\t100\tname\t500\t+\n";
    let reader = BedReader::new(&content[..]);
    let records: Vec<_> = reader.records().collect::<Result<_, _>>().unwrap();

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].interval.start, 100);
    assert_eq!(records[0].interval.end, 101); // Normalized!

    reset_config();
}

// =============================================================================
// Test 7: Config module functions
// =============================================================================

#[test]
#[serial]
fn test_config_normalize_end() {
    reset_config();

    // Strict mode: no normalization
    assert_eq!(config::normalize_end(100, 100), 100);
    assert_eq!(config::normalize_end(100, 200), 200);

    // Compatible mode: zero-length normalized
    config::set_bedtools_compatible(true);
    assert_eq!(config::normalize_end(100, 100), 101);
    assert_eq!(config::normalize_end(100, 200), 200); // Non-zero unchanged

    reset_config();
}

#[test]
#[serial]
fn test_config_is_bedtools_compatible() {
    reset_config();

    assert!(!config::is_bedtools_compatible());

    config::set_bedtools_compatible(true);
    assert!(config::is_bedtools_compatible());

    config::set_bedtools_compatible(false);
    assert!(!config::is_bedtools_compatible());
}
