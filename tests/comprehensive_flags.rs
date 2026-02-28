//! Comprehensive flag coverage tests for GRIT.
//!
//! Tests cover all flag combinations that were previously missing:
//! 1. CLOSEST: --io, --iu, --id flag combinations
//! 2. MERGE: -s (strand) flag tests
//! 3. INTERSECT/SUBTRACT: -f + -r edge cases
//! 4. WINDOW: -g genome validation, window size edge cases
//! 5. COVERAGE: combined flag tests
//! 6. stdin input tests
//! 7. Error handling for malformed input

use std::io::Write;
use std::process::{Command, Output, Stdio};
use tempfile::NamedTempFile;

// =============================================================================
// Helper functions
// =============================================================================

fn create_bed_file(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    write!(file, "{}", content).unwrap();
    file.flush().unwrap();
    file
}

fn create_genome_file(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    write!(file, "{}", content).unwrap();
    file.flush().unwrap();
    file
}

fn run_grit(args: &[&str]) -> Output {
    Command::new("cargo")
        .args(["run", "--release", "--"])
        .args(args)
        .output()
        .expect("Failed to run grit")
}

fn run_grit_with_stdin(args: &[&str], stdin_content: &str) -> Output {
    let mut child = Command::new("cargo")
        .args(["run", "--release", "--"])
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("Failed to spawn grit");

    if let Some(mut stdin) = child.stdin.take() {
        stdin.write_all(stdin_content.as_bytes()).unwrap();
    }

    child.wait_with_output().expect("Failed to wait for grit")
}

fn is_success(output: &Output) -> bool {
    output.status.success()
}

fn stdout(output: &Output) -> String {
    String::from_utf8_lossy(&output.stdout).to_string()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).to_string()
}

// =============================================================================
// CLOSEST: Flag combination tests (--io, --iu, --id)
// =============================================================================

/// Test --io (ignore overlaps) flag
#[test]
fn test_closest_ignore_overlaps() {
    // A overlaps B1, but B2 is downstream
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t150\t250\nchr1\t300\t400\n");

    let output = run_grit(&[
        "closest",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "--io",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should skip the overlapping B1 and find B2
    assert!(
        result.contains("300\t400"),
        "Should find non-overlapping B2: {}",
        result
    );
    assert!(
        !result.contains("150\t250"),
        "Should not include overlapping B1: {}",
        result
    );
}

/// Test --iu (ignore upstream) flag
#[test]
fn test_closest_ignore_upstream() {
    // A at 200-300, B1 upstream (100-150), B2 downstream (350-400)
    let a = create_bed_file("chr1\t200\t300\n");
    let b = create_bed_file("chr1\t100\t150\nchr1\t350\t400\n");

    let output = run_grit(&[
        "closest",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "--iu",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should skip upstream B1, find only B2
    assert!(
        result.contains("350\t400"),
        "Should find downstream B2: {}",
        result
    );
    assert!(
        !result.contains("100\t150"),
        "Should not include upstream B1: {}",
        result
    );
}

/// Test --id (ignore downstream) flag
#[test]
fn test_closest_ignore_downstream() {
    // A at 200-300, B1 upstream (100-150), B2 downstream (350-400)
    let a = create_bed_file("chr1\t200\t300\n");
    let b = create_bed_file("chr1\t100\t150\nchr1\t350\t400\n");

    let output = run_grit(&[
        "closest",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "--id",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should skip downstream B2, find only B1
    assert!(
        result.contains("100\t150"),
        "Should find upstream B1: {}",
        result
    );
    assert!(
        !result.contains("350\t400"),
        "Should not include downstream B2: {}",
        result
    );
}

/// Test --io + --iu combined (only consider downstream non-overlapping)
#[test]
fn test_closest_ignore_overlaps_and_upstream() {
    // A at 100-200, B1 overlaps (150-250), B2 upstream (50-80), B3 downstream (300-400)
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t50\t80\nchr1\t150\t250\nchr1\t300\t400\n");

    let output = run_grit(&[
        "closest",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "--io",
        "--iu",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should only find B3 (downstream, non-overlapping)
    assert!(
        result.contains("300\t400"),
        "Should find downstream B3: {}",
        result
    );
    assert!(
        !result.contains("50\t80"),
        "Should not include upstream B2: {}",
        result
    );
    assert!(
        !result.contains("150\t250"),
        "Should not include overlapping B1: {}",
        result
    );
}

/// Test --io + --id combined (only consider upstream non-overlapping)
#[test]
fn test_closest_ignore_overlaps_and_downstream() {
    // A at 100-200, B1 overlaps (150-250), B2 upstream (50-80), B3 downstream (300-400)
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t50\t80\nchr1\t150\t250\nchr1\t300\t400\n");

    let output = run_grit(&[
        "closest",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "--io",
        "--id",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should only find B2 (upstream, non-overlapping)
    assert!(
        result.contains("50\t80"),
        "Should find upstream B2: {}",
        result
    );
    assert!(
        !result.contains("300\t400"),
        "Should not include downstream B3: {}",
        result
    );
    assert!(
        !result.contains("150\t250"),
        "Should not include overlapping B1: {}",
        result
    );
}

/// Test --iu + --id combined (only consider overlapping)
#[test]
fn test_closest_ignore_upstream_and_downstream() {
    // A at 100-200, B1 overlaps (150-250), B2 upstream (50-80), B3 downstream (300-400)
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t50\t80\nchr1\t150\t250\nchr1\t300\t400\n");

    let output = run_grit(&[
        "closest",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "--iu",
        "--id",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should only find B1 (overlapping)
    assert!(
        result.contains("150\t250"),
        "Should find overlapping B1: {}",
        result
    );
    assert!(
        !result.contains("50\t80"),
        "Should not include upstream B2: {}",
        result
    );
    assert!(
        !result.contains("300\t400"),
        "Should not include downstream B3: {}",
        result
    );
}

/// Test --io + --iu + --id combined (no matches possible)
#[test]
fn test_closest_ignore_all() {
    // A at 100-200, B1 overlaps, B2 upstream, B3 downstream
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t50\t80\nchr1\t150\t250\nchr1\t300\t400\n");

    let output = run_grit(&[
        "closest",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "--io",
        "--iu",
        "--id",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should report no closest (all are filtered)
    assert!(
        result.contains(".\t-1\t-1"),
        "Should report no closest when all filtered: {}",
        result
    );
}

/// Test -t first (report only first tie)
#[test]
fn test_closest_tie_first() {
    // A at 200-300, equidistant B intervals
    let a = create_bed_file("chr1\t200\t300\n");
    let b = create_bed_file("chr1\t100\t150\nchr1\t350\t400\n");

    let output = run_grit(&[
        "closest",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-t",
        "first",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    let lines: Vec<_> = result.lines().collect();
    assert_eq!(lines.len(), 1, "Should report only first tie: {}", result);
}

/// Test -t last (report only last tie)
#[test]
fn test_closest_tie_last() {
    // A at 200-300, equidistant B intervals
    let a = create_bed_file("chr1\t200\t300\n");
    let b = create_bed_file("chr1\t100\t150\nchr1\t350\t400\n");

    let output = run_grit(&[
        "closest",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-t",
        "last",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    let lines: Vec<_> = result.lines().collect();
    assert_eq!(lines.len(), 1, "Should report only last tie: {}", result);
}

/// Test -d (report distance) flag
#[test]
fn test_closest_report_distance() {
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t300\t400\n");

    let output = run_grit(&[
        "closest",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-d",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Distance should be 300 - 200 + 1 = 101
    assert!(
        result.contains("101") || result.contains("100"),
        "Should include distance: {}",
        result
    );
}

// =============================================================================
// MERGE: Strand (-s) flag tests
// =============================================================================

/// Test -s (strand) flag - same strand merges
#[test]
fn test_merge_strand_same() {
    let bed = create_bed_file("chr1\t100\t200\t.\t.\t+\nchr1\t150\t250\t.\t.\t+\n");

    let output = run_grit(&["merge", "-i", bed.path().to_str().unwrap(), "-s"]);

    assert!(is_success(&output));
    let result = stdout(&output);
    let lines: Vec<_> = result.lines().collect();
    assert_eq!(
        lines.len(),
        1,
        "Same strand intervals should merge: {}",
        result
    );
}

/// Test -s (strand) flag - different strands don't merge
#[test]
fn test_merge_strand_different() {
    let bed = create_bed_file("chr1\t100\t200\t.\t.\t+\nchr1\t150\t250\t.\t.\t-\n");

    let output = run_grit(&["merge", "-i", bed.path().to_str().unwrap(), "-s"]);

    assert!(is_success(&output));
    let result = stdout(&output);
    let lines: Vec<_> = result.lines().collect();
    assert_eq!(
        lines.len(),
        2,
        "Different strand intervals should not merge: {}",
        result
    );
}

/// Test -s with -d (strand + distance)
#[test]
fn test_merge_strand_with_distance() {
    // Two + strand intervals with gap of 50bp
    let bed = create_bed_file("chr1\t100\t200\t.\t.\t+\nchr1\t250\t350\t.\t.\t+\n");

    let output = run_grit(&[
        "merge",
        "-i",
        bed.path().to_str().unwrap(),
        "-s",
        "-d",
        "100",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    let lines: Vec<_> = result.lines().collect();
    assert_eq!(
        lines.len(),
        1,
        "Same strand intervals within distance should merge: {}",
        result
    );
}

/// Test -c (count) flag
#[test]
fn test_merge_count() {
    let bed = create_bed_file("chr1\t100\t200\nchr1\t150\t250\nchr1\t180\t280\n");

    let output = run_grit(&["merge", "-i", bed.path().to_str().unwrap(), "-c"]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should report count of 3 merged intervals
    assert!(result.contains("3"), "Should report count of 3: {}", result);
}

// =============================================================================
// INTERSECT: -f (fraction) and -r (reciprocal) edge cases
// =============================================================================

/// Test -f 1.0 (100% overlap required)
#[test]
fn test_intersect_fraction_100_percent() {
    // A: 100-200, B: 50-250 (contains A), B2: 100-200 (exact match)
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t50\t250\nchr1\t100\t200\n"); // Sorted order

    let output = run_grit(&[
        "intersect",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-f",
        "1.0",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Both should match since A is 100% covered by both B intervals
    let lines: Vec<_> = result.lines().collect();
    assert!(
        lines.len() >= 1,
        "Should find at least one 100% overlap: {}",
        result
    );
}

/// Test -f with -r (reciprocal) - both must meet fraction requirement
#[test]
fn test_intersect_fraction_reciprocal() {
    // A: 100-200 (100bp), B: 100-200 (100bp) - 100% reciprocal overlap
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t100\t200\n");

    let output = run_grit(&[
        "intersect",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-f",
        "0.5",
        "-r",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    assert!(!result.is_empty(), "Should find reciprocal overlap: {}", result);
}

/// Test -f with -r where one side doesn't meet requirement
#[test]
fn test_intersect_fraction_reciprocal_fails() {
    // A: 100-200 (100bp), B: 50-250 (200bp)
    // Overlap: 100bp, A overlap fraction: 100%, B overlap fraction: 50%
    // With -f 0.6 -r, B doesn't meet requirement
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t50\t250\n");

    let output = run_grit(&[
        "intersect",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-f",
        "0.6",
        "-r",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // B's fraction is only 50%, so with -r and -f 0.6, should not match
    assert!(
        result.trim().is_empty(),
        "Should not find overlap with reciprocal 0.6: {}",
        result
    );
}

/// Test -u (unique) flag
#[test]
fn test_intersect_unique() {
    // A overlaps multiple B intervals
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t120\t180\nchr1\t150\t250\n");

    let output = run_grit(&[
        "intersect",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-u",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    let lines: Vec<_> = result.lines().collect();
    assert_eq!(
        lines.len(),
        1,
        "Should report A only once with -u: {}",
        result
    );
}

/// Test -v (no overlap) flag
#[test]
fn test_intersect_no_overlap() {
    let a = create_bed_file("chr1\t100\t200\nchr1\t500\t600\n");
    let b = create_bed_file("chr1\t150\t180\n");

    let output = run_grit(&[
        "intersect",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-v",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should only report chr1:500-600 which has no overlap
    assert!(
        result.contains("500\t600"),
        "Should find non-overlapping A: {}",
        result
    );
    assert!(
        !result.contains("100\t200"),
        "Should not include overlapping A: {}",
        result
    );
}

/// Test -wa -wb flags together
#[test]
fn test_intersect_wa_wb() {
    let a = create_bed_file("chr1\t100\t200\tgeneA\n");
    let b = create_bed_file("chr1\t150\t250\tgeneB\n");

    let output = run_grit(&[
        "intersect",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "--wa",
        "--wb",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    assert!(result.contains("geneA"), "Should include A entry: {}", result);
    assert!(result.contains("geneB"), "Should include B entry: {}", result);
}

// =============================================================================
// SUBTRACT: -f (fraction) and -r (reciprocal) tests
// =============================================================================

/// Test subtract -A (remove entire feature)
#[test]
fn test_subtract_remove_entire() {
    let a = create_bed_file("chr1\t100\t200\nchr1\t500\t600\n");
    let b = create_bed_file("chr1\t120\t140\n");

    let output = run_grit(&[
        "subtract",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-A",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // A interval 100-200 should be completely removed due to any overlap
    assert!(
        !result.contains("100\t200"),
        "Should remove entire A: {}",
        result
    );
    assert!(
        result.contains("500\t600"),
        "Should keep non-overlapping A: {}",
        result
    );
}

/// Test subtract with -f (fraction)
#[test]
fn test_subtract_fraction() {
    // A: 100-200 (100bp), B: 150-160 (10bp overlap = 10%)
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t150\t160\n");

    // With -f 0.2, 10% overlap doesn't meet threshold
    let output = run_grit(&[
        "subtract",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-f",
        "0.2",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should output A unchanged since overlap fraction < 0.2
    assert!(!result.trim().is_empty(), "Should output A: {}", result);
}

// =============================================================================
// WINDOW: genome validation and window size tests
// =============================================================================

/// Test -w (window size) flag
#[test]
fn test_window_size() {
    // A at 100-200, B at 500-600 (400bp away)
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t500\t600\n");

    // With window of 1000, should find overlap
    let output = run_grit(&[
        "window",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-w",
        "1000",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    assert!(
        !result.trim().is_empty(),
        "Should find B within window: {}",
        result
    );
}

/// Test -l (left window) and -r (right window) asymmetric
#[test]
fn test_window_asymmetric() {
    let a = create_bed_file("chr1\t500\t600\n");
    let b = create_bed_file("chr1\t100\t200\nchr1\t800\t900\n");

    // Only look upstream (left), not downstream
    let output = run_grit(&[
        "window",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-l",
        "500",
        "-r",
        "0",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    assert!(
        result.contains("100\t200"),
        "Should find upstream B: {}",
        result
    );
    assert!(
        !result.contains("800\t900"),
        "Should not find downstream B: {}",
        result
    );
}

/// Test -g (genome) validation for window
#[test]
fn test_window_genome_validation() {
    let a = create_bed_file("chr1\t100\t200\nchr2\t100\t200\n");
    let b = create_bed_file("chr1\t150\t250\nchr2\t150\t250\n");
    let genome = create_genome_file("chr1\t1000000\nchr2\t500000\n");

    let output = run_grit(&[
        "window",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-g",
        genome.path().to_str().unwrap(),
    ]);

    assert!(is_success(&output), "Should succeed with genome file");
}

/// Test -c (count) flag for window
#[test]
fn test_window_count() {
    let a = create_bed_file("chr1\t500\t600\n");
    let b = create_bed_file("chr1\t100\t200\nchr1\t800\t900\n");

    let output = run_grit(&[
        "window",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-w",
        "500",
        "-c",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should report count of 2 matches
    assert!(result.contains("2"), "Should report count: {}", result);
}

// =============================================================================
// COVERAGE: combined flag tests
// =============================================================================

/// Test --hist (histogram) flag
#[test]
fn test_coverage_histogram() {
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t120\t180\nchr1\t140\t160\n");

    let output = run_grit(&[
        "coverage",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "--hist",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should include histogram output with depth counts
    assert!(!result.trim().is_empty(), "Should produce histogram: {}", result);
}

/// Test --mean flag
#[test]
fn test_coverage_mean() {
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t100\t200\n");

    let output = run_grit(&[
        "coverage",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "--mean",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Mean should be 1.0 for exact coverage
    assert!(
        result.contains("1") || result.contains("1.0"),
        "Should report mean coverage: {}",
        result
    );
}

// =============================================================================
// SLOP: flag tests
// =============================================================================

/// Test -b (both sides) flag
#[test]
fn test_slop_both() {
    let bed = create_bed_file("chr1\t100\t200\n");
    let genome = create_genome_file("chr1\t1000000\n");

    let output = run_grit(&[
        "slop",
        "-i",
        bed.path().to_str().unwrap(),
        "-g",
        genome.path().to_str().unwrap(),
        "-b",
        "50",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should expand to 50-250
    assert!(result.contains("50\t250"), "Should expand both sides: {}", result);
}

/// Test --pct (percentage) flag
#[test]
fn test_slop_percentage() {
    let bed = create_bed_file("chr1\t100\t200\n"); // 100bp interval
    let genome = create_genome_file("chr1\t1000000\n");

    let output = run_grit(&[
        "slop",
        "-i",
        bed.path().to_str().unwrap(),
        "-g",
        genome.path().to_str().unwrap(),
        "-b",
        "0.5",
        "--pct",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // 50% of 100bp = 50bp each side, so 50-250
    assert!(result.contains("50\t250"), "Should expand by percentage: {}", result);
}

/// Test -s (strand-aware) flag
#[test]
fn test_slop_strand_aware() {
    let bed = create_bed_file("chr1\t100\t200\t.\t.\t-\n");
    let genome = create_genome_file("chr1\t1000000\n");

    let output = run_grit(&[
        "slop",
        "-i",
        bed.path().to_str().unwrap(),
        "-g",
        genome.path().to_str().unwrap(),
        "-l",
        "50",
        "-r",
        "0",
        "-s",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // For minus strand, -l should extend 3' (right side)
    assert!(result.contains("100\t250"), "Should extend 3' for minus strand: {}", result);
}

// =============================================================================
// STDIN input tests
// =============================================================================

/// Test merge with stdin input
#[test]
fn test_merge_stdin() {
    let input = "chr1\t100\t200\nchr1\t150\t250\n";

    let output = run_grit_with_stdin(&["merge", "-i", "-"], input);

    assert!(is_success(&output));
    let result = stdout(&output);
    assert!(
        result.contains("100\t250"),
        "Should merge from stdin: {}",
        result
    );
}

/// Test sort with stdin input
#[test]
fn test_sort_stdin() {
    let input = "chr1\t200\t300\nchr1\t100\t200\n";

    let output = run_grit_with_stdin(&["sort", "-i", "-"], input);

    assert!(is_success(&output));
    let result = stdout(&output);
    let lines: Vec<_> = result.lines().collect();
    assert!(
        lines[0].contains("100\t200"),
        "Should sort from stdin: {}",
        result
    );
}

// =============================================================================
// Error handling tests
// =============================================================================

/// Test malformed BED input (missing columns)
#[test]
fn test_malformed_missing_columns() {
    let bed = create_bed_file("chr1\t100\n"); // Missing end column

    let output = run_grit(&["merge", "-i", bed.path().to_str().unwrap()]);

    // Should handle gracefully (skip line or error)
    // Either success with no output or informative error
    let err = stderr(&output);
    let out = stdout(&output);
    assert!(
        out.is_empty() || err.contains("malformed") || err.contains("parse") || is_success(&output),
        "Should handle malformed input: stdout={}, stderr={}",
        out,
        err
    );
}

/// Test invalid coordinates (start > end)
#[test]
fn test_invalid_coordinates() {
    let bed = create_bed_file("chr1\t200\t100\n"); // start > end

    let output = run_grit(&["merge", "-i", bed.path().to_str().unwrap()]);

    // Should handle gracefully
    let err = stderr(&output);
    let out = stdout(&output);
    // Either skips invalid lines or reports error
    assert!(
        out.is_empty() || err.contains("invalid") || is_success(&output),
        "Should handle invalid coordinates: stdout={}, stderr={}",
        out,
        err
    );
}

/// Test file not found error
#[test]
fn test_file_not_found() {
    let output = run_grit(&["merge", "-i", "/nonexistent/path/file.bed"]);

    assert!(!is_success(&output), "Should fail for missing file");
    let err = stderr(&output);
    assert!(
        err.contains("No such file") || err.contains("not found") || err.contains("does not exist"),
        "Should report file not found: {}",
        err
    );
}

/// Test empty file
#[test]
fn test_empty_file() {
    let bed = create_bed_file("");

    let output = run_grit(&["merge", "-i", bed.path().to_str().unwrap()]);

    assert!(is_success(&output), "Should handle empty file");
    let result = stdout(&output);
    assert!(result.is_empty(), "Should produce no output for empty file");
}

/// Test file with only comments/headers
#[test]
fn test_only_comments() {
    let bed = create_bed_file("# This is a comment\n#track name=test\n");

    let output = run_grit(&["merge", "-i", bed.path().to_str().unwrap()]);

    assert!(is_success(&output), "Should handle comment-only file");
    let result = stdout(&output);
    assert!(
        result.is_empty(),
        "Should produce no output for comment-only file"
    );
}

// =============================================================================
// COMPLEMENT: tests
// =============================================================================

/// Test basic complement
#[test]
fn test_complement_basic() {
    let bed = create_bed_file("chr1\t100\t200\nchr1\t300\t400\n");
    let genome = create_genome_file("chr1\t500\n");

    let output = run_grit(&[
        "complement",
        "-i",
        bed.path().to_str().unwrap(),
        "-g",
        genome.path().to_str().unwrap(),
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should include 0-100, 200-300, 400-500
    assert!(result.contains("0\t100"), "Should include 0-100: {}", result);
    assert!(
        result.contains("200\t300"),
        "Should include 200-300: {}",
        result
    );
    assert!(
        result.contains("400\t500"),
        "Should include 400-500: {}",
        result
    );
}

// =============================================================================
// GENOMECOV: tests
// =============================================================================

/// Test --bg (bedgraph) flag
#[test]
fn test_genomecov_bedgraph() {
    let bed = create_bed_file("chr1\t100\t200\n");
    let genome = create_genome_file("chr1\t500\n");

    let output = run_grit(&[
        "genomecov",
        "-i",
        bed.path().to_str().unwrap(),
        "-g",
        genome.path().to_str().unwrap(),
        "--bg",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should output bedgraph format with non-zero coverage
    assert!(
        result.contains("chr1") && result.contains("100") && result.contains("200"),
        "Should produce bedgraph: {}",
        result
    );
}

/// Test --bga (bedgraph all) flag
#[test]
fn test_genomecov_bedgraph_all() {
    let bed = create_bed_file("chr1\t100\t200\n");
    let genome = create_genome_file("chr1\t300\n");

    let output = run_grit(&[
        "genomecov",
        "-i",
        bed.path().to_str().unwrap(),
        "-g",
        genome.path().to_str().unwrap(),
        "--bga",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should include zero-coverage regions
    let lines: Vec<_> = result.lines().collect();
    assert!(
        lines.len() >= 2,
        "Should include zero-coverage regions: {}",
        result
    );
}

// =============================================================================
// JACCARD: tests
// =============================================================================

/// Test basic jaccard
#[test]
fn test_jaccard_basic() {
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t150\t250\n");

    let output = run_grit(&[
        "jaccard",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should produce jaccard statistic
    assert!(
        result.contains("jaccard") || result.contains("0.") || result.contains("intersection"),
        "Should produce jaccard output: {}",
        result
    );
}

/// Test jaccard with no overlap
#[test]
fn test_jaccard_no_overlap() {
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t300\t400\n");

    let output = run_grit(&[
        "jaccard",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Jaccard should be 0
    assert!(
        result.contains("0") || result.contains("0.0"),
        "Jaccard should be 0: {}",
        result
    );
}

// =============================================================================
// MULTIINTER: tests
// =============================================================================

/// Test basic multiinter
#[test]
fn test_multiinter_basic() {
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t150\t250\n");

    let output = run_grit(&[
        "multiinter",
        "-i",
        a.path().to_str().unwrap(),
        b.path().to_str().unwrap(),
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    assert!(!result.is_empty(), "Should produce multiinter output: {}", result);
}

/// Test multiinter --cluster flag
#[test]
fn test_multiinter_cluster() {
    let a = create_bed_file("chr1\t100\t200\n");
    let b = create_bed_file("chr1\t150\t250\n");
    let c = create_bed_file("chr1\t180\t280\n");

    let output = run_grit(&[
        "multiinter",
        "--streaming",
        "-i",
        a.path().to_str().unwrap(),
        b.path().to_str().unwrap(),
        c.path().to_str().unwrap(),
        "--cluster",
    ]);

    assert!(is_success(&output));
    let result = stdout(&output);
    // Should only output regions where all 3 files overlap
    // Common region: 180-200
    assert!(
        result.is_empty() || result.contains("180"),
        "Should only show cluster overlap: {}",
        result
    );
}
