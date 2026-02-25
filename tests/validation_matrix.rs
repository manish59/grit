//! Comprehensive validation test matrix for GRIT.
//!
//! Tests cover:
//! 1. Unsorted input detection
//! 2. stdin validation
//! 3. Genome-order validation
//! 4. Flag combinations (--assume-sorted, --allow-unsorted)

use std::io::Write;
use std::process::{Command, Output};
use tempfile::NamedTempFile;

/// Helper to create a temporary BED file.
fn create_bed_file(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    write!(file, "{}", content).unwrap();
    file.flush().unwrap();
    file
}

/// Helper to create a temporary genome file.
fn create_genome_file(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    write!(file, "{}", content).unwrap();
    file.flush().unwrap();
    file
}

/// Helper to run grit command and return output.
fn run_grit(args: &[&str]) -> Output {
    Command::new("cargo")
        .args(["run", "--release", "--"])
        .args(args)
        .output()
        .expect("Failed to run grit")
}

/// Helper to check if command succeeded.
fn is_success(output: &Output) -> bool {
    output.status.success()
}

/// Helper to get stderr as string.
fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).to_string()
}

// =============================================================================
// Test fixtures
// =============================================================================

/// Lexicographically sorted BED file (chr1 < chr10 < chr2).
fn lex_sorted_bed() -> &'static str {
    "chr1\t100\t200\nchr1\t200\t300\nchr10\t100\t200\nchr2\t100\t200\n"
}

/// Genome-order sorted BED file (chr1 < chr2 < chr10).
fn genome_sorted_bed() -> &'static str {
    "chr1\t100\t200\nchr1\t200\t300\nchr2\t100\t200\nchr10\t100\t200\n"
}

/// Unsorted BED file (interleaved chromosomes).
fn unsorted_chrom_bed() -> &'static str {
    "chr1\t100\t200\nchr2\t100\t200\nchr1\t300\t400\n"
}

/// Unsorted BED file (out of order positions).
fn unsorted_position_bed() -> &'static str {
    "chr1\t200\t300\nchr1\t100\t200\n"
}

/// Genome file with natural order (chr1, chr2, chr10).
fn genome_natural() -> &'static str {
    "chr1\t1000000\nchr2\t500000\nchr10\t250000\n"
}

// =============================================================================
// Merge command tests
// =============================================================================

#[test]
fn test_merge_sorted_input_succeeds() {
    let bed = create_bed_file(lex_sorted_bed());
    let output = run_grit(&["merge", "-i", bed.path().to_str().unwrap()]);
    assert!(is_success(&output), "Sorted input should succeed");
}

#[test]
fn test_merge_unsorted_chrom_fails() {
    let bed = create_bed_file(unsorted_chrom_bed());
    let output = run_grit(&["merge", "-i", bed.path().to_str().unwrap()]);
    assert!(!is_success(&output), "Unsorted chromosomes should fail");
    assert!(
        stderr(&output).contains("not sorted"),
        "Error should mention sorting"
    );
}

#[test]
fn test_merge_unsorted_position_fails() {
    let bed = create_bed_file(unsorted_position_bed());
    let output = run_grit(&["merge", "-i", bed.path().to_str().unwrap()]);
    assert!(!is_success(&output), "Unsorted positions should fail");
    assert!(
        stderr(&output).contains("not sorted"),
        "Error should mention sorting"
    );
}

#[test]
fn test_merge_assume_sorted_bypasses_validation() {
    let bed = create_bed_file(unsorted_chrom_bed());
    let output = run_grit(&[
        "merge",
        "-i",
        bed.path().to_str().unwrap(),
        "--assume-sorted",
    ]);
    // May produce incorrect results but should not fail
    assert!(
        is_success(&output),
        "--assume-sorted should bypass validation"
    );
}

// =============================================================================
// Genome-order validation tests
// =============================================================================

#[test]
fn test_merge_genome_order_correct() {
    let bed = create_bed_file(genome_sorted_bed());
    let genome = create_genome_file(genome_natural());
    let output = run_grit(&[
        "merge",
        "-i",
        bed.path().to_str().unwrap(),
        "-g",
        genome.path().to_str().unwrap(),
    ]);
    assert!(
        is_success(&output),
        "Genome-order sorted should succeed with -g"
    );
}

#[test]
fn test_merge_genome_order_incorrect() {
    // lex_sorted has chr10 before chr2, which violates genome_natural order
    let bed = create_bed_file(lex_sorted_bed());
    let genome = create_genome_file(genome_natural());
    let output = run_grit(&[
        "merge",
        "-i",
        bed.path().to_str().unwrap(),
        "-g",
        genome.path().to_str().unwrap(),
    ]);
    assert!(
        !is_success(&output),
        "Wrong genome order should fail with -g"
    );
    assert!(
        stderr(&output).contains("genome order"),
        "Error should mention genome order"
    );
}

#[test]
fn test_merge_missing_chrom_in_genome() {
    let bed = create_bed_file("chr1\t100\t200\nchrX\t100\t200\n");
    let genome = create_genome_file(genome_natural()); // No chrX
    let output = run_grit(&[
        "merge",
        "-i",
        bed.path().to_str().unwrap(),
        "-g",
        genome.path().to_str().unwrap(),
    ]);
    assert!(
        !is_success(&output),
        "Missing chromosome should fail with -g"
    );
    assert!(
        stderr(&output).contains("not found in genome"),
        "Error should mention missing chromosome"
    );
}

// =============================================================================
// Intersect command tests
// =============================================================================

#[test]
fn test_intersect_streaming_sorted_succeeds() {
    let a = create_bed_file(lex_sorted_bed());
    let b = create_bed_file(lex_sorted_bed());
    let output = run_grit(&[
        "intersect",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
    ]);
    assert!(
        is_success(&output),
        "Streaming with sorted input should succeed"
    );
}

#[test]
fn test_intersect_streaming_unsorted_fails() {
    let a = create_bed_file(unsorted_chrom_bed());
    let b = create_bed_file(lex_sorted_bed());
    let output = run_grit(&[
        "intersect",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
    ]);
    assert!(
        !is_success(&output),
        "Streaming with unsorted A should fail"
    );
}

#[test]
fn test_intersect_allow_unsorted_processes() {
    let a = create_bed_file(unsorted_position_bed());
    let b = create_bed_file(lex_sorted_bed());
    let output = run_grit(&[
        "intersect",
        "--allow-unsorted",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
    ]);
    // Should succeed (silently re-sorts in memory)
    assert!(
        is_success(&output),
        "--allow-unsorted should allow unsorted input"
    );
}

#[test]
fn test_intersect_genome_order_validation() {
    let a = create_bed_file(genome_sorted_bed());
    let b = create_bed_file(genome_sorted_bed());
    let genome = create_genome_file(genome_natural());
    let output = run_grit(&[
        "intersect",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
        "-g",
        genome.path().to_str().unwrap(),
    ]);
    assert!(
        is_success(&output),
        "Genome-order sorted input should succeed with -g"
    );
}

// =============================================================================
// Subtract command tests
// =============================================================================

#[test]
fn test_subtract_streaming_sorted_succeeds() {
    let a = create_bed_file(lex_sorted_bed());
    let b = create_bed_file(lex_sorted_bed());
    let output = run_grit(&[
        "subtract",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
    ]);
    assert!(is_success(&output), "Streaming subtract should succeed");
}

#[test]
fn test_subtract_unsorted_fails_by_default() {
    let a = create_bed_file(unsorted_chrom_bed());
    let b = create_bed_file(lex_sorted_bed());
    let output = run_grit(&[
        "subtract",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
    ]);
    assert!(
        !is_success(&output),
        "Unsorted input should fail by default"
    );
}

// =============================================================================
// Closest command tests
// =============================================================================

#[test]
fn test_closest_streaming_sorted_succeeds() {
    let a = create_bed_file(lex_sorted_bed());
    let b = create_bed_file(lex_sorted_bed());
    let output = run_grit(&[
        "closest",
        "--streaming",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
    ]);
    assert!(is_success(&output), "Streaming closest should succeed");
}

#[test]
fn test_closest_unsorted_fails_by_default() {
    let a = create_bed_file(unsorted_position_bed());
    let b = create_bed_file(lex_sorted_bed());
    let output = run_grit(&[
        "closest",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
    ]);
    assert!(
        !is_success(&output),
        "Unsorted input should fail by default"
    );
}

// =============================================================================
// Window command tests
// =============================================================================

#[test]
fn test_window_sorted_succeeds() {
    let a = create_bed_file(lex_sorted_bed());
    let b = create_bed_file(lex_sorted_bed());
    let output = run_grit(&[
        "window",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
    ]);
    assert!(
        is_success(&output),
        "Window with sorted input should succeed"
    );
}

#[test]
fn test_window_unsorted_fails() {
    let a = create_bed_file(unsorted_chrom_bed());
    let b = create_bed_file(lex_sorted_bed());
    let output = run_grit(&[
        "window",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
    ]);
    assert!(
        !is_success(&output),
        "Window with unsorted input should fail"
    );
}

// =============================================================================
// Coverage command tests
// =============================================================================

#[test]
fn test_coverage_sorted_succeeds() {
    let a = create_bed_file(lex_sorted_bed());
    let b = create_bed_file(lex_sorted_bed());
    let output = run_grit(&[
        "coverage",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
    ]);
    assert!(
        is_success(&output),
        "Coverage with sorted input should succeed"
    );
}

#[test]
fn test_coverage_unsorted_fails() {
    let a = create_bed_file(unsorted_position_bed());
    let b = create_bed_file(lex_sorted_bed());
    let output = run_grit(&[
        "coverage",
        "-a",
        a.path().to_str().unwrap(),
        "-b",
        b.path().to_str().unwrap(),
    ]);
    assert!(
        !is_success(&output),
        "Coverage with unsorted input should fail"
    );
}

// =============================================================================
// Error message quality tests
// =============================================================================

#[test]
fn test_error_message_includes_fix_suggestion() {
    let bed = create_bed_file(unsorted_chrom_bed());
    let output = run_grit(&["merge", "-i", bed.path().to_str().unwrap()]);
    let err = stderr(&output);
    assert!(
        err.contains("Fix:") || err.contains("grit sort"),
        "Error should include fix suggestion: {}",
        err
    );
}

#[test]
fn test_error_message_includes_file_path() {
    let bed = create_bed_file(unsorted_chrom_bed());
    let output = run_grit(&[
        "intersect",
        "-a",
        bed.path().to_str().unwrap(),
        "-b",
        bed.path().to_str().unwrap(),
    ]);
    let err = stderr(&output);
    // Should mention "File A" or similar
    assert!(
        err.contains("File A") || err.contains("A is not"),
        "Error should identify which file is unsorted: {}",
        err
    );
}
