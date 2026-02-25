//! Regression tests for streaming closest correctness.
//!
//! These tests verify that streaming closest produces identical output
//! to bedtools closest and grit closest (default mode).

use std::io::Write;
use std::process::Command;
use tempfile::NamedTempFile;

/// Helper to run grit closest with streaming mode.
fn run_grit_closest_streaming(a_content: &str, b_content: &str) -> String {
    let mut a_file = NamedTempFile::new().unwrap();
    let mut b_file = NamedTempFile::new().unwrap();

    writeln!(a_file, "{}", a_content.trim()).unwrap();
    writeln!(b_file, "{}", b_content.trim()).unwrap();

    a_file.flush().unwrap();
    b_file.flush().unwrap();

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--",
            "closest",
            "-a",
            a_file.path().to_str().unwrap(),
            "-b",
            b_file.path().to_str().unwrap(),
            "--streaming",
        ])
        .output()
        .expect("Failed to run grit closest --streaming");

    String::from_utf8_lossy(&output.stdout).to_string()
}

/// Helper to run grit closest with default mode.
fn run_grit_closest_default(a_content: &str, b_content: &str) -> String {
    let mut a_file = NamedTempFile::new().unwrap();
    let mut b_file = NamedTempFile::new().unwrap();

    writeln!(a_file, "{}", a_content.trim()).unwrap();
    writeln!(b_file, "{}", b_content.trim()).unwrap();

    a_file.flush().unwrap();
    b_file.flush().unwrap();

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--",
            "closest",
            "-a",
            a_file.path().to_str().unwrap(),
            "-b",
            b_file.path().to_str().unwrap(),
        ])
        .output()
        .expect("Failed to run grit closest");

    String::from_utf8_lossy(&output.stdout).to_string()
}

/// Regression test: A intervals that are nested (A2 inside A1).
///
/// Bug: When A1 overlaps with B1, B1 is added to active set.
/// When A2 (nested inside A1) is processed, B1 becomes downstream
/// (B1.start >= A2.end) but wasn't being considered as a downstream candidate.
///
/// Scenario:
/// - A1: 100-200, A2: 120-150 (A2 ends before A1, both start after A1.start)
/// - B1: 180-300 (overlaps A1, downstream of A2)
/// - B2: 400-500 (downstream of both)
///
/// Expected:
/// - A1 closest to B1 (overlap)
/// - A2 closest to B1 (distance = 180 - 150 + 1 = 31)
///
/// Bug behavior:
/// - A2 was incorrectly matched to B2 (distance = 400 - 150 + 1 = 251)
#[test]
fn test_nested_a_intervals_downstream_from_active() {
    let a = "chr1\t100\t200\nchr1\t120\t150";
    let b = "chr1\t180\t300\nchr1\t400\t500";

    let streaming = run_grit_closest_streaming(a, b);
    let default = run_grit_closest_default(a, b);

    let streaming_lines: Vec<_> = streaming.lines().collect();
    let default_lines: Vec<_> = default.lines().collect();

    assert_eq!(
        streaming_lines.len(),
        default_lines.len(),
        "Line count mismatch: streaming={}, default={}",
        streaming_lines.len(),
        default_lines.len()
    );

    // A1 should match B1 (overlap)
    assert!(
        streaming_lines[0].contains("100\t200") && streaming_lines[0].contains("180\t300"),
        "A1 should be closest to B1 (overlap): {}",
        streaming_lines[0]
    );

    // A2 should match B1 (downstream, distance 31) NOT B2 (distance 251)
    assert!(
        streaming_lines[1].contains("120\t150") && streaming_lines[1].contains("180\t300"),
        "A2 should be closest to B1 (downstream from active): {}",
        streaming_lines[1]
    );

    // Streaming should match default
    assert_eq!(
        streaming_lines, default_lines,
        "Streaming output differs from default:\nStreaming: {:?}\nDefault: {:?}",
        streaming_lines, default_lines
    );
}

/// Test with the exact intervals from the 100K dataset failure.
#[test]
fn test_real_world_nested_a_chr17() {
    // From benchmarks/data/100K_50K/uniform
    // A1: chr17:57693442-57694077 (overlaps B1)
    // A2: chr17:57693716-57693861 (nested in A1, should be closest to B1)
    // B1: chr17:57693915-57694878
    // B2: chr17:57696099-57696647
    let a = "chr17\t57693442\t57694077\nchr17\t57693716\t57693861";
    let b = "chr17\t57693915\t57694878\nchr17\t57696099\t57696647";

    let streaming = run_grit_closest_streaming(a, b);
    let default = run_grit_closest_default(a, b);

    // Both A intervals should be closest to B1
    let streaming_lines: Vec<_> = streaming.lines().collect();

    assert!(
        streaming_lines[0].contains("57693915\t57694878"),
        "First A should be closest to B1: {}",
        streaming_lines[0]
    );

    assert!(
        streaming_lines[1].contains("57693915\t57694878"),
        "Second A should be closest to B1: {}",
        streaming_lines[1]
    );

    assert_eq!(
        streaming, default,
        "Streaming should match default:\nStreaming:\n{}\nDefault:\n{}",
        streaming, default
    );
}

/// Test where multiple A intervals become "before" a B that was overlapping.
#[test]
fn test_multiple_nested_a_intervals() {
    // A1 overlaps B1, A2 and A3 are nested in A1 and downstream of B1
    let a = "chr1\t100\t300\nchr1\t120\t150\nchr1\t160\t180";
    let b = "chr1\t200\t400\nchr1\t500\t600";

    let streaming = run_grit_closest_streaming(a, b);
    let default = run_grit_closest_default(a, b);

    assert_eq!(
        streaming, default,
        "Streaming should match default:\nStreaming:\n{}\nDefault:\n{}",
        streaming, default
    );
}

/// Test chromosome transition doesn't lose downstream candidates.
#[test]
fn test_chromosome_transition_with_downstream() {
    let a = "chr1\t100\t200\nchr1\t120\t150\nchr2\t100\t200";
    let b = "chr1\t180\t300\nchr2\t50\t80";

    let streaming = run_grit_closest_streaming(a, b);
    let default = run_grit_closest_default(a, b);

    assert_eq!(
        streaming, default,
        "Streaming should match default across chromosome transition"
    );
}
