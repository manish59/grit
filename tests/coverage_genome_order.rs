//! Regression test for coverage with genome-ordered chromosomes.
//!
//! Bug: streaming coverage used lexicographic chromosome comparison,
//! which failed when files were sorted in genome order (chr9 before chr10).
//!
//! Root cause: "chr10" < "chr9" lexicographically, so when B moved to chr10
//! while A was on chr9, the code thought chr10 was "earlier" and skipped it.
//!
//! Fix: Use equality comparison only, not < or >, to handle both
//! lexicographic and genome sort orders.

use std::io::Write;
use std::process::Command;
use tempfile::NamedTempFile;

/// Test coverage with genome-ordered chromosomes (chr9 before chr10).
///
/// This is a regression test for the bug where coverage missed overlaps
/// when chromosomes were sorted in genome order instead of lexicographic order.
#[test]
fn test_coverage_genome_order_chr9_chr10() {
    // Create test files with genome order (chr9 before chr10)
    let mut a_file = NamedTempFile::new().unwrap();
    let mut b_file = NamedTempFile::new().unwrap();

    writeln!(a_file, "chr9\t100\t200").unwrap();
    writeln!(a_file, "chr10\t100\t200").unwrap();

    writeln!(b_file, "chr9\t50\t150").unwrap();
    writeln!(b_file, "chr10\t50\t150").unwrap();

    a_file.flush().unwrap();
    b_file.flush().unwrap();

    // Run grit coverage
    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--",
            "coverage",
            "-a",
            a_file.path().to_str().unwrap(),
            "-b",
            b_file.path().to_str().unwrap(),
            "--assume-sorted",
        ])
        .output()
        .expect("Failed to run grit coverage");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = stdout.lines().collect();

    assert_eq!(lines.len(), 2, "Expected 2 output lines");

    // Parse coverage output: chrom start end count bases_covered length fraction
    // chr9 should have coverage: 1 overlap, 50 bases (150 - 100)
    assert!(
        lines[0].starts_with("chr9\t100\t200\t1\t50\t"),
        "chr9 should have 1 overlap with 50 bases covered, got: {}",
        lines[0]
    );

    // chr10 should also have coverage: 1 overlap, 50 bases (150 - 100)
    // This was the bug - chr10 was getting 0 coverage
    assert!(
        lines[1].starts_with("chr10\t100\t200\t1\t50\t"),
        "chr10 should have 1 overlap with 50 bases covered, got: {}",
        lines[1]
    );
}

/// Test coverage with multiple chromosome transitions in genome order.
///
/// Tests: chr1 -> chr2 -> chr10 -> chr11 -> chr20 -> chr9 (if genome order)
/// This exercises various chromosome comparison edge cases.
#[test]
fn test_coverage_genome_order_multiple_chroms() {
    let mut a_file = NamedTempFile::new().unwrap();
    let mut b_file = NamedTempFile::new().unwrap();

    // Genome order: chr1, chr2, chr10, chr11, chr20, chrX
    writeln!(a_file, "chr1\t100\t200").unwrap();
    writeln!(a_file, "chr2\t100\t200").unwrap();
    writeln!(a_file, "chr10\t100\t200").unwrap();
    writeln!(a_file, "chr11\t100\t200").unwrap();
    writeln!(a_file, "chr20\t100\t200").unwrap();
    writeln!(a_file, "chrX\t100\t200").unwrap();

    writeln!(b_file, "chr1\t50\t150").unwrap();
    writeln!(b_file, "chr2\t50\t150").unwrap();
    writeln!(b_file, "chr10\t50\t150").unwrap();
    writeln!(b_file, "chr11\t50\t150").unwrap();
    writeln!(b_file, "chr20\t50\t150").unwrap();
    writeln!(b_file, "chrX\t50\t150").unwrap();

    a_file.flush().unwrap();
    b_file.flush().unwrap();

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--",
            "coverage",
            "-a",
            a_file.path().to_str().unwrap(),
            "-b",
            b_file.path().to_str().unwrap(),
            "--assume-sorted",
        ])
        .output()
        .expect("Failed to run grit coverage");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = stdout.lines().collect();

    assert_eq!(lines.len(), 6, "Expected 6 output lines");

    // All chromosomes should have coverage
    for (i, chrom) in ["chr1", "chr2", "chr10", "chr11", "chr20", "chrX"]
        .iter()
        .enumerate()
    {
        assert!(
            lines[i].starts_with(&format!("{}\t100\t200\t1\t50\t", chrom)),
            "{} should have 1 overlap with 50 bases covered, got: {}",
            chrom,
            lines[i]
        );
    }
}

/// Test coverage when B has chromosomes that A doesn't have.
///
/// This tests that the algorithm correctly skips B chromosomes.
#[test]
fn test_coverage_b_has_extra_chroms() {
    let mut a_file = NamedTempFile::new().unwrap();
    let mut b_file = NamedTempFile::new().unwrap();

    // A only has chr2 and chr10
    writeln!(a_file, "chr2\t100\t200").unwrap();
    writeln!(a_file, "chr10\t100\t200").unwrap();

    // B has chr1, chr2, chr9, chr10 (genome order)
    writeln!(b_file, "chr1\t50\t150").unwrap();
    writeln!(b_file, "chr2\t50\t150").unwrap();
    writeln!(b_file, "chr9\t50\t150").unwrap();
    writeln!(b_file, "chr10\t50\t150").unwrap();

    a_file.flush().unwrap();
    b_file.flush().unwrap();

    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--",
            "coverage",
            "-a",
            a_file.path().to_str().unwrap(),
            "-b",
            b_file.path().to_str().unwrap(),
            "--assume-sorted",
        ])
        .output()
        .expect("Failed to run grit coverage");

    let stdout = String::from_utf8_lossy(&output.stdout);
    let lines: Vec<&str> = stdout.lines().collect();

    assert_eq!(lines.len(), 2, "Expected 2 output lines");

    // chr2 should have coverage
    assert!(
        lines[0].starts_with("chr2\t100\t200\t1\t50\t"),
        "chr2 should have coverage, got: {}",
        lines[0]
    );

    // chr10 should have coverage
    assert!(
        lines[1].starts_with("chr10\t100\t200\t1\t50\t"),
        "chr10 should have coverage, got: {}",
        lines[1]
    );
}
