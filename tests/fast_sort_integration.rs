//! Integration tests for sort command.
//!
//! Tests verify:
//! 1. Default mode and --fast mode produce identical output
//! 2. Lexicographic chromosome ordering matches `LC_ALL=C sort -k1,1 -k2,2n -k3,3n`
//! 3. Genome file ordering (-g) works correctly
//! 4. Stability: input order preserved for ties
//! 5. Determinism: multiple runs produce identical output
//! 6. GNU sort oracle: byte-for-byte match against `sort -k1,1 -k2,2n -k3,3n`
//!
//! Sort specification:
//! - Primary: chromosome (lexicographic or genome order)
//! - Secondary: start coordinate (ascending)
//! - Tertiary: end coordinate (ascending)
//! - Ties: input order preserved (stable sort)

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::process::Command;

/// Parsed BED record for comparison
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct BedKey {
    chrom: String,
    start: u64,
    end: u64,
}

/// Parse a BED line into its sort key
fn parse_bed_key(line: &str) -> Option<BedKey> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() >= 3 {
        Some(BedKey {
            chrom: fields[0].to_string(),
            start: fields[1].parse().ok()?,
            end: fields[2].parse().ok()?,
        })
    } else {
        None
    }
}

/// Verify that a file is sorted by (chrom, start, end) in lexicographic chrom order
fn verify_sorted(lines: &[String]) -> bool {
    let mut prev: Option<BedKey> = None;
    for line in lines {
        if let Some(key) = parse_bed_key(line) {
            if let Some(ref p) = prev {
                if key < *p {
                    return false;
                }
            }
            prev = Some(key);
        }
    }
    true
}

/// Generate a random BED file with specified number of records.
/// Uses unique (chrom, start) combinations to avoid tie-breaking issues.
fn generate_random_bed_unique(path: &str, num_records: usize) {
    let mut file = File::create(path).expect("Failed to create test file");

    // Use deterministic seed for reproducibility
    let mut rng_state: u64 = 12345;
    let next_rand = |state: &mut u64| -> u64 {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        *state >> 33
    };

    let chroms = [
        "chr1", "chr2", "chr3", "chr10", "chr11", "chr20", "chrX", "chrY", "chrM",
    ];

    // Use a set to track used (chrom, start) combinations
    let mut used = std::collections::HashSet::new();
    let mut count = 0;

    while count < num_records {
        let chrom_idx = (next_rand(&mut rng_state) % chroms.len() as u64) as usize;
        let chrom = chroms[chrom_idx];
        let start = (next_rand(&mut rng_state) % 100_000_000) as u32;

        let key = (chrom_idx, start);
        if used.contains(&key) {
            continue; // Skip duplicates
        }
        used.insert(key);

        let length = (next_rand(&mut rng_state) % 10000 + 1) as u32;
        let end = start + length;
        let name = format!("gene_{}", count);

        writeln!(file, "{}\t{}\t{}\t{}", chrom, start, end, name).expect("Failed to write record");
        count += 1;
    }
}

/// Generate a BED file with edge cases for sorting (unique chrom+start combinations).
fn generate_edge_case_bed(path: &str) {
    let mut file = File::create(path).expect("Failed to create test file");

    // Use unique (chrom, start) combinations to avoid tie-breaking issues
    let records = [
        // Lexicographic chromosome ordering test
        "chr2\t100\t200\tD",
        "chr10\t100\t200\tE",
        "chr1\t101\t200\tF", // Different start from others on chr1
        // Various chromosomes
        "chrX\t1000\t2000\tG",
        "chrY\t1000\t2000\tH",
        "chrM\t100\t200\tI",
        // Mixed positions on chr1
        "chr1\t50\t100\tJ",
        "chr1\t200\t300\tK",
        "chr1\t150\t250\tL",
    ];

    for record in &records {
        writeln!(file, "{}", record).expect("Failed to write record");
    }
}

/// Read lines from a file.
fn read_lines(path: &str) -> Vec<String> {
    let file = File::open(path).expect("Failed to open file");
    BufReader::new(file)
        .lines()
        .map(|l| l.expect("Failed to read line"))
        .collect()
}

#[test]
fn test_fast_sort_matches_bedtools_random_100k() {
    let input_path = "/tmp/grit_test_random_100k.bed";
    let bedtools_output = "/tmp/grit_test_bedtools_output.bed";
    let grit_output = "/tmp/grit_test_grit_output.bed";

    // Generate random BED file with unique (chrom, start) combinations
    // This avoids tie-breaking issues since bedtools sort is unstable
    generate_random_bed_unique(input_path, 100_000);

    // Run bedtools sort
    let bedtools_status = Command::new("bedtools")
        .args(["sort", "-i", input_path])
        .stdout(File::create(bedtools_output).unwrap())
        .status()
        .expect("Failed to run bedtools");
    assert!(bedtools_status.success(), "bedtools sort failed");

    // Run grit sort --fast
    let grit_status = Command::new("./target/release/grit")
        .args(["sort", "--fast", "-i", input_path])
        .stdout(File::create(grit_output).unwrap())
        .status()
        .expect("Failed to run grit");
    assert!(grit_status.success(), "grit sort failed");

    // Compare outputs
    let bedtools_lines = read_lines(bedtools_output);
    let grit_lines = read_lines(grit_output);

    assert_eq!(
        bedtools_lines.len(),
        grit_lines.len(),
        "Output line counts differ: bedtools={}, grit={}",
        bedtools_lines.len(),
        grit_lines.len()
    );

    // Verify both outputs are correctly sorted
    assert!(
        verify_sorted(&bedtools_lines),
        "bedtools output is not sorted"
    );
    assert!(verify_sorted(&grit_lines), "grit output is not sorted");

    // Compare line by line (since we have unique (chrom, start), order should match)
    for (i, (bt_line, gr_line)) in bedtools_lines.iter().zip(grit_lines.iter()).enumerate() {
        assert_eq!(
            bt_line,
            gr_line,
            "Line {} differs:\nbedtools: {}\ngrit:     {}",
            i + 1,
            bt_line,
            gr_line
        );
    }

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(bedtools_output).ok();
    std::fs::remove_file(grit_output).ok();
}

#[test]
fn test_fast_sort_matches_bedtools_edge_cases() {
    let input_path = "/tmp/grit_test_edge_cases.bed";
    let bedtools_output = "/tmp/grit_test_bedtools_edge.bed";
    let grit_output = "/tmp/grit_test_grit_edge.bed";

    // Generate edge case BED file (with unique chrom+start combinations)
    generate_edge_case_bed(input_path);

    // Run bedtools sort
    let bedtools_status = Command::new("bedtools")
        .args(["sort", "-i", input_path])
        .stdout(File::create(bedtools_output).unwrap())
        .status()
        .expect("Failed to run bedtools");
    assert!(bedtools_status.success(), "bedtools sort failed");

    // Run grit sort --fast
    let grit_status = Command::new("./target/release/grit")
        .args(["sort", "--fast", "-i", input_path])
        .stdout(File::create(grit_output).unwrap())
        .status()
        .expect("Failed to run grit");
    assert!(grit_status.success(), "grit sort failed");

    // Compare outputs
    let bedtools_lines = read_lines(bedtools_output);
    let grit_lines = read_lines(grit_output);

    assert_eq!(
        bedtools_lines.len(),
        grit_lines.len(),
        "Output line counts differ"
    );

    // Verify both are sorted correctly
    assert!(
        verify_sorted(&bedtools_lines),
        "bedtools output is not sorted"
    );
    assert!(verify_sorted(&grit_lines), "grit output is not sorted");

    // With unique (chrom, start), outputs should match exactly
    for (i, (bt_line, gr_line)) in bedtools_lines.iter().zip(grit_lines.iter()).enumerate() {
        assert_eq!(
            bt_line,
            gr_line,
            "Line {} differs:\nbedtools: {}\ngrit:     {}",
            i + 1,
            bt_line,
            gr_line
        );
    }

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(bedtools_output).ok();
    std::fs::remove_file(grit_output).ok();
}

#[test]
fn test_fast_sort_same_start_different_end() {
    let input_path = "/tmp/grit_test_same_start.bed";
    let grit_output = "/tmp/grit_test_gr_same_start.bed";

    // Create test file with same start, different end
    let mut file = File::create(input_path).unwrap();
    writeln!(file, "chr1\t100\t500").unwrap();
    writeln!(file, "chr1\t100\t200").unwrap();
    writeln!(file, "chr1\t100\t300").unwrap();
    writeln!(file, "chr1\t100\t400").unwrap();
    drop(file);

    // Run grit sort --fast
    let grit_status = Command::new("./target/release/grit")
        .args(["sort", "--fast", "-i", input_path])
        .stdout(File::create(grit_output).unwrap())
        .status()
        .expect("Failed to run grit");
    assert!(grit_status.success());

    let grit_lines = read_lines(grit_output);

    // Verify output is sorted by (chrom, start, end)
    assert!(verify_sorted(&grit_lines), "grit output is not sorted");

    // Same (chrom, start) — sorted by end ascending
    assert_eq!(grit_lines.len(), 4);
    assert_eq!(grit_lines[0], "chr1\t100\t200");
    assert_eq!(grit_lines[1], "chr1\t100\t300");
    assert_eq!(grit_lines[2], "chr1\t100\t400");
    assert_eq!(grit_lines[3], "chr1\t100\t500");

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(grit_output).ok();
}

#[test]
fn test_fast_sort_lexicographic_chromosomes() {
    let input_path = "/tmp/grit_test_lex_chrom.bed";
    let bedtools_output = "/tmp/grit_test_bt_lex.bed";
    let grit_output = "/tmp/grit_test_gr_lex.bed";

    // Create test file with chromosomes that differ in lexicographic vs natural order
    let mut file = File::create(input_path).unwrap();
    writeln!(file, "chr2\t100\t200").unwrap();
    writeln!(file, "chr10\t100\t200").unwrap();
    writeln!(file, "chr1\t100\t200").unwrap();
    writeln!(file, "chr20\t100\t200").unwrap();
    writeln!(file, "chr3\t100\t200").unwrap();
    drop(file);

    // Run bedtools sort (lexicographic by default)
    let bedtools_status = Command::new("bedtools")
        .args(["sort", "-i", input_path])
        .stdout(File::create(bedtools_output).unwrap())
        .status()
        .expect("Failed to run bedtools");
    assert!(bedtools_status.success());

    // Run grit sort --fast
    let grit_status = Command::new("./target/release/grit")
        .args(["sort", "--fast", "-i", input_path])
        .stdout(File::create(grit_output).unwrap())
        .status()
        .expect("Failed to run grit");
    assert!(grit_status.success());

    // Compare outputs
    let bedtools_lines = read_lines(bedtools_output);
    let grit_lines = read_lines(grit_output);

    // Lexicographic order: chr1 < chr10 < chr2 < chr20 < chr3
    assert_eq!(bedtools_lines, grit_lines);
    assert!(grit_lines[0].starts_with("chr1\t"));
    assert!(grit_lines[1].starts_with("chr10\t"));
    assert!(grit_lines[2].starts_with("chr2\t"));
    assert!(grit_lines[3].starts_with("chr20\t"));
    assert!(grit_lines[4].starts_with("chr3\t"));

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(bedtools_output).ok();
    std::fs::remove_file(grit_output).ok();
}

#[test]
fn test_fast_sort_identical_records() {
    let input_path = "/tmp/grit_test_identical.bed";
    let grit_output = "/tmp/grit_test_gr_identical.bed";

    // Create test file with identical (chrom, start, end) but different names
    let mut file = File::create(input_path).unwrap();
    writeln!(file, "chr1\t100\t200\tfirst").unwrap();
    writeln!(file, "chr1\t100\t200\tsecond").unwrap();
    writeln!(file, "chr1\t100\t200\tthird").unwrap();
    drop(file);

    // Run grit sort --fast
    let grit_status = Command::new("./target/release/grit")
        .args(["sort", "--fast", "-i", input_path])
        .stdout(File::create(grit_output).unwrap())
        .status()
        .expect("Failed to run grit");
    assert!(grit_status.success());

    let grit_lines = read_lines(grit_output);

    // Verify output is sorted
    assert!(verify_sorted(&grit_lines), "grit output is not sorted");

    // All 3 records should be present
    assert_eq!(grit_lines.len(), 3);

    // grit preserves input order (stable sort) - verify this
    assert!(grit_lines[0].ends_with("first"));
    assert!(grit_lines[1].ends_with("second"));
    assert!(grit_lines[2].ends_with("third"));

    // When (chrom, start, end) are all equal, grit preserves input order (stable sort)
    // GNU sort without -s uses full-line comparison instead

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(grit_output).ok();
}

#[test]
fn test_default_and_fast_mode_identical_output() {
    let input_path = "/tmp/grit_test_mode_compare.bed";
    let default_output = "/tmp/grit_test_default_mode.bed";
    let fast_output = "/tmp/grit_test_fast_mode.bed";

    // Create test file with various cases including ties
    let mut file = File::create(input_path).unwrap();
    // Mixed chromosomes
    writeln!(file, "chr2\t500\t600\tA").unwrap();
    writeln!(file, "chr1\t100\t200\tB").unwrap();
    writeln!(file, "chr10\t100\t150\tC").unwrap();
    // Same chrom, different start
    writeln!(file, "chr1\t50\t100\tD").unwrap();
    writeln!(file, "chr1\t200\t300\tE").unwrap();
    // Same chrom and start, different end (tests stability)
    writeln!(file, "chr1\t100\t300\tF").unwrap();
    writeln!(file, "chr1\t100\t150\tG").unwrap();
    // More records
    writeln!(file, "chr2\t100\t200\tH").unwrap();
    writeln!(file, "chr2\t100\t250\tI").unwrap();
    drop(file);

    // Run grit sort (default mode - no --fast)
    let default_status = Command::new("./target/release/grit")
        .args(["sort", "-i", input_path])
        .stdout(File::create(default_output).unwrap())
        .status()
        .expect("Failed to run grit default");
    assert!(default_status.success(), "grit sort (default) failed");

    // Run grit sort --fast
    let fast_status = Command::new("./target/release/grit")
        .args(["sort", "--fast", "-i", input_path])
        .stdout(File::create(fast_output).unwrap())
        .status()
        .expect("Failed to run grit fast");
    assert!(fast_status.success(), "grit sort --fast failed");

    // Compare outputs - they MUST be identical
    let default_lines = read_lines(default_output);
    let fast_lines = read_lines(fast_output);

    assert_eq!(
        default_lines.len(),
        fast_lines.len(),
        "Line counts differ: default={}, fast={}",
        default_lines.len(),
        fast_lines.len()
    );

    for (i, (def_line, fast_line)) in default_lines.iter().zip(fast_lines.iter()).enumerate() {
        assert_eq!(
            def_line,
            fast_line,
            "Line {} differs between default and --fast mode:\ndefault: {}\nfast:    {}",
            i + 1,
            def_line,
            fast_line
        );
    }

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(default_output).ok();
    std::fs::remove_file(fast_output).ok();
}

#[test]
fn test_determinism_multiple_runs() {
    let input_path = "/tmp/grit_test_determinism.bed";
    let output1 = "/tmp/grit_test_det1.bed";
    let output2 = "/tmp/grit_test_det2.bed";
    let output3 = "/tmp/grit_test_det3.bed";

    // Create test file with ties
    let mut file = File::create(input_path).unwrap();
    writeln!(file, "chr1\t100\t300\tA").unwrap();
    writeln!(file, "chr1\t100\t200\tB").unwrap();
    writeln!(file, "chr1\t100\t400\tC").unwrap();
    writeln!(file, "chr2\t50\t100\tD").unwrap();
    writeln!(file, "chr2\t50\t150\tE").unwrap();
    drop(file);

    // Run grit sort multiple times
    for (output, run) in [(output1, 1), (output2, 2), (output3, 3)] {
        let status = Command::new("./target/release/grit")
            .args(["sort", "--fast", "-i", input_path])
            .stdout(File::create(output).unwrap())
            .status()
            .unwrap_or_else(|_| panic!("Failed to run grit (run {})", run));
        assert!(status.success(), "grit sort failed on run {}", run);
    }

    // All outputs must be identical
    let lines1 = read_lines(output1);
    let lines2 = read_lines(output2);
    let lines3 = read_lines(output3);

    assert_eq!(lines1, lines2, "Run 1 and Run 2 produced different outputs");
    assert_eq!(lines2, lines3, "Run 2 and Run 3 produced different outputs");

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(output1).ok();
    std::fs::remove_file(output2).ok();
    std::fs::remove_file(output3).ok();
}

#[test]
fn test_adversarial_same_chrom_start_different_end() {
    let input_path = "/tmp/grit_test_adversarial.bed";
    let default_output = "/tmp/grit_test_adv_default.bed";
    let fast_output = "/tmp/grit_test_adv_fast.bed";

    // Adversarial test: many records with same (chrom, start) but different end
    // Ends are ascending (200, 201, ...) so sort-by-end matches input order
    let mut file = File::create(input_path).unwrap();
    for i in 0..100 {
        // All same chrom and start, different end
        writeln!(file, "chr1\t100\t{}\trecord_{}", 200 + i, i).unwrap();
    }
    drop(file);

    // Run default mode
    let default_status = Command::new("./target/release/grit")
        .args(["sort", "-i", input_path])
        .stdout(File::create(default_output).unwrap())
        .status()
        .expect("Failed to run grit default");
    assert!(default_status.success());

    // Run fast mode
    let fast_status = Command::new("./target/release/grit")
        .args(["sort", "--fast", "-i", input_path])
        .stdout(File::create(fast_output).unwrap())
        .status()
        .expect("Failed to run grit fast");
    assert!(fast_status.success());

    let default_lines = read_lines(default_output);
    let fast_lines = read_lines(fast_output);

    // Sorted by end ascending (200, 201, ..., 299)
    assert_eq!(default_lines.len(), 100);
    assert_eq!(fast_lines.len(), 100);

    // Verify sorted by end (which happens to match input order since ends are ascending)
    for (i, line) in default_lines.iter().enumerate() {
        assert!(
            line.ends_with(&format!("record_{}", i)),
            "Default mode: unexpected order at position {}: {}",
            i,
            line
        );
    }

    for (i, line) in fast_lines.iter().enumerate() {
        assert!(
            line.ends_with(&format!("record_{}", i)),
            "Fast mode: unexpected order at position {}: {}",
            i,
            line
        );
    }

    // Both modes must produce identical output
    assert_eq!(
        default_lines, fast_lines,
        "Default and fast mode outputs differ"
    );

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(default_output).ok();
    std::fs::remove_file(fast_output).ok();
}

#[test]
fn test_genome_mode_parity_with_bedtools() {
    let input_path = "/tmp/grit_test_genome_input.bed";
    let genome_path = "/tmp/grit_test_genome.txt";
    let bedtools_output = "/tmp/grit_test_bt_genome.bed";
    let grit_output = "/tmp/grit_test_gr_genome.bed";
    let grit_fast_output = "/tmp/grit_test_gr_genome_fast.bed";

    // Create genome file with custom order: chr2, chr1, chr10
    let mut genome_file = File::create(genome_path).unwrap();
    writeln!(genome_file, "chr2\t1000000").unwrap();
    writeln!(genome_file, "chr1\t1000000").unwrap();
    writeln!(genome_file, "chr10\t1000000").unwrap();
    drop(genome_file);

    // Create test BED file
    let mut bed_file = File::create(input_path).unwrap();
    writeln!(bed_file, "chr1\t100\t200\tA").unwrap();
    writeln!(bed_file, "chr10\t50\t100\tB").unwrap();
    writeln!(bed_file, "chr2\t100\t200\tC").unwrap();
    writeln!(bed_file, "chr1\t50\t100\tD").unwrap();
    writeln!(bed_file, "chr2\t50\t100\tE").unwrap();
    drop(bed_file);

    // Run bedtools sort -g
    let bedtools_status = Command::new("bedtools")
        .args(["sort", "-g", genome_path, "-i", input_path])
        .stdout(File::create(bedtools_output).unwrap())
        .status()
        .expect("Failed to run bedtools");
    assert!(bedtools_status.success(), "bedtools sort -g failed");

    // Run grit sort -g (default mode)
    let grit_status = Command::new("./target/release/grit")
        .args(["sort", "-g", genome_path, "-i", input_path])
        .stdout(File::create(grit_output).unwrap())
        .status()
        .expect("Failed to run grit");
    assert!(grit_status.success(), "grit sort -g failed");

    // Run grit sort --fast -g
    let grit_fast_status = Command::new("./target/release/grit")
        .args(["sort", "--fast", "-g", genome_path, "-i", input_path])
        .stdout(File::create(grit_fast_output).unwrap())
        .status()
        .expect("Failed to run grit fast");
    assert!(grit_fast_status.success(), "grit sort --fast -g failed");

    let bedtools_lines = read_lines(bedtools_output);
    let grit_lines = read_lines(grit_output);
    let grit_fast_lines = read_lines(grit_fast_output);

    // All outputs should have same number of lines
    assert_eq!(bedtools_lines.len(), grit_lines.len());
    assert_eq!(grit_lines.len(), grit_fast_lines.len());

    // Verify chromosome order matches genome file: chr2, chr1, chr10
    // chr2 records should come first
    assert!(grit_lines[0].starts_with("chr2\t"), "First should be chr2");
    assert!(grit_lines[1].starts_with("chr2\t"), "Second should be chr2");
    // Then chr1
    assert!(grit_lines[2].starts_with("chr1\t"), "Third should be chr1");
    assert!(grit_lines[3].starts_with("chr1\t"), "Fourth should be chr1");
    // Then chr10
    assert!(
        grit_lines[4].starts_with("chr10\t"),
        "Fifth should be chr10"
    );

    // Default and fast mode must be identical
    assert_eq!(
        grit_lines, grit_fast_lines,
        "Default and fast -g mode differ"
    );

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(genome_path).ok();
    std::fs::remove_file(bedtools_output).ok();
    std::fs::remove_file(grit_output).ok();
    std::fs::remove_file(grit_fast_output).ok();
}

#[test]
fn test_reverse_mode_consistency() {
    let input_path = "/tmp/grit_test_reverse_input.bed";
    let default_output = "/tmp/grit_test_rev_default.bed";
    let fast_output = "/tmp/grit_test_rev_fast.bed";
    let default_rev_output = "/tmp/grit_test_rev_default_r.bed";
    let fast_rev_output = "/tmp/grit_test_rev_fast_r.bed";

    // Create test file
    let mut file = File::create(input_path).unwrap();
    writeln!(file, "chr2\t100\t200\tA").unwrap();
    writeln!(file, "chr1\t100\t200\tB").unwrap();
    writeln!(file, "chr1\t50\t100\tC").unwrap();
    writeln!(file, "chr10\t100\t200\tD").unwrap();
    drop(file);

    // Run default sort
    Command::new("./target/release/grit")
        .args(["sort", "-i", input_path])
        .stdout(File::create(default_output).unwrap())
        .status()
        .expect("Failed");

    // Run fast sort
    Command::new("./target/release/grit")
        .args(["sort", "--fast", "-i", input_path])
        .stdout(File::create(fast_output).unwrap())
        .status()
        .expect("Failed");

    // Run default sort -r
    Command::new("./target/release/grit")
        .args(["sort", "-r", "-i", input_path])
        .stdout(File::create(default_rev_output).unwrap())
        .status()
        .expect("Failed");

    // Run fast sort -r
    Command::new("./target/release/grit")
        .args(["sort", "--fast", "-r", "-i", input_path])
        .stdout(File::create(fast_rev_output).unwrap())
        .status()
        .expect("Failed");

    let default_lines = read_lines(default_output);
    let fast_lines = read_lines(fast_output);
    let default_rev_lines = read_lines(default_rev_output);
    let fast_rev_lines = read_lines(fast_rev_output);

    // Default and fast must be identical
    assert_eq!(default_lines, fast_lines, "Default and fast differ");

    // Reverse outputs must be identical
    assert_eq!(
        default_rev_lines, fast_rev_lines,
        "Default -r and fast -r differ"
    );

    // Reverse must be exact reverse of forward
    let default_reversed: Vec<_> = default_lines.iter().rev().cloned().collect();
    assert_eq!(
        default_reversed, default_rev_lines,
        "Reverse is not exact reverse of forward"
    );

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(default_output).ok();
    std::fs::remove_file(fast_output).ok();
    std::fs::remove_file(default_rev_output).ok();
    std::fs::remove_file(fast_rev_output).ok();
}

/// Test sort idempotency: sorting an already-sorted file must produce identical output.
/// This is a critical correctness property: grit sort X > A; grit sort A > B; diff A B → identical
#[test]
fn test_sort_idempotency() {
    let input_path = "/tmp/grit_test_idempotency_input.bed";
    let sorted_a = "/tmp/grit_test_idempotency_a.bed";
    let sorted_b = "/tmp/grit_test_idempotency_b.bed";

    // Create test file with ties (same chrom+start, different end)
    let mut file = File::create(input_path).unwrap();
    writeln!(file, "chr2\t100\t300\tX").unwrap();
    writeln!(file, "chr1\t100\t200\tA").unwrap();
    writeln!(file, "chr1\t100\t150\tB").unwrap();
    writeln!(file, "chr1\t100\t250\tC").unwrap();
    writeln!(file, "chr1\t50\t100\tD").unwrap();
    writeln!(file, "chr10\t100\t200\tE").unwrap();
    writeln!(file, "chr10\t100\t300\tF").unwrap();
    writeln!(file, "chr2\t100\t200\tY").unwrap();
    drop(file);

    // First sort: X -> A
    let status_a = Command::new("./target/release/grit")
        .args(["sort", "-i", input_path])
        .stdout(File::create(sorted_a).unwrap())
        .status()
        .expect("Failed to run first sort");
    assert!(status_a.success(), "First sort failed");

    // Second sort: A -> B
    let status_b = Command::new("./target/release/grit")
        .args(["sort", "-i", sorted_a])
        .stdout(File::create(sorted_b).unwrap())
        .status()
        .expect("Failed to run second sort");
    assert!(status_b.success(), "Second sort failed");

    // A and B must be identical (idempotency)
    let lines_a = read_lines(sorted_a);
    let lines_b = read_lines(sorted_b);

    assert_eq!(
        lines_a.len(),
        lines_b.len(),
        "Idempotency failed: line counts differ"
    );

    for (i, (line_a, line_b)) in lines_a.iter().zip(lines_b.iter()).enumerate() {
        assert_eq!(
            line_a,
            line_b,
            "Idempotency failed at line {}:\nFirst sort:  {}\nSecond sort: {}",
            i + 1,
            line_a,
            line_b
        );
    }

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(sorted_a).ok();
    std::fs::remove_file(sorted_b).ok();
}

/// Test sort idempotency with fast mode
#[test]
fn test_sort_idempotency_fast_mode() {
    let input_path = "/tmp/grit_test_idempotency_fast_input.bed";
    let sorted_a = "/tmp/grit_test_idempotency_fast_a.bed";
    let sorted_b = "/tmp/grit_test_idempotency_fast_b.bed";

    // Create test file with many ties
    let mut file = File::create(input_path).unwrap();
    for i in 0..50 {
        writeln!(file, "chr1\t100\t{}\trecord_{}", 200 + i, i).unwrap();
    }
    for i in 0..30 {
        writeln!(file, "chr2\t500\t{}\trecord2_{}", 600 + i, i).unwrap();
    }
    writeln!(file, "chr10\t100\t200\tlast").unwrap();
    drop(file);

    // First sort: X -> A
    let status_a = Command::new("./target/release/grit")
        .args(["sort", "--fast", "-i", input_path])
        .stdout(File::create(sorted_a).unwrap())
        .status()
        .expect("Failed to run first sort");
    assert!(status_a.success(), "First sort failed");

    // Second sort: A -> B
    let status_b = Command::new("./target/release/grit")
        .args(["sort", "--fast", "-i", sorted_a])
        .stdout(File::create(sorted_b).unwrap())
        .status()
        .expect("Failed to run second sort");
    assert!(status_b.success(), "Second sort failed");

    // A and B must be identical
    let lines_a = read_lines(sorted_a);
    let lines_b = read_lines(sorted_b);

    assert_eq!(lines_a, lines_b, "Sort idempotency failed in fast mode");

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(sorted_a).ok();
    std::fs::remove_file(sorted_b).ok();
}

/// Test that grit generate output is idempotent with grit sort.
/// Generated files should already be sorted such that grit sort produces identical output.
#[test]
fn test_generate_sort_idempotency() {
    let gen_dir = "/tmp/grit_test_gen_idempotency";
    let sorted_output = "/tmp/grit_test_gen_sorted.bed";

    // Generate a small dataset with explicit sorting (--sorted yes)
    // Note: --sorted auto doesn't sort files < 1M intervals
    let gen_status = Command::new("./target/release/grit")
        .args([
            "generate", "--output", gen_dir, "--sizes", "10K", "--mode", "balanced", "--seed",
            "42", "--sorted", "yes", "--force",
        ])
        .status()
        .expect("Failed to run grit generate");
    assert!(gen_status.success(), "grit generate failed");

    let generated_file = format!("{}/balanced/10K/A.bed", gen_dir);

    // Sort the generated file
    let sort_status = Command::new("./target/release/grit")
        .args(["sort", "-i", &generated_file])
        .stdout(File::create(sorted_output).unwrap())
        .status()
        .expect("Failed to run grit sort");
    assert!(sort_status.success(), "grit sort failed");

    // Generated file and sorted output must be identical
    let gen_lines = read_lines(&generated_file);
    let sorted_lines = read_lines(sorted_output);

    assert_eq!(
        gen_lines.len(),
        sorted_lines.len(),
        "Generate-sort idempotency failed: line counts differ ({} vs {})",
        gen_lines.len(),
        sorted_lines.len()
    );

    for (i, (gen_line, sorted_line)) in gen_lines.iter().zip(sorted_lines.iter()).enumerate() {
        assert_eq!(
            gen_line,
            sorted_line,
            "Generate-sort idempotency failed at line {}:\nGenerated: {}\nSorted:    {}",
            i + 1,
            gen_line,
            sorted_line
        );
    }

    // Cleanup
    std::fs::remove_dir_all(gen_dir).ok();
    std::fs::remove_file(sorted_output).ok();
}

/// Test generate + sort idempotency with larger file (triggers different code paths)
#[test]
fn test_generate_sort_idempotency_1m() {
    let gen_dir = "/tmp/grit_test_gen_idempotency_1m";
    let sorted_output = "/tmp/grit_test_gen_sorted_1m.bed";

    // Generate 1M intervals (uses in-memory sort in generate)
    let gen_status = Command::new("./target/release/grit")
        .args([
            "generate", "--output", gen_dir, "--sizes", "1M", "--mode", "balanced", "--seed",
            "12345", "--force",
        ])
        .status()
        .expect("Failed to run grit generate");
    assert!(gen_status.success(), "grit generate failed");

    let generated_file = format!("{}/balanced/1M/A.bed", gen_dir);

    // Sort the generated file
    let sort_status = Command::new("./target/release/grit")
        .args(["sort", "-i", &generated_file])
        .stdout(File::create(sorted_output).unwrap())
        .status()
        .expect("Failed to run grit sort");
    assert!(sort_status.success(), "grit sort failed");

    // Compare files using diff for efficiency
    let diff_status = Command::new("diff")
        .args(["-q", &generated_file, sorted_output])
        .status()
        .expect("Failed to run diff");

    assert!(
        diff_status.success(),
        "Generate-sort idempotency failed for 1M intervals: files differ"
    );

    // Cleanup
    std::fs::remove_dir_all(gen_dir).ok();
    std::fs::remove_file(sorted_output).ok();
}

/// Test that grit sort output matches GNU sort -k1,1 -k2,2n -k3,3n byte-for-byte.
/// This is the oracle test: GNU sort is the ground truth.
#[test]
fn test_matches_gnu_sort() {
    let input_path = "/tmp/grit_test_gnu_sort.bed";
    let gnu_output = "/tmp/grit_test_gnu_sort_out.bed";
    let grit_output = "/tmp/grit_test_grit_gnu_out.bed";

    // Create test BED with ties on (chrom, start) to exercise end-sorting
    let mut file = File::create(input_path).unwrap();
    // Same chrom+start, different end
    writeln!(file, "chr1\t100\t500\tA").unwrap();
    writeln!(file, "chr1\t100\t200\tB").unwrap();
    writeln!(file, "chr1\t100\t300\tC").unwrap();
    // Different starts
    writeln!(file, "chr1\t50\t100\tD").unwrap();
    writeln!(file, "chr1\t200\t300\tE").unwrap();
    // Different chroms (lexicographic ordering)
    writeln!(file, "chr2\t100\t200\tF").unwrap();
    writeln!(file, "chr10\t100\t200\tG").unwrap();
    writeln!(file, "chr10\t100\t150\tH").unwrap();
    // Same chrom+start, different end
    writeln!(file, "chr2\t100\t300\tI").unwrap();
    writeln!(file, "chr2\t100\t200\tJ").unwrap();
    drop(file);

    // Run GNU sort: LC_ALL=C sort -k1,1 -k2,2n -k3,3n
    let gnu_status = Command::new("sort")
        .env("LC_ALL", "C")
        .args(["-k1,1", "-k2,2n", "-k3,3n", input_path])
        .stdout(File::create(gnu_output).unwrap())
        .status()
        .expect("Failed to run GNU sort");
    assert!(gnu_status.success(), "GNU sort failed");

    // Run grit sort
    let grit_status = Command::new("./target/release/grit")
        .args(["sort", "-i", input_path])
        .stdout(File::create(grit_output).unwrap())
        .status()
        .expect("Failed to run grit");
    assert!(grit_status.success(), "grit sort failed");

    let gnu_lines = read_lines(gnu_output);
    let grit_lines = read_lines(grit_output);

    assert_eq!(
        gnu_lines.len(),
        grit_lines.len(),
        "Line counts differ: GNU sort={}, grit={}",
        gnu_lines.len(),
        grit_lines.len()
    );

    // Compare line by line
    for (i, (gnu_line, grit_line)) in gnu_lines.iter().zip(grit_lines.iter()).enumerate() {
        assert_eq!(
            gnu_line,
            grit_line,
            "Line {} differs:\nGNU sort: {}\ngrit:     {}",
            i + 1,
            gnu_line,
            grit_line
        );
    }

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(gnu_output).ok();
    std::fs::remove_file(grit_output).ok();
}

/// Large-scale GNU sort oracle test with 100K records.
/// Uses unique (chrom, start, end) to avoid tie-breaking differences
/// (GNU sort without -s uses full-line comparison for ties, grit uses stable sort).
#[test]
fn test_matches_gnu_sort_large() {
    let input_path = "/tmp/grit_test_gnu_sort_100k.bed";
    let gnu_output = "/tmp/grit_test_gnu_sort_100k_out.bed";
    let grit_output = "/tmp/grit_test_grit_gnu_100k_out.bed";

    // Generate random BED file with unique (chrom, start) combinations
    // (unique start within a chrom means unique (chrom, start, end) tuples)
    generate_random_bed_unique(input_path, 100_000);

    // Run GNU sort
    let gnu_status = Command::new("sort")
        .env("LC_ALL", "C")
        .args(["-k1,1", "-k2,2n", "-k3,3n", input_path])
        .stdout(File::create(gnu_output).unwrap())
        .status()
        .expect("Failed to run GNU sort");
    assert!(gnu_status.success(), "GNU sort failed");

    // Run grit sort
    let grit_status = Command::new("./target/release/grit")
        .args(["sort", "-i", input_path])
        .stdout(File::create(grit_output).unwrap())
        .status()
        .expect("Failed to run grit");
    assert!(grit_status.success(), "grit sort failed");

    let gnu_lines = read_lines(gnu_output);
    let grit_lines = read_lines(grit_output);

    assert_eq!(
        gnu_lines.len(),
        grit_lines.len(),
        "Line counts differ: GNU sort={}, grit={}",
        gnu_lines.len(),
        grit_lines.len()
    );

    // Byte-for-byte comparison
    for (i, (gnu_line, grit_line)) in gnu_lines.iter().zip(grit_lines.iter()).enumerate() {
        assert_eq!(
            gnu_line,
            grit_line,
            "Line {} differs:\nGNU sort: {}\ngrit:     {}",
            i + 1,
            gnu_line,
            grit_line
        );
    }

    // Cleanup
    std::fs::remove_file(input_path).ok();
    std::fs::remove_file(gnu_output).ok();
    std::fs::remove_file(grit_output).ok();
}
