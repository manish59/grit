//! Sort validation for streaming operations.
//!
//! Streaming algorithms require sorted input. This module provides
//! utilities to verify sort order before processing.
//!
//! Sort validation checks that:
//! 1. All records for a chromosome are contiguous (no interleaving)
//! 2. Within a chromosome, positions are non-decreasing
//!
//! This supports both lexicographic order (chr1, chr10, chr2...) and
//! genome order (chr1, chr2, chr3...) - any consistent ordering works.

use crate::bed::{BedError, BedReader};
use std::collections::HashSet;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Verify that a BED file is sorted by chromosome and position.
///
/// Returns Ok(()) if sorted, Err with details if not.
///
/// # Example
///
/// ```rust,no_run
/// use grit_genomics::streaming::verify_sorted;
///
/// verify_sorted("input.bed").expect("File must be sorted");
/// ```
pub fn verify_sorted<P: AsRef<Path>>(path: P) -> Result<(), BedError> {
    let file = File::open(path.as_ref())?;
    let reader = BedReader::new(BufReader::new(file));

    let mut prev_chrom: Option<String> = None;
    let mut prev_start: u64 = 0;
    let mut seen_chroms: HashSet<String> = HashSet::new();
    let mut line_num = 0;

    for result in reader.records() {
        let rec = result?;
        line_num += 1;

        let chrom = rec.chrom();
        let start = rec.start();

        if let Some(ref pc) = prev_chrom {
            if chrom != pc {
                // Switching chromosomes - check we haven't seen this one before
                if seen_chroms.contains(chrom) {
                    return Err(BedError::InvalidFormat(format!(
                        "File not sorted: chromosome '{}' at line {} was seen earlier (chromosomes must be contiguous)",
                        chrom, line_num
                    )));
                }
                seen_chroms.insert(pc.clone());
            } else if start < prev_start {
                return Err(BedError::InvalidFormat(format!(
                    "File not sorted: position {} at line {} comes after {} on {}",
                    start, line_num, prev_start, chrom
                )));
            }
        }

        prev_chrom = Some(chrom.to_string());
        prev_start = start;
    }

    Ok(())
}

/// Inline sort validator for use within streaming loops.
///
/// This avoids the overhead of reading the file twice (once for validation,
/// once for processing) by validating sort order as records are processed.
///
/// Validates that:
/// 1. All records for a chromosome are contiguous (no interleaving)
/// 2. Within a chromosome, positions are non-decreasing
#[derive(Debug, Default)]
pub struct SortValidator {
    prev_chrom: Option<String>,
    prev_start: u64,
    seen_chroms: HashSet<String>,
    record_count: usize,
}

impl SortValidator {
    /// Create a new sort validator.
    pub fn new() -> Self {
        Self::default()
    }

    /// Validate that the given record maintains sort order.
    ///
    /// Returns Ok(()) if valid, Err if out of order.
    #[inline]
    pub fn validate(&mut self, chrom: &str, start: u64) -> Result<(), BedError> {
        self.record_count += 1;

        if let Some(ref pc) = self.prev_chrom {
            if chrom != pc {
                // Switching chromosomes - check we haven't seen this one before
                if self.seen_chroms.contains(chrom) {
                    return Err(BedError::InvalidFormat(format!(
                        "File not sorted: chromosome '{}' at record {} was seen earlier (chromosomes must be contiguous)",
                        chrom, self.record_count
                    )));
                }
                self.seen_chroms.insert(pc.clone());
            } else if start < self.prev_start {
                return Err(BedError::InvalidFormat(format!(
                    "File not sorted: position {} at record {} comes after {} on {}",
                    start, self.record_count, self.prev_start, chrom
                )));
            }
        }

        self.prev_chrom = Some(chrom.to_string());
        self.prev_start = start;

        Ok(())
    }

    /// Check if the record maintains sort order, returning a file-specific error.
    ///
    /// This variant includes the file identifier in the error message.
    #[inline]
    pub fn validate_with_file(
        &mut self,
        chrom: &str,
        start: u64,
        file_id: &str,
    ) -> Result<(), BedError> {
        self.record_count += 1;

        if let Some(ref pc) = self.prev_chrom {
            if chrom != pc {
                // Switching chromosomes - check we haven't seen this one before
                if self.seen_chroms.contains(chrom) {
                    return Err(BedError::InvalidFormat(format!(
                        "File {} not sorted: chromosome '{}' at record {} was seen earlier (chromosomes must be contiguous)",
                        file_id, chrom, self.record_count
                    )));
                }
                self.seen_chroms.insert(pc.clone());
            } else if start < self.prev_start {
                return Err(BedError::InvalidFormat(format!(
                    "File {} not sorted: position {} at record {} comes after {} on {}",
                    file_id, start, self.record_count, self.prev_start, chrom
                )));
            }
        }

        self.prev_chrom = Some(chrom.to_string());
        self.prev_start = start;

        Ok(())
    }

    /// Reset validator state (for a new chromosome or file).
    pub fn reset(&mut self) {
        self.prev_chrom = None;
        self.prev_start = 0;
        self.seen_chroms.clear();
        self.record_count = 0;
    }

    /// Get the number of records validated.
    pub fn record_count(&self) -> usize {
        self.record_count
    }
}

/// Verify that a BED file is sorted according to genome file order.
///
/// Validates that:
/// 1. Chromosomes appear in the order specified by the genome file
/// 2. Within a chromosome, positions are non-decreasing
/// 3. All chromosomes exist in the genome file
///
/// Returns Ok(()) if valid, Err with details if not.
pub fn verify_sorted_with_genome<P: AsRef<Path>>(
    path: P,
    genome: &crate::genome::Genome,
) -> Result<(), BedError> {
    let file = File::open(path.as_ref())?;
    let reader = BedReader::new(BufReader::new(file));

    // Build a map of chromosome -> position in genome order
    let chrom_order: std::collections::HashMap<&str, usize> = genome
        .chromosomes()
        .enumerate()
        .map(|(i, c)| (c.as_str(), i))
        .collect();

    let mut prev_chrom: Option<String> = None;
    let mut prev_chrom_order: Option<usize> = None;
    let mut prev_start: u64 = 0;
    let mut line_num = 0;

    for result in reader.records() {
        let rec = result?;
        line_num += 1;

        let chrom = rec.chrom();
        let start = rec.start();

        // Check if chromosome exists in genome file
        let current_order = match chrom_order.get(chrom) {
            Some(&order) => order,
            None => {
                return Err(BedError::InvalidFormat(format!(
                    "Chromosome '{}' at line {} not found in genome file",
                    chrom, line_num
                )));
            }
        };

        if let Some(ref pc) = prev_chrom {
            if chrom != pc {
                // Switching chromosomes - check it comes after previous in genome order
                if current_order < prev_chrom_order.unwrap() {
                    return Err(BedError::InvalidFormat(format!(
                        "File not sorted by genome order: chromosome '{}' at line {} should come before '{}'\n\n\
                         Fix: Run 'grit sort -i {} -g <genome.txt>' to sort by genome order.",
                        chrom, line_num, pc, path.as_ref().display()
                    )));
                }
            } else if start < prev_start {
                return Err(BedError::InvalidFormat(format!(
                    "File not sorted: position {} at line {} comes after {} on {}",
                    start, line_num, prev_start, chrom
                )));
            }
        }

        prev_chrom = Some(chrom.to_string());
        prev_chrom_order = Some(current_order);
        prev_start = start;
    }

    Ok(())
}

/// Inline genome-order validator for use within streaming loops.
///
/// Validates that:
/// 1. Chromosomes appear in the order specified by the genome file
/// 2. Within a chromosome, positions are non-decreasing
/// 3. All chromosomes exist in the genome file
#[derive(Debug)]
pub struct GenomeOrderValidator<'a> {
    genome: &'a crate::genome::Genome,
    chrom_order: std::collections::HashMap<String, usize>,
    prev_chrom: Option<String>,
    prev_chrom_order: Option<usize>,
    prev_start: u64,
    record_count: usize,
}

impl<'a> GenomeOrderValidator<'a> {
    /// Create a new genome-order validator.
    pub fn new(genome: &'a crate::genome::Genome) -> Self {
        let chrom_order: std::collections::HashMap<String, usize> = genome
            .chromosomes()
            .enumerate()
            .map(|(i, c)| (c.clone(), i))
            .collect();

        Self {
            genome,
            chrom_order,
            prev_chrom: None,
            prev_chrom_order: None,
            prev_start: 0,
            record_count: 0,
        }
    }

    /// Validate that the given record maintains genome order.
    ///
    /// Returns Ok(()) if valid, Err if out of order or chromosome not in genome.
    #[inline]
    pub fn validate(&mut self, chrom: &str, start: u64) -> Result<(), BedError> {
        self.record_count += 1;

        // Check if chromosome exists in genome file
        let current_order = match self.chrom_order.get(chrom) {
            Some(&order) => order,
            None => {
                return Err(BedError::InvalidFormat(format!(
                    "Chromosome '{}' at record {} not found in genome file",
                    chrom, self.record_count
                )));
            }
        };

        if let Some(ref pc) = self.prev_chrom {
            if chrom != pc {
                // Switching chromosomes - check it comes after previous in genome order
                if current_order < self.prev_chrom_order.unwrap() {
                    return Err(BedError::InvalidFormat(format!(
                        "File not sorted by genome order: chromosome '{}' at record {} should come before '{}'",
                        chrom, self.record_count, pc
                    )));
                }
            } else if start < self.prev_start {
                return Err(BedError::InvalidFormat(format!(
                    "File not sorted: position {} at record {} comes after {} on {}",
                    start, self.record_count, self.prev_start, chrom
                )));
            }
        }

        self.prev_chrom = Some(chrom.to_string());
        self.prev_chrom_order = Some(current_order);
        self.prev_start = start;

        Ok(())
    }

    /// Get the number of records validated.
    pub fn record_count(&self) -> usize {
        self.record_count
    }

    /// Get a reference to the genome.
    pub fn genome(&self) -> &crate::genome::Genome {
        self.genome
    }
}

/// Verify that a reader (e.g., stdin) contains sorted BED data.
///
/// Reads all data into a buffer, validates sort order, and returns the buffer
/// for subsequent processing. This is useful for stdin validation.
///
/// Returns Ok(buffer) if sorted, Err with details if not.
pub fn verify_sorted_reader<R: std::io::Read>(mut reader: R) -> Result<Vec<u8>, BedError> {
    let mut buffer = Vec::new();
    reader.read_to_end(&mut buffer)?;

    let cursor = std::io::Cursor::new(&buffer);
    let bed_reader = BedReader::new(std::io::BufReader::new(cursor));

    let mut prev_chrom: Option<String> = None;
    let mut prev_start: u64 = 0;
    let mut seen_chroms: HashSet<String> = HashSet::new();
    let mut line_num = 0;

    for result in bed_reader.records() {
        let rec = result?;
        line_num += 1;

        let chrom = rec.chrom();
        let start = rec.start();

        if let Some(ref pc) = prev_chrom {
            if chrom != pc {
                // Switching chromosomes - check we haven't seen this one before
                if seen_chroms.contains(chrom) {
                    return Err(BedError::InvalidFormat(format!(
                        "stdin not sorted: chromosome '{}' at line {} was seen earlier (chromosomes must be contiguous)",
                        chrom, line_num
                    )));
                }
                seen_chroms.insert(pc.clone());
            } else if start < prev_start {
                return Err(BedError::InvalidFormat(format!(
                    "stdin not sorted: position {} at line {} comes after {} on {}",
                    start, line_num, prev_start, chrom
                )));
            }
        }

        prev_chrom = Some(chrom.to_string());
        prev_start = start;
    }

    Ok(buffer)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_temp_bed(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_verify_sorted_valid() {
        let file = create_temp_bed("chr1\t100\t200\nchr1\t200\t300\nchr2\t100\t200\n");
        assert!(verify_sorted(file.path()).is_ok());
    }

    #[test]
    fn test_verify_sorted_interleaved_chrom() {
        // chr1 appears, then chr2, then chr1 again - this is invalid
        let file = create_temp_bed("chr1\t100\t200\nchr2\t100\t200\nchr1\t300\t400\n");
        let result = verify_sorted(file.path());
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("not sorted"));
    }

    #[test]
    fn test_verify_sorted_any_chrom_order() {
        // chr2 before chr1 is valid as long as chromosomes are contiguous
        let file = create_temp_bed("chr2\t100\t200\nchr1\t100\t200\n");
        assert!(verify_sorted(file.path()).is_ok());
    }

    #[test]
    fn test_verify_sorted_invalid_position() {
        let file = create_temp_bed("chr1\t200\t300\nchr1\t100\t200\n");
        let result = verify_sorted(file.path());
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("not sorted"));
    }

    #[test]
    fn test_sort_validator() {
        let mut validator = SortValidator::new();
        assert!(validator.validate("chr1", 100).is_ok());
        assert!(validator.validate("chr1", 200).is_ok());
        assert!(validator.validate("chr2", 100).is_ok());
        assert_eq!(validator.record_count(), 3);
    }

    #[test]
    fn test_sort_validator_invalid() {
        let mut validator = SortValidator::new();
        assert!(validator.validate("chr1", 200).is_ok());
        assert!(validator.validate("chr1", 100).is_err());
    }

    fn create_temp_genome(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_verify_sorted_with_genome_valid() {
        // Genome order: chr1, chr2, chr3
        let genome_file = create_temp_genome("chr1\t1000000\nchr2\t500000\nchr3\t250000\n");
        let genome = crate::genome::Genome::from_file(genome_file.path()).unwrap();

        // BED file in genome order
        let bed_file = create_temp_bed("chr1\t100\t200\nchr2\t100\t200\nchr3\t100\t200\n");
        assert!(verify_sorted_with_genome(bed_file.path(), &genome).is_ok());
    }

    #[test]
    fn test_verify_sorted_with_genome_wrong_order() {
        // Genome order: chr1, chr2, chr3
        let genome_file = create_temp_genome("chr1\t1000000\nchr2\t500000\nchr3\t250000\n");
        let genome = crate::genome::Genome::from_file(genome_file.path()).unwrap();

        // BED file NOT in genome order (chr2 before chr1)
        let bed_file = create_temp_bed("chr2\t100\t200\nchr1\t100\t200\n");
        let result = verify_sorted_with_genome(bed_file.path(), &genome);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("genome order"));
    }

    #[test]
    fn test_verify_sorted_with_genome_missing_chrom() {
        // Genome only has chr1, chr2
        let genome_file = create_temp_genome("chr1\t1000000\nchr2\t500000\n");
        let genome = crate::genome::Genome::from_file(genome_file.path()).unwrap();

        // BED file has chr3 which is not in genome
        let bed_file = create_temp_bed("chr1\t100\t200\nchr3\t100\t200\n");
        let result = verify_sorted_with_genome(bed_file.path(), &genome);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("not found in genome"));
    }

    #[test]
    fn test_genome_order_validator() {
        let genome_file = create_temp_genome("chr1\t1000000\nchr2\t500000\nchr3\t250000\n");
        let genome = crate::genome::Genome::from_file(genome_file.path()).unwrap();

        let mut validator = GenomeOrderValidator::new(&genome);
        assert!(validator.validate("chr1", 100).is_ok());
        assert!(validator.validate("chr1", 200).is_ok());
        assert!(validator.validate("chr2", 100).is_ok());
        assert!(validator.validate("chr3", 100).is_ok());
        assert_eq!(validator.record_count(), 4);
    }

    #[test]
    fn test_genome_order_validator_wrong_order() {
        let genome_file = create_temp_genome("chr1\t1000000\nchr2\t500000\nchr3\t250000\n");
        let genome = crate::genome::Genome::from_file(genome_file.path()).unwrap();

        let mut validator = GenomeOrderValidator::new(&genome);
        assert!(validator.validate("chr2", 100).is_ok());
        // chr1 comes after chr2 but should come before in genome order
        assert!(validator.validate("chr1", 100).is_err());
    }
}
