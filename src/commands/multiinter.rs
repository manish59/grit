//! Multiinter command implementation.
//!
//! Identifies common intervals across multiple BED files using sweep-line algorithm.
//! O(n log n) for sorting events, O(n) for sweep.

use crate::bed::{BedError, BedReader};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Multiinter command configuration.
#[derive(Debug, Clone)]
pub struct MultiinterCommand {
    /// Include header in output
    pub header: bool,
    /// Only report intervals present in all files
    pub cluster: bool,
    /// Empty placeholder for missing files
    pub empty: bool,
}

impl Default for MultiinterCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl MultiinterCommand {
    pub fn new() -> Self {
        Self {
            header: false,
            cluster: false,
            empty: false,
        }
    }

    /// Run multiinter on multiple files.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        inputs: &[P],
        output: &mut W,
    ) -> Result<(), BedError> {
        // Read all intervals from all files
        let mut all_intervals: Vec<Vec<(String, u64, u64)>> = Vec::with_capacity(inputs.len());

        for input in inputs {
            let file = File::open(input)?;
            let reader = BedReader::new(file);
            let mut intervals = Vec::new();

            for result in reader.records() {
                let record = result?;
                intervals.push((record.chrom().to_string(), record.start(), record.end()));
            }

            all_intervals.push(intervals);
        }

        self.multiinter_from_intervals(&all_intervals, output)
    }

    /// Run multiinter on pre-loaded intervals.
    pub fn multiinter_from_intervals<W: Write>(
        &self,
        all_intervals: &[Vec<(String, u64, u64)>],
        output: &mut W,
    ) -> Result<(), BedError> {
        let mut buf_output = BufWriter::with_capacity(256 * 1024, output);
        let n_files = all_intervals.len();

        // Group by chromosome
        let mut by_chrom: std::collections::HashMap<String, Vec<(u64, u64, usize)>> =
            std::collections::HashMap::new();

        for (file_idx, intervals) in all_intervals.iter().enumerate() {
            for (chrom, start, end) in intervals {
                by_chrom
                    .entry(chrom.clone())
                    .or_default()
                    .push((*start, *end, file_idx));
            }
        }

        // Get sorted chromosome list
        let mut chroms: Vec<_> = by_chrom.keys().cloned().collect();
        chroms.sort();

        // Process each chromosome
        for chrom in chroms {
            let intervals = by_chrom.get(&chrom).unwrap();
            self.process_chromosome(&chrom, intervals, n_files, &mut buf_output)?;
        }

        buf_output.flush().map_err(BedError::Io)?;
        Ok(())
    }

    /// Process a single chromosome using sweep-line algorithm.
    fn process_chromosome<W: Write>(
        &self,
        chrom: &str,
        intervals: &[(u64, u64, usize)],
        n_files: usize,
        output: &mut W,
    ) -> Result<(), BedError> {
        // Create events: (position, is_start, file_idx)
        // is_start: true for start, false for end
        let mut events: Vec<(u64, bool, usize)> = Vec::with_capacity(intervals.len() * 2);

        for &(start, end, file_idx) in intervals {
            events.push((start, true, file_idx));
            events.push((end, false, file_idx));
        }

        // Sort events: by position, then ends before starts at same position
        events.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        // Track which files have coverage at current position
        let mut file_depths: Vec<u32> = vec![0; n_files];
        let mut prev_pos: u64 = 0;
        let mut has_coverage = false;

        for (pos, is_start, file_idx) in events {
            // Output region if there was coverage
            if pos > prev_pos && has_coverage {
                self.output_region(chrom, prev_pos, pos, &file_depths, output)?;
            }

            // Update depth
            if is_start {
                file_depths[file_idx] += 1;
            } else {
                file_depths[file_idx] = file_depths[file_idx].saturating_sub(1);
            }

            // Check if any file has coverage
            has_coverage = file_depths.iter().any(|&d| d > 0);
            prev_pos = pos;
        }

        Ok(())
    }

    /// Output a region with coverage info.
    fn output_region<W: Write>(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        file_depths: &[u32],
        output: &mut W,
    ) -> Result<(), BedError> {
        // Count files with coverage
        let count: usize = file_depths.iter().filter(|&&d| d > 0).count();

        // Skip if cluster mode and not all files
        if self.cluster && count != file_depths.len() {
            return Ok(());
        }

        // Build list of file indices (1-based)
        let file_list: Vec<String> = file_depths
            .iter()
            .enumerate()
            .filter(|(_, &d)| d > 0)
            .map(|(i, _)| (i + 1).to_string())
            .collect();

        // Build presence flags
        let flags: Vec<&str> = file_depths
            .iter()
            .map(|&d| if d > 0 { "1" } else { "0" })
            .collect();

        writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{}",
            chrom,
            start,
            end,
            count,
            file_list.join(","),
            flags.join("\t")
        )
        .map_err(BedError::Io)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiinter_basic() {
        let cmd = MultiinterCommand::new();

        let file1 = vec![
            ("chr1".to_string(), 100u64, 200u64),
            ("chr1".to_string(), 300, 400),
        ];
        let file2 = vec![
            ("chr1".to_string(), 150, 250),
            ("chr1".to_string(), 350, 450),
        ];

        let all = vec![file1, file2];

        let mut output = Vec::new();
        cmd.multiinter_from_intervals(&all, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();

        // Should have multiple regions
        assert!(!lines.is_empty());

        // Check first region
        let first: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(first[0], "chr1");
        assert_eq!(first[1], "100");
        assert_eq!(first[2], "150");
        assert_eq!(first[3], "1"); // count
        assert_eq!(first[4], "1"); // file list
    }

    #[test]
    fn test_multiinter_three_files() {
        let cmd = MultiinterCommand::new();

        let file1 = vec![("chr1".to_string(), 100u64, 200u64)];
        let file2 = vec![("chr1".to_string(), 150, 250)];
        let file3 = vec![("chr1".to_string(), 180, 220)];

        let all = vec![file1, file2, file3];

        let mut output = Vec::new();
        cmd.multiinter_from_intervals(&all, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();

        // Find the region where all three overlap (180-200)
        let all_three = lines.iter().find(|line| {
            let parts: Vec<&str> = line.split('\t').collect();
            parts[3] == "3" // count == 3
        });

        assert!(all_three.is_some());
        let parts: Vec<&str> = all_three.unwrap().split('\t').collect();
        assert_eq!(parts[1], "180");
        assert_eq!(parts[2], "200");
    }

    #[test]
    fn test_multiinter_no_overlap() {
        let cmd = MultiinterCommand::new();

        let file1 = vec![("chr1".to_string(), 100u64, 200u64)];
        let file2 = vec![("chr1".to_string(), 300, 400)];

        let all = vec![file1, file2];

        let mut output = Vec::new();
        cmd.multiinter_from_intervals(&all, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();

        // Should have 2 regions with count=1 each
        assert_eq!(lines.len(), 2);

        for line in lines {
            let parts: Vec<&str> = line.split('\t').collect();
            assert_eq!(parts[3], "1"); // count == 1
        }
    }
}
