//! Streaming merge implementation with O(1) memory complexity.
//!
//! Merges overlapping intervals from a sorted BED file without loading
//! the entire file into memory.
//!
//! # Algorithm
//!
//! For sorted input:
//! 1. Read intervals one at a time
//! 2. Track current merge span (chrom, start, end)
//! 3. If next interval overlaps/touches current span, extend it
//! 4. If not, output current span and start new one
//!
//! # Memory Complexity
//!
//! O(1) - only tracks current merge span, regardless of input size.
//!
//! # Requirements
//!
//! Input file MUST be sorted by chromosome, then by start position.

use crate::bed::{BedError, BedReader};
use crate::interval::Strand;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::Path;

/// Streaming merge command configuration.
#[derive(Debug, Clone)]
pub struct StreamingMergeCommand {
    /// Maximum distance between intervals to merge (default: 0)
    pub distance: u64,
    /// Require strand to match for merging
    pub strand_specific: bool,
    /// Report count of merged intervals
    pub count: bool,
}

impl Default for StreamingMergeCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl StreamingMergeCommand {
    pub fn new() -> Self {
        Self {
            distance: 0,
            strand_specific: false,
            count: false,
        }
    }

    /// Set the maximum merge distance.
    pub fn with_distance(mut self, d: u64) -> Self {
        self.distance = d;
        self
    }

    /// Set strand-specific merging.
    pub fn with_strand(mut self, strand: bool) -> Self {
        self.strand_specific = strand;
        self
    }

    /// Execute streaming merge on a sorted BED file.
    ///
    /// Memory usage: O(1) - only tracks current merge span
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        input_path: P,
        output: &mut W,
    ) -> Result<StreamingMergeStats, BedError> {
        let file = File::open(input_path.as_ref())?;
        let reader = BedReader::new(BufReader::with_capacity(64 * 1024, file));
        self.run_streaming(reader, output)
    }

    /// Execute streaming merge from stdin.
    pub fn run_stdin<W: Write>(&self, output: &mut W) -> Result<StreamingMergeStats, BedError> {
        let stdin = io::stdin();
        let reader = BedReader::new(stdin.lock());
        self.run_streaming(reader, output)
    }

    /// Core streaming merge algorithm.
    ///
    /// Maintains only the current merge span in memory.
    pub fn run_streaming<R: io::Read, W: Write>(
        &self,
        reader: BedReader<R>,
        output: &mut W,
    ) -> Result<StreamingMergeStats, BedError> {
        let mut stats = StreamingMergeStats::default();
        let mut writer = BufWriter::with_capacity(64 * 1024, output);

        // Current merge span
        let mut current_chrom: Option<String> = None;
        let mut current_start: u64 = 0;
        let mut current_end: u64 = 0;
        let mut current_strand: Option<Strand> = None;
        let mut current_count: usize = 0;

        for result in reader.records() {
            let rec = result?;
            stats.intervals_read += 1;

            let rec_chrom = rec.chrom();
            let rec_start = rec.start();
            let rec_end = rec.end();
            let rec_strand = rec.strand;

            // Check if we should merge with current span
            let should_merge = if let Some(ref chrom) = current_chrom {
                let same_chrom = chrom == rec_chrom;
                let same_strand = !self.strand_specific || current_strand == rec_strand;
                let overlaps = rec_start <= current_end + self.distance;
                same_chrom && same_strand && overlaps
            } else {
                false
            };

            if should_merge {
                // Extend current span
                current_end = current_end.max(rec_end);
                current_count += 1;
            } else {
                // Output current span if exists
                if let Some(ref chrom) = current_chrom {
                    self.write_span(
                        &mut writer,
                        chrom,
                        current_start,
                        current_end,
                        current_strand,
                        current_count,
                    )?;
                    stats.intervals_written += 1;
                }

                // Start new span
                current_chrom = Some(rec_chrom.to_string());
                current_start = rec_start;
                current_end = rec_end;
                current_strand = rec_strand;
                current_count = 1;
            }
        }

        // Output final span
        if let Some(ref chrom) = current_chrom {
            self.write_span(
                &mut writer,
                chrom,
                current_start,
                current_end,
                current_strand,
                current_count,
            )?;
            stats.intervals_written += 1;
        }

        writer.flush().map_err(BedError::Io)?;
        Ok(stats)
    }

    #[inline]
    fn write_span<W: Write>(
        &self,
        writer: &mut W,
        chrom: &str,
        start: u64,
        end: u64,
        strand: Option<Strand>,
        count: usize,
    ) -> Result<(), BedError> {
        if self.strand_specific {
            if self.count {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}",
                    chrom,
                    start,
                    end,
                    strand
                        .map(|s| s.to_string())
                        .unwrap_or_else(|| ".".to_string()),
                    count
                )
                .map_err(BedError::Io)?;
            } else {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}",
                    chrom,
                    start,
                    end,
                    strand
                        .map(|s| s.to_string())
                        .unwrap_or_else(|| ".".to_string())
                )
                .map_err(BedError::Io)?;
            }
        } else if self.count {
            writeln!(writer, "{}\t{}\t{}\t{}", chrom, start, end, count).map_err(BedError::Io)?;
        } else {
            writeln!(writer, "{}\t{}\t{}", chrom, start, end).map_err(BedError::Io)?;
        }
        Ok(())
    }
}

/// Statistics from streaming merge operation.
#[derive(Debug, Default, Clone)]
pub struct StreamingMergeStats {
    /// Number of intervals read
    pub intervals_read: usize,
    /// Number of merged intervals written
    pub intervals_written: usize,
}

impl StreamingMergeStats {
    /// Compression ratio (how many input intervals per output interval)
    pub fn compression_ratio(&self) -> f64 {
        if self.intervals_written == 0 {
            0.0
        } else {
            self.intervals_read as f64 / self.intervals_written as f64
        }
    }
}

impl std::fmt::Display for StreamingMergeStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Read: {}, Written: {}, Compression: {:.2}x",
            self.intervals_read,
            self.intervals_written,
            self.compression_ratio()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_bed_content(intervals: &[(&str, u64, u64)]) -> String {
        intervals
            .iter()
            .map(|(c, s, e)| format!("{}\t{}\t{}", c, s, e))
            .collect::<Vec<_>>()
            .join("\n")
    }

    #[test]
    fn test_basic_streaming_merge() {
        let content =
            make_bed_content(&[("chr1", 100, 200), ("chr1", 150, 250), ("chr1", 300, 400)]);

        let cmd = StreamingMergeCommand::new();
        let reader = BedReader::new(content.as_bytes());
        let mut output = Vec::new();

        let stats = cmd.run_streaming(reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "chr1\t100\t250");
        assert_eq!(lines[1], "chr1\t300\t400");
        assert_eq!(stats.intervals_read, 3);
        assert_eq!(stats.intervals_written, 2);
    }

    #[test]
    fn test_streaming_merge_with_distance() {
        let content = make_bed_content(&[
            ("chr1", 100, 200),
            ("chr1", 250, 350), // Gap of 50
        ]);

        let cmd = StreamingMergeCommand::new().with_distance(50);
        let reader = BedReader::new(content.as_bytes());
        let mut output = Vec::new();

        cmd.run_streaming(reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0], "chr1\t100\t350");
    }

    #[test]
    fn test_streaming_merge_multiple_chroms() {
        let content = make_bed_content(&[
            ("chr1", 100, 200),
            ("chr1", 150, 250),
            ("chr2", 100, 200),
            ("chr2", 150, 250),
        ]);

        let cmd = StreamingMergeCommand::new();
        let reader = BedReader::new(content.as_bytes());
        let mut output = Vec::new();

        cmd.run_streaming(reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("chr1"));
        assert!(lines[1].starts_with("chr2"));
    }

    #[test]
    fn test_streaming_merge_with_count() {
        let content =
            make_bed_content(&[("chr1", 100, 200), ("chr1", 150, 250), ("chr1", 200, 300)]);

        let mut cmd = StreamingMergeCommand::new();
        cmd.count = true;

        let reader = BedReader::new(content.as_bytes());
        let mut output = Vec::new();

        cmd.run_streaming(reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 1);
        assert!(lines[0].ends_with("\t3")); // 3 intervals merged
    }

    #[test]
    fn test_streaming_merge_touching() {
        let content = make_bed_content(&[
            ("chr1", 100, 200),
            ("chr1", 200, 300), // Touching, should merge
        ]);

        let cmd = StreamingMergeCommand::new();
        let reader = BedReader::new(content.as_bytes());
        let mut output = Vec::new();

        cmd.run_streaming(reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0], "chr1\t100\t300");
    }

    // =============================================================================
    // Strand-specific merge tests
    // =============================================================================

    fn make_bed6_content(intervals: &[(&str, u64, u64, &str, &str, &str)]) -> String {
        intervals
            .iter()
            .map(|(c, s, e, name, score, strand)| {
                format!("{}\t{}\t{}\t{}\t{}\t{}", c, s, e, name, score, strand)
            })
            .collect::<Vec<_>>()
            .join("\n")
    }

    #[test]
    fn test_streaming_merge_strand_same() {
        // Two overlapping intervals with same strand should merge
        let content = make_bed6_content(&[
            ("chr1", 100, 200, ".", ".", "+"),
            ("chr1", 150, 250, ".", ".", "+"),
        ]);

        let mut cmd = StreamingMergeCommand::new();
        cmd.strand_specific = true;

        let reader = BedReader::new(content.as_bytes());
        let mut output = Vec::new();

        cmd.run_streaming(reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 1, "Same strand should merge: {}", result);
        assert!(lines[0].contains("100") && lines[0].contains("250"));
    }

    #[test]
    fn test_streaming_merge_strand_different() {
        // Two overlapping intervals with different strands should NOT merge
        let content = make_bed6_content(&[
            ("chr1", 100, 200, ".", ".", "+"),
            ("chr1", 150, 250, ".", ".", "-"),
        ]);

        let mut cmd = StreamingMergeCommand::new();
        cmd.strand_specific = true;

        let reader = BedReader::new(content.as_bytes());
        let mut output = Vec::new();

        cmd.run_streaming(reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(
            lines.len(),
            2,
            "Different strands should not merge: {}",
            result
        );
    }

    #[test]
    fn test_streaming_merge_strand_with_distance() {
        // Same strand intervals with gap, should merge with distance
        let content = make_bed6_content(&[
            ("chr1", 100, 200, ".", ".", "+"),
            ("chr1", 250, 350, ".", ".", "+"),
        ]);

        let mut cmd = StreamingMergeCommand::new();
        cmd.strand_specific = true;
        cmd.distance = 100;

        let reader = BedReader::new(content.as_bytes());
        let mut output = Vec::new();

        cmd.run_streaming(reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(
            lines.len(),
            1,
            "Same strand within distance should merge: {}",
            result
        );
    }

    #[test]
    fn test_streaming_merge_strand_mixed() {
        // Multiple intervals with mixed strands
        let content = make_bed6_content(&[
            ("chr1", 100, 200, ".", ".", "+"),
            ("chr1", 150, 250, ".", ".", "+"),
            ("chr1", 180, 280, ".", ".", "-"),
            ("chr1", 220, 320, ".", ".", "-"),
        ]);

        let mut cmd = StreamingMergeCommand::new();
        cmd.strand_specific = true;

        let reader = BedReader::new(content.as_bytes());
        let mut output = Vec::new();

        cmd.run_streaming(reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(
            lines.len(),
            2,
            "Should produce 2 merged intervals: {}",
            result
        );
    }

    #[test]
    fn test_streaming_merge_no_strand_column() {
        // Test with strand_specific but no strand column (should use default)
        let content = make_bed_content(&[("chr1", 100, 200), ("chr1", 150, 250)]);

        let mut cmd = StreamingMergeCommand::new();
        cmd.strand_specific = true;

        let reader = BedReader::new(content.as_bytes());
        let mut output = Vec::new();

        // Should still work - missing strand treated as same
        let result = cmd.run_streaming(reader, &mut output);
        assert!(result.is_ok(), "Should handle missing strand column");
    }
}
