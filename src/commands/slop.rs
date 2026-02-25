//! Slop command implementation.
//!
//! Extends intervals by a fixed number of bases upstream and/or downstream,
//! respecting chromosome boundaries.

use crate::bed::{BedError, BedReader};
use crate::genome::Genome;
use crate::interval::BedRecord;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

/// Slop command configuration.
#[derive(Debug, Clone)]
pub struct SlopCommand {
    /// Number of bases to extend on both sides (if left/right not specified)
    /// When pct=true, this is interpreted as a fraction (0.0-1.0)
    pub both: f64,
    /// Number of bases to extend on the left (upstream for + strand)
    pub left: Option<f64>,
    /// Number of bases to extend on the right (downstream for + strand)
    pub right: Option<f64>,
    /// Use strand information (left=upstream, right=downstream relative to strand)
    pub strand: bool,
    /// Use fraction of interval size instead of fixed bases
    pub pct: bool,
    /// Treat the slop values as header lines to skip
    pub header: bool,
}

impl Default for SlopCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl SlopCommand {
    pub fn new() -> Self {
        Self {
            both: 0.0,
            left: None,
            right: None,
            strand: false,
            pct: false,
            header: false,
        }
    }

    /// Get effective left extension.
    #[inline]
    fn get_left(&self, interval_len: u64) -> u64 {
        let base = self.left.unwrap_or(self.both);
        if self.pct {
            ((interval_len as f64) * base).round() as u64
        } else {
            base as u64
        }
    }

    /// Get effective right extension.
    #[inline]
    fn get_right(&self, interval_len: u64) -> u64 {
        let base = self.right.unwrap_or(self.both);
        if self.pct {
            ((interval_len as f64) * base).round() as u64
        } else {
            base as u64
        }
    }

    /// Apply slop to a single record.
    #[inline]
    pub fn slop_record(&self, record: &mut BedRecord, chrom_size: u64) {
        let interval_len = record.end() - record.start();
        let left_ext = self.get_left(interval_len);
        let right_ext = self.get_right(interval_len);

        // Handle strand-aware slop
        let (upstream, downstream) = if self.strand {
            match record.strand {
                Some(crate::interval::Strand::Minus) => (right_ext, left_ext),
                _ => (left_ext, right_ext),
            }
        } else {
            (left_ext, right_ext)
        };

        // Apply slop with boundary enforcement
        let new_start = record.start().saturating_sub(upstream);
        let new_end = (record.end() + downstream).min(chrom_size);

        record.interval.start = new_start;
        record.interval.end = new_end;
    }

    /// Run slop on a file with streaming output.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        input: P,
        genome: &Genome,
        output: &mut W,
    ) -> Result<(), BedError> {
        let file = std::fs::File::open(input)?;
        let reader = BedReader::new(file);
        self.slop_streaming(reader, genome, output)
    }

    /// Streaming slop processing.
    pub fn slop_streaming<R: Read, W: Write>(
        &self,
        reader: BedReader<R>,
        genome: &Genome,
        output: &mut W,
    ) -> Result<(), BedError> {
        let mut buf_output = BufWriter::with_capacity(256 * 1024, output);

        for result in reader.records() {
            let mut record = result?;

            // Get chromosome size, skip if not in genome
            let chrom_size = match genome.chrom_size(record.chrom()) {
                Some(size) => size,
                None => {
                    // bedtools skips intervals on unknown chromosomes
                    continue;
                }
            };

            self.slop_record(&mut record, chrom_size);

            // Only output if interval is valid (start < end)
            if record.start() < record.end() {
                writeln!(buf_output, "{}", record).map_err(BedError::Io)?;
            }
        }

        buf_output.flush().map_err(BedError::Io)?;
        Ok(())
    }

    /// Run slop from stdin to stdout.
    pub fn run_stdio(&self, genome: &Genome) -> Result<(), BedError> {
        let stdin = io::stdin();
        let reader = BedReader::new(stdin.lock());

        let stdout = io::stdout();
        let handle = stdout.lock();

        self.slop_streaming(reader, genome, &mut BufWriter::new(handle))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::interval::{BedRecord, Strand};

    fn make_record(chrom: &str, start: u64, end: u64) -> BedRecord {
        BedRecord::new(chrom, start, end)
    }

    fn make_stranded_record(chrom: &str, start: u64, end: u64, strand: Strand) -> BedRecord {
        let mut rec = BedRecord::new(chrom, start, end);
        rec.strand = Some(strand);
        rec
    }

    #[test]
    fn test_slop_both_sides() {
        let cmd = SlopCommand {
            both: 10.0,
            ..SlopCommand::new()
        };

        let mut rec = make_record("chr1", 100, 200);
        cmd.slop_record(&mut rec, 1000);

        assert_eq!(rec.start(), 90);
        assert_eq!(rec.end(), 210);
    }

    #[test]
    fn test_slop_left_right() {
        let cmd = SlopCommand {
            left: Some(5.0),
            right: Some(15.0),
            ..SlopCommand::new()
        };

        let mut rec = make_record("chr1", 100, 200);
        cmd.slop_record(&mut rec, 1000);

        assert_eq!(rec.start(), 95);
        assert_eq!(rec.end(), 215);
    }

    #[test]
    fn test_slop_boundary_left() {
        let cmd = SlopCommand {
            both: 100.0,
            ..SlopCommand::new()
        };

        let mut rec = make_record("chr1", 50, 150);
        cmd.slop_record(&mut rec, 1000);

        assert_eq!(rec.start(), 0); // Clamped at 0
        assert_eq!(rec.end(), 250);
    }

    #[test]
    fn test_slop_boundary_right() {
        let cmd = SlopCommand {
            both: 100.0,
            ..SlopCommand::new()
        };

        let mut rec = make_record("chr1", 900, 950);
        cmd.slop_record(&mut rec, 1000);

        assert_eq!(rec.start(), 800);
        assert_eq!(rec.end(), 1000); // Clamped at chrom size
    }

    #[test]
    fn test_slop_strand_plus() {
        let cmd = SlopCommand {
            left: Some(10.0),
            right: Some(20.0),
            strand: true,
            ..SlopCommand::new()
        };

        let mut rec = make_stranded_record("chr1", 100, 200, Strand::Plus);
        cmd.slop_record(&mut rec, 1000);

        // For + strand: left=upstream, right=downstream
        assert_eq!(rec.start(), 90); // -10 upstream
        assert_eq!(rec.end(), 220); // +20 downstream
    }

    #[test]
    fn test_slop_strand_minus() {
        let cmd = SlopCommand {
            left: Some(10.0),
            right: Some(20.0),
            strand: true,
            ..SlopCommand::new()
        };

        let mut rec = make_stranded_record("chr1", 100, 200, Strand::Minus);
        cmd.slop_record(&mut rec, 1000);

        // For - strand: left=downstream (at end), right=upstream (at start)
        assert_eq!(rec.start(), 80); // -20 (right becomes upstream)
        assert_eq!(rec.end(), 210); // +10 (left becomes downstream)
    }

    #[test]
    fn test_slop_percentage() {
        let cmd = SlopCommand {
            both: 0.0,
            left: Some(1.0), // 100% when pct=true
            right: Some(1.0),
            pct: true,
            ..SlopCommand::new()
        };

        let mut rec = make_record("chr1", 100, 200); // length = 100
        cmd.slop_record(&mut rec, 1000);

        assert_eq!(rec.start(), 0); // 100 - 100 = 0
        assert_eq!(rec.end(), 300); // 200 + 100 = 300
    }
}
