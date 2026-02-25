//! Streaming BED file parser.

use crate::config::normalize_end;
use crate::interval::{BedRecord, Interval, Strand};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;
use thiserror::Error;

/// Errors that can occur during BED parsing.
#[derive(Error, Debug)]
pub enum BedError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    #[error("Parse error at line {line}: {message}")]
    Parse { line: usize, message: String },

    #[error("Invalid BED format: {0}")]
    InvalidFormat(String),
}

pub type Result<T> = std::result::Result<T, BedError>;

/// A streaming BED file reader.
pub struct BedReader<R: Read> {
    reader: BufReader<R>,
    line_number: usize,
    buffer: String,
}

impl BedReader<File> {
    /// Open a BED file from a path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        Ok(Self::new(file))
    }
}

impl<R: Read> BedReader<R> {
    /// Create a new BED reader from any readable source.
    pub fn new(reader: R) -> Self {
        Self {
            reader: BufReader::new(reader),
            line_number: 0,
            buffer: String::with_capacity(1024),
        }
    }

    /// Create a BED reader with custom buffer capacity.
    pub fn with_capacity(reader: R, capacity: usize) -> Self {
        Self {
            reader: BufReader::with_capacity(capacity, reader),
            line_number: 0,
            buffer: String::with_capacity(1024),
        }
    }

    /// Read the next BED record.
    pub fn read_record(&mut self) -> Result<Option<BedRecord>> {
        loop {
            self.buffer.clear();
            let bytes_read = self.reader.read_line(&mut self.buffer)?;
            if bytes_read == 0 {
                return Ok(None);
            }
            self.line_number += 1;

            // Skip empty lines and comments
            let line = self.buffer.trim();
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track")
                || line.starts_with("browser")
            {
                continue;
            }

            return self.parse_line(line).map(Some);
        }
    }

    /// Parse a single BED line.
    fn parse_line(&self, line: &str) -> Result<BedRecord> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 3 {
            return Err(BedError::Parse {
                line: self.line_number,
                message: format!("Expected at least 3 fields, got {}", fields.len()),
            });
        }

        let chrom = fields[0].to_string();
        let start = self.parse_position(fields[1], "start")?;
        let end = self.parse_position(fields[2], "end")?;

        if start > end {
            return Err(BedError::Parse {
                line: self.line_number,
                message: format!("Start ({}) > end ({})", start, end),
            });
        }

        // Normalize zero-length intervals if bedtools-compatible mode is enabled
        let end = normalize_end(start, end);

        let mut record = BedRecord::new(chrom, start, end);

        // Parse optional fields
        if fields.len() > 3 {
            record.name = Some(fields[3].to_string());
        }
        if fields.len() > 4 {
            record.score = fields[4].parse().ok();
        }
        if fields.len() > 5 {
            record.strand = fields[5].chars().next().map(Strand::from_char);
        }
        if fields.len() > 6 {
            record.thick_start = fields[6].parse().ok();
        }
        if fields.len() > 7 {
            record.thick_end = fields[7].parse().ok();
        }
        if fields.len() > 8 {
            record.item_rgb = Some(fields[8].to_string());
        }
        if fields.len() > 9 {
            record.block_count = fields[9].parse().ok();
        }
        if fields.len() > 10 {
            record.block_sizes = Some(
                fields[10]
                    .split(',')
                    .filter(|s| !s.is_empty())
                    .filter_map(|s| s.parse().ok())
                    .collect(),
            );
        }
        if fields.len() > 11 {
            record.block_starts = Some(
                fields[11]
                    .split(',')
                    .filter(|s| !s.is_empty())
                    .filter_map(|s| s.parse().ok())
                    .collect(),
            );
        }
        if fields.len() > 12 {
            record.extra_fields = fields[12..].iter().map(|s| s.to_string()).collect();
        }

        Ok(record)
    }

    fn parse_position(&self, s: &str, field_name: &str) -> Result<u64> {
        s.parse().map_err(|_| BedError::Parse {
            line: self.line_number,
            message: format!("Invalid {} position: '{}'", field_name, s),
        })
    }

    /// Get an iterator over all records.
    pub fn records(self) -> BedRecordIter<R> {
        BedRecordIter { reader: self }
    }
}

/// Iterator over BED records.
pub struct BedRecordIter<R: Read> {
    reader: BedReader<R>,
}

impl<R: Read> Iterator for BedRecordIter<R> {
    type Item = Result<BedRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// Read all intervals from a BED file.
pub fn read_intervals<P: AsRef<Path>>(path: P) -> Result<Vec<Interval>> {
    let reader = BedReader::from_path(path)?;
    reader
        .records()
        .map(|r| r.map(|rec| rec.interval))
        .collect()
}

/// Read all BED records from a file.
pub fn read_records<P: AsRef<Path>>(path: P) -> Result<Vec<BedRecord>> {
    let reader = BedReader::from_path(path)?;
    reader.records().collect()
}

/// Parse intervals from a string (useful for testing).
pub fn parse_intervals(content: &str) -> Result<Vec<Interval>> {
    let reader = BedReader::new(content.as_bytes());
    reader
        .records()
        .map(|r| r.map(|rec| rec.interval))
        .collect()
}

/// Write intervals to a writer.
pub fn write_intervals<W: io::Write>(writer: &mut W, intervals: &[Interval]) -> io::Result<()> {
    for interval in intervals {
        writeln!(writer, "{}", interval)?;
    }
    Ok(())
}

/// Write BED records to a writer.
pub fn write_records<W: io::Write>(writer: &mut W, records: &[BedRecord]) -> io::Result<()> {
    for record in records {
        writeln!(writer, "{}", record)?;
    }
    Ok(())
}

/// Fast line parser using memchr for performance.
pub struct FastBedParser;

impl FastBedParser {
    pub fn new() -> Self {
        Self
    }

    /// Parse a line into an interval (BED3 only, for maximum speed).
    #[inline]
    pub fn parse_interval(&self, line: &[u8]) -> Option<Interval> {
        let mut fields = line.split(|&b| b == b'\t');

        let chrom = std::str::from_utf8(fields.next()?).ok()?;
        let start: u64 = std::str::from_utf8(fields.next()?).ok()?.parse().ok()?;
        let end: u64 = std::str::from_utf8(fields.next()?).ok()?.parse().ok()?;

        // Normalize zero-length intervals if bedtools-compatible mode is enabled
        let end = normalize_end(start, end);

        Some(Interval::new(chrom, start, end))
    }
}

impl Default for FastBedParser {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_bed3() {
        let content = "chr1\t100\t200\nchr1\t300\t400\n";
        let intervals = parse_intervals(content).unwrap();

        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].chrom, "chr1");
        assert_eq!(intervals[0].start, 100);
        assert_eq!(intervals[0].end, 200);
    }

    #[test]
    fn test_parse_bed6() {
        let content = "chr1\t100\t200\tgene1\t500\t+\n";
        let reader = BedReader::new(content.as_bytes());
        let records: Vec<_> = reader.records().collect::<Result<_>>().unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, Some("gene1".to_string()));
        assert_eq!(records[0].score, Some(500.0));
        assert_eq!(records[0].strand, Some(Strand::Plus));
    }

    #[test]
    fn test_skip_comments() {
        let content = "# comment\nchr1\t100\t200\n";
        let intervals = parse_intervals(content).unwrap();

        assert_eq!(intervals.len(), 1);
    }

    #[test]
    fn test_skip_track_lines() {
        let content = "track name=test\nbrowser position chr1:1-1000\nchr1\t100\t200\n";
        let intervals = parse_intervals(content).unwrap();

        assert_eq!(intervals.len(), 1);
    }

    #[test]
    fn test_invalid_bed() {
        let content = "chr1\t100\n"; // Only 2 fields
        let result = parse_intervals(content);
        assert!(result.is_err());
    }

    #[test]
    fn test_fast_parser() {
        let parser = FastBedParser::new();
        let line = b"chr1\t100\t200";
        let interval = parser.parse_interval(line).unwrap();

        assert_eq!(interval.chrom, "chr1");
        assert_eq!(interval.start, 100);
        assert_eq!(interval.end, 200);
    }
}
