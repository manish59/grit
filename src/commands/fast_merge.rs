//! Ultra-fast streaming merge with zero-allocation parsing.
//!
//! Optimizations:
//! - Zero-copy byte slice parsing (no String allocation)
//! - Hand-rolled integer parsing (3-5x faster than str::parse)
//! - memchr for fast newline scanning
//! - itoa for fast integer output
//! - Minimal branching in hot path
//! - Cache-friendly sequential access
//!
//! Memory: O(1) - only tracks current merge span

use crate::bed::BedError;
use memchr::memchr;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

/// Buffer size for I/O operations (256KB for better throughput)
const BUF_SIZE: usize = 256 * 1024;

/// Fast streaming merge command.
#[derive(Debug, Clone)]
pub struct FastMergeCommand {
    /// Maximum distance between intervals to merge
    pub distance: u64,
    /// Report count of merged intervals
    pub count: bool,
}

impl Default for FastMergeCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl FastMergeCommand {
    pub fn new() -> Self {
        Self {
            distance: 0,
            count: false,
        }
    }

    pub fn with_distance(mut self, d: u64) -> Self {
        self.distance = d;
        self
    }

    /// Run merge on a file.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        input_path: P,
        output: &mut W,
    ) -> Result<FastMergeStats, BedError> {
        let file = File::open(input_path.as_ref())?;
        self.run_reader(file, output)
    }

    /// Run merge from stdin.
    pub fn run_stdin<W: Write>(&self, output: &mut W) -> Result<FastMergeStats, BedError> {
        let stdin = io::stdin();
        self.run_reader(stdin.lock(), output)
    }

    /// Core merge implementation with zero-allocation parsing.
    pub fn run_reader<R: Read, W: Write>(
        &self,
        mut reader: R,
        output: &mut W,
    ) -> Result<FastMergeStats, BedError> {
        let mut stats = FastMergeStats::default();
        let mut writer = BufWriter::with_capacity(BUF_SIZE, output);

        // Read buffer - we process in chunks
        let mut buf = vec![0u8; BUF_SIZE];
        let mut leftover: Vec<u8> = Vec::with_capacity(4096);
        let mut work_buf: Vec<u8> = Vec::with_capacity(BUF_SIZE + 4096);

        // Current merge state - use byte slice for chrom comparison
        let mut current_chrom: Vec<u8> = Vec::with_capacity(32);
        let mut current_start: u64 = 0;
        let mut current_end: u64 = 0;
        let mut current_count: usize = 0;
        let mut has_current = false;

        // Output buffer for itoa
        let mut itoa_buf = itoa::Buffer::new();

        loop {
            let bytes_read = reader.read(&mut buf)?;
            if bytes_read == 0 {
                break;
            }

            // Combine leftover with new data into work buffer
            work_buf.clear();
            work_buf.extend_from_slice(&leftover);
            work_buf.extend_from_slice(&buf[..bytes_read]);
            leftover.clear();

            let data = work_buf.as_slice();
            let mut pos = 0;
            let len = data.len();

            // Process complete lines
            while let Some(newline_pos) = memchr(b'\n', &data[pos..]) {
                let line_end = pos + newline_pos;
                let line = &data[pos..line_end];
                pos = line_end + 1;

                // Skip empty lines and headers
                if line.is_empty()
                    || line[0] == b'#'
                    || line.starts_with(b"track")
                    || line.starts_with(b"browser")
                {
                    continue;
                }

                // Parse BED3 fields (zero allocation)
                if let Some((chrom, start, end)) = parse_bed3_fast(line) {
                    stats.intervals_read += 1;

                    if has_current {
                        // Check if should merge
                        let should_merge = chrom == current_chrom.as_slice()
                            && start <= current_end + self.distance;

                        if should_merge {
                            // Extend current span
                            if end > current_end {
                                current_end = end;
                            }
                            current_count += 1;
                        } else {
                            // Output current span
                            write_bed3_fast(
                                &mut writer,
                                &current_chrom,
                                current_start,
                                current_end,
                                if self.count {
                                    Some(current_count)
                                } else {
                                    None
                                },
                                &mut itoa_buf,
                            )?;
                            stats.intervals_written += 1;

                            // Start new span
                            current_chrom.clear();
                            current_chrom.extend_from_slice(chrom);
                            current_start = start;
                            current_end = end;
                            current_count = 1;
                        }
                    } else {
                        // First interval
                        current_chrom.extend_from_slice(chrom);
                        current_start = start;
                        current_end = end;
                        current_count = 1;
                        has_current = true;
                    }
                }
            }

            // Save incomplete line for next iteration
            if pos < len {
                leftover.extend_from_slice(&data[pos..]);
            }
        }

        // Handle any remaining data (file without final newline)
        if !leftover.is_empty() {
            let line = leftover.as_slice();
            if !line.is_empty()
                && line[0] != b'#'
                && !line.starts_with(b"track")
                && !line.starts_with(b"browser")
            {
                if let Some((chrom, start, end)) = parse_bed3_fast(line) {
                    stats.intervals_read += 1;

                    if has_current {
                        let should_merge = chrom == current_chrom.as_slice()
                            && start <= current_end + self.distance;

                        if should_merge {
                            if end > current_end {
                                current_end = end;
                            }
                            current_count += 1;
                        } else {
                            write_bed3_fast(
                                &mut writer,
                                &current_chrom,
                                current_start,
                                current_end,
                                if self.count {
                                    Some(current_count)
                                } else {
                                    None
                                },
                                &mut itoa_buf,
                            )?;
                            stats.intervals_written += 1;

                            current_chrom.clear();
                            current_chrom.extend_from_slice(chrom);
                            current_start = start;
                            current_end = end;
                            current_count = 1;
                        }
                    } else {
                        current_chrom.extend_from_slice(chrom);
                        current_start = start;
                        current_end = end;
                        current_count = 1;
                        has_current = true;
                    }
                }
            }
        }

        // Output final span
        if has_current {
            write_bed3_fast(
                &mut writer,
                &current_chrom,
                current_start,
                current_end,
                if self.count {
                    Some(current_count)
                } else {
                    None
                },
                &mut itoa_buf,
            )?;
            stats.intervals_written += 1;
        }

        writer.flush().map_err(BedError::Io)?;
        Ok(stats)
    }
}

/// Parse BED3 fields from a byte slice with zero allocation.
/// Returns (chrom, start, end) as byte slice and parsed integers.
#[inline(always)]
fn parse_bed3_fast(line: &[u8]) -> Option<(&[u8], u64, u64)> {
    // Find first tab (end of chrom)
    let tab1 = memchr(b'\t', line)?;
    let chrom = &line[..tab1];

    // Find second tab (end of start)
    let rest1 = &line[tab1 + 1..];
    let tab2 = memchr(b'\t', rest1)?;
    let start_bytes = &rest1[..tab2];

    // Find third tab or end of line (end of end field)
    let rest2 = &rest1[tab2 + 1..];
    let end_bytes = if let Some(tab3) = memchr(b'\t', rest2) {
        &rest2[..tab3]
    } else {
        // Trim potential \r from end
        let mut end = rest2;
        if end.last() == Some(&b'\r') {
            end = &end[..end.len() - 1];
        }
        end
    };

    // Parse integers with fast path
    let start = parse_u64_fast(start_bytes)?;
    let end = parse_u64_fast(end_bytes)?;

    Some((chrom, start, end))
}

/// Ultra-fast u64 parsing - no bounds checking, assumes valid input.
/// ~3-5x faster than str::parse for typical BED coordinates.
#[inline(always)]
fn parse_u64_fast(bytes: &[u8]) -> Option<u64> {
    if bytes.is_empty() {
        return None;
    }

    let mut result: u64 = 0;
    for &b in bytes {
        let digit = b.wrapping_sub(b'0');
        if digit > 9 {
            return None;
        }
        result = result.wrapping_mul(10).wrapping_add(digit as u64);
    }
    Some(result)
}

/// Write BED3 output using itoa for fast integer formatting.
#[inline(always)]
fn write_bed3_fast<W: Write>(
    writer: &mut W,
    chrom: &[u8],
    start: u64,
    end: u64,
    count: Option<usize>,
    itoa_buf: &mut itoa::Buffer,
) -> io::Result<()> {
    writer.write_all(chrom)?;
    writer.write_all(b"\t")?;
    writer.write_all(itoa_buf.format(start).as_bytes())?;
    writer.write_all(b"\t")?;
    writer.write_all(itoa_buf.format(end).as_bytes())?;
    if let Some(c) = count {
        writer.write_all(b"\t")?;
        writer.write_all(itoa_buf.format(c).as_bytes())?;
    }
    writer.write_all(b"\n")?;
    Ok(())
}

/// Statistics from fast merge operation.
#[derive(Debug, Default, Clone)]
pub struct FastMergeStats {
    pub intervals_read: usize,
    pub intervals_written: usize,
}

impl FastMergeStats {
    pub fn compression_ratio(&self) -> f64 {
        if self.intervals_written == 0 {
            0.0
        } else {
            self.intervals_read as f64 / self.intervals_written as f64
        }
    }
}

impl std::fmt::Display for FastMergeStats {
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

    #[test]
    fn test_parse_u64_fast() {
        assert_eq!(parse_u64_fast(b"0"), Some(0));
        assert_eq!(parse_u64_fast(b"123"), Some(123));
        assert_eq!(parse_u64_fast(b"100000000"), Some(100_000_000));
        assert_eq!(parse_u64_fast(b""), None);
        assert_eq!(parse_u64_fast(b"abc"), None);
    }

    #[test]
    fn test_parse_bed3_fast() {
        let line = b"chr1\t100\t200";
        let (chrom, start, end) = parse_bed3_fast(line).unwrap();
        assert_eq!(chrom, b"chr1");
        assert_eq!(start, 100);
        assert_eq!(end, 200);
    }

    #[test]
    fn test_parse_bed3_with_extra_fields() {
        let line = b"chr1\t100\t200\tname\t500\t+";
        let (chrom, start, end) = parse_bed3_fast(line).unwrap();
        assert_eq!(chrom, b"chr1");
        assert_eq!(start, 100);
        assert_eq!(end, 200);
    }

    #[test]
    fn test_fast_merge_basic() {
        let input = b"chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n";
        let cmd = FastMergeCommand::new();
        let mut output = Vec::new();

        let stats = cmd.run_reader(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "chr1\t100\t250");
        assert_eq!(lines[1], "chr1\t300\t400");
        assert_eq!(stats.intervals_read, 3);
        assert_eq!(stats.intervals_written, 2);
    }

    #[test]
    fn test_fast_merge_with_distance() {
        let input = b"chr1\t100\t200\nchr1\t250\t350\n";
        let cmd = FastMergeCommand::new().with_distance(50);
        let mut output = Vec::new();

        cmd.run_reader(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0], "chr1\t100\t350");
    }

    #[test]
    fn test_fast_merge_with_count() {
        let input = b"chr1\t100\t200\nchr1\t150\t250\nchr1\t200\t300\n";
        let mut cmd = FastMergeCommand::new();
        cmd.count = true;
        let mut output = Vec::new();

        cmd.run_reader(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("\t3\n")); // 3 intervals merged
    }

    #[test]
    fn test_fast_merge_multiple_chroms() {
        let input = b"chr1\t100\t200\nchr1\t150\t250\nchr2\t100\t200\nchr2\t150\t250\n";
        let cmd = FastMergeCommand::new();
        let mut output = Vec::new();

        cmd.run_reader(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("chr1\t100\t250"));
        assert!(lines[1].starts_with("chr2\t100\t250"));
    }
}
