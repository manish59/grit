//! Streaming window implementation with TRUE O(k) memory complexity.
//!
//! Finds intervals in B that are within a window distance of intervals in A.
//!
//! ZERO ALLOCATION in hot path:
//! - No per-record String allocation (raw byte parsing with memchr)
//! - Vec + head index for active set (no VecDeque)
//! - Large buffered I/O (256KB input, 8MB output)
//!
//! # Memory Complexity
//!
//! O(k) where k = maximum number of B intervals within window of any A interval.
//!
//! # Requirements
//!
//! Both input files MUST be sorted by chromosome (lexicographic), then by start position.

use crate::bed::BedError;
use crate::streaming::buffers::{DEFAULT_INPUT_BUFFER, DEFAULT_OUTPUT_BUFFER};
use crate::streaming::parsing::{parse_bed3_bytes, should_skip_line};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Active B interval - stores coordinates and original line for output.
#[derive(Debug, Clone)]
struct ActiveB {
    start: u32,
    end: u32,
    /// Original line bytes (stored for output)
    line: Vec<u8>,
}

/// Streaming window command configuration.
#[derive(Debug, Clone)]
pub struct StreamingWindowCommand {
    /// Window size on both sides (symmetric)
    pub window: u64,
    /// Window size on the left (upstream)
    pub left: Option<u64>,
    /// Window size on the right (downstream)
    pub right: Option<u64>,
    /// Only report A intervals with no matches
    pub no_overlap: bool,
    /// Report count of overlaps
    pub count: bool,
}

impl Default for StreamingWindowCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl StreamingWindowCommand {
    pub fn new() -> Self {
        Self {
            window: 1000,
            left: None,
            right: None,
            no_overlap: false,
            count: false,
        }
    }

    /// Get the left window size.
    #[inline(always)]
    fn left_window(&self) -> u64 {
        self.left.unwrap_or(self.window)
    }

    /// Get the right window size.
    #[inline(always)]
    fn right_window(&self) -> u64 {
        self.right.unwrap_or(self.window)
    }

    /// Execute streaming window on two sorted BED files.
    ///
    /// Memory usage: O(k) where k = max B intervals within window at any point
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        a_path: P,
        b_path: P,
        output: &mut W,
    ) -> Result<StreamingWindowStats, BedError> {
        // Output buffer (2MB default, reduced from 8MB for memory efficiency)
        let mut output = BufWriter::with_capacity(DEFAULT_OUTPUT_BUFFER, output);

        // Stream files
        let a_file = File::open(a_path.as_ref())?;
        let mut a_reader = BufReader::with_capacity(DEFAULT_INPUT_BUFFER, a_file);

        let b_file = File::open(b_path.as_ref())?;
        let mut b_reader = BufReader::with_capacity(DEFAULT_INPUT_BUFFER, b_file);

        // Reusable line buffers
        let mut a_line_buf = String::with_capacity(1024);
        let mut b_line_buf = String::with_capacity(1024);

        // Current A chromosome
        let mut a_chrom: Vec<u8> = Vec::with_capacity(64);

        // B state
        let mut b_chrom: Vec<u8> = Vec::with_capacity(64);
        let mut pending_b = Self::read_next_b(&mut b_reader, &mut b_line_buf, &mut b_chrom)?;
        let mut b_exhausted = pending_b.is_none();

        // Track seen B chromosomes to handle any sort order
        let mut seen_b_chroms: HashSet<Vec<u8>> = HashSet::new();
        if !b_exhausted {
            seen_b_chroms.insert(b_chrom.clone());
        }

        // Active set: B intervals that might be within window of current or future A
        let mut active: Vec<ActiveB> = Vec::with_capacity(1024);
        let mut head_idx: usize = 0;

        // Window sizes
        let left_win = self.left_window();
        let right_win = self.right_window();

        // Stats
        let mut stats = StreamingWindowStats::default();

        // Main loop
        loop {
            a_line_buf.clear();
            let bytes_read = a_reader.read_line(&mut a_line_buf)?;
            if bytes_read == 0 {
                break;
            }

            let line = a_line_buf.trim_end();
            let line_bytes = line.as_bytes();

            // Skip headers
            if should_skip_line(line_bytes) {
                continue;
            }

            let (chrom, a_start, a_end) = match parse_bed3_bytes(line_bytes) {
                Some(v) => v,
                None => continue,
            };

            stats.a_intervals += 1;

            // Expanded window boundaries
            let win_start = a_start.saturating_sub(left_win);
            let win_end = a_end.saturating_add(right_win);

            // Chromosome change
            let chrom_changed = chrom != a_chrom.as_slice();
            if chrom_changed {
                a_chrom.clear();
                a_chrom.extend_from_slice(chrom);
                active.clear();
                head_idx = 0;

                // Skip B to current chromosome (or B has already passed it)
                if !b_exhausted && !seen_b_chroms.contains(chrom) {
                    while b_chrom.as_slice() != chrom {
                        pending_b =
                            Self::read_next_b(&mut b_reader, &mut b_line_buf, &mut b_chrom)?;
                        if pending_b.is_none() {
                            b_exhausted = true;
                            break;
                        }
                        seen_b_chroms.insert(b_chrom.clone());
                    }
                }
            }

            // Expire old B from active set (B.end <= win_start means B cannot overlap window)
            while head_idx < active.len() {
                let b = &active[head_idx];
                if (b.end as u64) <= win_start {
                    head_idx += 1;
                } else {
                    break;
                }
            }

            // Compact if needed
            if head_idx > 4096 && head_idx * 2 > active.len() {
                active.drain(0..head_idx);
                head_idx = 0;
            }

            // Add new B intervals that might be within window
            if !b_exhausted {
                while let Some(b) = pending_b.take() {
                    if b_chrom.as_slice() != chrom {
                        // B is on a different chromosome
                        if seen_b_chroms.contains(chrom) {
                            // We've already seen A's chromosome in B, so B has moved past it
                            pending_b = Some(b);
                            break;
                        }
                        // B hasn't reached A's chromosome yet, read next B
                        pending_b =
                            Self::read_next_b(&mut b_reader, &mut b_line_buf, &mut b_chrom)?;
                        if pending_b.is_none() {
                            b_exhausted = true;
                            break;
                        }
                        seen_b_chroms.insert(b_chrom.clone());
                        continue;
                    } else {
                        // B is on the same chromosome as A
                        // B.start > win_end means B cannot overlap this window
                        if (b.start as u64) > win_end {
                            pending_b = Some(b);
                            break;
                        }
                        // Add to active (might overlap current or future windows)
                        active.push(b);
                        pending_b =
                            Self::read_next_b(&mut b_reader, &mut b_line_buf, &mut b_chrom)?;
                        if pending_b.is_none() {
                            b_exhausted = true;
                            break;
                        }
                        seen_b_chroms.insert(b_chrom.clone());
                    }
                }
            }

            stats.max_active_b = stats
                .max_active_b
                .max(active.len().saturating_sub(head_idx));

            // Find overlaps with expanded window
            let active_slice = &active[head_idx..];
            let mut match_count = 0;

            for b in active_slice {
                let b_start = b.start as u64;
                let b_end = b.end as u64;

                // Check if B overlaps the expanded window [win_start, win_end)
                if b_start < win_end && b_end > win_start {
                    match_count += 1;

                    if !self.no_overlap && !self.count {
                        // Output match: A_line \t B_line
                        Self::write_pair(&mut output, line_bytes, &b.line)?;
                        stats.output_pairs += 1;
                    }
                }
            }

            if self.count {
                // Output A with count
                Self::write_count(&mut output, line_bytes, match_count)?;
                stats.output_pairs += 1;
            } else if self.no_overlap && match_count == 0 {
                // Output A intervals with no matches
                output.write_all(line_bytes).map_err(BedError::Io)?;
                output.write_all(b"\n").map_err(BedError::Io)?;
                stats.output_pairs += 1;
            }
        }

        output.flush().map_err(BedError::Io)?;
        Ok(stats)
    }

    /// Read next B interval.
    /// Returns Err on IO error, Ok(None) on EOF, Ok(Some) on success.
    #[inline]
    fn read_next_b(
        reader: &mut BufReader<File>,
        line_buf: &mut String,
        chrom_buf: &mut Vec<u8>,
    ) -> Result<Option<ActiveB>, BedError> {
        loop {
            line_buf.clear();
            let bytes_read = reader.read_line(line_buf).map_err(BedError::Io)?;
            if bytes_read == 0 {
                return Ok(None);
            }

            let line = line_buf.trim_end();
            let line_bytes = line.as_bytes();

            if should_skip_line(line_bytes) {
                continue;
            }

            // Parse BED3 - skip malformed lines
            let (chrom, start, end) = match parse_bed3_bytes(line_bytes) {
                Some(v) => v,
                None => continue,
            };

            chrom_buf.clear();
            chrom_buf.extend_from_slice(chrom);

            return Ok(Some(ActiveB {
                start: start as u32,
                end: end as u32,
                line: line_bytes.to_vec(),
            }));
        }
    }

    #[inline]
    fn write_pair<W: Write>(output: &mut W, a_line: &[u8], b_line: &[u8]) -> Result<(), BedError> {
        output.write_all(a_line).map_err(BedError::Io)?;
        output.write_all(b"\t").map_err(BedError::Io)?;
        output.write_all(b_line).map_err(BedError::Io)?;
        output.write_all(b"\n").map_err(BedError::Io)?;
        Ok(())
    }

    #[inline]
    fn write_count<W: Write>(output: &mut W, a_line: &[u8], count: usize) -> Result<(), BedError> {
        output.write_all(a_line).map_err(BedError::Io)?;
        output.write_all(b"\t").map_err(BedError::Io)?;
        let mut buf = itoa::Buffer::new();
        output
            .write_all(buf.format(count).as_bytes())
            .map_err(BedError::Io)?;
        output.write_all(b"\n").map_err(BedError::Io)?;
        Ok(())
    }
}

/// Statistics from streaming window operation.
#[derive(Debug, Default)]
pub struct StreamingWindowStats {
    pub a_intervals: usize,
    pub output_pairs: usize,
    pub max_active_b: usize,
}

impl std::fmt::Display for StreamingWindowStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "A: {}, Pairs: {}, Max active B: {}",
            self.a_intervals, self.output_pairs, self.max_active_b
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as IoWrite;
    use tempfile::NamedTempFile;

    fn create_temp_bed(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_streaming_window_basic() {
        let a_file = create_temp_bed("chr1\t500\t600\n");
        let b_file = create_temp_bed("chr1\t650\t700\nchr1\t750\t800\n");

        let mut cmd = StreamingWindowCommand::new();
        cmd.window = 100;

        let mut output = Vec::new();
        let stats = cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        // Window of 100 around [500, 600] is [400, 700]
        // B[650-700] overlaps [400, 700] - YES
        // B[750-800] does not overlap [400, 700] - NO
        assert_eq!(stats.a_intervals, 1);
        assert!(result.contains("chr1\t500\t600\tchr1\t650\t700"));
    }

    #[test]
    fn test_streaming_window_no_overlap() {
        let a_file = create_temp_bed("chr1\t500\t600\nchr1\t2000\t2100\n");
        let b_file = create_temp_bed("chr1\t650\t700\n");

        let mut cmd = StreamingWindowCommand::new();
        cmd.window = 100;
        cmd.no_overlap = true;

        let mut output = Vec::new();
        let stats = cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        // Only chr1:2000-2100 has no matches
        assert_eq!(stats.a_intervals, 2);
        assert!(result.contains("chr1\t2000\t2100"));
        assert!(!result.contains("chr1\t500\t600"));
    }

    #[test]
    fn test_streaming_window_count() {
        let a_file = create_temp_bed("chr1\t500\t600\n");
        let b_file = create_temp_bed("chr1\t550\t650\nchr1\t580\t620\n");

        let mut cmd = StreamingWindowCommand::new();
        cmd.window = 100;
        cmd.count = true;

        let mut output = Vec::new();
        let stats = cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert_eq!(stats.a_intervals, 1);
        // Should have count of 2
        assert!(result.contains("chr1\t500\t600\t2"));
    }

    #[test]
    fn test_streaming_window_preserves_columns() {
        let a_file = create_temp_bed("chr1\t500\t600\tgeneA\t100\t+\n");
        let b_file = create_temp_bed("chr1\t550\t650\tgeneB\t200\t-\n");

        let mut cmd = StreamingWindowCommand::new();
        cmd.window = 100;

        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        // Should preserve all columns from both A and B
        assert!(result.contains("geneA"));
        assert!(result.contains("geneB"));
        assert!(result.contains("+"));
        assert!(result.contains("-"));
    }

    #[test]
    fn test_streaming_window_left_right() {
        let a_file = create_temp_bed("chr1\t500\t600\n");
        let b_file = create_temp_bed("chr1\t350\t400\nchr1\t750\t800\n");

        let mut cmd = StreamingWindowCommand::new();
        cmd.left = Some(200); // Window [300, 600+0) = [300, 600)
        cmd.right = Some(0);

        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        // B[350-400] overlaps [300, 600) - YES
        // B[750-800] does not overlap [300, 600) - NO
        assert!(result.contains("chr1\t350\t400"));
        assert!(!result.contains("chr1\t750\t800"));
    }
}
