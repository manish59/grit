//! Streaming coverage command - TRUE O(k) memory coverage computation.
//!
//! Memory complexity: O(k) where k = max overlapping B intervals at any position.
//! Both A and B are streamed - neither file is fully loaded into memory.
//!
//! ZERO ALLOCATION in hot path:
//! - No per-B chromosome allocation (chrom buffer reused)
//! - No temporary vectors in basic coverage mode
//! - No sorting (active set is naturally start-sorted)
//! - Float formatting uses ryu (stack-allocated buffer)
//!
//! REQUIREMENT: Both A and B must be sorted by (chrom, start) in same order.
//! Use `--assume-sorted` flag or pre-sort with `grit sort`.

use crate::bed::BedError;
use crate::streaming::parsing::{parse_bed3_bytes, should_skip_line};
use crate::streaming::ActiveInterval;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Streaming coverage command configuration.
#[derive(Debug, Clone)]
pub struct StreamingCoverageCommand {
    pub histogram: bool,
    pub per_base: bool,
    pub mean: bool,
}

impl Default for StreamingCoverageCommand {
    fn default() -> Self {
        Self::new()
    }
}

/// Pending B interval - stores coordinates only, chrom tracked separately.
#[derive(Debug, Clone, Copy)]
struct PendingB {
    start: u32,
    end: u32,
}

impl StreamingCoverageCommand {
    pub fn new() -> Self {
        Self {
            histogram: false,
            per_base: false,
            mean: false,
        }
    }

    /// Execute TRUE O(k) streaming coverage.
    ///
    /// Memory: O(k) where k = max overlapping B intervals.
    /// Both A and B files are streamed - never fully loaded.
    ///
    /// REQUIREMENT: Both files must be sorted by (chrom, start) in lexicographic order.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        a_path: P,
        b_path: P,
        output: &mut W,
    ) -> Result<(), BedError> {
        // Large output buffer (8MB)
        let mut output = BufWriter::with_capacity(8 * 1024 * 1024, output);

        // Stream A file
        let a_file = File::open(&a_path)?;
        let mut a_reader = BufReader::with_capacity(256 * 1024, a_file);

        // Stream B file
        let b_file = File::open(&b_path)?;
        let mut b_reader = BufReader::with_capacity(256 * 1024, b_file);

        // Reusable line buffers (no per-line allocation)
        let mut a_line_buf = String::with_capacity(1024);
        let mut b_line_buf = String::with_capacity(1024);

        // Current A chromosome (reused buffer)
        let mut a_chrom: Vec<u8> = Vec::with_capacity(64);

        // Pending B record: chrom stored separately, only (start, end) in struct
        let mut b_chrom: Vec<u8> = Vec::with_capacity(64);
        let mut pending_b = Self::read_next_b(&mut b_reader, &mut b_line_buf, &mut b_chrom)?;
        let mut b_exhausted = pending_b.is_none();

        // Active set: Vec with head index (no VecDeque, no make_contiguous)
        let mut active: Vec<ActiveInterval> = Vec::with_capacity(1024);
        let mut head_idx: usize = 0;

        // itoa buffer for fast integer formatting
        let mut itoa_buf = itoa::Buffer::new();

        // Reusable event buffer for mean/histogram modes
        let mut events_buf: Vec<(u64, i32)> = Vec::with_capacity(2048);

        // Main loop: stream A records
        loop {
            a_line_buf.clear();
            let bytes_read = a_reader.read_line(&mut a_line_buf)?;
            if bytes_read == 0 {
                break;
            }

            let line = a_line_buf.trim_end();
            let line_bytes = line.as_bytes();

            // Skip empty lines and headers
            if should_skip_line(line_bytes) {
                continue;
            }

            // Parse A record (zero allocation)
            let (chrom, a_start, a_end) = match parse_bed3_bytes(line_bytes) {
                Some(v) => v,
                None => continue,
            };

            let a_len = a_end.saturating_sub(a_start);

            // Check chromosome change
            let chrom_changed = chrom != a_chrom.as_slice();
            if chrom_changed {
                // Update A chromosome (reuses buffer)
                a_chrom.clear();
                a_chrom.extend_from_slice(chrom);

                // Clear active set on chromosome change
                active.clear();
                head_idx = 0;

                // Skip B records until we reach this chromosome or exhaust B.
                // NOTE: We use != instead of < to handle both lexicographic and genome sort orders.
                // Both A and B must be sorted in the SAME order, but that order can be either.
                if !b_exhausted {
                    while b_chrom.as_slice() != chrom {
                        pending_b =
                            Self::read_next_b(&mut b_reader, &mut b_line_buf, &mut b_chrom)?;
                        if pending_b.is_none() {
                            b_exhausted = true;
                            break;
                        }
                    }
                }
            }

            // Handle zero-length intervals
            if a_len == 0 {
                self.write_zero_coverage(&mut output, line, 0, &mut itoa_buf)?;
                continue;
            }

            // Step 1: Remove expired B intervals (head index advancement)
            while head_idx < active.len() && (active[head_idx].end as u64) <= a_start {
                head_idx += 1;
            }

            // Periodic compaction: avoid unbounded head_idx growth
            if head_idx > 4096 && head_idx * 2 > active.len() {
                active.drain(0..head_idx);
                head_idx = 0;
            }

            // Step 2: Add new B intervals to active set
            // B intervals are added in START-SORTED order (B file is sorted)
            // NOTE: We only use == comparison, not < or >, to handle both
            // lexicographic and genome sort orders.
            if !b_exhausted {
                while let Some(b) = pending_b {
                    // Check chromosome: B must match A's chromosome
                    if b_chrom.as_slice() == chrom {
                        // Same chromosome - check position
                        if (b.start as u64) >= a_end {
                            // B starts at or after A ends - defer to future A
                            break;
                        }
                        // Add to active set (no chrom stored, just start/end)
                        active.push(ActiveInterval {
                            start: b.start,
                            end: b.end,
                        });
                        // Read next B
                        pending_b =
                            Self::read_next_b(&mut b_reader, &mut b_line_buf, &mut b_chrom)?;
                        if pending_b.is_none() {
                            b_exhausted = true;
                            break;
                        }
                    } else {
                        // Different chromosome - stop adding for this A interval.
                        // B will be processed when A reaches B's chromosome.
                        break;
                    }
                }
            }

            // Step 3: Compute coverage from active slice (ZERO ALLOCATION)
            let active_slice = &active[head_idx..];

            // Step 4: Output based on mode
            if self.per_base {
                self.write_per_base_coverage(
                    &mut output,
                    line,
                    a_start,
                    a_end,
                    active_slice,
                    &mut events_buf,
                )?;
            } else if self.histogram {
                self.write_histogram_coverage(
                    &mut output,
                    line,
                    a_start,
                    a_end,
                    a_len,
                    active_slice,
                    &mut events_buf,
                )?;
            } else if self.mean {
                self.write_mean_coverage(
                    &mut output,
                    line,
                    a_start,
                    a_end,
                    a_len,
                    active_slice,
                    &mut events_buf,
                )?;
            } else {
                // Basic coverage - most common path, ZERO ALLOCATION
                let (num_overlaps, bases_covered) =
                    Self::compute_coverage_inline(active_slice, a_start, a_end);

                self.write_basic_coverage_fast(
                    &mut output,
                    line,
                    num_overlaps,
                    bases_covered,
                    a_len,
                    &mut itoa_buf,
                )?;
            }
        }

        output.flush()?;
        Ok(())
    }

    /// Read next B interval. Updates b_chrom buffer in place.
    /// Returns Err on IO error, Ok(None) on EOF, Ok(Some) on success.
    /// ZERO ALLOCATION per call (reuses buffers).
    #[inline]
    fn read_next_b(
        reader: &mut BufReader<File>,
        line_buf: &mut String,
        chrom_buf: &mut Vec<u8>,
    ) -> Result<Option<PendingB>, BedError> {
        loop {
            line_buf.clear();
            let bytes_read = reader.read_line(line_buf).map_err(BedError::Io)?;
            if bytes_read == 0 {
                return Ok(None);
            }

            let line = line_buf.trim_end().as_bytes();

            // Skip empty lines and headers
            if should_skip_line(line) {
                continue;
            }

            // Parse BED3 - skip malformed lines
            let (chrom, start, end) = match parse_bed3_bytes(line) {
                Some(v) => v,
                None => continue,
            };

            // Update chromosome buffer (reuses allocation)
            chrom_buf.clear();
            chrom_buf.extend_from_slice(chrom);

            return Ok(Some(PendingB {
                start: start as u32,
                end: end as u32,
            }));
        }
    }

    /// Compute coverage INLINE without any allocation.
    ///
    /// INVARIANT: Active intervals are sorted by START because:
    /// 1. B file is sorted by (chrom, start)
    /// 2. Intervals are added to active set in file order
    /// 3. Therefore active set is naturally start-sorted
    ///
    /// This enables O(n) single-pass union computation.
    #[inline]
    fn compute_coverage_inline(
        active: &[ActiveInterval],
        a_start: u64,
        a_end: u64,
    ) -> (usize, u64) {
        let mut num_overlaps: usize = 0;
        let mut bases_covered: u64 = 0;

        // Union tracking - stack variables only, no allocation
        let mut union_start: u64 = 0;
        let mut union_end: u64 = 0;
        let mut in_union = false;

        for b in active {
            let b_start = b.start as u64;
            let b_end = b.end as u64;

            // Check overlap: B.end > A.start AND B.start < A.end
            if b_end > a_start && b_start < a_end {
                num_overlaps += 1;

                // Clip to A bounds
                let clip_start = b_start.max(a_start);
                let clip_end = b_end.min(a_end);

                if !in_union {
                    union_start = clip_start;
                    union_end = clip_end;
                    in_union = true;
                } else if clip_start > union_end {
                    // Gap - finalize previous span
                    bases_covered += union_end - union_start;
                    union_start = clip_start;
                    union_end = clip_end;
                } else {
                    // Extend union
                    union_end = union_end.max(clip_end);
                }
            }
        }

        // Finalize last span
        if in_union {
            bases_covered += union_end - union_start;
        }

        (num_overlaps, bases_covered)
    }

    /// Fast basic coverage output using itoa (ZERO ALLOCATION for integers).
    #[inline]
    fn write_basic_coverage_fast<W: Write>(
        &self,
        output: &mut W,
        original_line: &str,
        num_overlaps: usize,
        bases_covered: u64,
        a_len: u64,
        itoa_buf: &mut itoa::Buffer,
    ) -> Result<(), BedError> {
        // Use f32 to match bedtools precision (bedtools uses float internally)
        let fraction: f32 = if a_len > 0 {
            bases_covered as f32 / a_len as f32
        } else {
            0.0
        };

        output
            .write_all(original_line.as_bytes())
            .map_err(BedError::Io)?;
        output.write_all(b"\t").map_err(BedError::Io)?;
        output
            .write_all(itoa_buf.format(num_overlaps).as_bytes())
            .map_err(BedError::Io)?;
        output.write_all(b"\t").map_err(BedError::Io)?;
        output
            .write_all(itoa_buf.format(bases_covered).as_bytes())
            .map_err(BedError::Io)?;
        output.write_all(b"\t").map_err(BedError::Io)?;
        output
            .write_all(itoa_buf.format(a_len).as_bytes())
            .map_err(BedError::Io)?;

        // Format with {:.7} to match bedtools (uses C printf rounding)
        writeln!(output, "\t{:.7}", fraction).map_err(BedError::Io)?;

        Ok(())
    }

    /// Write zero coverage output.
    #[inline]
    fn write_zero_coverage<W: Write>(
        &self,
        output: &mut W,
        original_line: &str,
        a_len: u64,
        itoa_buf: &mut itoa::Buffer,
    ) -> Result<(), BedError> {
        if self.per_base {
            for pos in 1..=a_len {
                output
                    .write_all(original_line.as_bytes())
                    .map_err(BedError::Io)?;
                output.write_all(b"\t").map_err(BedError::Io)?;
                output
                    .write_all(itoa_buf.format(pos).as_bytes())
                    .map_err(BedError::Io)?;
                output.write_all(b"\t0\n").map_err(BedError::Io)?;
            }
        } else if self.histogram {
            output
                .write_all(original_line.as_bytes())
                .map_err(BedError::Io)?;
            output.write_all(b"\t0\t").map_err(BedError::Io)?;
            output
                .write_all(itoa_buf.format(a_len).as_bytes())
                .map_err(BedError::Io)?;
            output.write_all(b"\t").map_err(BedError::Io)?;
            output
                .write_all(itoa_buf.format(a_len).as_bytes())
                .map_err(BedError::Io)?;
            output.write_all(b"\t1.0000000\n").map_err(BedError::Io)?;
        } else if self.mean {
            output
                .write_all(original_line.as_bytes())
                .map_err(BedError::Io)?;
            output.write_all(b"\t0.0000000\n").map_err(BedError::Io)?;
        } else {
            output
                .write_all(original_line.as_bytes())
                .map_err(BedError::Io)?;
            output.write_all(b"\t0\t0\t").map_err(BedError::Io)?;
            output
                .write_all(itoa_buf.format(a_len).as_bytes())
                .map_err(BedError::Io)?;
            output.write_all(b"\t0.0000000\n").map_err(BedError::Io)?;
        }
        Ok(())
    }

    /// Mean coverage using reusable event buffer.
    #[inline]
    fn write_mean_coverage<W: Write>(
        &self,
        output: &mut W,
        original_line: &str,
        a_start: u64,
        a_end: u64,
        a_len: u64,
        active: &[ActiveInterval],
        events: &mut Vec<(u64, i32)>,
    ) -> Result<(), BedError> {
        let total_depth = Self::compute_total_depth(active, a_start, a_end, events);
        // Use f32 to match bedtools precision (bedtools uses float internally)
        let mean: f32 = if a_len > 0 {
            total_depth as f32 / a_len as f32
        } else {
            0.0
        };

        output
            .write_all(original_line.as_bytes())
            .map_err(BedError::Io)?;

        // Format mean with {:.7} to match bedtools (uses C printf rounding)
        writeln!(output, "\t{:.7}", mean).map_err(BedError::Io)?;

        Ok(())
    }

    /// Compute total depth using reusable event buffer.
    #[inline]
    fn compute_total_depth(
        active: &[ActiveInterval],
        a_start: u64,
        a_end: u64,
        events: &mut Vec<(u64, i32)>,
    ) -> u64 {
        events.clear();

        for b in active {
            let b_start = b.start as u64;
            let b_end = b.end as u64;

            if b_end > a_start && b_start < a_end {
                let clip_start = b_start.max(a_start);
                let clip_end = b_end.min(a_end);
                events.push((clip_start, 1));
                events.push((clip_end, -1));
            }
        }

        if events.is_empty() {
            return 0;
        }

        events.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        let mut depth: i32 = 0;
        let mut prev_pos = a_start;
        let mut total_depth: u64 = 0;

        for &(pos, delta) in events.iter() {
            if pos > prev_pos && depth > 0 {
                total_depth += (pos - prev_pos) * (depth as u64);
            }
            depth += delta;
            prev_pos = pos;
        }

        total_depth
    }

    /// Histogram coverage using reusable event buffer.
    fn write_histogram_coverage<W: Write>(
        &self,
        output: &mut W,
        original_line: &str,
        a_start: u64,
        a_end: u64,
        a_len: u64,
        active: &[ActiveInterval],
        events: &mut Vec<(u64, i32)>,
    ) -> Result<(), BedError> {
        use std::collections::BTreeMap;

        events.clear();
        events.push((a_start, 0));
        events.push((a_end, 0));

        for b in active {
            let b_start = b.start as u64;
            let b_end = b.end as u64;

            if b_end > a_start && b_start < a_end {
                let clip_start = b_start.max(a_start);
                let clip_end = b_end.min(a_end);
                events.push((clip_start, 1));
                events.push((clip_end, -1));
            }
        }

        events.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        let mut histogram: BTreeMap<u32, u64> = BTreeMap::new();
        let mut depth: i32 = 0;
        let mut prev_pos = a_start;

        for &(pos, delta) in events.iter() {
            if pos > prev_pos && pos <= a_end && prev_pos >= a_start {
                let span = pos - prev_pos;
                *histogram.entry(depth as u32).or_insert(0) += span;
            }
            depth += delta;
            prev_pos = pos;
        }

        for (d, count) in histogram {
            // Use f32 to match bedtools precision (bedtools uses float internally)
            let fraction: f32 = count as f32 / a_len as f32;
            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{:.7}",
                original_line, d, count, a_len, fraction
            )
            .map_err(BedError::Io)?;
        }

        Ok(())
    }

    /// Per-base coverage using reusable event buffer.
    fn write_per_base_coverage<W: Write>(
        &self,
        output: &mut W,
        original_line: &str,
        a_start: u64,
        a_end: u64,
        active: &[ActiveInterval],
        events: &mut Vec<(u64, i32)>,
    ) -> Result<(), BedError> {
        events.clear();

        for b in active {
            let b_start = b.start as u64;
            let b_end = b.end as u64;

            if b_end > a_start && b_start < a_end {
                let clip_start = b_start.max(a_start);
                let clip_end = b_end.min(a_end);
                events.push((clip_start, 1));
                events.push((clip_end, -1));
            }
        }

        events.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        let mut depth: i32 = 0;
        let mut event_idx = 0;

        for pos in a_start..a_end {
            while event_idx < events.len() && events[event_idx].0 <= pos {
                depth += events[event_idx].1;
                event_idx += 1;
            }
            let one_based = pos - a_start + 1;
            writeln!(output, "{}\t{}\t{}", original_line, one_based, depth)
                .map_err(BedError::Io)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_u64_fast() {
        use crate::streaming::parsing::parse_u64_fast;
        assert_eq!(parse_u64_fast(b"12345"), Some(12345));
        assert_eq!(parse_u64_fast(b"0"), Some(0));
        assert_eq!(parse_u64_fast(b""), None);
        assert_eq!(parse_u64_fast(b"abc"), None);
    }

    #[test]
    fn test_compute_coverage_inline() {
        let active = vec![
            ActiveInterval {
                start: 100,
                end: 150,
            },
            ActiveInterval {
                start: 125,
                end: 175,
            },
        ];
        let (num, bases) = StreamingCoverageCommand::compute_coverage_inline(&active, 100, 200);
        assert_eq!(num, 2);
        assert_eq!(bases, 75);
    }

    #[test]
    fn test_compute_coverage_inline_disjoint() {
        let active = vec![
            ActiveInterval {
                start: 100,
                end: 120,
            },
            ActiveInterval {
                start: 150,
                end: 180,
            },
        ];
        let (num, bases) = StreamingCoverageCommand::compute_coverage_inline(&active, 100, 200);
        assert_eq!(num, 2);
        assert_eq!(bases, 50);
    }

    #[test]
    fn test_active_interval_size() {
        assert_eq!(std::mem::size_of::<ActiveInterval>(), 8);
    }

    #[test]
    fn test_streaming_basic_coverage() {
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        let mut a_file = NamedTempFile::new().unwrap();
        let mut b_file = NamedTempFile::new().unwrap();

        writeln!(a_file, "chr1\t100\t200").unwrap();
        writeln!(b_file, "chr1\t100\t150").unwrap();
        writeln!(b_file, "chr1\t125\t175").unwrap();

        a_file.flush().unwrap();
        b_file.flush().unwrap();

        let cmd = StreamingCoverageCommand::new();
        let mut output = Vec::new();

        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("2\t75\t100"));
    }

    #[test]
    fn test_streaming_multiple_chromosomes() {
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        let mut a_file = NamedTempFile::new().unwrap();
        let mut b_file = NamedTempFile::new().unwrap();

        // A has chr1 and chr2
        writeln!(a_file, "chr1\t100\t200").unwrap();
        writeln!(a_file, "chr2\t100\t200").unwrap();

        // B has only chr1
        writeln!(b_file, "chr1\t100\t150").unwrap();

        a_file.flush().unwrap();
        b_file.flush().unwrap();

        let cmd = StreamingCoverageCommand::new();
        let mut output = Vec::new();

        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();
        assert_eq!(lines.len(), 2);
        // chr1 should have coverage
        assert!(lines[0].contains("chr1") && lines[0].contains("1\t50\t100"));
        // chr2 should have zero coverage
        assert!(lines[1].contains("chr2") && lines[1].contains("0\t0\t100"));
    }

    #[test]
    fn test_b_before_a_chromosome() {
        use std::io::Write as IoWrite;
        use tempfile::NamedTempFile;

        let mut a_file = NamedTempFile::new().unwrap();
        let mut b_file = NamedTempFile::new().unwrap();

        // A has only chr2
        writeln!(a_file, "chr2\t100\t200").unwrap();

        // B has chr1 (before chr2 lexicographically)
        writeln!(b_file, "chr1\t100\t150").unwrap();
        writeln!(b_file, "chr2\t100\t150").unwrap();

        a_file.flush().unwrap();
        b_file.flush().unwrap();

        let cmd = StreamingCoverageCommand::new();
        let mut output = Vec::new();

        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        // chr2 should have coverage from B's chr2 interval
        assert!(result.contains("1\t50\t100"));
    }
}
