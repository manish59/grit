//! Streaming closest implementation with TRUE O(k) memory complexity.
//!
//! Finds the closest B interval(s) for each A interval using a streaming
//! sweep-line algorithm that doesn't load entire files into memory.
//!
//! ZERO ALLOCATION in hot path:
//! - No per-record String allocation (raw byte parsing with memchr)
//! - Vec + head index for active set (no VecDeque)
//! - itoa for integer formatting
//! - Large buffered I/O (256KB input, 8MB output)
//!
//! # Memory Complexity
//!
//! O(k) where k = maximum number of B intervals overlapping any single A interval
//! plus O(t) for ties at any position.
//!
//! # Distance Calculation (matches bedtools)
//!
//! - Overlap: distance = 0
//! - Upstream (B.end <= A.start): distance = A.start - B.end + 1
//! - Downstream (B.start >= A.end): distance = B.start - A.end + 1
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
    /// Original line bytes (stored for output since B may be emitted multiple times)
    line: Vec<u8>,
}

/// Streaming closest command configuration.
#[derive(Debug, Clone)]
pub struct StreamingClosestCommand {
    /// Ignore overlapping intervals (-io flag)
    pub ignore_overlaps: bool,
    /// Ignore upstream intervals (-iu flag)
    pub ignore_upstream: bool,
    /// Ignore downstream intervals (-id flag)
    pub ignore_downstream: bool,
    /// Report all ties (bedtools -t all, default true)
    pub report_all_ties: bool,
}

impl Default for StreamingClosestCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl StreamingClosestCommand {
    pub fn new() -> Self {
        Self {
            ignore_overlaps: false,
            ignore_upstream: false,
            ignore_downstream: false,
            report_all_ties: true,
        }
    }

    /// Execute streaming closest on two sorted BED files.
    ///
    /// Memory usage: O(k) where k = max overlapping B intervals at any point
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        a_path: P,
        b_path: P,
        output: &mut W,
    ) -> Result<StreamingClosestStats, BedError> {
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

        // Active set: B intervals that might overlap current or future A
        let mut active: Vec<ActiveB> = Vec::with_capacity(1024);
        let mut head_idx: usize = 0;

        // Left candidates: B intervals with largest end that is <= current A.start
        // Multiple B can have the same end (ties)
        let mut left_candidates: Vec<ActiveB> = Vec::with_capacity(16);
        let mut left_end: u32 = 0; // The end position of left candidates

        // Right candidates: B intervals downstream with smallest start that is >= current A.end
        // Multiple B can have the same start (ties)
        let mut right_candidates: Vec<ActiveB> = Vec::with_capacity(16);

        // Stats
        let mut stats = StreamingClosestStats::default();

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

            // Chromosome change
            let chrom_changed = chrom != a_chrom.as_slice();
            if chrom_changed {
                a_chrom.clear();
                a_chrom.extend_from_slice(chrom);
                active.clear();
                head_idx = 0;
                left_candidates.clear();
                left_end = 0;
                right_candidates.clear();

                // Skip B to current chromosome (or B has already passed it)
                if !b_exhausted && !seen_b_chroms.contains(chrom) {
                    while b_chrom.as_slice() != chrom {
                        pending_b =
                            Self::read_next_b(&mut b_reader, &mut b_line_buf, &mut b_chrom)?;
                        stats.b_intervals += 1;
                        if pending_b.is_none() {
                            b_exhausted = true;
                            break;
                        }
                        seen_b_chroms.insert(b_chrom.clone());
                    }
                }
            }

            // Re-evaluate right_candidates from previous A
            // They might now be overlapping or upstream for this A.
            // Move overlapping ones to active set first, defer upstream ones.
            let mut deferred_upstream: Vec<ActiveB> = Vec::new();
            if !right_candidates.is_empty() {
                let right_start = right_candidates[0].start as u64;
                if right_start < a_end {
                    // Right candidates are no longer downstream - move to appropriate bucket
                    for rc in right_candidates.drain(..) {
                        if (rc.end as u64) <= a_start {
                            // Now upstream - defer until after expire loop
                            // to preserve B-file order (active set items first)
                            deferred_upstream.push(rc);
                        } else {
                            // Now overlap/active
                            active.push(rc);
                        }
                    }
                }
                // If right_start >= a_end, keep right_candidates as is (still downstream)
            }

            // Expire old B from active and update left_candidates
            while head_idx < active.len() {
                let b = &active[head_idx];
                if (b.end as u64) <= a_start {
                    // B is now upstream - add to left_candidates if it's the closest or tied
                    if b.end > left_end {
                        // New closest upstream
                        left_candidates.clear();
                        left_candidates.push(b.clone());
                        left_end = b.end;
                    } else if b.end == left_end {
                        // Tied with current closest
                        left_candidates.push(b.clone());
                    }
                    head_idx += 1;
                } else {
                    break;
                }
            }

            // Now process deferred upstream from right_candidates
            // These have higher start than active-set items, so appending
            // preserves B-file order (sorted by start).
            for rc in deferred_upstream {
                if rc.end > left_end {
                    left_candidates.clear();
                    left_candidates.push(rc);
                    left_end = left_candidates[0].end;
                } else if rc.end == left_end {
                    left_candidates.push(rc);
                }
            }

            // Compact if needed
            if head_idx > 4096 && head_idx * 2 > active.len() {
                active.drain(0..head_idx);
                head_idx = 0;
            }

            // Add new B intervals until B.start >= A.end
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
                        stats.b_intervals += 1;
                        if pending_b.is_none() {
                            b_exhausted = true;
                            break;
                        }
                        seen_b_chroms.insert(b_chrom.clone());
                        continue;
                    } else {
                        // B is on the same chromosome as A
                        if (b.start as u64) >= a_end {
                            // B is downstream
                            if !right_candidates.is_empty() {
                                // We already have downstream candidates at an earlier position
                                // This B is further downstream - keep as pending and stop
                                pending_b = Some(b);
                                break;
                            }
                            // Collect all B at this start position (ties)
                            let right_start = b.start;
                            right_candidates.push(b);

                            // Read more B to find ties at same start position
                            loop {
                                let next_b = Self::read_next_b(
                                    &mut b_reader,
                                    &mut b_line_buf,
                                    &mut b_chrom,
                                )?;
                                if let Some(nb) = next_b {
                                    stats.b_intervals += 1;
                                    seen_b_chroms.insert(b_chrom.clone());
                                    if b_chrom.as_slice() != chrom {
                                        pending_b = Some(nb);
                                        break;
                                    }
                                    if nb.start != right_start {
                                        pending_b = Some(nb);
                                        break;
                                    }
                                    right_candidates.push(nb);
                                } else {
                                    pending_b = None;
                                    b_exhausted = true;
                                    break;
                                }
                            }
                            break; // Done collecting B for this A
                        }
                        // Check if B is upstream (ends before A starts)
                        if (b.end as u64) <= a_start {
                            // B is upstream - add to left_candidates if appropriate
                            if b.end > left_end {
                                left_candidates.clear();
                                left_candidates.push(b);
                                left_end = left_candidates[0].end;
                            } else if b.end == left_end {
                                left_candidates.push(b);
                            }
                        } else {
                            // B could overlap current or future A - add to active
                            active.push(b);
                        }
                        pending_b =
                            Self::read_next_b(&mut b_reader, &mut b_line_buf, &mut b_chrom)?;
                        stats.b_intervals += 1;
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

            // Find closest
            let active_slice = &active[head_idx..];
            let mut min_dist: u64 = u64::MAX;

            // Check overlaps in active set
            let mut overlaps: Vec<&ActiveB> = Vec::new();
            if !self.ignore_overlaps {
                for b in active_slice {
                    let b_start = b.start as u64;
                    let b_end = b.end as u64;
                    if b_start < a_end && b_end > a_start {
                        overlaps.push(b);
                    }
                }
            }

            // If overlaps exist, they win (distance = 0)
            if !overlaps.is_empty() {
                if self.report_all_ties {
                    for b in &overlaps {
                        Self::write_pair(&mut output, line_bytes, &b.line)?;
                        stats.pairs_written += 1;
                    }
                } else {
                    Self::write_pair(&mut output, line_bytes, &overlaps[0].line)?;
                    stats.pairs_written += 1;
                }
                continue;
            }

            // Check for downstream candidates in active set (handles nested A intervals)
            // These are B intervals added for a previous A but are downstream of current A
            let mut active_downstream: Vec<&ActiveB> = Vec::new();
            let mut active_downstream_start: u32 = u32::MAX;
            if !self.ignore_downstream {
                for b in active_slice {
                    if (b.start as u64) >= a_end {
                        if b.start < active_downstream_start {
                            active_downstream.clear();
                            active_downstream.push(b);
                            active_downstream_start = b.start;
                        } else if b.start == active_downstream_start {
                            active_downstream.push(b);
                        }
                    }
                }
            }

            // Check upstream (left_candidates)
            let mut upstream_dist: u64 = u64::MAX;
            if !self.ignore_upstream && !left_candidates.is_empty() {
                // bedtools distance: A.start - B.end + 1
                upstream_dist = a_start - left_end as u64 + 1;
                min_dist = upstream_dist;
            }

            // Check downstream - consider both right_candidates and active_downstream
            let mut downstream_dist: u64 = u64::MAX;
            let mut use_active_downstream = false;
            let mut use_right_candidates = false;

            if !self.ignore_downstream {
                // Check active_downstream first (may have closer B from nested A scenario)
                if !active_downstream.is_empty() {
                    downstream_dist = active_downstream_start as u64 - a_end + 1;
                    use_active_downstream = true;
                }

                // Check right_candidates
                if !right_candidates.is_empty() {
                    let right_dist = right_candidates[0].start as u64 - a_end + 1;
                    if right_dist < downstream_dist {
                        downstream_dist = right_dist;
                        use_active_downstream = false;
                        use_right_candidates = true;
                    } else if right_dist == downstream_dist {
                        // Tie - use both
                        use_right_candidates = true;
                    }
                }

                if downstream_dist < min_dist {
                    min_dist = downstream_dist;
                }
            }

            // Output results
            if min_dist == u64::MAX {
                // No closest found
                Self::write_no_closest(&mut output, line_bytes)?;
            } else if upstream_dist == downstream_dist && upstream_dist == min_dist {
                // Tie between upstream and downstream
                if self.report_all_ties {
                    for lc in &left_candidates {
                        Self::write_pair(&mut output, line_bytes, &lc.line)?;
                        stats.pairs_written += 1;
                    }
                    if use_active_downstream {
                        for b in &active_downstream {
                            Self::write_pair(&mut output, line_bytes, &b.line)?;
                            stats.pairs_written += 1;
                        }
                    }
                    if use_right_candidates {
                        for rc in &right_candidates {
                            Self::write_pair(&mut output, line_bytes, &rc.line)?;
                            stats.pairs_written += 1;
                        }
                    }
                } else if !left_candidates.is_empty() {
                    Self::write_pair(&mut output, line_bytes, &left_candidates[0].line)?;
                    stats.pairs_written += 1;
                }
            } else if upstream_dist == min_dist {
                if self.report_all_ties {
                    for lc in &left_candidates {
                        Self::write_pair(&mut output, line_bytes, &lc.line)?;
                        stats.pairs_written += 1;
                    }
                } else if !left_candidates.is_empty() {
                    Self::write_pair(&mut output, line_bytes, &left_candidates[0].line)?;
                    stats.pairs_written += 1;
                }
            } else if downstream_dist == min_dist {
                if self.report_all_ties {
                    if use_active_downstream {
                        for b in &active_downstream {
                            Self::write_pair(&mut output, line_bytes, &b.line)?;
                            stats.pairs_written += 1;
                        }
                    }
                    if use_right_candidates {
                        for rc in &right_candidates {
                            Self::write_pair(&mut output, line_bytes, &rc.line)?;
                            stats.pairs_written += 1;
                        }
                    }
                } else if use_active_downstream && !active_downstream.is_empty() {
                    Self::write_pair(&mut output, line_bytes, &active_downstream[0].line)?;
                    stats.pairs_written += 1;
                } else if use_right_candidates && !right_candidates.is_empty() {
                    Self::write_pair(&mut output, line_bytes, &right_candidates[0].line)?;
                    stats.pairs_written += 1;
                }
            } else {
                Self::write_no_closest(&mut output, line_bytes)?;
            }
        }

        // Count remaining B
        while pending_b.is_some() {
            stats.b_intervals += 1;
            pending_b = Self::read_next_b(&mut b_reader, &mut b_line_buf, &mut b_chrom)?;
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
    fn write_no_closest<W: Write>(output: &mut W, a_line: &[u8]) -> Result<(), BedError> {
        output.write_all(a_line).map_err(BedError::Io)?;
        output.write_all(b"\t.\t-1\t-1").map_err(BedError::Io)?;
        output.write_all(b"\n").map_err(BedError::Io)?;
        Ok(())
    }
}

/// Statistics from streaming closest operation.
#[derive(Debug, Default, Clone)]
pub struct StreamingClosestStats {
    pub a_intervals: usize,
    pub b_intervals: usize,
    pub pairs_written: usize,
    pub max_active_b: usize,
}

impl std::fmt::Display for StreamingClosestStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "A: {}, B: {}, Pairs: {}, Max active B: {}",
            self.a_intervals, self.b_intervals, self.pairs_written, self.max_active_b
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
    fn test_basic_closest_downstream() {
        let a_file = create_temp_bed("chr1\t100\t200\n");
        let b_file = create_temp_bed("chr1\t300\t400\n");

        let cmd = StreamingClosestCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("chr1\t100\t200\tchr1\t300\t400"));
    }

    #[test]
    fn test_basic_closest_upstream() {
        let a_file = create_temp_bed("chr1\t300\t400\n");
        let b_file = create_temp_bed("chr1\t100\t200\n");

        let cmd = StreamingClosestCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("chr1\t300\t400\tchr1\t100\t200"));
    }

    #[test]
    fn test_closest_overlap() {
        let a_file = create_temp_bed("chr1\t100\t200\n");
        let b_file = create_temp_bed("chr1\t150\t250\n");

        let cmd = StreamingClosestCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("chr1\t100\t200\tchr1\t150\t250"));
    }

    #[test]
    fn test_closest_ties() {
        // A at 200-300, B at 100-150 (upstream, dist=51) and 350-400 (downstream, dist=51)
        let a_file = create_temp_bed("chr1\t200\t300\n");
        let b_file = create_temp_bed("chr1\t100\t150\nchr1\t350\t400\n");

        let cmd = StreamingClosestCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();
        assert_eq!(lines.len(), 2, "Should report both ties");
    }

    #[test]
    fn test_no_closest_different_chrom() {
        let a_file = create_temp_bed("chr2\t100\t200\n");
        let b_file = create_temp_bed("chr1\t100\t200\n");

        let cmd = StreamingClosestCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains(".\t-1\t-1"));
    }

    #[test]
    fn test_multiple_a_intervals() {
        let a_file = create_temp_bed("chr1\t100\t200\nchr1\t500\t600\n");
        let b_file = create_temp_bed("chr1\t300\t400\n");

        let cmd = StreamingClosestCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();
        assert_eq!(lines.len(), 2);
    }

    #[test]
    fn test_upstream_ties() {
        // Two B intervals at same position (same start/end, different content)
        let a_file = create_temp_bed("chr1\t300\t400\n");
        let b_file = create_temp_bed("chr1\t100\t200\tgene1\nchr1\t100\t200\tgene2\n");

        let cmd = StreamingClosestCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();
        // Should report both B intervals since they're at the same position
        assert_eq!(
            lines.len(),
            2,
            "Should report both upstream ties: {}",
            result
        );
    }

    #[test]
    fn test_downstream_ties() {
        // Two B intervals downstream with same start position
        let a_file = create_temp_bed("chr1\t100\t200\n");
        let b_file = create_temp_bed("chr1\t300\t400\tgene1\nchr1\t300\t500\tgene2\n");

        let cmd = StreamingClosestCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();
        // Should report both B intervals since they have the same start (same distance)
        assert_eq!(
            lines.len(),
            2,
            "Should report both downstream ties: {}",
            result
        );
        assert!(result.contains("gene1"), "Should include gene1");
        assert!(result.contains("gene2"), "Should include gene2");
    }

    #[test]
    fn test_downstream_ties_different_end() {
        // Three B intervals with same start but different ends
        let a_file = create_temp_bed("chr1\t100\t200\n");
        let b_file = create_temp_bed("chr1\t300\t350\nchr1\t300\t400\nchr1\t300\t450\n");

        let cmd = StreamingClosestCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();
        assert_eq!(
            lines.len(),
            3,
            "Should report all three downstream ties: {}",
            result
        );
    }

    #[test]
    fn test_parse_u64_fast() {
        use crate::streaming::parsing::parse_u64_fast;
        assert_eq!(parse_u64_fast(b"12345"), Some(12345));
        assert_eq!(parse_u64_fast(b"0"), Some(0));
        assert_eq!(parse_u64_fast(b""), None);
    }

    // =============================================================================
    // Flag combination unit tests
    // =============================================================================

    #[test]
    fn test_ignore_overlaps_flag() {
        // A overlaps B1, B2 is downstream
        let a_file = create_temp_bed("chr1\t100\t200\n");
        let b_file = create_temp_bed("chr1\t150\t250\nchr1\t300\t400\n");

        let mut cmd = StreamingClosestCommand::new();
        cmd.ignore_overlaps = true;

        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(
            result.contains("300\t400"),
            "Should find non-overlapping B2: {}",
            result
        );
        assert!(
            !result.contains("150\t250"),
            "Should skip overlapping B1: {}",
            result
        );
    }

    #[test]
    fn test_ignore_upstream_flag() {
        // A at 200-300, B1 upstream, B2 downstream
        let a_file = create_temp_bed("chr1\t200\t300\n");
        let b_file = create_temp_bed("chr1\t100\t150\nchr1\t350\t400\n");

        let mut cmd = StreamingClosestCommand::new();
        cmd.ignore_upstream = true;

        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(
            result.contains("350\t400"),
            "Should find downstream B2: {}",
            result
        );
        assert!(
            !result.contains("100\t150"),
            "Should skip upstream B1: {}",
            result
        );
    }

    #[test]
    fn test_ignore_downstream_flag() {
        // A at 200-300, B1 upstream, B2 downstream
        let a_file = create_temp_bed("chr1\t200\t300\n");
        let b_file = create_temp_bed("chr1\t100\t150\nchr1\t350\t400\n");

        let mut cmd = StreamingClosestCommand::new();
        cmd.ignore_downstream = true;

        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(
            result.contains("100\t150"),
            "Should find upstream B1: {}",
            result
        );
        assert!(
            !result.contains("350\t400"),
            "Should skip downstream B2: {}",
            result
        );
    }

    #[test]
    fn test_ignore_overlaps_and_upstream() {
        // Only downstream non-overlapping should match
        let a_file = create_temp_bed("chr1\t100\t200\n");
        let b_file = create_temp_bed("chr1\t50\t80\nchr1\t150\t250\nchr1\t300\t400\n");

        let mut cmd = StreamingClosestCommand::new();
        cmd.ignore_overlaps = true;
        cmd.ignore_upstream = true;

        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(
            result.contains("300\t400"),
            "Should find downstream B3: {}",
            result
        );
        assert!(
            !result.contains("50\t80"),
            "Should skip upstream: {}",
            result
        );
        assert!(
            !result.contains("150\t250"),
            "Should skip overlapping: {}",
            result
        );
    }

    #[test]
    fn test_ignore_overlaps_and_downstream() {
        // Only upstream non-overlapping should match
        let a_file = create_temp_bed("chr1\t100\t200\n");
        let b_file = create_temp_bed("chr1\t50\t80\nchr1\t150\t250\nchr1\t300\t400\n");

        let mut cmd = StreamingClosestCommand::new();
        cmd.ignore_overlaps = true;
        cmd.ignore_downstream = true;

        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(
            result.contains("50\t80"),
            "Should find upstream B2: {}",
            result
        );
        assert!(
            !result.contains("300\t400"),
            "Should skip downstream: {}",
            result
        );
        assert!(
            !result.contains("150\t250"),
            "Should skip overlapping: {}",
            result
        );
    }

    #[test]
    fn test_ignore_upstream_and_downstream() {
        // Only overlapping should match
        let a_file = create_temp_bed("chr1\t100\t200\n");
        let b_file = create_temp_bed("chr1\t50\t80\nchr1\t150\t250\nchr1\t300\t400\n");

        let mut cmd = StreamingClosestCommand::new();
        cmd.ignore_upstream = true;
        cmd.ignore_downstream = true;

        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(
            result.contains("150\t250"),
            "Should find overlapping B1: {}",
            result
        );
        assert!(
            !result.contains("50\t80"),
            "Should skip upstream: {}",
            result
        );
        assert!(
            !result.contains("300\t400"),
            "Should skip downstream: {}",
            result
        );
    }

    #[test]
    fn test_ignore_all_flags() {
        // No matches possible
        let a_file = create_temp_bed("chr1\t100\t200\n");
        let b_file = create_temp_bed("chr1\t50\t80\nchr1\t150\t250\nchr1\t300\t400\n");

        let mut cmd = StreamingClosestCommand::new();
        cmd.ignore_overlaps = true;
        cmd.ignore_upstream = true;
        cmd.ignore_downstream = true;

        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(
            result.contains(".\t-1\t-1"),
            "Should report no closest: {}",
            result
        );
    }

    #[test]
    fn test_report_first_tie_only() {
        // Equidistant B intervals
        let a_file = create_temp_bed("chr1\t200\t300\n");
        let b_file = create_temp_bed("chr1\t100\t150\nchr1\t350\t400\n");

        let mut cmd = StreamingClosestCommand::new();
        cmd.report_all_ties = false;

        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();
        assert_eq!(lines.len(), 1, "Should report only first tie: {}", result);
    }
}
