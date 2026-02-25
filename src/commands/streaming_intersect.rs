//! Streaming intersect implementation with O(k) memory complexity.
//!
//! This module implements a true streaming sweep-line intersection algorithm
//! that processes sorted BED files without loading them entirely into memory.
//!
//! # Algorithm
//!
//! For sorted inputs A and B:
//! 1. Read intervals from A one at a time
//! 2. Maintain a sliding window of "active" B intervals that could overlap current A
//! 3. Advance B reader, adding intervals to active set when B.start < A.end
//! 4. Remove from active set when B.end <= A.start (can't overlap anymore)
//! 5. Report overlaps between current A and active B intervals
//!
//! # Memory Complexity
//!
//! - O(k) where k = maximum number of B intervals overlapping any single A interval
//! - For typical genomic data, k << n, often k < 100
//! - Worst case (all intervals overlap): O(m) where m = total B intervals
//!
//! # Worst Case Behavior
//!
//! The pathological case occurs when ALL B intervals overlap a single A interval.
//! Example: A = [(chr1, 0, 1_000_000_000)] and B = [(chr1, i, i+100) for i in 0..10_000_000]
//! In this case, active_b grows to contain all 10M B intervals.
//!
//! To detect this, the implementation tracks `max_active_b` and emits a warning
//! to stderr if it exceeds `ACTIVE_WINDOW_WARNING_THRESHOLD`.
//!
//! # Requirements
//!
//! - Both input files MUST be sorted by chromosome, then by start position
//! - Use `grit sort` or `sort -k1,1 -k2,2n` to pre-sort if needed
//! - The implementation validates sorted order inline and fails fast on violation
//!
//! # bedtools CLI Compatibility
//!
//! Output modes exactly match bedtools intersect:
//!
//! | Flags     | Output                                    |
//! |-----------|-------------------------------------------|
//! | (default) | Overlap region (intersection of A and B)  |
//! | -wa       | A record (once per overlapping B)         |
//! | -wb       | Overlap region + B record                 |
//! | -wa -wb   | A record + B record (tab-separated)       |
//! | -c        | A record + overlap count                  |
//! | -u        | A record (once if ANY overlap)            |
//! | -v        | A record (only if NO overlaps)            |

use crate::bed::{BedError, BedReader};
use crate::interval::BedRecord;
use crate::streaming::parsing::{parse_bed3_bytes, parse_bed3_bytes_with_rest, should_skip_line};
use std::collections::{HashSet, VecDeque};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Warning threshold for active window size (potential pathological case)
const ACTIVE_WINDOW_WARNING_THRESHOLD: usize = 100_000;

/// Compaction threshold for active set - trigger when head_idx exceeds this.
const COMPACTION_THRESHOLD: usize = 4096;

/// Active B interval - stores coordinates and original line for output.
/// Coordinates use u32 (4GB max position) for memory efficiency.
#[derive(Debug, Clone)]
struct ActiveB {
    start: u32,
    end: u32,
    /// Original line bytes (stored for output)
    line: Vec<u8>,
}

/// Output mode computed once before processing to reduce branch entropy.
/// This replaces repeated flag checks in the hot loop.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputMode {
    /// Default: print overlap region (intersection of A and B)
    Default,
    /// -wa: print A record once per overlap
    WriteA,
    /// -wb: print overlap region + B record
    WriteB,
    /// -wa -wb: print A + B records
    WriteBoth,
    /// -u: print A once if any overlap
    Unique,
    /// -c: print A + overlap count
    Count,
    /// -v: print A only if no overlaps
    NoOverlap,
}

/// Streaming intersect command configuration.
#[derive(Debug, Clone)]
pub struct StreamingIntersectCommand {
    /// Write original A entry (-wa)
    pub write_a: bool,
    /// Write original B entry (-wb)
    pub write_b: bool,
    /// Only report unique A intervals (first overlap only) (-u)
    pub unique: bool,
    /// Only report A intervals with NO overlap (-v)
    pub no_overlap: bool,
    /// Minimum overlap fraction for A (-f)
    pub fraction_a: Option<f64>,
    /// Minimum overlap fraction for B (-F)
    pub fraction_b: Option<f64>,
    /// Require reciprocal fraction overlap (-r)
    pub reciprocal: bool,
    /// Report the number of overlaps (-c)
    pub count: bool,
    /// Require same strand (-s)
    pub same_strand: bool,
    /// Require opposite strand (-S)
    pub opposite_strand: bool,
    /// Skip sorted validation (use --assume-sorted)
    pub assume_sorted: bool,
    /// Warn if active window exceeds threshold
    pub warn_large_window: bool,
}

impl Default for StreamingIntersectCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl StreamingIntersectCommand {
    pub fn new() -> Self {
        Self {
            write_a: false,
            write_b: false,
            unique: false,
            no_overlap: false,
            fraction_a: None,
            fraction_b: None,
            reciprocal: false,
            count: false,
            same_strand: false,
            opposite_strand: false,
            assume_sorted: false,
            warn_large_window: true,
        }
    }

    /// Compute output mode once before processing.
    /// This eliminates repeated flag checks in the hot loop.
    #[inline]
    fn compute_output_mode(&self) -> OutputMode {
        if self.no_overlap {
            OutputMode::NoOverlap
        } else if self.count {
            OutputMode::Count
        } else if self.unique {
            OutputMode::Unique
        } else if self.write_a && self.write_b {
            OutputMode::WriteBoth
        } else if self.write_b {
            OutputMode::WriteB
        } else if self.write_a {
            OutputMode::WriteA
        } else {
            // Default bedtools behavior: print A record
            OutputMode::Default
        }
    }

    /// Check if any filters are active (fraction or strand)
    #[inline]
    fn has_filters(&self) -> bool {
        self.fraction_a.is_some()
            || self.fraction_b.is_some()
            || self.reciprocal
            || self.same_strand
            || self.opposite_strand
    }

    /// Execute streaming intersect on two sorted BED files.
    ///
    /// Memory usage: O(k) where k = max overlapping B intervals at any point
    ///
    /// This uses the optimized zero-allocation parsing path for maximum performance,
    /// falling back to the original path when strand filters are used.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        a_path: P,
        b_path: P,
        output: &mut W,
    ) -> Result<StreamingStats, BedError> {
        // Fall back to original path for strand filtering (not supported in optimized path)
        if self.same_strand || self.opposite_strand {
            let a_file = File::open(a_path.as_ref())?;
            let b_file = File::open(b_path.as_ref())?;
            let a_reader = BedReader::new(BufReader::with_capacity(256 * 1024, a_file));
            let b_reader = BedReader::new(BufReader::with_capacity(256 * 1024, b_file));
            return self.run_streaming(a_reader, b_reader, output);
        }

        // Use optimized path with raw line parsing
        self.run_optimized(a_path, b_path, output)
    }

    /// Optimized streaming intersect with zero-allocation parsing.
    ///
    /// Uses:
    /// - Raw line parsing with memchr (no String allocation per record)
    /// - Vec + head_idx instead of VecDeque (better cache locality)
    /// - Stores raw line bytes for output (avoids formatting overhead)
    fn run_optimized<P: AsRef<Path>, W: Write>(
        &self,
        a_path: P,
        b_path: P,
        output: &mut W,
    ) -> Result<StreamingStats, BedError> {
        let mut stats = StreamingStats::default();

        // Large output buffer (8MB)
        let mut writer = BufWriter::with_capacity(8 * 1024 * 1024, output);

        // Stream A file with large buffer
        let a_file = File::open(a_path.as_ref())?;
        let mut a_reader = BufReader::with_capacity(256 * 1024, a_file);

        // Stream B file with large buffer
        let b_file = File::open(b_path.as_ref())?;
        let mut b_reader = BufReader::with_capacity(256 * 1024, b_file);

        // Reusable line buffers
        let mut a_line_buf = String::with_capacity(1024);
        let mut b_line_buf = String::with_capacity(1024);

        // Current A chromosome (reused buffer)
        let mut a_chrom: Vec<u8> = Vec::with_capacity(64);

        // Pending B: chrom stored separately
        let mut b_chrom: Vec<u8> = Vec::with_capacity(64);
        let mut pending_b =
            Self::read_next_b_optimized(&mut b_reader, &mut b_line_buf, &mut b_chrom)?;
        let mut b_exhausted = pending_b.is_none();

        // Track seen chromosomes for sort validation
        let mut seen_a_chroms: HashSet<Vec<u8>> = HashSet::new();
        let mut seen_b_chroms: HashSet<Vec<u8>> = HashSet::new();
        if !b_exhausted {
            seen_b_chroms.insert(b_chrom.clone());
        }

        // Active set: Vec with head index (better cache locality than VecDeque)
        let mut active: Vec<ActiveB> = Vec::with_capacity(1024);
        let mut head_idx: usize = 0;

        // Sorted validation state
        let mut prev_a_start: u64 = 0;
        let mut prev_b_start: u64 = 0;
        let mut warned_large_window = false;

        // Compute output mode once
        let output_mode = self.compute_output_mode();
        let has_filters = self.has_filters();

        // itoa buffer for fast integer formatting
        let mut itoa_buf = itoa::Buffer::new();

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
            let (chrom, a_start, a_end, rest_start) = match parse_bed3_bytes_with_rest(line_bytes) {
                Some(v) => v,
                None => continue,
            };

            stats.a_intervals += 1;

            // Sorted validation for A
            if !self.assume_sorted {
                let chrom_changed = chrom != a_chrom.as_slice();
                if chrom_changed {
                    if seen_a_chroms.contains(chrom) {
                        return Err(BedError::InvalidFormat(format!(
                            "File A not sorted: chromosome '{}' at record {} seen before",
                            String::from_utf8_lossy(chrom),
                            stats.a_intervals
                        )));
                    }
                } else if a_start < prev_a_start {
                    return Err(BedError::InvalidFormat(format!(
                        "File A not sorted: position {} at record {} comes after {} on {}",
                        a_start,
                        stats.a_intervals,
                        prev_a_start,
                        String::from_utf8_lossy(chrom)
                    )));
                }
                prev_a_start = a_start;
            }

            // Check chromosome change
            let chrom_changed = chrom != a_chrom.as_slice();
            if chrom_changed {
                // Update A chromosome (reuses buffer)
                seen_a_chroms.insert(a_chrom.clone());
                a_chrom.clear();
                a_chrom.extend_from_slice(chrom);

                // Clear active set
                active.clear();
                head_idx = 0;
                prev_b_start = 0;

                // Skip B records until we reach this chromosome (or B has already passed it)
                if !b_exhausted && !seen_b_chroms.contains(chrom) {
                    while b_chrom.as_slice() != chrom {
                        pending_b = Self::read_next_b_optimized(
                            &mut b_reader,
                            &mut b_line_buf,
                            &mut b_chrom,
                        )?;
                        stats.b_intervals += 1;
                        if pending_b.is_none() {
                            b_exhausted = true;
                            break;
                        }
                        seen_b_chroms.insert(b_chrom.clone());
                    }
                }
            }

            // Step 1: Remove expired B intervals (head index advancement)
            while head_idx < active.len() && (active[head_idx].end as u64) <= a_start {
                head_idx += 1;
            }

            // Periodic compaction to prevent memory growth
            if head_idx > COMPACTION_THRESHOLD && head_idx * 2 > active.len() {
                active.drain(0..head_idx);
                head_idx = 0;
            }

            // Step 2: Add new B intervals to active set
            if !b_exhausted {
                while let Some(b) = pending_b.take() {
                    // Sorted validation for B
                    if !self.assume_sorted && b_chrom.as_slice() == chrom {
                        if (b.start as u64) < prev_b_start {
                            return Err(BedError::InvalidFormat(format!(
                                "File B not sorted: position {} comes after {} on {}",
                                b.start,
                                prev_b_start,
                                String::from_utf8_lossy(&b_chrom)
                            )));
                        }
                        prev_b_start = b.start as u64;
                    }

                    if b_chrom.as_slice() != chrom {
                        // B is on a different chromosome
                        if seen_b_chroms.contains(chrom) {
                            // We've already seen A's chromosome in B, so B has moved past it
                            pending_b = Some(b);
                            break;
                        }
                        // B hasn't reached A's chromosome yet, skip it
                        stats.b_intervals += 1;
                        pending_b = Self::read_next_b_optimized(
                            &mut b_reader,
                            &mut b_line_buf,
                            &mut b_chrom,
                        )?;
                        if pending_b.is_none() {
                            b_exhausted = true;
                            break;
                        }
                        // Check for sort violation
                        if !self.assume_sorted && seen_b_chroms.contains(&b_chrom) {
                            return Err(BedError::InvalidFormat(format!(
                                "File B not sorted: chromosome '{}' seen before",
                                String::from_utf8_lossy(&b_chrom)
                            )));
                        }
                        seen_b_chroms.insert(b_chrom.clone());
                        continue;
                    }

                    // B is on the same chromosome as A
                    if (b.start as u64) >= a_end {
                        // B starts after A ends, put back for later
                        pending_b = Some(b);
                        break;
                    }

                    // Only add if B could overlap current A (B.end > A.start)
                    if (b.end as u64) > a_start {
                        active.push(b);
                    }

                    // Read next B
                    stats.b_intervals += 1;
                    pending_b =
                        Self::read_next_b_optimized(&mut b_reader, &mut b_line_buf, &mut b_chrom)?;
                    if pending_b.is_none() {
                        b_exhausted = true;
                        break;
                    }
                    if !self.assume_sorted
                        && b_chrom.as_slice() != chrom
                        && seen_b_chroms.contains(&b_chrom)
                    {
                        return Err(BedError::InvalidFormat(format!(
                            "File B not sorted: chromosome '{}' seen before",
                            String::from_utf8_lossy(&b_chrom)
                        )));
                    }
                    seen_b_chroms.insert(b_chrom.clone());
                }
            }

            // Track max active set size
            let active_size = active.len() - head_idx;
            stats.max_active_b = stats.max_active_b.max(active_size);

            // Warn on pathological case (only once)
            if self.warn_large_window
                && !warned_large_window
                && active_size > ACTIVE_WINDOW_WARNING_THRESHOLD
            {
                eprintln!(
                    "Warning: Large active window detected ({} intervals). Memory usage: O({})",
                    active_size, active_size
                );
                warned_large_window = true;
            }

            // Step 3: Process overlaps based on output mode
            let active_slice = &active[head_idx..];

            match output_mode {
                OutputMode::NoOverlap => {
                    // -v mode: output A if no overlaps found
                    let has_overlap = active_slice.iter().any(|b| {
                        let b_start = b.start as u64;
                        let b_end = b.end as u64;
                        b_end > a_start
                            && b_start < a_end
                            && (!has_filters
                                || self.passes_filters_raw(a_start, a_end, b_start, b_end))
                    });

                    if !has_overlap {
                        writer.write_all(line_bytes)?;
                        writer.write_all(b"\n")?;
                    }
                }

                OutputMode::Count => {
                    // -c mode: output A with overlap count
                    let count = active_slice
                        .iter()
                        .filter(|b| {
                            let b_start = b.start as u64;
                            let b_end = b.end as u64;
                            b_end > a_start
                                && b_start < a_end
                                && (!has_filters
                                    || self.passes_filters_raw(a_start, a_end, b_start, b_end))
                        })
                        .count();

                    writer.write_all(line_bytes)?;
                    writer.write_all(b"\t")?;
                    writer.write_all(itoa_buf.format(count).as_bytes())?;
                    writer.write_all(b"\n")?;
                }

                OutputMode::Unique => {
                    // -u mode: output A once if any overlap exists
                    let has_overlap = active_slice.iter().any(|b| {
                        let b_start = b.start as u64;
                        let b_end = b.end as u64;
                        b_end > a_start
                            && b_start < a_end
                            && (!has_filters
                                || self.passes_filters_raw(a_start, a_end, b_start, b_end))
                    });

                    if has_overlap {
                        writer.write_all(line_bytes)?;
                        writer.write_all(b"\n")?;
                        stats.overlaps_found += 1;
                    }
                }

                OutputMode::Default => {
                    // Default: output overlap region with A's extra fields
                    for b in active_slice {
                        let b_start = b.start as u64;
                        let b_end = b.end as u64;

                        if b_end <= a_start || b_start >= a_end {
                            continue;
                        }

                        if has_filters && !self.passes_filters_raw(a_start, a_end, b_start, b_end) {
                            continue;
                        }

                        let overlap_start = a_start.max(b_start);
                        let overlap_end = a_end.min(b_end);

                        // Write overlap region with A's extra fields
                        writer.write_all(chrom)?;
                        writer.write_all(b"\t")?;
                        writer.write_all(itoa_buf.format(overlap_start).as_bytes())?;
                        writer.write_all(b"\t")?;
                        writer.write_all(itoa_buf.format(overlap_end).as_bytes())?;
                        // Write A's extra fields if present
                        if rest_start < line_bytes.len() {
                            writer.write_all(&line_bytes[rest_start..])?;
                        }
                        writer.write_all(b"\n")?;
                        stats.overlaps_found += 1;
                    }
                }

                OutputMode::WriteA => {
                    // -wa: output A record once per overlap
                    for b in active_slice {
                        let b_start = b.start as u64;
                        let b_end = b.end as u64;

                        if b_end <= a_start || b_start >= a_end {
                            continue;
                        }

                        if has_filters && !self.passes_filters_raw(a_start, a_end, b_start, b_end) {
                            continue;
                        }

                        writer.write_all(line_bytes)?;
                        writer.write_all(b"\n")?;
                        stats.overlaps_found += 1;
                    }
                }

                OutputMode::WriteB => {
                    // -wb: output overlap region + B record
                    for b in active_slice {
                        let b_start = b.start as u64;
                        let b_end = b.end as u64;

                        if b_end <= a_start || b_start >= a_end {
                            continue;
                        }

                        if has_filters && !self.passes_filters_raw(a_start, a_end, b_start, b_end) {
                            continue;
                        }

                        let overlap_start = a_start.max(b_start);
                        let overlap_end = a_end.min(b_end);

                        // Write overlap region with A's extra fields
                        writer.write_all(chrom)?;
                        writer.write_all(b"\t")?;
                        writer.write_all(itoa_buf.format(overlap_start).as_bytes())?;
                        writer.write_all(b"\t")?;
                        writer.write_all(itoa_buf.format(overlap_end).as_bytes())?;
                        if rest_start < line_bytes.len() {
                            writer.write_all(&line_bytes[rest_start..])?;
                        }
                        // Tab separator + B record
                        writer.write_all(b"\t")?;
                        // Write B's raw line (already trimmed)
                        writer.write_all(&b.line)?;
                        writer.write_all(b"\n")?;
                        stats.overlaps_found += 1;
                    }
                }

                OutputMode::WriteBoth => {
                    // -wa -wb: output A + B for each overlap
                    for b in active_slice {
                        let b_start = b.start as u64;
                        let b_end = b.end as u64;

                        if b_end <= a_start || b_start >= a_end {
                            continue;
                        }

                        if has_filters && !self.passes_filters_raw(a_start, a_end, b_start, b_end) {
                            continue;
                        }

                        // Write A record
                        writer.write_all(line_bytes)?;
                        // Tab separator + B record
                        writer.write_all(b"\t")?;
                        writer.write_all(&b.line)?;
                        writer.write_all(b"\n")?;
                        stats.overlaps_found += 1;
                    }
                }
            }
        }

        // Count remaining B intervals for stats
        while pending_b.is_some() {
            stats.b_intervals += 1;
            pending_b = Self::read_next_b_optimized(&mut b_reader, &mut b_line_buf, &mut b_chrom)?;
        }

        writer.flush().map_err(BedError::Io)?;
        Ok(stats)
    }

    /// Read next B interval with zero-allocation parsing.
    #[inline]
    fn read_next_b_optimized(
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

            return Ok(Some(ActiveB {
                start: start as u32,
                end: end as u32,
                line: line.to_vec(),
            }));
        }
    }

    /// Check fraction filters without BedRecord (raw coordinates).
    #[inline]
    fn passes_filters_raw(&self, a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> bool {
        // Note: strand filtering is not supported in optimized path (no strand info stored)
        // For -s/-S flags, the old path should be used

        if let Some(frac) = self.fraction_a {
            let overlap_start = a_start.max(b_start);
            let overlap_end = a_end.min(b_end);
            let overlap_len = overlap_end.saturating_sub(overlap_start);
            let a_len = a_end.saturating_sub(a_start);
            if a_len == 0 || (overlap_len as f64 / a_len as f64) < frac {
                return false;
            }
        }

        if let Some(frac) = self.fraction_b {
            let overlap_start = a_start.max(b_start);
            let overlap_end = a_end.min(b_end);
            let overlap_len = overlap_end.saturating_sub(overlap_start);
            let b_len = b_end.saturating_sub(b_start);
            if b_len == 0 || (overlap_len as f64 / b_len as f64) < frac {
                return false;
            }
        }

        if self.reciprocal {
            if let Some(frac) = self.fraction_a.or(self.fraction_b) {
                let overlap_start = a_start.max(b_start);
                let overlap_end = a_end.min(b_end);
                let overlap_len = overlap_end.saturating_sub(overlap_start);
                let a_len = a_end.saturating_sub(a_start);
                let b_len = b_end.saturating_sub(b_start);
                if a_len == 0 || b_len == 0 {
                    return false;
                }
                let a_frac = overlap_len as f64 / a_len as f64;
                let b_frac = overlap_len as f64 / b_len as f64;
                if a_frac < frac || b_frac < frac {
                    return false;
                }
            }
        }

        true
    }

    /// Core streaming intersect algorithm.
    ///
    /// Uses a sliding window approach:
    /// - `active_b`: VecDeque of B intervals that could still overlap future A intervals
    /// - Only keeps B intervals where B.end > current_a.start
    ///
    /// # Zero-Clone Design
    ///
    /// B records are moved (not cloned) from pending_b into active_b using Option::take().
    /// This eliminates heap allocations in the hot path.
    pub fn run_streaming<R1: io::Read, R2: io::Read, W: Write>(
        &self,
        a_reader: BedReader<R1>,
        mut b_reader: BedReader<R2>,
        output: &mut W,
    ) -> Result<StreamingStats, BedError> {
        let mut stats = StreamingStats::default();
        let mut writer = BufWriter::with_capacity(256 * 1024, output);

        // Compute output mode once to avoid repeated flag checks
        let output_mode = self.compute_output_mode();
        let has_filters = self.has_filters();

        // Active B intervals that could still overlap current/future A intervals
        // Zero-clone: records are MOVED into this deque, never cloned
        let mut active_b: VecDeque<BedRecord> = VecDeque::with_capacity(256);

        // Next B record to potentially add to active set
        let mut pending_b: Option<BedRecord> = b_reader.read_record()?;

        // Current chromosome we're processing (stored as bytes to avoid allocation)
        let mut current_chrom: Option<String> = None;

        // Track if we've exhausted B for current chromosome
        let mut b_exhausted_for_chrom = false;

        // Sorted validation state
        let mut prev_a_chrom: Option<String> = None;
        let mut prev_a_start: u64 = 0;
        let mut prev_b_chrom: Option<String> = None;
        let mut prev_b_start: u64 = 0;
        let mut warned_large_window = false;

        // Track seen chromosomes to handle any sort order (genome or lexicographic)
        let mut seen_a_chroms: std::collections::HashSet<String> = std::collections::HashSet::new();
        let mut seen_b_chroms: std::collections::HashSet<String> = std::collections::HashSet::new();
        // Note: Don't add pending_b's chrom to seen_b_chroms here
        // It should only be added after validation when we process the record

        // Pre-allocated output buffer (16KB) to avoid allocations in hot path
        let mut output_buf: Vec<u8> = Vec::with_capacity(16 * 1024);

        // Cached itoa buffer for fast integer formatting (reused across all writes)
        let mut itoa_buf = itoa::Buffer::new();

        for a_result in a_reader.records() {
            let a_rec = a_result?;
            stats.a_intervals += 1;

            let a_chrom = a_rec.chrom();
            let a_start = a_rec.start();
            let a_end = a_rec.end();

            // Inline sorted validation for A
            if !self.assume_sorted {
                if let Some(ref pc) = prev_a_chrom {
                    // Detect unsorted: chromosome changed but we've seen this one before
                    if a_chrom != pc.as_str() && seen_a_chroms.contains(a_chrom) {
                        return Err(BedError::InvalidFormat(format!(
                            "File A not sorted: chromosome '{}' at record {} comes after '{}'",
                            a_chrom, stats.a_intervals, pc
                        )));
                    }
                    if a_chrom == pc && a_start < prev_a_start {
                        return Err(BedError::InvalidFormat(format!(
                            "File A not sorted: position {} at record {} comes after {} on {}",
                            a_start, stats.a_intervals, prev_a_start, a_chrom
                        )));
                    }
                }
                seen_a_chroms.insert(a_chrom.to_string());
                prev_a_chrom = Some(a_chrom.to_string());
                prev_a_start = a_start;
            }

            // Check if we've moved to a new chromosome
            let chrom_changed = current_chrom.as_ref().is_none_or(|c| c != a_chrom);

            if chrom_changed {
                // Clear active set - different chromosome
                active_b.clear();
                current_chrom = Some(a_chrom.to_string());
                b_exhausted_for_chrom = false;

                // Skip B intervals until we reach a_chrom or B has passed it
                while let Some(ref b_rec) = pending_b {
                    if b_rec.chrom() == a_chrom {
                        break; // Found matching chromosome
                    }
                    if seen_b_chroms.contains(a_chrom) {
                        // B has already passed a_chrom
                        break;
                    }
                    // Inline sorted validation for B
                    if !self.assume_sorted {
                        if let Some(ref pc) = prev_b_chrom {
                            if b_rec.chrom() != pc.as_str() && seen_b_chroms.contains(b_rec.chrom())
                            {
                                return Err(BedError::InvalidFormat(format!(
                                    "File B not sorted: chromosome '{}' comes after '{}'",
                                    b_rec.chrom(),
                                    pc
                                )));
                            }
                        }
                        prev_b_chrom = Some(b_rec.chrom().to_string());
                        prev_b_start = b_rec.start();
                    }
                    seen_b_chroms.insert(b_rec.chrom().to_string());
                    stats.b_intervals += 1;
                    pending_b = b_reader.read_record()?;
                    // Note: Don't add pending_b's chrom to seen_b_chroms here
                    // It should only be added after validation in Step 2
                }

                // Skip B intervals that end before this A starts (on same chrom)
                while let Some(ref b_rec) = pending_b {
                    if b_rec.chrom() != a_chrom {
                        break;
                    }
                    if b_rec.end() > a_start {
                        break;
                    }
                    // Inline sorted validation for B
                    if !self.assume_sorted {
                        if let Some(ref pc) = prev_b_chrom {
                            if b_rec.chrom() == pc && b_rec.start() < prev_b_start {
                                return Err(BedError::InvalidFormat(format!(
                                    "File B not sorted: position {} comes after {} on {}",
                                    b_rec.start(),
                                    prev_b_start,
                                    b_rec.chrom()
                                )));
                            }
                        }
                        prev_b_chrom = Some(b_rec.chrom().to_string());
                        prev_b_start = b_rec.start();
                    }
                    stats.b_intervals += 1;
                    pending_b = b_reader.read_record()?;
                }
            }

            // Step 1: Remove B intervals from active set that can't overlap anymore
            // Condition: B.end <= A.start (no overlap possible)
            while let Some(front) = active_b.front() {
                if front.end() <= a_start {
                    active_b.pop_front();
                } else {
                    break;
                }
            }

            // Step 2: Add new B intervals to active set (ZERO-CLONE using take())
            // Add all B intervals where B.start < A.end (could overlap A)
            if !b_exhausted_for_chrom {
                while let Some(b_rec) = pending_b.take() {
                    // Inline sorted validation for B
                    if !self.assume_sorted {
                        if let Some(ref pc) = prev_b_chrom {
                            // Detect unsorted: chromosome changed but we've seen this one before
                            if b_rec.chrom() != pc.as_str() && seen_b_chroms.contains(b_rec.chrom())
                            {
                                return Err(BedError::InvalidFormat(format!(
                                    "File B not sorted: chromosome '{}' comes after '{}'",
                                    b_rec.chrom(),
                                    pc
                                )));
                            }
                            if b_rec.chrom() == pc.as_str() && b_rec.start() < prev_b_start {
                                return Err(BedError::InvalidFormat(format!(
                                    "File B not sorted: position {} comes after {} on {}",
                                    b_rec.start(),
                                    prev_b_start,
                                    b_rec.chrom()
                                )));
                            }
                        }
                        prev_b_chrom = Some(b_rec.chrom().to_string());
                        prev_b_start = b_rec.start();
                    }
                    seen_b_chroms.insert(b_rec.chrom().to_string());

                    // Different chromosome - stop adding, put record back
                    if b_rec.chrom() != a_chrom {
                        pending_b = Some(b_rec);
                        b_exhausted_for_chrom = true;
                        break;
                    }

                    // B starts at or after A ends - can't overlap current A
                    // But might overlap future A, so put back and break
                    if b_rec.start() >= a_end {
                        pending_b = Some(b_rec);
                        break;
                    }

                    // B.start < A.end: potential overlap
                    // Only add if B.end > A.start (actual overlap region exists)
                    // Note: we check B.end > a_start, not B.end <= a_start
                    // This is the ONLY overlap condition needed since B.start < A.end is guaranteed
                    stats.b_intervals += 1;

                    if b_rec.end() > a_start {
                        // MOVE record into active set (no clone!)
                        active_b.push_back(b_rec);
                    }
                    // If B.end <= a_start, the record is simply dropped (no overlap)

                    // Read next B record
                    pending_b = b_reader.read_record()?;
                }
            }

            // Track max active set size for memory stats
            let active_size = active_b.len();
            stats.max_active_b = stats.max_active_b.max(active_size);

            // Warn on pathological case (only once)
            if self.warn_large_window
                && !warned_large_window
                && active_size > ACTIVE_WINDOW_WARNING_THRESHOLD
            {
                eprintln!(
                    "Warning: Large active window detected ({} intervals). \
                     This may indicate pathological input where many B intervals \
                     overlap a single A interval. Memory usage: O({})",
                    active_size, active_size
                );
                warned_large_window = true;
            }

            // Step 3: Process overlaps based on output mode
            // Helper closure to check if B overlaps A
            let overlaps = |b: &BedRecord| b.end() > a_start && b.start() < a_end;

            match output_mode {
                OutputMode::NoOverlap => {
                    // -v mode: output A if no overlaps found
                    let has_overlap = if has_filters {
                        active_b
                            .iter()
                            .any(|b| overlaps(b) && self.passes_filters(&a_rec, b))
                    } else {
                        active_b.iter().any(overlaps)
                    };

                    if !has_overlap {
                        output_buf.clear();
                        self.write_record(&mut output_buf, &a_rec, &mut itoa_buf);
                        writer.write_all(&output_buf)?;
                    }
                }

                OutputMode::Count => {
                    // -c mode: output A with overlap count
                    let count = if has_filters {
                        active_b
                            .iter()
                            .filter(|b| overlaps(b) && self.passes_filters(&a_rec, b))
                            .count()
                    } else {
                        active_b.iter().filter(|b| overlaps(b)).count()
                    };

                    output_buf.clear();
                    self.write_record_with_count(&mut output_buf, &a_rec, count, &mut itoa_buf);
                    writer.write_all(&output_buf)?;
                }

                OutputMode::Unique => {
                    // -u mode: output A once if any overlap exists
                    let has_overlap = if has_filters {
                        active_b
                            .iter()
                            .any(|b| overlaps(b) && self.passes_filters(&a_rec, b))
                    } else {
                        active_b.iter().any(overlaps)
                    };

                    if has_overlap {
                        output_buf.clear();
                        self.write_record(&mut output_buf, &a_rec, &mut itoa_buf);
                        writer.write_all(&output_buf)?;
                        stats.overlaps_found += 1;
                    }
                }

                OutputMode::Default => {
                    // Default: output overlap region (intersection of A and B)
                    for b_rec in active_b.iter() {
                        // Check both overlap conditions: B.end > A.start AND B.start < A.end
                        if b_rec.end() <= a_start || b_rec.start() >= a_end {
                            continue;
                        }

                        if has_filters && !self.passes_filters(&a_rec, b_rec) {
                            continue;
                        }

                        output_buf.clear();
                        self.write_overlap_region(&mut output_buf, &a_rec, b_rec, &mut itoa_buf);
                        writer.write_all(&output_buf)?;
                        stats.overlaps_found += 1;
                    }
                }

                OutputMode::WriteA => {
                    // -wa: output A record once per overlap
                    for b_rec in active_b.iter() {
                        // Check both overlap conditions
                        if b_rec.end() <= a_start || b_rec.start() >= a_end {
                            continue;
                        }

                        if has_filters && !self.passes_filters(&a_rec, b_rec) {
                            continue;
                        }

                        output_buf.clear();
                        self.write_record(&mut output_buf, &a_rec, &mut itoa_buf);
                        writer.write_all(&output_buf)?;
                        stats.overlaps_found += 1;
                    }
                }

                OutputMode::WriteB => {
                    // -wb: output overlap region + B record
                    for b_rec in active_b.iter() {
                        // Check both overlap conditions
                        if b_rec.end() <= a_start || b_rec.start() >= a_end {
                            continue;
                        }

                        if has_filters && !self.passes_filters(&a_rec, b_rec) {
                            continue;
                        }

                        output_buf.clear();
                        self.write_overlap_with_b(&mut output_buf, &a_rec, b_rec, &mut itoa_buf);
                        writer.write_all(&output_buf)?;
                        stats.overlaps_found += 1;
                    }
                }

                OutputMode::WriteBoth => {
                    // -wa -wb: output A + B for each overlap
                    for b_rec in active_b.iter() {
                        // Check both overlap conditions
                        if b_rec.end() <= a_start || b_rec.start() >= a_end {
                            continue;
                        }

                        if has_filters && !self.passes_filters(&a_rec, b_rec) {
                            continue;
                        }

                        output_buf.clear();
                        self.write_both_records(&mut output_buf, &a_rec, b_rec, &mut itoa_buf);
                        writer.write_all(&output_buf)?;
                        stats.overlaps_found += 1;
                    }
                }
            }
        }

        // Count remaining B intervals for stats
        while pending_b.is_some() {
            stats.b_intervals += 1;
            pending_b = b_reader.read_record()?;
        }

        writer.flush().map_err(BedError::Io)?;
        Ok(stats)
    }

    /// Check if overlap passes fraction and strand filters.
    #[inline]
    fn passes_filters(&self, a: &BedRecord, b: &BedRecord) -> bool {
        // Strand filtering
        if self.same_strand {
            match (a.strand, b.strand) {
                (Some(a_s), Some(b_s)) if a_s != b_s => return false,
                _ => {}
            }
        }
        if self.opposite_strand {
            match (a.strand, b.strand) {
                (Some(a_s), Some(b_s)) if a_s == b_s => return false,
                _ => {}
            }
        }

        // Fraction filtering
        if let Some(frac) = self.fraction_a {
            if !a.interval.overlaps_by_fraction(&b.interval, frac) {
                return false;
            }
        }

        if let Some(frac) = self.fraction_b {
            if !b.interval.overlaps_by_fraction(&a.interval, frac) {
                return false;
            }
        }

        if self.reciprocal {
            if let Some(frac) = self.fraction_a.or(self.fraction_b) {
                if !a.interval.overlaps_reciprocal(&b.interval, frac) {
                    return false;
                }
            }
        }

        true
    }

    /// Write chrom, start, end using fast itoa formatting.
    #[inline]
    fn write_bed3(
        &self,
        buf: &mut Vec<u8>,
        chrom: &str,
        start: u64,
        end: u64,
        itoa_buf: &mut itoa::Buffer,
    ) {
        buf.extend_from_slice(chrom.as_bytes());
        buf.push(b'\t');
        buf.extend_from_slice(itoa_buf.format(start).as_bytes());
        buf.push(b'\t');
        buf.extend_from_slice(itoa_buf.format(end).as_bytes());
    }

    /// Write a single BED record.
    #[inline]
    fn write_record(&self, buf: &mut Vec<u8>, rec: &BedRecord, itoa_buf: &mut itoa::Buffer) {
        self.write_bed3(buf, rec.chrom(), rec.start(), rec.end(), itoa_buf);
        self.write_optional_fields(buf, rec, itoa_buf);
        buf.push(b'\n');
    }

    /// Write overlap region (intersection of A and B) with A's extra fields.
    #[inline]
    fn write_overlap_region(
        &self,
        buf: &mut Vec<u8>,
        a: &BedRecord,
        b: &BedRecord,
        itoa_buf: &mut itoa::Buffer,
    ) {
        let overlap_start = a.start().max(b.start());
        let overlap_end = a.end().min(b.end());
        self.write_bed3(buf, a.chrom(), overlap_start, overlap_end, itoa_buf);
        self.write_optional_fields(buf, a, itoa_buf);
        buf.push(b'\n');
    }

    /// Write overlap region + B record (for -wb mode).
    #[inline]
    fn write_overlap_with_b(
        &self,
        buf: &mut Vec<u8>,
        a: &BedRecord,
        b: &BedRecord,
        itoa_buf: &mut itoa::Buffer,
    ) {
        let overlap_start = a.start().max(b.start());
        let overlap_end = a.end().min(b.end());
        // Write overlap region with A's extra fields
        self.write_bed3(buf, a.chrom(), overlap_start, overlap_end, itoa_buf);
        self.write_optional_fields(buf, a, itoa_buf);
        // Tab separator
        buf.push(b'\t');
        // Write B record
        self.write_bed3(buf, b.chrom(), b.start(), b.end(), itoa_buf);
        self.write_optional_fields(buf, b, itoa_buf);
        buf.push(b'\n');
    }

    /// Write A record with overlap count.
    #[inline]
    fn write_record_with_count(
        &self,
        buf: &mut Vec<u8>,
        rec: &BedRecord,
        count: usize,
        itoa_buf: &mut itoa::Buffer,
    ) {
        self.write_bed3(buf, rec.chrom(), rec.start(), rec.end(), itoa_buf);
        self.write_optional_fields(buf, rec, itoa_buf);
        buf.push(b'\t');
        buf.extend_from_slice(itoa_buf.format(count).as_bytes());
        buf.push(b'\n');
    }

    /// Write both A and B records (for -wa -wb mode).
    #[inline]
    fn write_both_records(
        &self,
        buf: &mut Vec<u8>,
        a: &BedRecord,
        b: &BedRecord,
        itoa_buf: &mut itoa::Buffer,
    ) {
        // Write A record
        self.write_bed3(buf, a.chrom(), a.start(), a.end(), itoa_buf);
        self.write_optional_fields(buf, a, itoa_buf);
        // Tab separator
        buf.push(b'\t');
        // Write B record
        self.write_bed3(buf, b.chrom(), b.start(), b.end(), itoa_buf);
        self.write_optional_fields(buf, b, itoa_buf);
        buf.push(b'\n');
    }

    #[inline]
    fn write_optional_fields(
        &self,
        buf: &mut Vec<u8>,
        rec: &BedRecord,
        itoa_buf: &mut itoa::Buffer,
    ) {
        use crate::interval::Strand;
        if let Some(ref name) = rec.name {
            buf.push(b'\t');
            buf.extend_from_slice(name.as_bytes());
            if let Some(score) = rec.score {
                buf.push(b'\t');
                buf.extend_from_slice(itoa_buf.format(score as i64).as_bytes());
                if let Some(strand) = rec.strand {
                    buf.push(b'\t');
                    match strand {
                        Strand::Plus => buf.push(b'+'),
                        Strand::Minus => buf.push(b'-'),
                        Strand::Unknown => buf.push(b'.'),
                    }
                }
            }
        }
    }
}

/// Statistics from streaming intersect operation.
#[derive(Debug, Default, Clone)]
pub struct StreamingStats {
    /// Number of A intervals processed
    pub a_intervals: usize,
    /// Number of B intervals processed
    pub b_intervals: usize,
    /// Number of overlaps found
    pub overlaps_found: usize,
    /// Maximum size of active B set (memory high-water mark)
    pub max_active_b: usize,
}

impl std::fmt::Display for StreamingStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "A intervals: {}, B intervals: {}, Overlaps: {}, Max active B: {}",
            self.a_intervals, self.b_intervals, self.overlaps_found, self.max_active_b
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

    // ==================== bedtools CLI Compatibility Tests ====================

    #[test]
    fn test_default_output_prints_overlap_region() {
        // Default bedtools behavior: print overlap region (intersection)
        let a_content = make_bed_content(&[("chr1", 100, 200)]);
        let b_content = make_bed_content(&[("chr1", 150, 250), ("chr1", 175, 225)]);

        let cmd = StreamingIntersectCommand::new();
        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        cmd.run_streaming(a_reader, b_reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        // Overlap region printed for each B overlap
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "chr1\t150\t200"); // intersection of 100-200 and 150-250
        assert_eq!(lines[1], "chr1\t175\t200"); // intersection of 100-200 and 175-225
    }

    #[test]
    fn test_wa_flag_prints_a_record() {
        let a_content = make_bed_content(&[("chr1", 100, 200)]);
        let b_content = make_bed_content(&[("chr1", 150, 250)]);

        let mut cmd = StreamingIntersectCommand::new();
        cmd.write_a = true;

        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        cmd.run_streaming(a_reader, b_reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert_eq!(result.trim(), "chr1\t100\t200");
    }

    #[test]
    fn test_wb_flag_prints_overlap_plus_b_record() {
        let a_content = make_bed_content(&[("chr1", 100, 200)]);
        let b_content = make_bed_content(&[("chr1", 150, 250)]);

        let mut cmd = StreamingIntersectCommand::new();
        cmd.write_b = true;

        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        cmd.run_streaming(a_reader, b_reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        // -wb outputs: overlap_region + B_record
        assert_eq!(result.trim(), "chr1\t150\t200\tchr1\t150\t250");
    }

    #[test]
    fn test_wa_wb_flags_print_both_records() {
        let a_content = make_bed_content(&[("chr1", 100, 200)]);
        let b_content = make_bed_content(&[("chr1", 150, 250)]);

        let mut cmd = StreamingIntersectCommand::new();
        cmd.write_a = true;
        cmd.write_b = true;

        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        cmd.run_streaming(a_reader, b_reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert_eq!(result.trim(), "chr1\t100\t200\tchr1\t150\t250");
    }

    #[test]
    fn test_c_flag_prints_count() {
        let a_content = make_bed_content(&[("chr1", 100, 500), ("chr1", 600, 700)]);
        let b_content =
            make_bed_content(&[("chr1", 150, 200), ("chr1", 250, 300), ("chr1", 350, 400)]);

        let mut cmd = StreamingIntersectCommand::new();
        cmd.count = true;

        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        cmd.run_streaming(a_reader, b_reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 2);
        assert!(lines[0].ends_with("\t3")); // 3 overlaps
        assert!(lines[1].ends_with("\t0")); // 0 overlaps
    }

    #[test]
    fn test_u_flag_prints_unique() {
        let a_content = make_bed_content(&[("chr1", 100, 500)]);
        let b_content = make_bed_content(&[("chr1", 150, 200), ("chr1", 250, 300)]);

        let mut cmd = StreamingIntersectCommand::new();
        cmd.unique = true;

        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        let stats = cmd.run_streaming(a_reader, b_reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        // Only one output even though 2 overlaps
        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0], "chr1\t100\t500");
        assert_eq!(stats.overlaps_found, 1);
    }

    #[test]
    fn test_v_flag_prints_no_overlap() {
        let a_content = make_bed_content(&[("chr1", 100, 200), ("chr1", 500, 600)]);
        let b_content = make_bed_content(&[("chr1", 300, 400)]);

        let mut cmd = StreamingIntersectCommand::new();
        cmd.no_overlap = true;

        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        cmd.run_streaming(a_reader, b_reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        // Both A intervals have no overlap with B
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "chr1\t100\t200");
        assert_eq!(lines[1], "chr1\t500\t600");
    }

    // ==================== Sorted Validation Tests ====================

    #[test]
    fn test_unsorted_a_fails() {
        // Test actual interleaving: chr1, chr2, chr1 (chr1 appears again after chr2)
        // Note: chr2 then chr1 is valid (non-interleaved, just different order)
        let a_content =
            make_bed_content(&[("chr1", 100, 200), ("chr2", 100, 200), ("chr1", 300, 400)]);
        let b_content = make_bed_content(&[("chr1", 150, 250)]);

        let cmd = StreamingIntersectCommand::new();
        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        let result = cmd.run_streaming(a_reader, b_reader, &mut output);

        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("not sorted"));
    }

    #[test]
    fn test_unsorted_a_position_fails() {
        let a_content = make_bed_content(&[("chr1", 200, 300), ("chr1", 100, 200)]);
        let b_content = make_bed_content(&[("chr1", 150, 250)]);

        let cmd = StreamingIntersectCommand::new();
        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        let result = cmd.run_streaming(a_reader, b_reader, &mut output);

        assert!(result.is_err());
    }

    #[test]
    fn test_assume_sorted_skips_validation() {
        // Even with unsorted input, --assume-sorted skips validation
        // (This may produce incorrect results, but that's the user's choice)
        let a_content = make_bed_content(&[("chr1", 100, 200)]);
        let b_content = make_bed_content(&[("chr1", 150, 250)]);

        let mut cmd = StreamingIntersectCommand::new();
        cmd.assume_sorted = true;

        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        let result = cmd.run_streaming(a_reader, b_reader, &mut output);

        assert!(result.is_ok());
    }

    // ==================== Core Streaming Tests ====================

    #[test]
    fn test_basic_streaming_intersect() {
        let a_content =
            make_bed_content(&[("chr1", 100, 200), ("chr1", 300, 400), ("chr1", 500, 600)]);
        let b_content = make_bed_content(&[("chr1", 150, 250), ("chr1", 350, 450)]);

        let cmd = StreamingIntersectCommand::new();
        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        let stats = cmd.run_streaming(a_reader, b_reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 2);
        // Default outputs overlap region
        assert!(lines[0].contains("150\t200")); // Overlap of 100-200 and 150-250
        assert!(lines[1].contains("350\t400")); // Overlap of 300-400 and 350-450
        assert_eq!(stats.max_active_b, 1);
    }

    #[test]
    fn test_streaming_multiple_chromosomes() {
        let a_content =
            make_bed_content(&[("chr1", 100, 200), ("chr2", 100, 200), ("chr3", 100, 200)]);
        let b_content = make_bed_content(&[("chr1", 150, 250), ("chr2", 150, 250)]);

        let cmd = StreamingIntersectCommand::new();
        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        let stats = cmd.run_streaming(a_reader, b_reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 2);
        assert_eq!(stats.a_intervals, 3);
    }

    // ==================== Pathological Case Tests ====================

    #[test]
    fn test_pathological_all_b_overlap_one_a() {
        // Worst case: all B intervals overlap a single A interval
        let a_content = make_bed_content(&[("chr1", 100, 1000)]);
        let b_content = make_bed_content(&[
            ("chr1", 100, 200),
            ("chr1", 150, 250),
            ("chr1", 200, 300),
            ("chr1", 250, 350),
            ("chr1", 300, 400),
            ("chr1", 350, 450),
            ("chr1", 400, 500),
        ]);

        let mut cmd = StreamingIntersectCommand::new();
        cmd.warn_large_window = false; // Disable warning for test

        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        let stats = cmd.run_streaming(a_reader, b_reader, &mut output).unwrap();

        assert_eq!(stats.overlaps_found, 7);
        assert_eq!(stats.max_active_b, 7);
    }

    #[test]
    fn test_many_b_per_a_stress() {
        // Generate 1000 B intervals all overlapping one A
        let a_content = "chr1\t0\t1000000\n".to_string();
        let b_content: String = (0..1000)
            .map(|i| format!("chr1\t{}\t{}", i * 100, i * 100 + 150))
            .collect::<Vec<_>>()
            .join("\n");

        let mut cmd = StreamingIntersectCommand::new();
        cmd.count = true;
        cmd.warn_large_window = false;

        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        let stats = cmd.run_streaming(a_reader, b_reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("\t1000")); // All 1000 B intervals overlap
        assert!(stats.max_active_b <= 1000);
    }

    // ==================== Chromosome Boundary Tests ====================

    #[test]
    fn test_chromosome_boundary_transition() {
        // Test proper handling when chromosomes change
        let a_content = make_bed_content(&[
            ("chr1", 100, 200),
            ("chr1", 900, 1000),
            ("chr2", 100, 200),
            ("chr3", 100, 200),
        ]);
        let b_content = make_bed_content(&[
            ("chr1", 50, 150),
            ("chr1", 850, 950),
            ("chr2", 150, 250),
            ("chr4", 100, 200), // No A on chr4
        ]);

        let cmd = StreamingIntersectCommand::new();
        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        let stats = cmd.run_streaming(a_reader, b_reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 3); // chr1:100-200, chr1:900-1000, chr2:100-200
        assert_eq!(stats.a_intervals, 4);
        assert_eq!(stats.b_intervals, 4);
    }

    #[test]
    fn test_b_chromosome_before_all_a() {
        // B has chromosomes that come before any A chromosomes
        let a_content = make_bed_content(&[("chr2", 100, 200)]);
        let b_content = make_bed_content(&[("chr1", 100, 200), ("chr2", 150, 250)]);

        let cmd = StreamingIntersectCommand::new();
        let a_reader = BedReader::new(a_content.as_bytes());
        let b_reader = BedReader::new(b_content.as_bytes());

        let mut output = Vec::new();
        let stats = cmd.run_streaming(a_reader, b_reader, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 1);
        assert_eq!(stats.b_intervals, 2);
    }

    // ==================== Zero-Clone Verification ====================

    #[test]
    fn test_output_mode_computation() {
        let mut cmd = StreamingIntersectCommand::new();
        assert_eq!(cmd.compute_output_mode(), OutputMode::Default);

        cmd.write_a = true;
        assert_eq!(cmd.compute_output_mode(), OutputMode::WriteA);

        cmd.write_b = true;
        assert_eq!(cmd.compute_output_mode(), OutputMode::WriteBoth);

        cmd.write_a = false;
        assert_eq!(cmd.compute_output_mode(), OutputMode::WriteB);

        cmd.write_b = false;
        cmd.count = true;
        assert_eq!(cmd.compute_output_mode(), OutputMode::Count);

        cmd.count = false;
        cmd.unique = true;
        assert_eq!(cmd.compute_output_mode(), OutputMode::Unique);

        cmd.unique = false;
        cmd.no_overlap = true;
        assert_eq!(cmd.compute_output_mode(), OutputMode::NoOverlap);
    }
}
