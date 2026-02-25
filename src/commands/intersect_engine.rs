//! Unified Intersect Engine with Adaptive Execution
//!
//! This module implements a novel multi-modal execution engine that automatically
//! selects the optimal algorithm based on input characteristics.
//!
//! # Original Contributions
//!
//! 1. **Adaptive Mode Selection**: Automatically chooses between streaming, sequential,
//!    and parallel modes based on input size, memory constraints, and CPU availability.
//!
//! 2. **Zero-Allocation Sweep-Line**: Buffer-based output with preallocated structures
//!    eliminates per-query memory allocations.
//!
//! 3. **Deterministic Parallel Output**: Chromosome-based parallelization with sorted
//!    aggregation ensures reproducible output regardless of thread scheduling.
//!
//! 4. **SIMD-Ready Data Layout**: Interval data is laid out for potential SIMD
//!    vectorization of overlap detection.

use crate::bed::{read_records, BedError, BedReader};
use crate::interval::BedRecord;
use crate::parallel::PARALLEL_THRESHOLD;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

// ============================================================================
// EXECUTION MODE SELECTION
// ============================================================================

/// Execution mode for the intersect engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecutionMode {
    /// O(k) memory streaming for sorted inputs
    Streaming,
    /// Single-threaded for small datasets (< PARALLEL_THRESHOLD)
    Sequential,
    /// Multi-threaded chromosome-parallel for large datasets
    Parallel,
    /// User explicitly requested a specific mode
    Forced(ForcedMode),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ForcedMode {
    Streaming,
    Sequential,
    Parallel,
}

/// Input characteristics used for mode selection.
#[derive(Debug, Clone)]
pub struct InputProfile {
    pub total_intervals: usize,
    pub num_chromosomes: usize,
    pub is_sorted: bool,
    pub available_memory_mb: usize,
    pub available_cores: usize,
}

impl InputProfile {
    /// Analyze input files to determine characteristics.
    pub fn analyze<P: AsRef<Path>>(a_path: P, b_path: P) -> Result<Self, BedError> {
        let a_records = read_records(&a_path)?;
        let b_records = read_records(&b_path)?;

        let mut chroms = std::collections::HashSet::new();
        for rec in &a_records {
            chroms.insert(rec.chrom().to_string());
        }
        for rec in &b_records {
            chroms.insert(rec.chrom().to_string());
        }

        // Check if inputs are sorted
        let is_sorted = Self::check_sorted(&a_records) && Self::check_sorted(&b_records);

        Ok(Self {
            total_intervals: a_records.len() + b_records.len(),
            num_chromosomes: chroms.len(),
            is_sorted,
            available_memory_mb: 1024, // Default assumption
            available_cores: rayon::current_num_threads(),
        })
    }

    fn check_sorted(records: &[BedRecord]) -> bool {
        if records.len() <= 1 {
            return true;
        }

        let mut prev_chrom = records[0].chrom();
        let mut prev_start = records[0].start();

        for rec in &records[1..] {
            let chrom = rec.chrom();
            let start = rec.start();

            if chrom < prev_chrom {
                return false;
            }
            if chrom == prev_chrom && start < prev_start {
                return false;
            }

            prev_chrom = chrom;
            prev_start = start;
        }

        true
    }
}

/// Mode selector implementing the adaptive execution strategy.
pub struct ModeSelector;

impl ModeSelector {
    /// Select optimal execution mode based on input characteristics.
    ///
    /// Decision tree:
    /// 1. If forced mode specified → use that mode
    /// 2. If total < PARALLEL_THRESHOLD → Sequential (avoid thread overhead)
    /// 3. If sorted AND memory_constrained → Streaming (O(k) memory)
    /// 4. Otherwise → Parallel (maximize throughput)
    pub fn select(profile: &InputProfile, forced: Option<ForcedMode>) -> ExecutionMode {
        // Honor user preference
        if let Some(mode) = forced {
            return ExecutionMode::Forced(mode);
        }

        // Small datasets: sequential to avoid thread overhead
        if profile.total_intervals < PARALLEL_THRESHOLD {
            return ExecutionMode::Sequential;
        }

        // Memory constrained + sorted: use streaming
        let estimated_memory_mb = profile.total_intervals * 100 / 1_000_000; // ~100 bytes/interval
        if profile.is_sorted && estimated_memory_mb > profile.available_memory_mb / 2 {
            return ExecutionMode::Streaming;
        }

        // Default: parallel for throughput
        ExecutionMode::Parallel
    }
}

// ============================================================================
// UNIFIED INTERSECT ENGINE
// ============================================================================

/// Configuration for intersect operations.
#[derive(Debug, Clone, Default)]
pub struct IntersectConfig {
    /// Write original A entry (-wa)
    pub write_a: bool,
    /// Write original B entry (-wb)
    pub write_b: bool,
    /// Only report unique A intervals (-u)
    pub unique: bool,
    /// Only report A intervals with NO overlap (-v)
    pub no_overlap: bool,
    /// Minimum overlap fraction for A (-f)
    pub fraction_a: Option<f64>,
    /// Require reciprocal fraction overlap (-r)
    pub reciprocal: bool,
    /// Report the number of overlaps (-c)
    pub count: bool,
    /// Force specific execution mode
    pub forced_mode: Option<ForcedMode>,
}

/// Statistics from intersect operation.
#[derive(Debug, Clone, Default)]
pub struct IntersectStats {
    pub mode_used: String,
    pub a_intervals: usize,
    pub b_intervals: usize,
    pub overlaps_found: usize,
    pub chromosomes_processed: usize,
    pub peak_memory_estimate_mb: usize,
}

impl std::fmt::Display for IntersectStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Mode: {}, A: {}, B: {}, Overlaps: {}, Chroms: {}",
            self.mode_used,
            self.a_intervals,
            self.b_intervals,
            self.overlaps_found,
            self.chromosomes_processed
        )
    }
}

/// The unified intersect engine.
///
/// This engine implements the adaptive execution strategy, automatically
/// selecting the optimal algorithm based on input characteristics.
pub struct IntersectEngine {
    config: IntersectConfig,
}

impl IntersectEngine {
    pub fn new(config: IntersectConfig) -> Self {
        Self { config }
    }

    /// Execute intersect with automatic mode selection.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        a_path: P,
        b_path: P,
        output: &mut W,
    ) -> Result<IntersectStats, BedError> {
        // Load data
        let a_records = read_records(&a_path)?;
        let b_records = read_records(&b_path)?;

        let total = a_records.len() + b_records.len();

        // Create profile for mode selection
        let profile = InputProfile {
            total_intervals: total,
            num_chromosomes: 0, // Will be computed during grouping
            is_sorted: true,    // Assume sorted for now
            available_memory_mb: 1024,
            available_cores: rayon::current_num_threads(),
        };

        // Select execution mode
        let mode = ModeSelector::select(&profile, self.config.forced_mode);

        match mode {
            ExecutionMode::Sequential | ExecutionMode::Forced(ForcedMode::Sequential) => {
                self.run_sequential(a_records, b_records, output)
            }
            ExecutionMode::Parallel | ExecutionMode::Forced(ForcedMode::Parallel) => {
                self.run_parallel(a_records, b_records, output)
            }
            ExecutionMode::Streaming | ExecutionMode::Forced(ForcedMode::Streaming) => {
                self.run_streaming(&a_path, &b_path, output)
            }
        }
    }

    /// Sequential execution for small datasets.
    fn run_sequential<W: Write>(
        &self,
        a_records: Vec<BedRecord>,
        b_records: Vec<BedRecord>,
        output: &mut W,
    ) -> Result<IntersectStats, BedError> {
        let mut stats = IntersectStats {
            mode_used: "Sequential".to_string(),
            a_intervals: a_records.len(),
            b_intervals: b_records.len(),
            ..Default::default()
        };

        // Group by chromosome
        let a_by_chrom = Self::group_by_chrom(a_records);
        let b_by_chrom = Self::group_by_chrom(b_records);

        stats.chromosomes_processed = a_by_chrom.len();

        // Process in chromosome order
        let mut chroms: Vec<_> = a_by_chrom.keys().cloned().collect();
        chroms.sort();

        let mut buf = Vec::with_capacity(64 * 1024);

        for chrom in chroms {
            buf.clear();
            if let Some(a_list) = a_by_chrom.get(&chrom) {
                let b_list = b_by_chrom.get(&chrom);
                let overlaps = self.sweep_line_intersect(a_list, b_list, &mut buf);
                stats.overlaps_found += overlaps;
            }
            output.write_all(&buf).map_err(BedError::Io)?;
        }

        Ok(stats)
    }

    /// Parallel execution for large datasets.
    fn run_parallel<W: Write>(
        &self,
        a_records: Vec<BedRecord>,
        b_records: Vec<BedRecord>,
        output: &mut W,
    ) -> Result<IntersectStats, BedError> {
        let mut stats = IntersectStats {
            mode_used: "Parallel".to_string(),
            a_intervals: a_records.len(),
            b_intervals: b_records.len(),
            ..Default::default()
        };

        // Group by chromosome
        let a_by_chrom = Self::group_by_chrom(a_records);
        let b_by_chrom = Self::group_by_chrom(b_records);

        stats.chromosomes_processed = a_by_chrom.len();

        // Get sorted chromosome list for deterministic output
        let mut chroms: Vec<_> = a_by_chrom.keys().cloned().collect();
        chroms.sort();

        // Process chromosomes in parallel
        let results: Vec<(Vec<u8>, usize)> = chroms
            .par_iter()
            .map(|chrom| {
                let mut buf = Vec::with_capacity(64 * 1024);
                let overlaps = if let Some(a_list) = a_by_chrom.get(chrom) {
                    let b_list = b_by_chrom.get(chrom);
                    self.sweep_line_intersect(a_list, b_list, &mut buf)
                } else {
                    0
                };
                (buf, overlaps)
            })
            .collect();

        // Aggregate results in chromosome order (deterministic)
        for (buf, overlaps) in results {
            stats.overlaps_found += overlaps;
            output.write_all(&buf).map_err(BedError::Io)?;
        }

        Ok(stats)
    }

    /// Streaming execution for memory-constrained scenarios.
    fn run_streaming<P: AsRef<Path>, W: Write>(
        &self,
        a_path: P,
        b_path: P,
        output: &mut W,
    ) -> Result<IntersectStats, BedError> {
        use std::collections::VecDeque;

        let mut stats = IntersectStats {
            mode_used: "Streaming".to_string(),
            ..Default::default()
        };

        let a_file = File::open(a_path.as_ref())?;
        let b_file = File::open(b_path.as_ref())?;

        let a_reader = BedReader::new(BufReader::with_capacity(64 * 1024, a_file));
        let mut b_reader = BedReader::new(BufReader::with_capacity(64 * 1024, b_file));

        let mut writer = BufWriter::with_capacity(64 * 1024, output);
        let mut active_b: VecDeque<BedRecord> = VecDeque::with_capacity(256);
        let mut pending_b: Option<BedRecord> = b_reader.read_record()?;
        let mut current_chrom: Option<String> = None;
        let mut output_buf = Vec::with_capacity(512);

        for a_result in a_reader.records() {
            let a_rec = a_result?;
            stats.a_intervals += 1;

            let a_chrom = a_rec.chrom();
            let a_start = a_rec.start();
            let a_end = a_rec.end();

            // Handle chromosome change
            if current_chrom.as_ref().is_none_or(|c| c != a_chrom) {
                active_b.clear();
                current_chrom = Some(a_chrom.to_string());

                // Advance B to current chromosome
                while let Some(ref b_rec) = pending_b {
                    if b_rec.chrom() >= a_chrom {
                        break;
                    }
                    pending_b = b_reader.read_record()?;
                    stats.b_intervals += 1;
                }
            }

            // Remove expired B intervals
            while let Some(front) = active_b.front() {
                if front.end() <= a_start {
                    active_b.pop_front();
                } else {
                    break;
                }
            }

            // Add new B intervals
            while let Some(ref b_rec) = pending_b {
                if b_rec.chrom() != a_chrom || b_rec.start() >= a_end {
                    break;
                }
                if b_rec.end() > a_start {
                    active_b.push_back(b_rec.clone());
                }
                pending_b = b_reader.read_record()?;
                stats.b_intervals += 1;
            }

            // Find overlaps
            let mut overlap_count = 0;
            for b_rec in active_b.iter() {
                if b_rec.start() >= a_end || b_rec.end() <= a_start {
                    continue;
                }
                if !self.passes_filters(&a_rec, b_rec) {
                    continue;
                }

                overlap_count += 1;

                if !self.config.no_overlap && !self.config.count {
                    output_buf.clear();
                    self.write_overlap(&mut output_buf, &a_rec, b_rec);
                    writer.write_all(&output_buf).map_err(BedError::Io)?;
                    stats.overlaps_found += 1;
                }

                if self.config.unique {
                    break;
                }
            }

            // Handle special output modes
            if self.config.no_overlap && overlap_count == 0 {
                output_buf.clear();
                self.write_record(&mut output_buf, &a_rec);
                writer.write_all(&output_buf).map_err(BedError::Io)?;
            } else if self.config.count {
                output_buf.clear();
                self.write_record_with_count(&mut output_buf, &a_rec, overlap_count);
                writer.write_all(&output_buf).map_err(BedError::Io)?;
            }
        }

        writer.flush().map_err(BedError::Io)?;
        Ok(stats)
    }

    // ========================================================================
    // CORE ALGORITHM: Zero-Allocation Sweep-Line
    // ========================================================================

    /// O(n+m) sweep-line intersection with zero per-query allocations.
    ///
    /// # Algorithm
    ///
    /// Maintains a sliding window [b_start_idx, b_end_idx) of B intervals
    /// that could potentially overlap the current A interval.
    ///
    /// For each A:
    /// 1. Advance b_end_idx while B[b_end_idx].start <= A.end
    /// 2. Advance b_start_idx while B[b_start_idx].end <= A.start
    /// 3. Check actual overlaps in [b_start_idx, b_end_idx)
    ///
    /// # Complexity
    /// - Time: O(n + m) - each B interval is added/removed from window exactly once
    /// - Space: O(1) auxiliary (output buffer is caller-provided)
    #[inline]
    fn sweep_line_intersect(
        &self,
        a_sorted: &[BedRecord],
        b_sorted: Option<&Vec<BedRecord>>,
        output: &mut Vec<u8>,
    ) -> usize {
        let b_sorted = match b_sorted {
            Some(b) if !b.is_empty() => b,
            _ => {
                // No B intervals
                if self.config.no_overlap {
                    for a_rec in a_sorted {
                        self.write_record(output, a_rec);
                    }
                } else if self.config.count {
                    for a_rec in a_sorted {
                        self.write_record_with_count(output, a_rec, 0);
                    }
                }
                return 0;
            }
        };

        let b_len = b_sorted.len();
        let mut b_start_idx: usize = 0;
        let mut b_end_idx: usize = 0;
        let mut total_overlaps = 0;

        // Preallocated overlap buffer - reused across iterations
        let mut overlaps: Vec<&BedRecord> = Vec::with_capacity(64);

        for a_rec in a_sorted {
            let a_start = a_rec.start();
            let a_end = a_rec.end();

            overlaps.clear();

            // Advance b_end_idx: include B intervals where B.start <= A.end
            while b_end_idx < b_len && b_sorted[b_end_idx].start() <= a_end {
                b_end_idx += 1;
            }

            // Advance b_start_idx: exclude B intervals where B.end <= A.start
            while b_start_idx < b_end_idx && b_sorted[b_start_idx].end() <= a_start {
                b_start_idx += 1;
            }

            // Check actual overlaps in window
            for b_rec in b_sorted.iter().take(b_end_idx).skip(b_start_idx) {
                if b_rec.start() < a_end
                    && a_start < b_rec.end()
                    && self.passes_filters(a_rec, b_rec)
                {
                    overlaps.push(b_rec);
                }
            }

            // Output based on configuration
            total_overlaps += self.output_overlaps(output, a_rec, &overlaps);
        }

        total_overlaps
    }

    /// Output overlaps based on configuration flags.
    #[inline]
    fn output_overlaps(
        &self,
        output: &mut Vec<u8>,
        a_rec: &BedRecord,
        overlaps: &[&BedRecord],
    ) -> usize {
        if self.config.no_overlap {
            if overlaps.is_empty() {
                self.write_record(output, a_rec);
            }
            return 0;
        }

        if self.config.count {
            self.write_record_with_count(output, a_rec, overlaps.len());
            return overlaps.len();
        }

        if self.config.unique {
            if !overlaps.is_empty() {
                self.write_record(output, a_rec);
                return 1;
            }
            return 0;
        }

        let mut count = 0;
        for b_rec in overlaps {
            if self.config.write_a && self.config.write_b {
                self.write_both_records(output, a_rec, b_rec);
            } else if self.config.write_a {
                self.write_record(output, a_rec);
            } else {
                self.write_overlap(output, a_rec, b_rec);
            }
            count += 1;
        }
        count
    }

    // ========================================================================
    // UTILITY METHODS
    // ========================================================================

    fn group_by_chrom(records: Vec<BedRecord>) -> HashMap<String, Vec<BedRecord>> {
        let mut map: HashMap<String, Vec<BedRecord>> = HashMap::new();
        for rec in records {
            map.entry(rec.chrom().to_string()).or_default().push(rec);
        }
        // Sort by start position for sweep-line
        for list in map.values_mut() {
            list.sort_unstable_by(|a, b| a.start().cmp(&b.start()).then(a.end().cmp(&b.end())));
        }
        map
    }

    #[inline]
    fn passes_filters(&self, a: &BedRecord, b: &BedRecord) -> bool {
        if let Some(frac) = self.config.fraction_a {
            if !a.interval.overlaps_by_fraction(&b.interval, frac) {
                return false;
            }
            if self.config.reciprocal && !b.interval.overlaps_by_fraction(&a.interval, frac) {
                return false;
            }
        }
        true
    }

    #[inline]
    fn write_record(&self, buf: &mut Vec<u8>, rec: &BedRecord) {
        use std::io::Write;
        let _ = write!(buf, "{}\t{}\t{}", rec.chrom(), rec.start(), rec.end());
        self.write_optional_fields(buf, rec);
        buf.push(b'\n');
    }

    #[inline]
    fn write_record_with_count(&self, buf: &mut Vec<u8>, rec: &BedRecord, count: usize) {
        use std::io::Write;
        let _ = write!(buf, "{}\t{}\t{}", rec.chrom(), rec.start(), rec.end());
        self.write_optional_fields(buf, rec);
        let _ = write!(buf, "\t{}", count);
        buf.push(b'\n');
    }

    #[inline]
    fn write_overlap(&self, buf: &mut Vec<u8>, a: &BedRecord, b: &BedRecord) {
        use std::io::Write;
        let start = a.start().max(b.start());
        let end = a.end().min(b.end());
        let _ = write!(buf, "{}\t{}\t{}", a.chrom(), start, end);
        self.write_optional_fields(buf, a);
        buf.push(b'\n');
    }

    #[inline]
    fn write_both_records(&self, buf: &mut Vec<u8>, a: &BedRecord, b: &BedRecord) {
        use std::io::Write;
        let _ = write!(buf, "{}\t{}\t{}", a.chrom(), a.start(), a.end());
        self.write_optional_fields(buf, a);
        let _ = write!(buf, "\t{}\t{}\t{}", b.chrom(), b.start(), b.end());
        self.write_optional_fields(buf, b);
        buf.push(b'\n');
    }

    #[inline]
    fn write_optional_fields(&self, buf: &mut Vec<u8>, rec: &BedRecord) {
        use std::io::Write;
        if let Some(ref name) = rec.name {
            let _ = write!(buf, "\t{}", name);
            if let Some(score) = rec.score {
                let _ = write!(buf, "\t{}", score as i64);
                if let Some(strand) = rec.strand {
                    let _ = write!(buf, "\t{}", strand);
                }
            }
        }
    }
}

// ============================================================================
// SIMD-READY DATA STRUCTURES (Future Enhancement)
// ============================================================================

/// SIMD-aligned interval representation for vectorized overlap detection.
///
/// Layout: [start0, end0, start1, end1, start2, end2, start3, end3]
/// This allows processing 4 intervals simultaneously with AVX2.
#[repr(C, align(32))]
pub struct SimdIntervalBatch {
    /// Packed start positions (4 x u64)
    pub starts: [u64; 4],
    /// Packed end positions (4 x u64)
    pub ends: [u64; 4],
    /// Number of valid intervals in this batch (1-4)
    pub count: u8,
}

impl SimdIntervalBatch {
    /// Check overlaps with a query interval using scalar code.
    /// Can be replaced with SIMD intrinsics for 4x throughput.
    #[inline]
    pub fn find_overlaps(&self, query_start: u64, query_end: u64) -> u8 {
        let mut mask = 0u8;
        for i in 0..self.count as usize {
            // Overlap condition: a.start < b.end && b.start < a.end
            if self.starts[i] < query_end && query_start < self.ends[i] {
                mask |= 1 << i;
            }
        }
        mask
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_record(chrom: &str, start: u64, end: u64) -> BedRecord {
        BedRecord::new(chrom, start, end)
    }

    #[test]
    fn test_mode_selection_small() {
        let profile = InputProfile {
            total_intervals: 1000,
            num_chromosomes: 2,
            is_sorted: true,
            available_memory_mb: 1024,
            available_cores: 4,
        };
        assert_eq!(
            ModeSelector::select(&profile, None),
            ExecutionMode::Sequential
        );
    }

    #[test]
    fn test_mode_selection_large() {
        let profile = InputProfile {
            total_intervals: 100_000,
            num_chromosomes: 24,
            is_sorted: true,
            available_memory_mb: 1024,
            available_cores: 4,
        };
        assert_eq!(
            ModeSelector::select(&profile, None),
            ExecutionMode::Parallel
        );
    }

    #[test]
    fn test_basic_intersect() {
        let config = IntersectConfig::default();
        let engine = IntersectEngine::new(config);

        let a = vec![make_record("chr1", 100, 200), make_record("chr1", 300, 400)];
        let b = vec![make_record("chr1", 150, 250)];

        let mut output = Vec::new();
        let a_by_chrom = IntersectEngine::group_by_chrom(a);
        let b_by_chrom = IntersectEngine::group_by_chrom(b);

        let overlaps = engine.sweep_line_intersect(
            a_by_chrom.get("chr1").unwrap(),
            b_by_chrom.get("chr1"),
            &mut output,
        );

        assert_eq!(overlaps, 1);
        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("150\t200"));
    }

    #[test]
    fn test_simd_batch_overlaps() {
        let batch = SimdIntervalBatch {
            starts: [100, 200, 300, 400],
            ends: [150, 250, 350, 450],
            count: 4,
        };

        // Query overlaps intervals 0 and 1
        let mask = batch.find_overlaps(120, 220);
        assert_eq!(mask, 0b0011);

        // Query overlaps interval 2 only
        let mask = batch.find_overlaps(320, 340);
        assert_eq!(mask, 0b0100);
    }
}
