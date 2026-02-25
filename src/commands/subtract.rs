//! Subtract command implementation.
//!
//! Uses O(n + m) sweep-line algorithm per chromosome for optimal performance.

use crate::bed::{read_records, BedError};
use crate::index::IntervalIndex;
use crate::interval::{BedRecord, Interval};
use crate::parallel::{group_by_chromosome, PARALLEL_THRESHOLD};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// Subtract command configuration.
#[derive(Debug, Clone)]
pub struct SubtractCommand {
    /// Remove entire A feature if any overlap (like bedtools -A)
    pub remove_entire: bool,
    /// Minimum overlap fraction required to subtract
    pub fraction: Option<f64>,
    /// Require reciprocal fraction overlap
    pub reciprocal: bool,
    /// Require same strand
    pub same_strand: bool,
    /// Process in parallel by chromosome
    pub parallel: bool,
}

impl Default for SubtractCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl SubtractCommand {
    pub fn new() -> Self {
        Self {
            remove_entire: false,
            fraction: None,
            reciprocal: false,
            same_strand: false,
            parallel: true,
        }
    }

    /// Subtract B intervals from A intervals.
    pub fn subtract(&self, a_intervals: &[Interval], b_index: &IntervalIndex) -> Vec<Interval> {
        let mut results = Vec::new();

        for a in a_intervals {
            let overlaps = b_index.find_overlaps(a);

            if overlaps.is_empty() {
                results.push(a.clone());
                continue;
            }

            // Filter by fraction if specified
            let filtered: Vec<&Interval> = overlaps
                .into_iter()
                .filter(|b| self.passes_filters(a, b))
                .collect();

            if filtered.is_empty() {
                results.push(a.clone());
                continue;
            }

            if self.remove_entire {
                // -A flag: remove entire feature if any overlap
                continue;
            }

            // Subtract each overlapping region
            let mut remaining = vec![a.clone()];
            for b in filtered {
                let mut new_remaining = Vec::new();
                for r in remaining {
                    new_remaining.extend(r.subtract(b));
                }
                remaining = new_remaining;
            }

            results.extend(remaining);
        }

        results
    }

    /// Subtract in parallel by chromosome.
    pub fn subtract_parallel(
        &self,
        a_intervals: Vec<Interval>,
        b_intervals: Vec<Interval>,
    ) -> Vec<Interval> {
        let b_index = IntervalIndex::from_intervals(b_intervals);
        let a_groups = group_by_chromosome(a_intervals);

        let results: Vec<Vec<Interval>> = a_groups
            .into_par_iter()
            .map(|(_, intervals)| self.subtract(&intervals, &b_index))
            .collect();

        let mut final_results: Vec<Interval> = results.into_iter().flatten().collect();
        final_results.sort();
        final_results
    }

    /// Check if an overlap passes the fraction filter.
    #[inline]
    fn passes_filters(&self, a: &Interval, b: &Interval) -> bool {
        if let Some(frac) = self.fraction {
            if self.reciprocal {
                if !a.overlaps_reciprocal(b, frac) {
                    return false;
                }
            } else if !a.overlaps_by_fraction(b, frac) {
                return false;
            }
        }
        true
    }

    /// Check if an overlap passes the fraction filter for BedRecords.
    #[inline]
    fn passes_record_filters(&self, a: &BedRecord, b: &BedRecord) -> bool {
        self.passes_filters(&a.interval, &b.interval)
    }

    /// Execute subtract command on files using O(n+m) sweep-line algorithm.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        a_path: P,
        b_path: P,
        output: &mut W,
    ) -> Result<(), BedError> {
        let a_records = read_records(a_path)?;
        let b_records = read_records(b_path)?;

        if a_records.is_empty() {
            return Ok(());
        }

        // Group by chromosome
        let a_by_chrom = Self::group_records_by_chrom_owned(a_records);
        let b_by_chrom = Self::group_records_by_chrom_owned(b_records);

        // Get sorted chromosome list for deterministic output
        let mut chroms: Vec<_> = a_by_chrom.keys().cloned().collect();
        chroms.sort();

        // Calculate total intervals for threshold check
        let total: usize = a_by_chrom.values().map(|v| v.len()).sum::<usize>()
            + b_by_chrom.values().map(|v| v.len()).sum::<usize>();

        if total < PARALLEL_THRESHOLD {
            // Sequential processing for small datasets
            for chrom in &chroms {
                let mut buf = Vec::with_capacity(64 * 1024);
                if let Some(a_list) = a_by_chrom.get(chrom) {
                    let b_list = b_by_chrom.get(chrom);
                    self.subtract_chromosome_sweepline(a_list, b_list, &mut buf);
                }
                output.write_all(&buf).map_err(BedError::Io)?;
            }
        } else {
            // Parallel processing for large datasets
            let results: Vec<Vec<u8>> = chroms
                .par_iter()
                .map(|chrom| {
                    let mut buf = Vec::with_capacity(64 * 1024);
                    if let Some(a_list) = a_by_chrom.get(chrom) {
                        let b_list = b_by_chrom.get(chrom);
                        self.subtract_chromosome_sweepline(a_list, b_list, &mut buf);
                    }
                    buf
                })
                .collect();

            // Write results in chromosome order
            for buf in results {
                output.write_all(&buf).map_err(BedError::Io)?;
            }
        }

        Ok(())
    }

    /// O(n+m) sweep-line subtract for a single chromosome.
    ///
    /// Algorithm:
    /// 1. A and B are sorted by start position
    /// 2. Maintain pointer j into B
    /// 3. For each A:
    ///    a. Advance j while B[j].end <= A.start (skip non-overlapping)
    ///    b. Collect overlapping B intervals (B.start < A.end)
    ///    c. Subtract and emit fragments directly
    fn subtract_chromosome_sweepline(
        &self,
        a_sorted: &[BedRecord],
        b_sorted: Option<&Vec<BedRecord>>,
        output: &mut Vec<u8>,
    ) {
        let b_sorted = match b_sorted {
            Some(b) if !b.is_empty() => b,
            _ => {
                // No B intervals - output all A unchanged
                for a_rec in a_sorted {
                    self.write_record_to_buf(output, a_rec);
                }
                return;
            }
        };

        let b_len = b_sorted.len();

        // Sweep-line pointer into B
        let mut b_start: usize = 0;

        for a_rec in a_sorted {
            let a_start = a_rec.start();
            let a_end = a_rec.end();

            // Advance b_start: skip B intervals that end before A starts
            while b_start < b_len && b_sorted[b_start].end() <= a_start {
                b_start += 1;
            }

            // Find overlapping B intervals: B.start < A.end AND B.end > A.start
            // Since B is sorted by start, scan from b_start until B.start >= A.end
            let mut has_valid_overlap = false;
            let mut overlap_start = b_start;
            let mut overlap_end = b_start;

            #[allow(clippy::needless_range_loop)] // j used for overlap_start/overlap_end tracking
            for j in b_start..b_len {
                let b_rec = &b_sorted[j];
                if b_rec.start() >= a_end {
                    break; // No more overlaps possible
                }
                // B.start < A.end (checked above) AND B.end > A.start (from b_start advancement)
                if b_rec.end() > a_start {
                    // Check fraction filter
                    if self.passes_record_filters(a_rec, b_rec) {
                        if !has_valid_overlap {
                            overlap_start = j;
                            has_valid_overlap = true;
                        }
                        overlap_end = j + 1;
                    }
                }
            }

            if !has_valid_overlap {
                // No valid overlaps - output A unchanged
                self.write_record_to_buf(output, a_rec);
                continue;
            }

            if self.remove_entire {
                // -A flag: remove entire feature if any overlap
                continue;
            }

            // Subtract overlapping B intervals from A
            // Use in-place subtraction to avoid allocations
            self.subtract_and_emit(output, a_rec, &b_sorted[overlap_start..overlap_end]);
        }
    }

    /// Subtract B intervals from A and emit results directly to buffer.
    /// Avoids intermediate Vec allocations.
    #[inline]
    fn subtract_and_emit(
        &self,
        output: &mut Vec<u8>,
        a_rec: &BedRecord,
        b_intervals: &[BedRecord],
    ) {
        // Sort B intervals by start for correct subtraction order
        // (They should already be sorted, but ensure correctness)
        let mut b_sorted: Vec<&BedRecord> = b_intervals.iter().collect();
        b_sorted.sort_unstable_by_key(|b| b.start());

        let a_start = a_rec.start();
        let a_end = a_rec.end();

        // Current position in A that we're processing
        let mut current_pos = a_start;

        for b_rec in b_sorted {
            let b_start = b_rec.start();
            let b_end = b_rec.end();

            // Skip B intervals that don't overlap current remaining A
            if b_end <= current_pos {
                continue;
            }
            if b_start >= a_end {
                break;
            }

            // Emit fragment before B (if any)
            if b_start > current_pos {
                let frag_end = b_start.min(a_end);
                if frag_end > current_pos {
                    self.write_fragment_to_buf(output, a_rec, current_pos, frag_end);
                }
            }

            // Advance current position past B
            current_pos = current_pos.max(b_end);
        }

        // Emit remaining fragment after all B intervals
        if current_pos < a_end {
            self.write_fragment_to_buf(output, a_rec, current_pos, a_end);
        }
    }

    /// Write a full record to buffer.
    #[inline]
    fn write_record_to_buf(&self, buf: &mut Vec<u8>, rec: &BedRecord) {
        use std::io::Write as IoWrite;
        let _ = write!(buf, "{}\t{}\t{}", rec.chrom(), rec.start(), rec.end());
        if let Some(ref name) = rec.name {
            let _ = write!(buf, "\t{}", name);
            if let Some(score) = rec.score {
                let _ = write!(buf, "\t{}", score as i64);
                if let Some(strand) = rec.strand {
                    let _ = write!(buf, "\t{}", strand);
                }
            }
        }
        buf.push(b'\n');
    }

    /// Write a fragment (with modified start/end) to buffer.
    #[inline]
    fn write_fragment_to_buf(&self, buf: &mut Vec<u8>, rec: &BedRecord, start: u64, end: u64) {
        use std::io::Write as IoWrite;
        let _ = write!(buf, "{}\t{}\t{}", rec.chrom(), start, end);
        if let Some(ref name) = rec.name {
            let _ = write!(buf, "\t{}", name);
            if let Some(score) = rec.score {
                let _ = write!(buf, "\t{}", score as i64);
                if let Some(strand) = rec.strand {
                    let _ = write!(buf, "\t{}", strand);
                }
            }
        }
        buf.push(b'\n');
    }

    /// Group records by chromosome, returning owned sorted records.
    fn group_records_by_chrom_owned(records: Vec<BedRecord>) -> HashMap<String, Vec<BedRecord>> {
        let mut map: HashMap<String, Vec<BedRecord>> = HashMap::new();
        for rec in records {
            map.entry(rec.chrom().to_string()).or_default().push(rec);
        }
        // Sort each chromosome's records by start position
        map.par_iter_mut().for_each(|(_, records)| {
            records.sort_unstable_by(|a, b| a.start().cmp(&b.start()).then(a.end().cmp(&b.end())));
        });
        map
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_subtract() {
        let cmd = SubtractCommand::new();
        let a = vec![Interval::new("chr1", 100, 300)];
        let b = vec![Interval::new("chr1", 150, 200)];
        let b_index = IntervalIndex::from_intervals(b);

        let results = cmd.subtract(&a, &b_index);

        assert_eq!(results.len(), 2);
        assert_eq!(results[0].start, 100);
        assert_eq!(results[0].end, 150);
        assert_eq!(results[1].start, 200);
        assert_eq!(results[1].end, 300);
    }

    #[test]
    fn test_subtract_no_overlap() {
        let cmd = SubtractCommand::new();
        let a = vec![Interval::new("chr1", 100, 200)];
        let b = vec![Interval::new("chr1", 300, 400)];
        let b_index = IntervalIndex::from_intervals(b);

        let results = cmd.subtract(&a, &b_index);

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].start, 100);
        assert_eq!(results[0].end, 200);
    }

    #[test]
    fn test_subtract_complete_overlap() {
        let cmd = SubtractCommand::new();
        let a = vec![Interval::new("chr1", 100, 200)];
        let b = vec![Interval::new("chr1", 50, 250)];
        let b_index = IntervalIndex::from_intervals(b);

        let results = cmd.subtract(&a, &b_index);

        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_subtract_remove_entire() {
        let mut cmd = SubtractCommand::new();
        cmd.remove_entire = true;

        let a = vec![Interval::new("chr1", 100, 300)];
        let b = vec![Interval::new("chr1", 150, 200)];
        let b_index = IntervalIndex::from_intervals(b);

        let results = cmd.subtract(&a, &b_index);

        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_subtract_multiple_overlaps() {
        let cmd = SubtractCommand::new();
        let a = vec![Interval::new("chr1", 100, 500)];
        let b = vec![
            Interval::new("chr1", 150, 200),
            Interval::new("chr1", 300, 350),
        ];
        let b_index = IntervalIndex::from_intervals(b);

        let results = cmd.subtract(&a, &b_index);

        assert_eq!(results.len(), 3);
        // 100-150, 200-300, 350-500
        assert_eq!(results[0].start, 100);
        assert_eq!(results[0].end, 150);
        assert_eq!(results[1].start, 200);
        assert_eq!(results[1].end, 300);
        assert_eq!(results[2].start, 350);
        assert_eq!(results[2].end, 500);
    }

    #[test]
    fn test_subtract_with_fraction() {
        let mut cmd = SubtractCommand::new();
        cmd.fraction = Some(0.5);

        let a = vec![Interval::new("chr1", 100, 200)];
        let b = vec![Interval::new("chr1", 175, 225)]; // Only 25% overlap
        let b_index = IntervalIndex::from_intervals(b);

        let results = cmd.subtract(&a, &b_index);

        // Overlap doesn't meet threshold, so A is returned unchanged
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].start, 100);
        assert_eq!(results[0].end, 200);
    }

    #[test]
    fn test_parallel_subtract() {
        let cmd = SubtractCommand::new();
        let a = vec![
            Interval::new("chr1", 100, 300),
            Interval::new("chr2", 100, 300),
        ];
        let b = vec![
            Interval::new("chr1", 150, 200),
            Interval::new("chr2", 150, 200),
        ];

        let results = cmd.subtract_parallel(a, b);

        assert_eq!(results.len(), 4); // 2 pieces from each chromosome
    }
}
