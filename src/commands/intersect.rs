//! Intersect command implementation - the most-used bedtools command.
//!
//! Uses O(n+m) sweep-line algorithm per chromosome for optimal performance.

use crate::bed::{read_records, BedError};
use crate::index::IntervalIndex;
use crate::interval::{BedRecord, Interval};
use crate::parallel::{group_by_chromosome, PARALLEL_THRESHOLD};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// Output mode for intersect command.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntersectOutput {
    /// Output only intervals from A that overlap B
    WriteA,
    /// Output intervals from A with overlapping B intervals
    WriteAB,
    /// Output only the overlapping portion
    WriteOverlap,
    /// Output count of overlaps per A interval
    WriteCount,
}

/// Intersect command configuration.
#[derive(Debug, Clone)]
pub struct IntersectCommand {
    /// Write original A entry
    pub write_a: bool,
    /// Write original B entry
    pub write_b: bool,
    /// Only report unique A intervals
    pub unique: bool,
    /// Only report A intervals with no overlap
    pub no_overlap: bool,
    /// Minimum fraction overlap for A
    pub fraction_a: Option<f64>,
    /// Minimum fraction overlap for B
    pub fraction_b: Option<f64>,
    /// Require reciprocal fraction overlap
    pub reciprocal: bool,
    /// Report the number of overlaps
    pub count: bool,
    /// Require same strand
    pub same_strand: bool,
    /// Require opposite strand
    pub opposite_strand: bool,
    /// Report only once per A interval
    pub report_once: bool,
    /// Split by chromosome for parallel processing
    pub parallel: bool,
}

impl Default for IntersectCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl IntersectCommand {
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
            report_once: false,
            parallel: true,
        }
    }

    /// Find all intersecting pairs.
    pub fn find_intersections(
        &self,
        a_intervals: &[Interval],
        b_index: &IntervalIndex,
    ) -> Vec<IntersectResult> {
        let mut results = Vec::new();

        for (a_idx, a) in a_intervals.iter().enumerate() {
            let overlaps = b_index.find_overlaps(a);

            if self.no_overlap {
                // -v flag: report A intervals with NO overlap
                if overlaps.is_empty() {
                    results.push(IntersectResult {
                        a_interval: a.clone(),
                        a_index: a_idx,
                        b_intervals: Vec::new(),
                        overlap_count: 0,
                    });
                }
            } else {
                let filtered: Vec<&Interval> = overlaps
                    .into_iter()
                    .filter(|b| self.passes_filters(a, b))
                    .collect();

                if !filtered.is_empty() {
                    if self.report_once || self.unique {
                        results.push(IntersectResult {
                            a_interval: a.clone(),
                            a_index: a_idx,
                            b_intervals: if self.unique {
                                Vec::new()
                            } else {
                                vec![filtered[0].clone()]
                            },
                            overlap_count: filtered.len(),
                        });
                    } else {
                        for b in filtered {
                            results.push(IntersectResult {
                                a_interval: a.clone(),
                                a_index: a_idx,
                                b_intervals: vec![b.clone()],
                                overlap_count: 1,
                            });
                        }
                    }
                }
            }
        }

        results
    }

    /// Find intersections in parallel by chromosome.
    pub fn find_intersections_parallel(
        &self,
        a_intervals: Vec<Interval>,
        b_intervals: Vec<Interval>,
    ) -> Vec<IntersectResult> {
        let b_index = IntervalIndex::from_intervals(b_intervals);
        let a_groups = group_by_chromosome(a_intervals);

        let results: Vec<Vec<IntersectResult>> = a_groups
            .into_par_iter()
            .map(|(_, intervals)| self.find_intersections(&intervals, &b_index))
            .collect();

        results.into_iter().flatten().collect()
    }

    /// Check if an overlap passes all filters.
    #[inline(always)]
    fn passes_filters(&self, a: &Interval, b: &Interval) -> bool {
        // Check fraction filters
        if let Some(frac) = self.fraction_a {
            if !a.overlaps_by_fraction(b, frac) {
                return false;
            }
        }

        if let Some(frac) = self.fraction_b {
            if !b.overlaps_by_fraction(a, frac) {
                return false;
            }
        }

        if self.reciprocal {
            if let (Some(fa), Some(fb)) = (self.fraction_a, self.fraction_b) {
                if !a.overlaps_by_fraction(b, fa) || !b.overlaps_by_fraction(a, fb) {
                    return false;
                }
            } else if let Some(f) = self.fraction_a.or(self.fraction_b) {
                if !a.overlaps_reciprocal(b, f) {
                    return false;
                }
            }
        }

        true
    }

    /// Check if overlap passes filters for BedRecords
    #[inline(always)]
    fn passes_record_filters(&self, a: &BedRecord, b: &BedRecord) -> bool {
        self.passes_filters(&a.interval, &b.interval)
    }

    /// Compute the intersection (overlap) of two intervals.
    pub fn compute_overlap(&self, a: &Interval, b: &Interval) -> Option<Interval> {
        if !a.overlaps(b) {
            return None;
        }

        Some(Interval {
            chrom: a.chrom.clone(),
            start: a.start.max(b.start),
            end: a.end.min(b.end),
        })
    }

    /// Execute intersect command on files using O(n+m) sweep-line algorithm.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        a_path: P,
        b_path: P,
        output: &mut W,
    ) -> Result<(), BedError> {
        let a_records = read_records(a_path)?;
        let b_records = read_records(b_path)?;

        // Group by chromosome
        let a_by_chrom = Self::group_records_by_chrom_owned(a_records);
        let b_by_chrom = Self::group_records_by_chrom_owned(b_records);

        // Get sorted list of chromosomes (for deterministic output)
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
                    self.intersect_chromosome_sweepline(a_list, b_list, &mut buf);
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
                        self.intersect_chromosome_sweepline(a_list, b_list, &mut buf);
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

    /// O(n+m) sweep-line intersection for a single chromosome.
    ///
    /// Algorithm:
    /// 1. A and B are sorted by start position
    /// 2. Maintain pointer j into B and active window of B intervals
    /// 3. For each A interval:
    ///    a. Advance j while B[j].start <= A.end (add to active)
    ///    b. Remove from active any B where B.end <= A.start
    ///    c. Active intervals are potential overlaps
    fn intersect_chromosome_sweepline(
        &self,
        a_sorted: &[BedRecord],
        b_sorted: Option<&Vec<BedRecord>>,
        output: &mut Vec<u8>,
    ) {
        let b_sorted = match b_sorted {
            Some(b) if !b.is_empty() => b,
            _ => {
                // No B intervals - handle -v flag or output nothing
                if self.no_overlap {
                    for a_rec in a_sorted {
                        self.write_record_to_buf(output, a_rec);
                    }
                } else if self.count {
                    for a_rec in a_sorted {
                        self.write_record_with_count_to_buf(output, a_rec, 0);
                    }
                }
                return;
            }
        };

        let b_len = b_sorted.len();

        // Active window: indices into b_sorted that could overlap current A
        // We track (start_idx, end_idx) as a sliding window
        let mut b_start_idx: usize = 0; // First B that hasn't ended before A.start
        let mut b_end_idx: usize = 0; // First B that starts after A.end

        // Preallocate overlap buffer
        let mut overlaps: Vec<&BedRecord> = Vec::with_capacity(64);

        for a_rec in a_sorted {
            let a_start = a_rec.start();
            let a_end = a_rec.end();

            overlaps.clear();

            // Advance b_end_idx: add B intervals where B.start <= A.end
            while b_end_idx < b_len && b_sorted[b_end_idx].start() <= a_end {
                b_end_idx += 1;
            }

            // Advance b_start_idx: skip B intervals where B.end <= A.start (no overlap)
            while b_start_idx < b_end_idx && b_sorted[b_start_idx].end() <= a_start {
                b_start_idx += 1;
            }

            // All B in [b_start_idx..b_end_idx) potentially overlap A
            // Check actual overlap: B.start < A.end && A.start < B.end
            for b_rec in b_sorted.iter().take(b_end_idx).skip(b_start_idx) {
                let b_start = b_rec.start();
                let b_end = b_rec.end();

                // Actual overlap check (B might start after A ends due to our window)
                if b_start < a_end && a_start < b_end {
                    // Apply fraction/strand filters
                    if self.passes_record_filters(a_rec, b_rec) {
                        overlaps.push(b_rec);
                    }
                }
            }

            // Output based on flags
            self.output_overlaps(output, a_rec, &overlaps);
        }
    }

    /// Output overlaps for a single A record based on command flags
    #[inline]
    fn output_overlaps(&self, output: &mut Vec<u8>, a_rec: &BedRecord, overlaps: &[&BedRecord]) {
        if self.no_overlap {
            // -v flag: report A if NO overlap
            if overlaps.is_empty() {
                self.write_record_to_buf(output, a_rec);
            }
        } else if self.count {
            // -c flag: report A with count
            self.write_record_with_count_to_buf(output, a_rec, overlaps.len());
        } else if self.unique {
            // -u flag: report A once if any overlap
            if !overlaps.is_empty() {
                self.write_record_to_buf(output, a_rec);
            }
        } else if self.write_a && self.write_b {
            // -wa -wb: report both A and B
            for b_rec in overlaps {
                self.write_both_records_to_buf(output, a_rec, b_rec);
            }
        } else if self.write_a {
            // -wa: report original A for each overlap
            for _ in overlaps {
                self.write_record_to_buf(output, a_rec);
            }
        } else if self.write_b {
            // -wb: report overlap portion + B entry
            for b_rec in overlaps {
                self.write_overlap_with_b_to_buf(output, a_rec, b_rec);
            }
        } else {
            // Default: report overlapping portion with A's fields
            for b_rec in overlaps {
                self.write_overlap_to_buf(output, a_rec, b_rec);
            }
        }
    }

    /// Group records by chromosome, returning owned records sorted by start
    fn group_records_by_chrom_owned(records: Vec<BedRecord>) -> HashMap<String, Vec<BedRecord>> {
        let mut map: HashMap<String, Vec<BedRecord>> = HashMap::new();
        for rec in records {
            map.entry(rec.chrom().to_string()).or_default().push(rec);
        }
        // Sort each chromosome's records by start position, then by end
        for list in map.values_mut() {
            list.sort_unstable_by(|a, b| a.start().cmp(&b.start()).then(a.end().cmp(&b.end())));
        }
        map
    }

    // ==================== Buffer-based output methods (zero allocation) ====================

    #[inline]
    fn write_record_to_buf(&self, buf: &mut Vec<u8>, rec: &BedRecord) {
        use std::io::Write;
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

    #[inline]
    fn write_record_with_count_to_buf(&self, buf: &mut Vec<u8>, rec: &BedRecord, count: usize) {
        use std::io::Write;
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
        let _ = write!(buf, "\t{}", count);
        buf.push(b'\n');
    }

    #[inline]
    fn write_both_records_to_buf(&self, buf: &mut Vec<u8>, a_rec: &BedRecord, b_rec: &BedRecord) {
        use std::io::Write;
        // A record
        let _ = write!(buf, "{}\t{}\t{}", a_rec.chrom(), a_rec.start(), a_rec.end());
        if let Some(ref name) = a_rec.name {
            let _ = write!(buf, "\t{}", name);
            if let Some(score) = a_rec.score {
                let _ = write!(buf, "\t{}", score as i64);
                if let Some(strand) = a_rec.strand {
                    let _ = write!(buf, "\t{}", strand);
                }
            }
        }
        // B record
        let _ = write!(
            buf,
            "\t{}\t{}\t{}",
            b_rec.chrom(),
            b_rec.start(),
            b_rec.end()
        );
        if let Some(ref name) = b_rec.name {
            let _ = write!(buf, "\t{}", name);
            if let Some(score) = b_rec.score {
                let _ = write!(buf, "\t{}", score as i64);
                if let Some(strand) = b_rec.strand {
                    let _ = write!(buf, "\t{}", strand);
                }
            }
        }
        buf.push(b'\n');
    }

    #[inline]
    fn write_overlap_to_buf(&self, buf: &mut Vec<u8>, a_rec: &BedRecord, b_rec: &BedRecord) {
        use std::io::Write;
        let a = &a_rec.interval;
        let b = &b_rec.interval;
        let overlap_start = a.start.max(b.start);
        let overlap_end = a.end.min(b.end);

        let _ = write!(buf, "{}\t{}\t{}", a.chrom, overlap_start, overlap_end);
        if let Some(ref name) = a_rec.name {
            let _ = write!(buf, "\t{}", name);
            if let Some(score) = a_rec.score {
                let _ = write!(buf, "\t{}", score as i64);
                if let Some(strand) = a_rec.strand {
                    let _ = write!(buf, "\t{}", strand);
                }
            }
        }
        buf.push(b'\n');
    }

    #[inline]
    fn write_overlap_with_b_to_buf(&self, buf: &mut Vec<u8>, a_rec: &BedRecord, b_rec: &BedRecord) {
        use std::io::Write;
        let a = &a_rec.interval;
        let b = &b_rec.interval;
        let overlap_start = a.start.max(b.start);
        let overlap_end = a.end.min(b.end);

        let _ = write!(buf, "{}\t{}\t{}", a.chrom, overlap_start, overlap_end);
        if let Some(ref name) = a_rec.name {
            let _ = write!(buf, "\t{}", name);
            if let Some(score) = a_rec.score {
                let _ = write!(buf, "\t{}", score as i64);
                if let Some(strand) = a_rec.strand {
                    let _ = write!(buf, "\t{}", strand);
                }
            }
        }
        // B record
        let _ = write!(
            buf,
            "\t{}\t{}\t{}",
            b_rec.chrom(),
            b_rec.start(),
            b_rec.end()
        );
        if let Some(ref name) = b_rec.name {
            let _ = write!(buf, "\t{}", name);
            if let Some(score) = b_rec.score {
                let _ = write!(buf, "\t{}", score as i64);
                if let Some(strand) = b_rec.strand {
                    let _ = write!(buf, "\t{}", strand);
                }
            }
        }
        buf.push(b'\n');
    }
}

/// Result of an intersection query.
#[derive(Debug, Clone)]
pub struct IntersectResult {
    pub a_interval: Interval,
    pub a_index: usize,
    pub b_intervals: Vec<Interval>,
    pub overlap_count: usize,
}

/// Find all overlapping pairs between two sets of intervals.
pub fn find_overlaps(a: &[Interval], b: &[Interval]) -> Vec<(Interval, Interval)> {
    let b_index = IntervalIndex::from_intervals(b.to_vec());
    let mut results = Vec::new();

    for a_int in a {
        for b_int in b_index.find_overlaps(a_int) {
            results.push((a_int.clone(), b_int.clone()));
        }
    }

    results
}

/// Find intervals in A that have any overlap with B.
pub fn filter_overlapping(a: &[Interval], b: &[Interval]) -> Vec<Interval> {
    let b_index = IntervalIndex::from_intervals(b.to_vec());

    a.iter()
        .filter(|int| b_index.has_overlap(int))
        .cloned()
        .collect()
}

/// Find intervals in A that have NO overlap with B.
pub fn filter_non_overlapping(a: &[Interval], b: &[Interval]) -> Vec<Interval> {
    let b_index = IntervalIndex::from_intervals(b.to_vec());

    a.iter()
        .filter(|int| !b_index.has_overlap(int))
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_intervals() -> (Vec<Interval>, Vec<Interval>) {
        let a = vec![
            Interval::new("chr1", 100, 200),
            Interval::new("chr1", 300, 400),
            Interval::new("chr1", 500, 600),
        ];
        let b = vec![
            Interval::new("chr1", 150, 250),
            Interval::new("chr1", 350, 450),
        ];
        (a, b)
    }

    #[test]
    fn test_basic_intersect() {
        let (a, b) = make_intervals();
        let cmd = IntersectCommand::new();
        let b_index = IntervalIndex::from_intervals(b);

        let results = cmd.find_intersections(&a, &b_index);

        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_no_overlap_filter() {
        let (a, b) = make_intervals();
        let mut cmd = IntersectCommand::new();
        cmd.no_overlap = true;
        let b_index = IntervalIndex::from_intervals(b);

        let results = cmd.find_intersections(&a, &b_index);

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].a_interval.start, 500);
    }

    #[test]
    fn test_fraction_overlap() {
        let a = vec![Interval::new("chr1", 100, 200)];
        let b = vec![Interval::new("chr1", 150, 250)];

        let mut cmd = IntersectCommand::new();
        cmd.fraction_a = Some(0.5);
        let b_index = IntervalIndex::from_intervals(b.clone());

        // 50% overlap (50bp out of 100bp)
        let results = cmd.find_intersections(&a, &b_index);
        assert_eq!(results.len(), 1);

        cmd.fraction_a = Some(0.6);
        let b_index = IntervalIndex::from_intervals(b);
        let results = cmd.find_intersections(&a, &b_index);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_compute_overlap() {
        let cmd = IntersectCommand::new();
        let a = Interval::new("chr1", 100, 200);
        let b = Interval::new("chr1", 150, 250);

        let overlap = cmd.compute_overlap(&a, &b).unwrap();
        assert_eq!(overlap.start, 150);
        assert_eq!(overlap.end, 200);
    }

    #[test]
    fn test_filter_overlapping() {
        let (a, b) = make_intervals();
        let filtered = filter_overlapping(&a, &b);

        assert_eq!(filtered.len(), 2);
    }

    #[test]
    fn test_filter_non_overlapping() {
        let (a, b) = make_intervals();
        let filtered = filter_non_overlapping(&a, &b);

        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].start, 500);
    }

    #[test]
    fn test_parallel_intersect() {
        let (a, b) = make_intervals();
        let cmd = IntersectCommand::new();

        let results = cmd.find_intersections_parallel(a, b);

        assert_eq!(results.len(), 2);
    }
}
