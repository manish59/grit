//! Closest command implementation - find nearest non-overlapping intervals.
//!
//! Uses O(n log m) algorithm per chromosome with binary search and limited scans.

use crate::bed::{read_records, BedError};
use crate::interval::{BedRecord, Interval};
use crate::parallel::{group_by_chromosome, PARALLEL_THRESHOLD};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// How to handle ties (multiple equally close intervals).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TieHandling {
    /// Report all ties
    All,
    /// Report first tie only
    First,
    /// Report last tie only
    Last,
}

/// Closest command configuration.
#[derive(Debug, Clone)]
pub struct ClosestCommand {
    /// Report distance in output
    pub report_distance: bool,
    /// How to handle ties
    pub tie_handling: TieHandling,
    /// Ignore overlapping intervals
    pub ignore_overlaps: bool,
    /// Ignore upstream intervals
    pub ignore_upstream: bool,
    /// Ignore downstream intervals
    pub ignore_downstream: bool,
    /// Require same strand
    pub same_strand: bool,
    /// Require opposite strand
    pub opposite_strand: bool,
    /// Maximum distance to report
    pub max_distance: Option<u64>,
    /// Process in parallel by chromosome
    pub parallel: bool,
}

impl Default for ClosestCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl ClosestCommand {
    pub fn new() -> Self {
        Self {
            report_distance: false,
            tie_handling: TieHandling::All,
            ignore_overlaps: false,
            ignore_upstream: false,
            ignore_downstream: false,
            same_strand: false,
            opposite_strand: false,
            max_distance: None,
            parallel: true,
        }
    }

    /// Find the closest B intervals for each A interval.
    pub fn find_closest(
        &self,
        a_intervals: &[Interval],
        b_intervals: &[Interval],
    ) -> Vec<ClosestResult> {
        // Sort B intervals by start position for binary search
        let mut b_sorted: Vec<&Interval> = b_intervals.iter().collect();
        b_sorted.sort_by(|a, b| a.start.cmp(&b.start).then(a.end.cmp(&b.end)));

        let mut results = Vec::new();

        for a in a_intervals {
            let closest = self.find_closest_for_interval(a, &b_sorted);
            results.push(ClosestResult {
                a_interval: a.clone(),
                closest_intervals: closest,
            });
        }

        results
    }

    /// Find closest intervals for a single query interval.
    fn find_closest_for_interval(
        &self,
        a: &Interval,
        b_sorted: &[&Interval],
    ) -> Vec<(Interval, i64)> {
        let mut candidates: Vec<(Interval, i64)> = Vec::new();
        let mut min_distance = i64::MAX;

        for b in b_sorted.iter() {
            if a.chrom != b.chrom {
                continue;
            }

            // Check if overlapping
            let overlaps = a.overlaps(b);
            if overlaps && self.ignore_overlaps {
                continue;
            }

            // Calculate signed distance (bedtools semantics: +1 for non-overlapping)
            let distance = if overlaps {
                0i64
            } else if b.end <= a.start {
                // B is upstream (left)
                if self.ignore_upstream {
                    continue;
                }
                -((a.start - b.end + 1) as i64)
            } else {
                // B is downstream (right)
                if self.ignore_downstream {
                    continue;
                }
                (b.start - a.end + 1) as i64
            };

            let abs_distance = distance.abs();

            // Check max distance
            if let Some(max_d) = self.max_distance {
                if abs_distance > max_d as i64 {
                    continue;
                }
            }

            // Update closest
            match abs_distance.cmp(&min_distance.abs()) {
                Ordering::Less => {
                    min_distance = distance;
                    candidates.clear();
                    candidates.push(((*b).clone(), distance));
                }
                Ordering::Equal => {
                    candidates.push(((*b).clone(), distance));
                }
                Ordering::Greater => {}
            }
        }

        // Handle ties
        match self.tie_handling {
            TieHandling::All => candidates,
            TieHandling::First => candidates.into_iter().take(1).collect(),
            TieHandling::Last => candidates.into_iter().last().into_iter().collect(),
        }
    }

    /// Find closest in parallel by chromosome.
    pub fn find_closest_parallel(
        &self,
        a_intervals: Vec<Interval>,
        b_intervals: Vec<Interval>,
    ) -> Vec<ClosestResult> {
        let a_groups = group_by_chromosome(a_intervals);
        let b_groups = group_by_chromosome(b_intervals);

        let results: Vec<Vec<ClosestResult>> = a_groups
            .into_par_iter()
            .map(|(chrom, a_ints)| {
                let b_ints = b_groups.get(&chrom).map(|v| v.as_slice()).unwrap_or(&[]);
                self.find_closest(&a_ints, b_ints)
            })
            .collect();

        results.into_iter().flatten().collect()
    }

    /// Execute closest command on files using optimized algorithm.
    ///
    /// Algorithm complexity: O(n log m) per chromosome where n = |A|, m = |B|
    /// - Binary search for downstream: O(log m)
    /// - Binary search for upstream: O(log m)
    /// - Limited scan for overlaps: O(k) where k = number of overlapping intervals
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

        // Get sorted list of chromosomes
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
                    self.closest_chromosome_optimized(a_list, b_list, &mut buf);
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
                        self.closest_chromosome_optimized(a_list, b_list, &mut buf);
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

    /// Optimized closest search for a single chromosome.
    ///
    /// Algorithm:
    /// 1. Sort A and B by start position
    /// 2. For each A:
    ///    a. Binary search for downstream B: O(log m)
    ///    b. Binary search for upstream B using secondary end-sorted index: O(log m)
    ///    c. Scan limited window for overlaps: O(k) where k = overlapping count
    fn closest_chromosome_optimized(
        &self,
        a_sorted: &[BedRecord],
        b_sorted: Option<&Vec<BedRecord>>,
        output: &mut Vec<u8>,
    ) {
        let b_sorted = match b_sorted {
            Some(b) if !b.is_empty() => b,
            _ => {
                // No B intervals - output "no closest" for all A
                for a_rec in a_sorted {
                    self.write_no_closest_to_buf(output, a_rec);
                }
                return;
            }
        };

        let b_len = b_sorted.len();

        // Create index of B sorted by end position (for upstream queries)
        let b_by_end: Vec<usize> = {
            let mut idx: Vec<usize> = (0..b_len).collect();
            idx.sort_unstable_by_key(|&i| b_sorted[i].end());
            idx
        };

        // Precompute: for each index in b_sorted, its position in b_by_end
        // This allows O(1) check if a B interval has end > threshold
        let _b_end_rank: Vec<usize> = {
            let mut rank = vec![0usize; b_len];
            for (pos, &idx) in b_by_end.iter().enumerate() {
                rank[idx] = pos;
            }
            rank
        };

        // Precompute: max B.end for each prefix of b_sorted
        // max_end_prefix[i] = max(b_sorted[0..=i].end())
        // This allows O(1) check: is there any B in [0..k] with end > threshold?
        let max_end_prefix: Vec<u64> = {
            let mut v = Vec::with_capacity(b_len);
            let mut max_end = 0u64;
            for b in b_sorted.iter() {
                max_end = max_end.max(b.end());
                v.push(max_end);
            }
            v
        };

        // Sweep-line pointer for upstream queries
        let mut end_ptr: usize = 0;

        // Reusable candidates buffer
        let mut candidates: Vec<(&BedRecord, i64)> = Vec::with_capacity(16);

        for a_rec in a_sorted {
            let a_start = a_rec.start();
            let a_end = a_rec.end();

            candidates.clear();
            let mut min_dist = i64::MAX;

            // ========== Step 1: Find closest downstream B ==========
            // Downstream: B.start >= A.end
            // Binary search for first B where B.start >= A.end
            let ds_start = b_sorted.partition_point(|b| b.start() < a_end);

            if !self.ignore_downstream && ds_start < b_len {
                let b_rec = &b_sorted[ds_start];
                let dist = (b_rec.start() - a_end + 1) as i64;

                // Check max distance
                if self.max_distance.is_none_or(|max_d| dist <= max_d as i64) {
                    min_dist = dist;
                    candidates.push((b_rec, dist));

                    // Check for ties (consecutive B at same distance)
                    for b_rec2 in &b_sorted[ds_start + 1..] {
                        let dist2 = (b_rec2.start() - a_end + 1) as i64;
                        if dist2 > dist {
                            break;
                        }
                        candidates.push((b_rec2, dist2));
                    }
                }
            }

            // ========== Step 2: Find closest upstream B ==========
            // Upstream: B.end <= A.start
            // Advance sweep pointer
            while end_ptr < b_len && b_sorted[b_by_end[end_ptr]].end() <= a_start {
                end_ptr += 1;
            }

            // Best upstream is at end_ptr - 1
            if !self.ignore_upstream && end_ptr > 0 {
                let b_idx = b_by_end[end_ptr - 1];
                let b_rec = &b_sorted[b_idx];
                let dist = -((a_start - b_rec.end() + 1) as i64);
                let abs_dist = dist.abs();

                if self
                    .max_distance
                    .is_none_or(|max_d| abs_dist <= max_d as i64)
                {
                    if abs_dist < min_dist.abs() {
                        min_dist = dist;
                        candidates.clear();
                        candidates.push((b_rec, dist));
                    } else if abs_dist == min_dist.abs() {
                        candidates.push((b_rec, dist));
                    }

                    // Check for ties (other B at same end position)
                    let target_end = b_rec.end();
                    let mut tie_ptr = end_ptr - 1;
                    while tie_ptr > 0 {
                        tie_ptr -= 1;
                        let tie_idx = b_by_end[tie_ptr];
                        let tie_rec = &b_sorted[tie_idx];
                        if tie_rec.end() < target_end {
                            break;
                        }
                        if abs_dist == min_dist.abs() {
                            candidates.push((tie_rec, dist));
                        }
                    }
                }
            }

            // ========== Step 3: Find overlapping B intervals ==========
            // Overlap: B.start < A.end AND B.end > A.start
            // OPTIMIZED: Iterate over b_by_end[end_ptr..] which have end > A.start
            // Then filter for start < A.end (index in b_sorted < ds_start)
            if !self.ignore_overlaps && ds_start > 0 {
                // Quick check: is there ANY B in [0..ds_start) with end > A.start?
                // If max_end_prefix[ds_start-1] <= A.start, no overlaps exist
                if max_end_prefix[ds_start - 1] > a_start {
                    // Iterate only over B intervals with end > A.start
                    // These are at positions b_by_end[end_ptr..]
                    let mut found_overlaps = false;

                    for &idx in &b_by_end[end_ptr..] {
                        // idx is index in b_sorted
                        // Check if B.start < A.end (i.e., idx < ds_start)
                        if idx < ds_start {
                            let b_rec = &b_sorted[idx];
                            // This B overlaps A
                            if !found_overlaps {
                                if min_dist != 0 {
                                    candidates.clear();
                                }
                                found_overlaps = true;
                            }
                            candidates.push((b_rec, 0));
                            // Early termination for First tie handling
                            if self.tie_handling == TieHandling::First {
                                break;
                            }
                        }
                    }
                }
            }

            // ========== Apply tie handling ==========
            self.apply_tie_handling(&mut candidates);

            // ========== Output results ==========
            self.output_closest_to_buf(output, a_rec, &candidates);
        }
    }

    /// Apply tie handling policy to candidates.
    #[inline]
    fn apply_tie_handling(&self, candidates: &mut Vec<(&BedRecord, i64)>) {
        if candidates.is_empty() {
            return;
        }

        match self.tie_handling {
            TieHandling::All => {} // Keep all
            TieHandling::First => {
                candidates.truncate(1);
            }
            TieHandling::Last => {
                if let Some(last) = candidates.pop() {
                    candidates.clear();
                    candidates.push(last);
                }
            }
        }
    }

    /// Output closest results to buffer.
    #[inline]
    fn output_closest_to_buf(
        &self,
        output: &mut Vec<u8>,
        a_rec: &BedRecord,
        candidates: &[(&BedRecord, i64)],
    ) {
        if candidates.is_empty() {
            self.write_no_closest_to_buf(output, a_rec);
        } else {
            for (b_rec, distance) in candidates {
                self.write_closest_pair_to_buf(output, a_rec, b_rec, *distance);
            }
        }
    }

    /// Write "no closest found" output.
    #[inline]
    fn write_no_closest_to_buf(&self, buf: &mut Vec<u8>, a_rec: &BedRecord) {
        use std::io::Write as IoWrite;
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
        if self.report_distance {
            let _ = write!(buf, "\t.\t-1\t-1\t-1");
        } else {
            let _ = write!(buf, "\t.\t-1\t-1");
        }
        buf.push(b'\n');
    }

    /// Write a closest pair to buffer.
    #[inline]
    fn write_closest_pair_to_buf(
        &self,
        buf: &mut Vec<u8>,
        a_rec: &BedRecord,
        b_rec: &BedRecord,
        distance: i64,
    ) {
        use std::io::Write as IoWrite;

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

        // Distance (if requested)
        if self.report_distance {
            let _ = write!(buf, "\t{}", distance.abs());
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
        for list in map.values_mut() {
            list.sort_unstable_by(|a, b| a.start().cmp(&b.start()).then(a.end().cmp(&b.end())));
        }
        map
    }
}

/// Result of a closest query.
#[derive(Debug, Clone)]
pub struct ClosestResult {
    pub a_interval: Interval,
    /// Closest intervals with their signed distances (negative = upstream)
    pub closest_intervals: Vec<(Interval, i64)>,
}

impl ClosestResult {
    /// Get the minimum distance.
    pub fn min_distance(&self) -> Option<i64> {
        self.closest_intervals.first().map(|(_, d)| *d)
    }

    /// Check if there are any closest intervals.
    pub fn has_closest(&self) -> bool {
        !self.closest_intervals.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_closest() {
        let cmd = ClosestCommand::new();
        let a = vec![Interval::new("chr1", 100, 200)];
        let b = vec![
            Interval::new("chr1", 300, 400),
            Interval::new("chr1", 500, 600),
        ];

        let results = cmd.find_closest(&a, &b);

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].closest_intervals.len(), 1);
        assert_eq!(results[0].closest_intervals[0].0.start, 300);
        assert_eq!(results[0].closest_intervals[0].1, 101); // Distance: 300 - 200 + 1 = 101
    }

    #[test]
    fn test_closest_with_ties() {
        let cmd = ClosestCommand::new();
        let a = vec![Interval::new("chr1", 200, 300)];
        let b = vec![
            Interval::new("chr1", 100, 150), // 51bp upstream (300 - 150 + 1)
            Interval::new("chr1", 350, 400), // 51bp downstream (350 - 300 + 1)
        ];

        let results = cmd.find_closest(&a, &b);

        assert_eq!(results[0].closest_intervals.len(), 2);
    }

    #[test]
    fn test_closest_tie_handling_first() {
        let mut cmd = ClosestCommand::new();
        cmd.tie_handling = TieHandling::First;

        let a = vec![Interval::new("chr1", 200, 300)];
        let b = vec![
            Interval::new("chr1", 100, 150),
            Interval::new("chr1", 350, 400),
        ];

        let results = cmd.find_closest(&a, &b);

        assert_eq!(results[0].closest_intervals.len(), 1);
    }

    #[test]
    fn test_closest_ignore_overlaps() {
        let mut cmd = ClosestCommand::new();
        cmd.ignore_overlaps = true;

        let a = vec![Interval::new("chr1", 100, 200)];
        let b = vec![
            Interval::new("chr1", 150, 250), // Overlaps
            Interval::new("chr1", 300, 400), // Doesn't overlap
        ];

        let results = cmd.find_closest(&a, &b);

        assert_eq!(results[0].closest_intervals.len(), 1);
        assert_eq!(results[0].closest_intervals[0].0.start, 300);
    }

    #[test]
    fn test_closest_ignore_upstream() {
        let mut cmd = ClosestCommand::new();
        cmd.ignore_upstream = true;

        let a = vec![Interval::new("chr1", 200, 300)];
        let b = vec![
            Interval::new("chr1", 50, 100),  // Upstream
            Interval::new("chr1", 400, 500), // Downstream
        ];

        let results = cmd.find_closest(&a, &b);

        assert_eq!(results[0].closest_intervals.len(), 1);
        assert_eq!(results[0].closest_intervals[0].0.start, 400);
    }

    #[test]
    fn test_closest_with_max_distance() {
        let mut cmd = ClosestCommand::new();
        cmd.max_distance = Some(50);

        let a = vec![Interval::new("chr1", 100, 200)];
        let b = vec![Interval::new("chr1", 300, 400)]; // 100bp away

        let results = cmd.find_closest(&a, &b);

        assert!(results[0].closest_intervals.is_empty());
    }

    #[test]
    fn test_closest_overlapping() {
        let cmd = ClosestCommand::new();
        let a = vec![Interval::new("chr1", 100, 200)];
        let b = vec![Interval::new("chr1", 150, 250)];

        let results = cmd.find_closest(&a, &b);

        assert_eq!(results[0].closest_intervals.len(), 1);
        assert_eq!(results[0].closest_intervals[0].1, 0); // Distance is 0
    }

    #[test]
    fn test_parallel_closest() {
        let cmd = ClosestCommand::new();
        let a = vec![
            Interval::new("chr1", 100, 200),
            Interval::new("chr2", 100, 200),
        ];
        let b = vec![
            Interval::new("chr1", 300, 400),
            Interval::new("chr2", 300, 400),
        ];

        let results = cmd.find_closest_parallel(a, b);

        assert_eq!(results.len(), 2);
    }
}
