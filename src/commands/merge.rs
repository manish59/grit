//! Merge command implementation.
//!
//! Uses O(n log n) sort + O(n) single-pass sweep-line merge.

use crate::bed::{read_records, BedError, BedReader};
use crate::interval::{BedRecord, Interval};
use crate::parallel::{group_by_chromosome, parallel_sort_records, PARALLEL_THRESHOLD};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{self, Write};
use std::path::Path;

/// Operations for merging column values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MergeOperation {
    Sum,
    Min,
    Max,
    Mean,
    Median,
    Count,
    CountDistinct,
    Collapse,
    Distinct,
    First,
    Last,
}

impl MergeOperation {
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "sum" => Some(MergeOperation::Sum),
            "min" => Some(MergeOperation::Min),
            "max" => Some(MergeOperation::Max),
            "mean" => Some(MergeOperation::Mean),
            "median" => Some(MergeOperation::Median),
            "count" => Some(MergeOperation::Count),
            "count_distinct" => Some(MergeOperation::CountDistinct),
            "collapse" => Some(MergeOperation::Collapse),
            "distinct" => Some(MergeOperation::Distinct),
            "first" => Some(MergeOperation::First),
            "last" => Some(MergeOperation::Last),
            _ => None,
        }
    }
}

/// Merge command configuration.
#[derive(Debug, Clone)]
pub struct MergeCommand {
    /// Maximum distance between intervals to merge (default: 0)
    pub distance: u64,
    /// Require strand to match for merging
    pub strand_specific: bool,
    /// Column indices to operate on (1-based)
    pub columns: Vec<usize>,
    /// Operations to perform on columns
    pub operations: Vec<MergeOperation>,
    /// Delimiter for collapsed values
    pub delimiter: String,
}

impl Default for MergeCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl MergeCommand {
    pub fn new() -> Self {
        Self {
            distance: 0,
            strand_specific: false,
            columns: Vec::new(),
            operations: Vec::new(),
            delimiter: ",".to_string(),
        }
    }

    /// Set the maximum merge distance.
    pub fn with_distance(mut self, d: u64) -> Self {
        self.distance = d;
        self
    }

    /// Set strand-specific merging.
    pub fn with_strand(mut self, strand: bool) -> Self {
        self.strand_specific = strand;
        self
    }

    /// Merge intervals, returning merged intervals.
    pub fn merge(&self, intervals: Vec<Interval>) -> Vec<Interval> {
        if intervals.is_empty() {
            return Vec::new();
        }

        // Group by chromosome
        let groups = group_by_chromosome(intervals);

        // Merge each chromosome in parallel
        let merged: Vec<Vec<Interval>> = groups
            .into_par_iter()
            .map(|(_, mut intervals)| {
                intervals.sort_by(|a, b| a.start.cmp(&b.start).then(a.end.cmp(&b.end)));
                self.merge_sorted(&intervals)
            })
            .collect();

        // Flatten and sort results
        let mut result: Vec<Interval> = merged.into_iter().flatten().collect();
        result.sort();
        result
    }

    /// Merge already-sorted intervals for a single chromosome.
    fn merge_sorted(&self, intervals: &[Interval]) -> Vec<Interval> {
        if intervals.is_empty() {
            return Vec::new();
        }

        let mut result = Vec::new();
        let mut current = intervals[0].clone();

        for interval in intervals.iter().skip(1) {
            // Check if intervals should be merged
            if self.should_merge(&current, interval) {
                // Extend current interval
                current.end = current.end.max(interval.end);
            } else {
                result.push(current);
                current = interval.clone();
            }
        }

        result.push(current);
        result
    }

    /// Check if two intervals should be merged.
    #[inline]
    fn should_merge(&self, a: &Interval, b: &Interval) -> bool {
        a.chrom == b.chrom && b.start <= a.end + self.distance
    }

    /// Merge BED records, preserving and combining extra fields.
    pub fn merge_records(&self, records: Vec<BedRecord>) -> Vec<MergedRecord> {
        if records.is_empty() {
            return Vec::new();
        }

        // Sort records by chrom, strand (if strand-specific), start, end
        let mut records = if self.strand_specific {
            let mut recs = records;
            recs.sort_by(|a, b| {
                a.chrom()
                    .cmp(b.chrom())
                    .then_with(|| {
                        let a_strand = a.strand.map(|s| s.to_string()).unwrap_or_default();
                        let b_strand = b.strand.map(|s| s.to_string()).unwrap_or_default();
                        a_strand.cmp(&b_strand)
                    })
                    .then(a.start().cmp(&b.start()))
                    .then(a.end().cmp(&b.end()))
            });
            recs
        } else {
            parallel_sort_records(records)
        };

        // Group and merge
        let mut result = Vec::new();
        let first_record = records.remove(0);
        let mut current_chrom = first_record.chrom().to_string();
        let mut current_strand = first_record.strand;
        let mut current_end = first_record.end();
        let mut current_group: Vec<BedRecord> = vec![first_record];

        for record in records {
            let rec_chrom = record.chrom();
            let rec_start = record.start();
            let rec_end = record.end();
            let rec_strand = record.strand;

            // Check if this record should be merged with current group
            let same_chrom = rec_chrom == current_chrom;
            let same_strand = !self.strand_specific || rec_strand == current_strand;
            let overlaps = rec_start <= current_end + self.distance;

            if same_chrom && same_strand && overlaps {
                // Extend the current group's end position
                current_end = current_end.max(rec_end);
                current_group.push(record);
            } else {
                // Emit merged record
                result.push(self.emit_merged_record(&current_group));
                current_chrom = rec_chrom.to_string();
                current_strand = rec_strand;
                current_end = rec_end;
                current_group = vec![record];
            }
        }

        if !current_group.is_empty() {
            result.push(self.emit_merged_record(&current_group));
        }

        result
    }

    /// Create a merged record from a group of overlapping records.
    fn emit_merged_record(&self, group: &[BedRecord]) -> MergedRecord {
        let chrom = group[0].chrom().to_string();
        let start = group.iter().map(|r| r.start()).min().unwrap();
        let end = group.iter().map(|r| r.end()).max().unwrap();

        MergedRecord {
            interval: Interval::new(chrom, start, end),
            count: group.len(),
            names: group.iter().filter_map(|r| r.name.clone()).collect(),
            scores: group.iter().filter_map(|r| r.score).collect(),
        }
    }

    /// Execute merge command on a file using optimized sweep-line algorithm.
    ///
    /// Algorithm:
    /// 1. Read all records
    /// 2. Group by chromosome
    /// 3. Sort each chromosome's records in parallel
    /// 4. Single-pass merge per chromosome with direct buffer output
    /// 5. Write results in chromosome order
    pub fn run<P: AsRef<Path>, W: Write>(&self, input: P, output: &mut W) -> Result<(), BedError> {
        let records = read_records(input)?;

        if records.is_empty() {
            return Ok(());
        }

        // Group by chromosome
        let grouped = Self::group_records_by_chrom_owned(records, self.strand_specific);

        // Get sorted chromosome list for deterministic output
        let mut chroms: Vec<_> = grouped.keys().cloned().collect();
        chroms.sort();

        // Calculate total intervals for threshold check
        let total: usize = grouped.values().map(|v| v.len()).sum();

        if total < PARALLEL_THRESHOLD {
            // Sequential processing for small datasets
            for chrom in &chroms {
                let mut buf = Vec::with_capacity(32 * 1024);
                if let Some(records) = grouped.get(chrom) {
                    self.merge_chromosome_sweepline(records, &mut buf);
                }
                output.write_all(&buf).map_err(BedError::Io)?;
            }
        } else {
            // Parallel processing for large datasets
            let results: Vec<Vec<u8>> = chroms
                .par_iter()
                .map(|chrom| {
                    let mut buf = Vec::with_capacity(32 * 1024);
                    if let Some(records) = grouped.get(chrom) {
                        self.merge_chromosome_sweepline(records, &mut buf);
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

    /// O(n) single-pass sweep-line merge for a single chromosome.
    ///
    /// Records must be pre-sorted by start position.
    /// Outputs directly to buffer to avoid allocations.
    #[inline]
    fn merge_chromosome_sweepline(&self, records: &[BedRecord], output: &mut Vec<u8>) {
        if records.is_empty() {
            return;
        }

        // Current merge span
        let mut current_chrom: &str = records[0].chrom();
        let mut current_start: u64 = records[0].start();
        let mut current_end: u64 = records[0].end();
        let mut current_strand = records[0].strand;

        for rec in &records[1..] {
            let rec_start = rec.start();
            let rec_end = rec.end();

            // Check if should merge with current span
            let should_merge = if self.strand_specific {
                rec.strand == current_strand && rec_start <= current_end + self.distance
            } else {
                rec_start <= current_end + self.distance
            };

            if should_merge {
                // Extend current span
                current_end = current_end.max(rec_end);
            } else {
                // Output current span
                self.write_interval_to_buf(output, current_chrom, current_start, current_end);

                // Start new span
                current_chrom = rec.chrom();
                current_start = rec_start;
                current_end = rec_end;
                current_strand = rec.strand;
            }
        }

        // Output final span
        self.write_interval_to_buf(output, current_chrom, current_start, current_end);
    }

    /// Write interval directly to buffer (zero allocation).
    #[inline]
    fn write_interval_to_buf(&self, buf: &mut Vec<u8>, chrom: &str, start: u64, end: u64) {
        use std::io::Write as IoWrite;
        let _ = writeln!(buf, "{}\t{}\t{}", chrom, start, end);
    }

    /// Group records by chromosome (and strand if strand-specific), returning sorted records.
    fn group_records_by_chrom_owned(
        records: Vec<BedRecord>,
        strand_specific: bool,
    ) -> HashMap<String, Vec<BedRecord>> {
        let mut map: HashMap<String, Vec<BedRecord>> = HashMap::new();

        if strand_specific {
            // Group by chrom+strand
            for rec in records {
                let strand_char = rec.strand.map(|s| s.to_string()).unwrap_or_default();
                let key = format!("{}_{}", rec.chrom(), strand_char);
                map.entry(key).or_default().push(rec);
            }
        } else {
            // Group by chrom only
            for rec in records {
                map.entry(rec.chrom().to_string()).or_default().push(rec);
            }
        }

        // Sort each group by start position
        map.par_iter_mut().for_each(|(_, records)| {
            records.sort_unstable_by(|a, b| a.start().cmp(&b.start()).then(a.end().cmp(&b.end())));
        });

        map
    }

    /// Streaming merge for sorted input.
    pub fn merge_streaming<R: io::BufRead, W: Write>(
        &self,
        reader: BedReader<R>,
        output: &mut W,
    ) -> Result<(), BedError> {
        let mut current: Option<Interval> = None;

        for result in reader.records() {
            let record = result?;
            let interval = record.interval;

            if let Some(ref mut curr) = current {
                if curr.chrom == interval.chrom && interval.start <= curr.end + self.distance {
                    curr.end = curr.end.max(interval.end);
                } else {
                    writeln!(output, "{}", curr).map_err(BedError::Io)?;
                    *curr = interval;
                }
            } else {
                current = Some(interval);
            }
        }

        if let Some(curr) = current {
            writeln!(output, "{}", curr).map_err(BedError::Io)?;
        }

        Ok(())
    }
}

/// A merged record with aggregated information.
#[derive(Debug, Clone)]
pub struct MergedRecord {
    pub interval: Interval,
    pub count: usize,
    pub names: Vec<String>,
    pub scores: Vec<f64>,
}

impl MergedRecord {
    pub fn mean_score(&self) -> Option<f64> {
        if self.scores.is_empty() {
            return None;
        }
        Some(self.scores.iter().sum::<f64>() / self.scores.len() as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_merge() {
        let cmd = MergeCommand::new();
        let intervals = vec![
            Interval::new("chr1", 100, 200),
            Interval::new("chr1", 150, 250),
            Interval::new("chr1", 300, 400),
        ];

        let merged = cmd.merge(intervals);

        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0].start, 100);
        assert_eq!(merged[0].end, 250);
        assert_eq!(merged[1].start, 300);
        assert_eq!(merged[1].end, 400);
    }

    #[test]
    fn test_merge_with_distance() {
        let cmd = MergeCommand::new().with_distance(50);
        let intervals = vec![
            Interval::new("chr1", 100, 200),
            Interval::new("chr1", 250, 350), // Gap of 50, should merge
        ];

        let merged = cmd.merge(intervals);

        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start, 100);
        assert_eq!(merged[0].end, 350);
    }

    #[test]
    fn test_merge_different_chroms() {
        let cmd = MergeCommand::new();
        let intervals = vec![
            Interval::new("chr1", 100, 200),
            Interval::new("chr2", 100, 200),
        ];

        let merged = cmd.merge(intervals);

        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn test_merge_adjacent() {
        let cmd = MergeCommand::new();
        let intervals = vec![
            Interval::new("chr1", 100, 200),
            Interval::new("chr1", 200, 300), // Adjacent
        ];

        let merged = cmd.merge(intervals);

        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start, 100);
        assert_eq!(merged[0].end, 300);
    }

    #[test]
    fn test_merge_contained() {
        let cmd = MergeCommand::new();
        let intervals = vec![
            Interval::new("chr1", 100, 400),
            Interval::new("chr1", 150, 250), // Contained within first
        ];

        let merged = cmd.merge(intervals);

        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start, 100);
        assert_eq!(merged[0].end, 400);
    }

    #[test]
    fn test_merge_unsorted_input() {
        let cmd = MergeCommand::new();
        let intervals = vec![
            Interval::new("chr1", 300, 400),
            Interval::new("chr1", 100, 200),
            Interval::new("chr1", 150, 250),
        ];

        let merged = cmd.merge(intervals);

        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0].start, 100);
        assert_eq!(merged[0].end, 250);
    }
}
