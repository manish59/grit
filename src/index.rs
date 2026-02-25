//! Interval indexing for fast overlap queries.

use crate::interval::{BedRecord, Interval};
use std::collections::HashMap;

/// An indexed collection of intervals organized by chromosome.
/// Uses a sorted list with binary search for efficient queries.
pub struct IntervalIndex {
    intervals_by_chrom: HashMap<String, Vec<(Interval, usize)>>,
    intervals: Vec<Interval>,
}

impl IntervalIndex {
    /// Create a new empty index.
    pub fn new() -> Self {
        Self {
            intervals_by_chrom: HashMap::new(),
            intervals: Vec::new(),
        }
    }

    /// Build an index from a collection of intervals.
    pub fn from_intervals(intervals: Vec<Interval>) -> Self {
        let mut by_chrom: HashMap<String, Vec<(Interval, usize)>> = HashMap::new();

        for (idx, interval) in intervals.iter().enumerate() {
            by_chrom
                .entry(interval.chrom.clone())
                .or_default()
                .push((interval.clone(), idx));
        }

        // Sort each chromosome's intervals by start position
        for chrom_intervals in by_chrom.values_mut() {
            chrom_intervals.sort_by(|a, b| a.0.start.cmp(&b.0.start).then(a.0.end.cmp(&b.0.end)));
        }

        Self {
            intervals_by_chrom: by_chrom,
            intervals,
        }
    }

    /// Build an index from BED records.
    pub fn from_records(records: &[BedRecord]) -> Self {
        let intervals: Vec<Interval> = records.iter().map(|r| r.interval.clone()).collect();
        Self::from_intervals(intervals)
    }

    /// Find all intervals overlapping a query interval.
    pub fn find_overlaps(&self, query: &Interval) -> Vec<&Interval> {
        let mut results = Vec::new();

        if let Some(chrom_intervals) = self.intervals_by_chrom.get(&query.chrom) {
            // Binary search to find starting point
            let start_idx = chrom_intervals
                .binary_search_by(|(i, _)| {
                    if i.end <= query.start {
                        std::cmp::Ordering::Less
                    } else {
                        std::cmp::Ordering::Greater
                    }
                })
                .unwrap_or_else(|i| i);

            // Scan forward to find all overlaps
            for (interval, idx) in chrom_intervals.iter().skip(start_idx) {
                if interval.start >= query.end {
                    break;
                }
                if query.overlaps(interval) {
                    results.push(&self.intervals[*idx]);
                }
            }
        }

        results
    }

    /// Find all intervals overlapping a query, returning their indices.
    pub fn find_overlap_indices(&self, query: &Interval) -> Vec<usize> {
        let mut results = Vec::new();

        if let Some(chrom_intervals) = self.intervals_by_chrom.get(&query.chrom) {
            let start_idx = chrom_intervals
                .binary_search_by(|(i, _)| {
                    if i.end <= query.start {
                        std::cmp::Ordering::Less
                    } else {
                        std::cmp::Ordering::Greater
                    }
                })
                .unwrap_or_else(|i| i);

            for (interval, idx) in chrom_intervals.iter().skip(start_idx) {
                if interval.start >= query.end {
                    break;
                }
                if query.overlaps(interval) {
                    results.push(*idx);
                }
            }
        }

        results
    }

    /// Count overlapping intervals.
    pub fn count_overlaps(&self, query: &Interval) -> usize {
        self.find_overlap_indices(query).len()
    }

    /// Check if any interval overlaps the query.
    pub fn has_overlap(&self, query: &Interval) -> bool {
        if let Some(chrom_intervals) = self.intervals_by_chrom.get(&query.chrom) {
            let start_idx = chrom_intervals
                .binary_search_by(|(i, _)| {
                    if i.end <= query.start {
                        std::cmp::Ordering::Less
                    } else {
                        std::cmp::Ordering::Greater
                    }
                })
                .unwrap_or_else(|i| i);

            for (interval, _) in chrom_intervals.iter().skip(start_idx) {
                if interval.start >= query.end {
                    break;
                }
                if query.overlaps(interval) {
                    return true;
                }
            }
        }
        false
    }

    /// Get all chromosomes in the index.
    pub fn chromosomes(&self) -> impl Iterator<Item = &String> {
        self.intervals_by_chrom.keys()
    }

    /// Get all intervals.
    pub fn intervals(&self) -> &[Interval] {
        &self.intervals
    }

    /// Get an interval by index.
    pub fn get(&self, index: usize) -> Option<&Interval> {
        self.intervals.get(index)
    }

    /// Get the total number of intervals.
    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    /// Check if the index is empty.
    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }
}

impl Default for IntervalIndex {
    fn default() -> Self {
        Self::new()
    }
}

/// A simple interval index using a sorted list with binary search.
/// Alternative implementation for compatibility.
pub struct SimpleIndex {
    intervals_by_chrom: HashMap<String, Vec<Interval>>,
}

impl SimpleIndex {
    /// Create a new index from intervals.
    pub fn new(intervals: Vec<Interval>) -> Self {
        let mut by_chrom: HashMap<String, Vec<Interval>> = HashMap::new();

        for interval in intervals {
            by_chrom
                .entry(interval.chrom.clone())
                .or_default()
                .push(interval);
        }

        // Sort each chromosome's intervals
        for intervals in by_chrom.values_mut() {
            intervals.sort_by(|a, b| a.start.cmp(&b.start).then(a.end.cmp(&b.end)));
        }

        Self {
            intervals_by_chrom: by_chrom,
        }
    }

    /// Find overlapping intervals.
    pub fn find_overlaps(&self, query: &Interval) -> Vec<&Interval> {
        let mut results = Vec::new();

        if let Some(intervals) = self.intervals_by_chrom.get(&query.chrom) {
            let start_idx = intervals
                .binary_search_by(|i| {
                    if i.end <= query.start {
                        std::cmp::Ordering::Less
                    } else {
                        std::cmp::Ordering::Greater
                    }
                })
                .unwrap_or_else(|i| i);

            for interval in intervals.iter().skip(start_idx) {
                if interval.start >= query.end {
                    break;
                }
                if query.overlaps(interval) {
                    results.push(interval);
                }
            }
        }

        results
    }

    /// Get intervals for a specific chromosome.
    pub fn get_chrom(&self, chrom: &str) -> Option<&[Interval]> {
        self.intervals_by_chrom.get(chrom).map(|v| v.as_slice())
    }

    /// Get all chromosomes.
    pub fn chromosomes(&self) -> impl Iterator<Item = &String> {
        self.intervals_by_chrom.keys()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_intervals() -> Vec<Interval> {
        vec![
            Interval::new("chr1", 100, 200),
            Interval::new("chr1", 150, 250),
            Interval::new("chr1", 300, 400),
            Interval::new("chr2", 100, 200),
        ]
    }

    #[test]
    fn test_build_index() {
        let intervals = sample_intervals();
        let index = IntervalIndex::from_intervals(intervals);

        assert_eq!(index.len(), 4);
    }

    #[test]
    fn test_find_overlaps() {
        let intervals = sample_intervals();
        let index = IntervalIndex::from_intervals(intervals);

        let query = Interval::new("chr1", 175, 225);
        let overlaps = index.find_overlaps(&query);

        assert_eq!(overlaps.len(), 2);
    }

    #[test]
    fn test_count_overlaps() {
        let intervals = sample_intervals();
        let index = IntervalIndex::from_intervals(intervals);

        let query = Interval::new("chr1", 175, 225);
        assert_eq!(index.count_overlaps(&query), 2);
    }

    #[test]
    fn test_no_overlap() {
        let intervals = sample_intervals();
        let index = IntervalIndex::from_intervals(intervals);

        let query = Interval::new("chr1", 500, 600);
        assert_eq!(index.count_overlaps(&query), 0);
        assert!(!index.has_overlap(&query));
    }

    #[test]
    fn test_different_chrom() {
        let intervals = sample_intervals();
        let index = IntervalIndex::from_intervals(intervals);

        let query = Interval::new("chr3", 100, 200);
        assert_eq!(index.count_overlaps(&query), 0);
    }

    #[test]
    fn test_simple_index() {
        let intervals = sample_intervals();
        let index = SimpleIndex::new(intervals);

        let query = Interval::new("chr1", 175, 225);
        let overlaps = index.find_overlaps(&query);

        assert_eq!(overlaps.len(), 2);
    }
}
