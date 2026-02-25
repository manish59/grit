//! Parallel processing utilities using Rayon.

use crate::bed::Result as BedResult;

/// Minimum number of intervals before enabling parallelization.
/// Below this threshold, sequential processing is faster due to
/// thread spawn overhead.
pub const PARALLEL_THRESHOLD: usize = 10_000;
use crate::interval::{BedRecord, Interval};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Group intervals by chromosome for parallel processing.
pub fn group_by_chromosome(intervals: Vec<Interval>) -> HashMap<String, Vec<Interval>> {
    let mut groups: HashMap<String, Vec<Interval>> = HashMap::new();

    for interval in intervals {
        groups
            .entry(interval.chrom.clone())
            .or_default()
            .push(interval);
    }

    groups
}

/// Group BED records by chromosome.
pub fn group_records_by_chromosome(records: Vec<BedRecord>) -> HashMap<String, Vec<BedRecord>> {
    let mut groups: HashMap<String, Vec<BedRecord>> = HashMap::new();

    for record in records {
        groups
            .entry(record.chrom().to_string())
            .or_default()
            .push(record);
    }

    groups
}

/// Process chromosomes in parallel.
pub fn process_chromosomes<F, T>(groups: HashMap<String, Vec<Interval>>, f: F) -> HashMap<String, T>
where
    F: Fn(&str, Vec<Interval>) -> T + Sync + Send,
    T: Send,
{
    groups
        .into_par_iter()
        .map(|(chrom, intervals)| {
            let result = f(&chrom, intervals);
            (chrom, result)
        })
        .collect()
}

/// Parallel sort of intervals.
pub fn parallel_sort(mut intervals: Vec<Interval>) -> Vec<Interval> {
    intervals.par_sort_unstable_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });
    intervals
}

/// Parallel sort of BED records.
pub fn parallel_sort_records(mut records: Vec<BedRecord>) -> Vec<BedRecord> {
    records.par_sort_unstable_by(|a, b| {
        a.chrom()
            .cmp(b.chrom())
            .then(a.start().cmp(&b.start()))
            .then(a.end().cmp(&b.end()))
    });
    records
}

/// Read a BED file in parallel chunks.
/// Note: This reads the entire file into memory first, then processes in parallel.
pub fn parallel_read_bed<P: AsRef<Path>>(path: P) -> BedResult<Vec<Interval>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read all lines (skip lines that fail to read - intentional behavior)
    #[allow(clippy::lines_filter_map_ok)]
    let lines: Vec<String> = reader
        .lines()
        .filter_map(|l| l.ok())
        .filter(|l| {
            let trimmed = l.trim();
            !trimmed.is_empty()
                && !trimmed.starts_with('#')
                && !trimmed.starts_with("track")
                && !trimmed.starts_with("browser")
        })
        .collect();

    // Parse in parallel
    let intervals: Vec<Option<Interval>> = lines
        .par_iter()
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 3 {
                let chrom = fields[0].to_string();
                let start: u64 = fields[1].parse().ok()?;
                let end: u64 = fields[2].parse().ok()?;
                Some(Interval::new(chrom, start, end))
            } else {
                None
            }
        })
        .collect();

    Ok(intervals.into_iter().flatten().collect())
}

/// Statistics for parallel work distribution.
#[derive(Debug, Clone)]
pub struct ParallelStats {
    pub total_intervals: usize,
    pub num_chromosomes: usize,
    pub intervals_per_chrom: Vec<(String, usize)>,
}

impl ParallelStats {
    pub fn from_groups(groups: &HashMap<String, Vec<Interval>>) -> Self {
        let mut intervals_per_chrom: Vec<(String, usize)> = groups
            .iter()
            .map(|(chrom, intervals)| (chrom.clone(), intervals.len()))
            .collect();
        intervals_per_chrom.sort_by(|a, b| b.1.cmp(&a.1));

        Self {
            total_intervals: groups.values().map(|v| v.len()).sum(),
            num_chromosomes: groups.len(),
            intervals_per_chrom,
        }
    }
}

/// A trait for operations that can be parallelized by chromosome.
pub trait ChromosomeParallel {
    type Output: Send;

    fn process_chromosome(&self, chrom: &str, intervals: &[Interval]) -> Self::Output;

    fn merge_results(&self, results: Vec<(String, Self::Output)>) -> Self::Output;
}

/// Execute a chromosome-parallel operation.
pub fn execute_parallel<T: ChromosomeParallel + Sync>(
    operation: &T,
    groups: &HashMap<String, Vec<Interval>>,
) -> T::Output {
    let results: Vec<(String, T::Output)> = groups
        .par_iter()
        .map(|(chrom, intervals)| {
            (
                chrom.clone(),
                operation.process_chromosome(chrom, intervals),
            )
        })
        .collect();

    operation.merge_results(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_group_by_chromosome() {
        let intervals = vec![
            Interval::new("chr1", 100, 200),
            Interval::new("chr2", 100, 200),
            Interval::new("chr1", 300, 400),
        ];

        let groups = group_by_chromosome(intervals);

        assert_eq!(groups.len(), 2);
        assert_eq!(groups.get("chr1").unwrap().len(), 2);
        assert_eq!(groups.get("chr2").unwrap().len(), 1);
    }

    #[test]
    fn test_parallel_sort() {
        let intervals = vec![
            Interval::new("chr2", 100, 200),
            Interval::new("chr1", 300, 400),
            Interval::new("chr1", 100, 200),
        ];

        let sorted = parallel_sort(intervals);

        assert_eq!(sorted[0].chrom, "chr1");
        assert_eq!(sorted[0].start, 100);
        assert_eq!(sorted[1].chrom, "chr1");
        assert_eq!(sorted[1].start, 300);
        assert_eq!(sorted[2].chrom, "chr2");
    }
}
