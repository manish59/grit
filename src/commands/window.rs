//! Window command implementation - proximity-based matching.

use crate::bed::{read_records, BedError};
use crate::index::IntervalIndex;
use crate::interval::Interval;
use crate::parallel::group_by_chromosome;
use rayon::prelude::*;
use std::io::Write;
use std::path::Path;

/// Window command configuration.
#[derive(Debug, Clone)]
pub struct WindowCommand {
    /// Window size on both sides (symmetric)
    pub window: u64,
    /// Window size on the left (upstream)
    pub left: Option<u64>,
    /// Window size on the right (downstream)
    pub right: Option<u64>,
    /// Require same strand
    pub same_strand: bool,
    /// Require opposite strand
    pub opposite_strand: bool,
    /// Report original A entry
    pub write_a: bool,
    /// Report original B entry
    pub write_b: bool,
    /// Report number of overlaps
    pub count: bool,
    /// Only report A intervals with no matches
    pub no_overlap: bool,
    /// Process in parallel by chromosome
    pub parallel: bool,
}

impl Default for WindowCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl WindowCommand {
    pub fn new() -> Self {
        Self {
            window: 1000,
            left: None,
            right: None,
            same_strand: false,
            opposite_strand: false,
            write_a: true,
            write_b: true,
            count: false,
            no_overlap: false,
            parallel: true,
        }
    }

    /// Get the left window size.
    fn left_window(&self) -> u64 {
        self.left.unwrap_or(self.window)
    }

    /// Get the right window size.
    fn right_window(&self) -> u64 {
        self.right.unwrap_or(self.window)
    }

    /// Expand an interval by the window sizes.
    fn expand_interval(&self, interval: &Interval) -> Interval {
        Interval {
            chrom: interval.chrom.clone(),
            start: interval.start.saturating_sub(self.left_window()),
            end: interval.end.saturating_add(self.right_window()),
        }
    }

    /// Find all B intervals within the window of each A interval.
    pub fn find_window_matches(
        &self,
        a_intervals: &[Interval],
        b_index: &IntervalIndex,
    ) -> Vec<WindowResult> {
        let mut results = Vec::new();

        for a in a_intervals {
            // Expand A by window size
            let expanded = self.expand_interval(a);

            // Find overlaps with expanded window
            let matches: Vec<Interval> = b_index
                .find_overlaps(&expanded)
                .into_iter()
                .cloned()
                .collect();

            if self.no_overlap {
                if matches.is_empty() {
                    results.push(WindowResult {
                        a_interval: a.clone(),
                        b_intervals: Vec::new(),
                    });
                }
            } else if !matches.is_empty() {
                results.push(WindowResult {
                    a_interval: a.clone(),
                    b_intervals: matches,
                });
            }
        }

        results
    }

    /// Find window matches in parallel by chromosome.
    pub fn find_window_matches_parallel(
        &self,
        a_intervals: Vec<Interval>,
        b_intervals: Vec<Interval>,
    ) -> Vec<WindowResult> {
        let b_index = IntervalIndex::from_intervals(b_intervals);
        let a_groups = group_by_chromosome(a_intervals);

        let results: Vec<Vec<WindowResult>> = a_groups
            .into_par_iter()
            .map(|(_, intervals)| self.find_window_matches(&intervals, &b_index))
            .collect();

        results.into_iter().flatten().collect()
    }

    /// Execute window command on files.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        a_path: P,
        b_path: P,
        output: &mut W,
    ) -> Result<(), BedError> {
        let a_records = read_records(a_path)?;
        let b_records = read_records(b_path)?;

        let a_intervals: Vec<Interval> = a_records.iter().map(|r| r.interval.clone()).collect();
        let b_intervals: Vec<Interval> = b_records.iter().map(|r| r.interval.clone()).collect();

        let results = if self.parallel {
            self.find_window_matches_parallel(a_intervals, b_intervals)
        } else {
            let b_index = IntervalIndex::from_intervals(b_intervals);
            self.find_window_matches(&a_intervals, &b_index)
        };

        for result in results {
            if self.count {
                writeln!(
                    output,
                    "{}\t{}",
                    result.a_interval,
                    result.b_intervals.len()
                )
                .map_err(BedError::Io)?;
            } else if result.b_intervals.is_empty() {
                writeln!(output, "{}", result.a_interval).map_err(BedError::Io)?;
            } else {
                for b in &result.b_intervals {
                    if self.write_a && self.write_b {
                        writeln!(output, "{}\t{}", result.a_interval, b).map_err(BedError::Io)?;
                    } else if self.write_a {
                        writeln!(output, "{}", result.a_interval).map_err(BedError::Io)?;
                    } else {
                        writeln!(output, "{}", b).map_err(BedError::Io)?;
                    }
                }
            }
        }

        Ok(())
    }
}

/// Result of a window query.
#[derive(Debug, Clone)]
pub struct WindowResult {
    pub a_interval: Interval,
    pub b_intervals: Vec<Interval>,
}

impl WindowResult {
    /// Get the number of matches.
    pub fn count(&self) -> usize {
        self.b_intervals.len()
    }

    /// Check if there are any matches.
    pub fn has_matches(&self) -> bool {
        !self.b_intervals.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_window() {
        let mut cmd = WindowCommand::new();
        cmd.window = 100;

        // A: 500-600, with 100bp window -> expanded to 400-700
        let a = vec![Interval::new("chr1", 500, 600)];
        let b = vec![
            Interval::new("chr1", 350, 450), // Overlaps with expanded (400-700)
            Interval::new("chr1", 650, 750), // Overlaps with expanded (400-700)
            Interval::new("chr1", 100, 150), // Too far (doesn't overlap 400-700)
        ];
        let b_index = IntervalIndex::from_intervals(b);

        let results = cmd.find_window_matches(&a, &b_index);

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].b_intervals.len(), 2);
    }

    #[test]
    fn test_asymmetric_window() {
        let mut cmd = WindowCommand::new();
        cmd.left = Some(50);
        cmd.right = Some(200);

        // A: 500-600, with left=50, right=200 -> expanded to 450-800
        let a = vec![Interval::new("chr1", 500, 600)];
        let b = vec![
            Interval::new("chr1", 420, 470), // Overlaps with expanded (450-800)
            Interval::new("chr1", 750, 850), // Overlaps with expanded (450-800)
            Interval::new("chr1", 200, 300), // Outside window (doesn't overlap 450-800)
        ];
        let b_index = IntervalIndex::from_intervals(b);

        let results = cmd.find_window_matches(&a, &b_index);

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].b_intervals.len(), 2);
    }

    #[test]
    fn test_window_no_overlap_flag() {
        let mut cmd = WindowCommand::new();
        cmd.window = 100;
        cmd.no_overlap = true;

        let a = vec![
            Interval::new("chr1", 500, 600),   // Has matches
            Interval::new("chr1", 1000, 1100), // No matches
        ];
        let b = vec![Interval::new("chr1", 650, 700)];
        let b_index = IntervalIndex::from_intervals(b);

        let results = cmd.find_window_matches(&a, &b_index);

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].a_interval.start, 1000);
    }

    #[test]
    fn test_window_direct_overlap() {
        let mut cmd = WindowCommand::new();
        cmd.window = 100;

        let a = vec![Interval::new("chr1", 500, 600)];
        let b = vec![Interval::new("chr1", 550, 650)]; // Direct overlap
        let b_index = IntervalIndex::from_intervals(b);

        let results = cmd.find_window_matches(&a, &b_index);

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].b_intervals.len(), 1);
    }

    #[test]
    fn test_expand_interval() {
        let mut cmd = WindowCommand::new();
        cmd.left = Some(100);
        cmd.right = Some(200);

        let interval = Interval::new("chr1", 500, 600);
        let expanded = cmd.expand_interval(&interval);

        assert_eq!(expanded.start, 400);
        assert_eq!(expanded.end, 800);
    }

    #[test]
    fn test_expand_interval_no_underflow() {
        let mut cmd = WindowCommand::new();
        cmd.window = 1000;

        let interval = Interval::new("chr1", 500, 600);
        let expanded = cmd.expand_interval(&interval);

        assert_eq!(expanded.start, 0); // Clamped at 0
        assert_eq!(expanded.end, 1600);
    }

    #[test]
    fn test_parallel_window() {
        let mut cmd = WindowCommand::new();
        cmd.window = 100;

        let a = vec![
            Interval::new("chr1", 500, 600),
            Interval::new("chr2", 500, 600),
        ];
        let b = vec![
            Interval::new("chr1", 650, 700),
            Interval::new("chr2", 650, 700),
        ];

        let results = cmd.find_window_matches_parallel(a, b);

        assert_eq!(results.len(), 2);
    }
}
