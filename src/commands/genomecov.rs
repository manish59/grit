//! Genomecov command implementation.
//!
//! Computes genome-wide coverage using event-based sweep-line algorithm.
//! O(n log n) for sorting events, O(n) for sweep.

use crate::bed::{BedError, BedReader};
use crate::genome::Genome;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

/// Format fraction like bedtools (uses %g style formatting with 6 significant digits).
fn format_fraction(f: f64) -> String {
    if f == 0.0 {
        "0".to_string()
    } else if f == 1.0 {
        "1".to_string()
    } else {
        // Implement %g-style formatting: 6 significant figures, scientific if exponent < -4
        let abs = f.abs();
        let exp = abs.log10().floor() as i32;

        if (-4..6).contains(&exp) {
            // Fixed point notation with up to 6 significant figures
            let precision = (5 - exp).max(0) as usize;
            let s = format!("{:.prec$}", f, prec = precision);
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            // Scientific notation with 6 significant figures
            let mantissa = f / 10_f64.powi(exp);
            let mantissa_str = format!("{:.5}", mantissa);
            let mantissa_trimmed = mantissa_str.trim_end_matches('0').trim_end_matches('.');
            format!("{}e{:+03}", mantissa_trimmed, exp)
        }
    }
}

/// Genomecov output mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputMode {
    /// Default histogram output (chrom, depth, bases, chrom_size, fraction)
    Histogram,
    /// Per-base depth (-d): chrom, pos, depth
    PerBase,
    /// BedGraph (-bg): chrom, start, end, depth (non-zero only)
    BedGraph,
    /// BedGraph all (-bga): chrom, start, end, depth (including zero)
    BedGraphAll,
}

/// Genomecov command configuration.
#[derive(Debug, Clone)]
pub struct GenomecovCommand {
    /// Output mode
    pub mode: OutputMode,
    /// Scale factor for depth
    pub scale: f64,
    /// Maximum depth to report (for histogram)
    pub max_depth: Option<u32>,
    /// Report depth >= 1 (not zero coverage)
    pub report_zero: bool,
    /// Split by strand
    pub strand: bool,
    /// 5' end only
    pub five_prime: bool,
    /// 3' end only
    pub three_prime: bool,
}

impl Default for GenomecovCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl GenomecovCommand {
    pub fn new() -> Self {
        Self {
            mode: OutputMode::Histogram,
            scale: 1.0,
            max_depth: None,
            report_zero: true,
            strand: false,
            five_prime: false,
            three_prime: false,
        }
    }

    /// Process a chromosome's intervals using sweep-line algorithm.
    /// Returns sorted (start, end, depth) tuples with adjacent same-depth regions merged.
    fn sweep_chromosome(&self, intervals: &[(u64, u64)], chrom_size: u64) -> Vec<(u64, u64, u32)> {
        if intervals.is_empty() {
            return vec![(0, chrom_size, 0)];
        }

        // Create events: (position, delta) where delta is +1 for start, -1 for end
        let mut events: Vec<(u64, i32)> = Vec::with_capacity(intervals.len() * 2);
        for &(start, end) in intervals {
            events.push((start, 1));
            events.push((end, -1));
        }

        // Sort events by position, with starts before ends at same position
        events.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));

        // Sweep and collect depth regions
        let mut result: Vec<(u64, u64, u32)> = Vec::new();
        let mut depth: i32 = 0;
        let mut prev_pos: u64 = 0;

        for (pos, delta) in events {
            // Clamp to chromosome size
            let pos = pos.min(chrom_size);

            if pos > prev_pos {
                // Merge with previous region if same depth
                let cur_depth = depth as u32;
                if let Some(last) = result.last_mut() {
                    if last.2 == cur_depth && last.1 == prev_pos {
                        last.1 = pos;
                    } else {
                        result.push((prev_pos, pos, cur_depth));
                    }
                } else {
                    result.push((prev_pos, pos, cur_depth));
                }
            }

            depth += delta;
            prev_pos = pos;
        }

        // Handle trailing region
        if prev_pos < chrom_size {
            let trailing_depth = depth as u32;
            if let Some(last) = result.last_mut() {
                if last.2 == trailing_depth && last.1 == prev_pos {
                    last.1 = chrom_size;
                } else {
                    result.push((prev_pos, chrom_size, trailing_depth));
                }
            } else {
                result.push((prev_pos, chrom_size, trailing_depth));
            }
        }

        result
    }

    /// Build histogram from depth regions.
    fn build_histogram(&self, regions: &[(u64, u64, u32)]) -> HashMap<u32, u64> {
        let mut hist: HashMap<u32, u64> = HashMap::new();
        for &(start, end, depth) in regions {
            *hist.entry(depth).or_insert(0) += end - start;
        }
        hist
    }

    /// Run genomecov with streaming output.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        input: P,
        genome: &Genome,
        output: &mut W,
    ) -> Result<(), BedError> {
        let file = File::open(input)?;
        let reader = BedReader::new(file);
        self.genomecov_streaming(reader, genome, output)
    }

    /// Streaming genomecov implementation.
    pub fn genomecov_streaming<R: Read, W: Write>(
        &self,
        reader: BedReader<R>,
        genome: &Genome,
        output: &mut W,
    ) -> Result<(), BedError> {
        let mut buf_output = BufWriter::with_capacity(256 * 1024, output);

        // Group intervals by chromosome
        let mut by_chrom: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        for result in reader.records() {
            let record = result?;
            let chrom = record.chrom().to_string();

            // Skip chromosomes not in genome
            if !genome.has_chrom(&chrom) {
                continue;
            }

            by_chrom
                .entry(chrom)
                .or_default()
                .push((record.start(), record.end()));
        }

        // Process chromosomes in genome order
        let mut genome_hist: HashMap<u32, u64> = HashMap::new();
        let mut total_bases: u64 = 0;

        for chrom in genome.chromosomes() {
            let chrom_size = genome.chrom_size(chrom).unwrap();
            total_bases += chrom_size;

            let intervals = by_chrom.get(chrom).map(|v| v.as_slice()).unwrap_or(&[]);
            let regions = self.sweep_chromosome(intervals, chrom_size);

            match self.mode {
                OutputMode::Histogram => {
                    let chrom_hist = self.build_histogram(&regions);

                    // Output per-chromosome histogram
                    let mut depths: Vec<_> = chrom_hist.keys().copied().collect();
                    depths.sort_unstable();

                    for depth in depths {
                        let bases = chrom_hist[&depth];
                        if depth == 0 && !self.report_zero {
                            continue;
                        }
                        let fraction = bases as f64 / chrom_size as f64;
                        writeln!(
                            buf_output,
                            "{}\t{}\t{}\t{}\t{}",
                            chrom,
                            depth,
                            bases,
                            chrom_size,
                            format_fraction(fraction)
                        )
                        .map_err(BedError::Io)?;

                        // Accumulate for genome-wide
                        *genome_hist.entry(depth).or_insert(0) += bases;
                    }
                }

                OutputMode::BedGraph | OutputMode::BedGraphAll => {
                    for (start, end, depth) in &regions {
                        if self.mode == OutputMode::BedGraph && *depth == 0 {
                            continue;
                        }
                        let scaled_depth = (*depth as f64 * self.scale) as u32;
                        writeln!(
                            buf_output,
                            "{}\t{}\t{}\t{}",
                            chrom, start, end, scaled_depth
                        )
                        .map_err(BedError::Io)?;
                    }
                }

                OutputMode::PerBase => {
                    for (start, end, depth) in &regions {
                        let scaled_depth = (*depth as f64 * self.scale) as u32;
                        for pos in *start..*end {
                            // 1-based position for output
                            writeln!(buf_output, "{}\t{}\t{}", chrom, pos + 1, scaled_depth)
                                .map_err(BedError::Io)?;
                        }
                    }
                }
            }
        }

        // Output genome-wide histogram
        if self.mode == OutputMode::Histogram {
            let mut depths: Vec<_> = genome_hist.keys().copied().collect();
            depths.sort_unstable();

            for depth in depths {
                let bases = genome_hist[&depth];
                if depth == 0 && !self.report_zero {
                    continue;
                }
                let fraction = bases as f64 / total_bases as f64;
                writeln!(
                    buf_output,
                    "genome\t{}\t{}\t{}\t{}",
                    depth,
                    bases,
                    total_bases,
                    format_fraction(fraction)
                )
                .map_err(BedError::Io)?;
            }
        }

        buf_output.flush().map_err(BedError::Io)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sweep_single_interval() {
        let cmd = GenomecovCommand::new();
        let intervals = vec![(100, 200)];
        let regions = cmd.sweep_chromosome(&intervals, 1000);

        assert_eq!(regions.len(), 3);
        assert_eq!(regions[0], (0, 100, 0));
        assert_eq!(regions[1], (100, 200, 1));
        assert_eq!(regions[2], (200, 1000, 0));
    }

    #[test]
    fn test_sweep_overlapping() {
        let cmd = GenomecovCommand::new();
        let intervals = vec![(100, 200), (150, 250)];
        let regions = cmd.sweep_chromosome(&intervals, 1000);

        // Regions: [0-100)=0, [100-150)=1, [150-200)=2, [200-250)=1, [250-1000)=0
        assert_eq!(regions.len(), 5);
        assert_eq!(regions[0], (0, 100, 0));
        assert_eq!(regions[1], (100, 150, 1));
        assert_eq!(regions[2], (150, 200, 2));
        assert_eq!(regions[3], (200, 250, 1));
        assert_eq!(regions[4], (250, 1000, 0));
    }

    #[test]
    fn test_sweep_adjacent() {
        let cmd = GenomecovCommand::new();
        let intervals = vec![(100, 200), (200, 300)];
        let regions = cmd.sweep_chromosome(&intervals, 1000);

        // Adjacent intervals with same depth are merged
        assert_eq!(regions.len(), 3);
        assert_eq!(regions[0], (0, 100, 0));
        assert_eq!(regions[1], (100, 300, 1)); // Merged
        assert_eq!(regions[2], (300, 1000, 0));
    }

    #[test]
    fn test_histogram() {
        let cmd = GenomecovCommand::new();
        let regions = vec![(0, 100, 0), (100, 150, 1), (150, 200, 2), (200, 1000, 0)];
        let hist = cmd.build_histogram(&regions);

        assert_eq!(hist.get(&0), Some(&900));
        assert_eq!(hist.get(&1), Some(&50));
        assert_eq!(hist.get(&2), Some(&50));
    }

    #[test]
    fn test_empty_chromosome() {
        let cmd = GenomecovCommand::new();
        let intervals: Vec<(u64, u64)> = vec![];
        let regions = cmd.sweep_chromosome(&intervals, 1000);

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0], (0, 1000, 0));
    }
}
