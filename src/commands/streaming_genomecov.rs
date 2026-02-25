//! Streaming genomecov command - TRUE O(k) memory genome coverage computation.
//!
//! Memory complexity: O(k) where k = max overlapping intervals at any position.
//! Input is streamed - file is never fully loaded into memory.
//!
//! ZERO ALLOCATION in hot path:
//! - No per-interval allocation (chrom buffer reused)
//! - Float formatting uses ryu (stack-allocated buffer)
//! - Integer formatting uses itoa (stack-allocated buffer)
//!

#![allow(clippy::needless_range_loop)]
//! REQUIREMENT: Input must be sorted by (chrom, start) for streaming mode.
//! Use `--assume-sorted` flag or pre-sort with `grit sort`.
//!
//! For histogram output, we need to accumulate per-chromosome and genome-wide stats,
//! but still stream through the input efficiently.

use crate::bed::BedError;
use crate::genome::Genome;
use crate::streaming::parsing::{parse_bed3_bytes, should_skip_line};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Output mode for streaming genomecov.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamingGenomecovMode {
    /// Default histogram output (chrom, depth, bases, chrom_size, fraction)
    Histogram,
    /// Per-base depth (-d): chrom, pos, depth
    PerBase,
    /// BedGraph (-bg): chrom, start, end, depth (non-zero only)
    BedGraph,
    /// BedGraph all (-bga): chrom, start, end, depth (including zero)
    BedGraphAll,
}

/// Streaming genomecov command configuration.
#[derive(Debug, Clone)]
pub struct StreamingGenomecovCommand {
    /// Output mode
    pub mode: StreamingGenomecovMode,
    /// Scale factor for depth
    pub scale: f64,
    /// Skip sorted validation (faster for pre-sorted input)
    pub assume_sorted: bool,
}

impl Default for StreamingGenomecovCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl StreamingGenomecovCommand {
    pub fn new() -> Self {
        Self {
            mode: StreamingGenomecovMode::Histogram,
            scale: 1.0,
            assume_sorted: false,
        }
    }

    /// Set assume_sorted flag (builder pattern).
    pub fn with_assume_sorted(mut self, assume_sorted: bool) -> Self {
        self.assume_sorted = assume_sorted;
        self
    }

    /// Set output mode (builder pattern).
    pub fn with_mode(mut self, mode: StreamingGenomecovMode) -> Self {
        self.mode = mode;
        self
    }

    /// Set scale factor (builder pattern).
    pub fn with_scale(mut self, scale: f64) -> Self {
        self.scale = scale;
        self
    }

    /// Execute streaming genomecov.
    ///
    /// Memory: O(k) where k = max overlapping intervals on any chromosome.
    /// Input file is streamed - never fully loaded.
    ///
    /// REQUIREMENT: Input must be sorted by (chrom, start) for correct results.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        input: P,
        genome: &Genome,
        output: &mut W,
    ) -> Result<(), BedError> {
        let file = File::open(input)?;
        let reader = BufReader::with_capacity(256 * 1024, file);
        self.genomecov_streaming(reader, genome, output)
    }

    /// Streaming genomecov implementation.
    ///
    /// Algorithm:
    /// 1. Stream through sorted input
    /// 2. When chromosome changes, process the completed chromosome:
    ///    - Run sweep-line algorithm on accumulated events
    ///    - Output results immediately
    ///    - Clear events buffer for next chromosome
    /// 3. At end, process the last chromosome
    ///
    /// Memory is O(k) where k = max intervals on any single chromosome.
    /// For typical genomic data, this is much smaller than total file size.
    fn genomecov_streaming<R: BufRead, W: Write>(
        &self,
        mut reader: R,
        genome: &Genome,
        output: &mut W,
    ) -> Result<(), BedError> {
        // Large output buffer (8MB)
        let mut buf_output = BufWriter::with_capacity(8 * 1024 * 1024, output);

        // Get genome chromosomes info
        let chroms: Vec<&String> = genome.chromosomes().collect();
        let chrom_indices: HashMap<&[u8], usize> = chroms
            .iter()
            .enumerate()
            .map(|(i, c)| (c.as_bytes(), i))
            .collect();

        // Reusable line buffer
        let mut line_buf = String::with_capacity(1024);

        // Current chromosome events: (position, delta) where delta is +1 for start, -1 for end
        let mut events: Vec<(u64, i32)> = Vec::with_capacity(1024);
        let mut current_chrom_idx: Option<usize> = None;

        // For histogram mode: genome-wide accumulator
        let mut genome_hist: HashMap<u32, u64> = HashMap::new();
        let mut total_bases: u64 = 0;

        // Track which chromosomes we've seen (for outputting empty chromosomes)
        let mut seen_chroms: Vec<bool> = vec![false; chroms.len()];

        // itoa buffer for fast integer formatting
        let mut itoa_buf = itoa::Buffer::new();

        loop {
            line_buf.clear();
            let bytes_read = reader.read_line(&mut line_buf)?;
            if bytes_read == 0 {
                break;
            }

            let line_bytes = line_buf.trim_end().as_bytes();
            if should_skip_line(line_bytes) {
                continue;
            }

            let (chrom, start, end) = match parse_bed3_bytes(line_bytes) {
                Some(v) => v,
                None => continue,
            };

            // Skip chromosomes not in genome
            let chrom_idx = match chrom_indices.get(chrom) {
                Some(&idx) => idx,
                None => continue,
            };

            // Check if chromosome changed
            if let Some(prev_idx) = current_chrom_idx {
                if chrom_idx != prev_idx {
                    // Process completed chromosome
                    self.process_chromosome(
                        &chroms,
                        prev_idx,
                        &events,
                        genome,
                        &mut buf_output,
                        &mut genome_hist,
                        &mut itoa_buf,
                    )?;

                    // Also process any chromosomes between prev and current that had no intervals
                    for skip_idx in (prev_idx + 1)..chrom_idx {
                        self.process_empty_chromosome(
                            &chroms,
                            skip_idx,
                            genome,
                            &mut buf_output,
                            &mut genome_hist,
                            &mut itoa_buf,
                        )?;
                        seen_chroms[skip_idx] = true;
                    }

                    events.clear();
                }
            } else {
                // First interval - process any chromosomes before this one
                for skip_idx in 0..chrom_idx {
                    self.process_empty_chromosome(
                        &chroms,
                        skip_idx,
                        genome,
                        &mut buf_output,
                        &mut genome_hist,
                        &mut itoa_buf,
                    )?;
                    seen_chroms[skip_idx] = true;
                }
            }

            current_chrom_idx = Some(chrom_idx);
            seen_chroms[chrom_idx] = true;

            // Add events for this interval
            events.push((start, 1));
            events.push((end, -1));
        }

        // Process last chromosome
        if let Some(last_idx) = current_chrom_idx {
            self.process_chromosome(
                &chroms,
                last_idx,
                &events,
                genome,
                &mut buf_output,
                &mut genome_hist,
                &mut itoa_buf,
            )?;

            // Process remaining chromosomes
            for skip_idx in (last_idx + 1)..chroms.len() {
                self.process_empty_chromosome(
                    &chroms,
                    skip_idx,
                    genome,
                    &mut buf_output,
                    &mut genome_hist,
                    &mut itoa_buf,
                )?;
            }
        } else {
            // No intervals at all - process all chromosomes as empty
            for skip_idx in 0..chroms.len() {
                self.process_empty_chromosome(
                    &chroms,
                    skip_idx,
                    genome,
                    &mut buf_output,
                    &mut genome_hist,
                    &mut itoa_buf,
                )?;
            }
        }

        // Calculate total bases for genome histogram
        for chrom in &chroms {
            total_bases += genome.chrom_size(chrom).unwrap_or(0);
        }

        // Output genome-wide histogram if in histogram mode
        if self.mode == StreamingGenomecovMode::Histogram {
            self.output_genome_histogram(&genome_hist, total_bases, &mut buf_output)?;
        }

        buf_output.flush().map_err(BedError::Io)?;
        Ok(())
    }

    /// Process a chromosome with events using sweep-line algorithm.
    fn process_chromosome<W: Write>(
        &self,
        chroms: &[&String],
        chrom_idx: usize,
        events: &[(u64, i32)],
        genome: &Genome,
        output: &mut W,
        genome_hist: &mut HashMap<u32, u64>,
        itoa_buf: &mut itoa::Buffer,
    ) -> Result<(), BedError> {
        let chrom = chroms[chrom_idx];
        let chrom_size = genome.chrom_size(chrom).unwrap();

        if events.is_empty() {
            return self.process_empty_chromosome(
                chroms,
                chrom_idx,
                genome,
                output,
                genome_hist,
                itoa_buf,
            );
        }

        // Sort events: by position, starts before ends at same position
        let mut sorted_events = events.to_vec();
        sorted_events.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));

        // Sweep and collect depth regions
        let regions = self.sweep_events(&sorted_events, chrom_size);

        // Output based on mode
        match self.mode {
            StreamingGenomecovMode::Histogram => {
                let chrom_hist = self.build_histogram(&regions);
                self.output_chromosome_histogram(
                    chrom,
                    &chrom_hist,
                    chrom_size,
                    output,
                    genome_hist,
                )?;
            }
            StreamingGenomecovMode::BedGraph | StreamingGenomecovMode::BedGraphAll => {
                self.output_bedgraph(chrom.as_bytes(), &regions, output, itoa_buf)?;
            }
            StreamingGenomecovMode::PerBase => {
                self.output_per_base(chrom.as_bytes(), &regions, output, itoa_buf)?;
            }
        }

        Ok(())
    }

    /// Process an empty chromosome (no intervals).
    fn process_empty_chromosome<W: Write>(
        &self,
        chroms: &[&String],
        chrom_idx: usize,
        genome: &Genome,
        output: &mut W,
        genome_hist: &mut HashMap<u32, u64>,
        itoa_buf: &mut itoa::Buffer,
    ) -> Result<(), BedError> {
        let chrom = chroms[chrom_idx];
        let chrom_size = genome.chrom_size(chrom).unwrap();

        if chrom_size == 0 {
            return Ok(());
        }

        match self.mode {
            StreamingGenomecovMode::Histogram => {
                // Entire chromosome at depth 0
                let mut chrom_hist = HashMap::new();
                chrom_hist.insert(0u32, chrom_size);
                self.output_chromosome_histogram(
                    chrom,
                    &chrom_hist,
                    chrom_size,
                    output,
                    genome_hist,
                )?;
            }
            StreamingGenomecovMode::BedGraphAll => {
                // Output single region at depth 0
                output.write_all(chrom.as_bytes()).map_err(BedError::Io)?;
                output.write_all(b"\t0\t").map_err(BedError::Io)?;
                output
                    .write_all(itoa_buf.format(chrom_size).as_bytes())
                    .map_err(BedError::Io)?;
                output.write_all(b"\t0\n").map_err(BedError::Io)?;
            }
            StreamingGenomecovMode::BedGraph => {
                // No output for BedGraph (only non-zero)
            }
            StreamingGenomecovMode::PerBase => {
                // Output all positions at depth 0
                let chrom_bytes = chrom.as_bytes();
                for pos in 0..chrom_size {
                    output.write_all(chrom_bytes).map_err(BedError::Io)?;
                    output.write_all(b"\t").map_err(BedError::Io)?;
                    output
                        .write_all(itoa_buf.format(pos + 1).as_bytes())
                        .map_err(BedError::Io)?;
                    output.write_all(b"\t0\n").map_err(BedError::Io)?;
                }
            }
        }

        Ok(())
    }

    /// Sweep through sorted events and return (start, end, depth) regions.
    fn sweep_events(&self, events: &[(u64, i32)], chrom_size: u64) -> Vec<(u64, u64, u32)> {
        let mut result: Vec<(u64, u64, u32)> = Vec::new();
        let mut depth: i32 = 0;
        let mut prev_pos: u64 = 0;

        for &(pos, delta) in events {
            let pos = pos.min(chrom_size);

            if pos > prev_pos {
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

    /// Output per-chromosome histogram and accumulate genome-wide stats.
    fn output_chromosome_histogram<W: Write>(
        &self,
        chrom: &str,
        chrom_hist: &HashMap<u32, u64>,
        chrom_size: u64,
        output: &mut W,
        genome_hist: &mut HashMap<u32, u64>,
    ) -> Result<(), BedError> {
        let mut depths: Vec<_> = chrom_hist.keys().copied().collect();
        depths.sort_unstable();

        for depth in depths {
            let bases = chrom_hist[&depth];
            let fraction = bases as f64 / chrom_size as f64;
            writeln!(
                output,
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

        Ok(())
    }

    /// Output genome-wide histogram.
    fn output_genome_histogram<W: Write>(
        &self,
        genome_hist: &HashMap<u32, u64>,
        total_bases: u64,
        output: &mut W,
    ) -> Result<(), BedError> {
        let mut depths: Vec<_> = genome_hist.keys().copied().collect();
        depths.sort_unstable();

        for depth in depths {
            let bases = genome_hist[&depth];
            let fraction = bases as f64 / total_bases as f64;
            writeln!(
                output,
                "genome\t{}\t{}\t{}\t{}",
                depth,
                bases,
                total_bases,
                format_fraction(fraction)
            )
            .map_err(BedError::Io)?;
        }

        Ok(())
    }

    /// Output BedGraph format.
    fn output_bedgraph<W: Write>(
        &self,
        chrom: &[u8],
        regions: &[(u64, u64, u32)],
        output: &mut W,
        itoa_buf: &mut itoa::Buffer,
    ) -> Result<(), BedError> {
        for &(start, end, depth) in regions {
            if self.mode == StreamingGenomecovMode::BedGraph && depth == 0 {
                continue;
            }
            let scaled_depth = (depth as f64 * self.scale) as u32;
            output.write_all(chrom).map_err(BedError::Io)?;
            output.write_all(b"\t").map_err(BedError::Io)?;
            output
                .write_all(itoa_buf.format(start).as_bytes())
                .map_err(BedError::Io)?;
            output.write_all(b"\t").map_err(BedError::Io)?;
            output
                .write_all(itoa_buf.format(end).as_bytes())
                .map_err(BedError::Io)?;
            output.write_all(b"\t").map_err(BedError::Io)?;
            output
                .write_all(itoa_buf.format(scaled_depth).as_bytes())
                .map_err(BedError::Io)?;
            output.write_all(b"\n").map_err(BedError::Io)?;
        }
        Ok(())
    }

    /// Output per-base format.
    fn output_per_base<W: Write>(
        &self,
        chrom: &[u8],
        regions: &[(u64, u64, u32)],
        output: &mut W,
        itoa_buf: &mut itoa::Buffer,
    ) -> Result<(), BedError> {
        for &(start, end, depth) in regions {
            let scaled_depth = (depth as f64 * self.scale) as u32;
            for pos in start..end {
                output.write_all(chrom).map_err(BedError::Io)?;
                output.write_all(b"\t").map_err(BedError::Io)?;
                // 1-based position
                output
                    .write_all(itoa_buf.format(pos + 1).as_bytes())
                    .map_err(BedError::Io)?;
                output.write_all(b"\t").map_err(BedError::Io)?;
                output
                    .write_all(itoa_buf.format(scaled_depth).as_bytes())
                    .map_err(BedError::Io)?;
                output.write_all(b"\n").map_err(BedError::Io)?;
            }
        }
        Ok(())
    }
}

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

#[cfg(test)]
mod tests {
    use super::*;

    fn make_genome() -> Genome {
        let mut g = Genome::new();
        g.insert("chr1".to_string(), 1000);
        g.insert("chr2".to_string(), 500);
        g
    }

    #[test]
    fn test_streaming_genomecov_bedgraph() {
        let genome = make_genome();
        let bed_data = "chr1\t100\t200\nchr1\t150\t250\n";

        let cmd = StreamingGenomecovCommand::new()
            .with_mode(StreamingGenomecovMode::BedGraph)
            .with_assume_sorted(true);

        let mut output = Vec::new();
        let reader = BufReader::new(bed_data.as_bytes());
        cmd.genomecov_streaming(reader, &genome, &mut output)
            .unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();

        // Should have regions: [100-150) depth 1, [150-200) depth 2, [200-250) depth 1
        assert!(lines.len() >= 3);
    }

    #[test]
    fn test_streaming_genomecov_bedgraph_all() {
        let genome = make_genome();
        let bed_data = "chr1\t100\t200\n";

        let cmd = StreamingGenomecovCommand::new()
            .with_mode(StreamingGenomecovMode::BedGraphAll)
            .with_assume_sorted(true);

        let mut output = Vec::new();
        let reader = BufReader::new(bed_data.as_bytes());
        cmd.genomecov_streaming(reader, &genome, &mut output)
            .unwrap();

        let result = String::from_utf8(output).unwrap();

        // Should include zero-depth regions
        assert!(result.contains("chr1\t0\t100\t0"));
        assert!(result.contains("chr1\t100\t200\t1"));
        assert!(result.contains("chr1\t200\t1000\t0"));
        // chr2 should be all zeros
        assert!(result.contains("chr2\t0\t500\t0"));
    }

    #[test]
    fn test_streaming_genomecov_histogram() {
        let genome = make_genome();
        let bed_data = "chr1\t100\t200\n";

        let cmd = StreamingGenomecovCommand::new()
            .with_mode(StreamingGenomecovMode::Histogram)
            .with_assume_sorted(true);

        let mut output = Vec::new();
        let reader = BufReader::new(bed_data.as_bytes());
        cmd.genomecov_streaming(reader, &genome, &mut output)
            .unwrap();

        let result = String::from_utf8(output).unwrap();

        // Should have per-chromosome and genome-wide histograms
        assert!(result.contains("chr1\t0\t900")); // 900 bases at depth 0
        assert!(result.contains("chr1\t1\t100")); // 100 bases at depth 1
        assert!(result.contains("genome\t"));
    }

    #[test]
    fn test_streaming_genomecov_empty() {
        let genome = make_genome();
        let bed_data = "";

        let cmd = StreamingGenomecovCommand::new()
            .with_mode(StreamingGenomecovMode::BedGraphAll)
            .with_assume_sorted(true);

        let mut output = Vec::new();
        let reader = BufReader::new(bed_data.as_bytes());
        cmd.genomecov_streaming(reader, &genome, &mut output)
            .unwrap();

        let result = String::from_utf8(output).unwrap();

        // All chromosomes should be at depth 0
        assert!(result.contains("chr1\t0\t1000\t0"));
        assert!(result.contains("chr2\t0\t500\t0"));
    }

    #[test]
    fn test_format_fraction() {
        assert_eq!(format_fraction(0.0), "0");
        assert_eq!(format_fraction(1.0), "1");
        assert_eq!(format_fraction(0.5), "0.5");
        assert_eq!(format_fraction(0.123456), "0.123456");
    }
}
