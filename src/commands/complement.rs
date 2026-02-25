//! Complement command implementation.
//!
//! Returns genomic regions NOT covered by intervals.
//! Single-pass O(n) streaming algorithm.

#![allow(clippy::needless_range_loop)]

use crate::bed::{BedError, BedReader};
use crate::genome::Genome;
use crate::interval::Interval;
use crate::streaming::parsing::{parse_bed3_bytes, should_skip_line};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Complement command configuration.
#[derive(Debug, Clone)]
pub struct ComplementCommand {
    /// Only compute complement for chromosomes in the genome file
    pub genome_only: bool,
    /// Assume input is sorted in genome order (enables O(1) memory streaming)
    pub assume_sorted: bool,
}

impl Default for ComplementCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl ComplementCommand {
    pub fn new() -> Self {
        Self {
            genome_only: true,
            assume_sorted: false,
        }
    }

    /// Set assume_sorted flag (builder pattern).
    pub fn with_assume_sorted(mut self, assume_sorted: bool) -> Self {
        self.assume_sorted = assume_sorted;
        self
    }

    /// Streaming complement - assumes input is sorted by chrom, start, end.
    /// O(n) single pass through input, outputs in genome file order.
    pub fn complement_streaming<R: Read, W: Write>(
        &self,
        reader: BedReader<R>,
        genome: &Genome,
        output: &mut W,
    ) -> Result<(), BedError> {
        use std::collections::HashMap;

        let mut buf_output = BufWriter::with_capacity(256 * 1024, output);

        // Accumulate gaps per chromosome
        let mut gaps: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let mut chrom_last_end: HashMap<String, u64> = HashMap::new();

        for result in reader.records() {
            let record = result?;
            let chrom = record.chrom();

            // Skip chromosomes not in genome
            let chrom_size = match genome.chrom_size(chrom) {
                Some(size) => size,
                None => continue,
            };

            let prev_end = chrom_last_end.entry(chrom.to_string()).or_insert(0);

            // Record gap if there's space before this interval
            if record.start() > *prev_end {
                gaps.entry(chrom.to_string())
                    .or_default()
                    .push((*prev_end, record.start()));
            }

            // Update last end (handle overlaps)
            *prev_end = (*prev_end).max(record.end().min(chrom_size));
        }

        // Output in genome order
        for chrom in genome.chromosomes() {
            let chrom_size = genome.chrom_size(chrom).unwrap();
            let last_end = chrom_last_end.get(chrom).copied();

            if let Some(chrom_gaps) = gaps.get(chrom) {
                // Output accumulated gaps
                for (start, end) in chrom_gaps {
                    writeln!(buf_output, "{}\t{}\t{}", chrom, start, end).map_err(BedError::Io)?;
                }
                // Output trailing gap
                if let Some(end_pos) = last_end {
                    if end_pos < chrom_size {
                        writeln!(buf_output, "{}\t{}\t{}", chrom, end_pos, chrom_size)
                            .map_err(BedError::Io)?;
                    }
                }
            } else if let Some(end_pos) = last_end {
                // We saw intervals on this chromosome but no gaps (e.g., interval starts at 0)
                // Just output trailing gap if any
                if end_pos < chrom_size {
                    writeln!(buf_output, "{}\t{}\t{}", chrom, end_pos, chrom_size)
                        .map_err(BedError::Io)?;
                }
            } else {
                // No intervals on this chromosome - entire chromosome is complement
                if chrom_size > 0 {
                    writeln!(buf_output, "{}\t{}\t{}", chrom, 0, chrom_size)
                        .map_err(BedError::Io)?;
                }
            }
        }

        buf_output.flush().map_err(BedError::Io)?;
        Ok(())
    }

    /// True streaming complement for sorted input - O(1) memory.
    ///
    /// When input is sorted in genome file order, we can stream output directly
    /// without accumulating gaps in memory. This reduces memory from ~96MB to ~2MB.
    ///
    /// Algorithm:
    /// - Track current chromosome and last_end position
    /// - Output gaps immediately as we process each interval
    /// - Handle chromosome transitions and missing chromosomes
    pub fn complement_streaming_sorted<R: Read, W: Write>(
        &self,
        reader: BedReader<R>,
        genome: &Genome,
        output: &mut W,
    ) -> Result<(), BedError> {
        let mut buf_output = BufWriter::with_capacity(256 * 1024, output);

        // Get genome chromosomes as ordered list
        let chroms: Vec<&String> = genome.chromosomes().collect();
        let chrom_indices: std::collections::HashMap<&str, usize> = chroms
            .iter()
            .enumerate()
            .map(|(i, c)| (c.as_str(), i))
            .collect();

        // State: current chromosome index and last end position
        let mut current_chrom_idx: Option<usize> = None;
        let mut last_end: u64 = 0;

        for result in reader.records() {
            let record = result?;
            let chrom = record.chrom();

            // Skip chromosomes not in genome
            let chrom_idx = match chrom_indices.get(chrom) {
                Some(&idx) => idx,
                None => continue,
            };

            let chrom_size = genome.chrom_size(chrom).unwrap();

            match current_chrom_idx {
                None => {
                    // First interval - output full chromosomes before this one
                    for i in 0..chrom_idx {
                        let c = chroms[i];
                        let size = genome.chrom_size(c).unwrap();
                        if size > 0 {
                            writeln!(buf_output, "{}\t0\t{}", c, size).map_err(BedError::Io)?;
                        }
                    }
                    // Output leading gap on current chromosome
                    if record.start() > 0 {
                        writeln!(buf_output, "{}\t0\t{}", chrom, record.start())
                            .map_err(BedError::Io)?;
                    }
                    current_chrom_idx = Some(chrom_idx);
                    last_end = record.end().min(chrom_size);
                }
                Some(prev_idx) if chrom_idx != prev_idx => {
                    // Chromosome changed
                    // 1. Output trailing gap for previous chromosome
                    let prev_chrom = chroms[prev_idx];
                    let prev_size = genome.chrom_size(prev_chrom).unwrap();
                    if last_end < prev_size {
                        writeln!(buf_output, "{}\t{}\t{}", prev_chrom, last_end, prev_size)
                            .map_err(BedError::Io)?;
                    }

                    // 2. Output full chromosomes between prev and current
                    for i in (prev_idx + 1)..chrom_idx {
                        let c = chroms[i];
                        let size = genome.chrom_size(c).unwrap();
                        if size > 0 {
                            writeln!(buf_output, "{}\t0\t{}", c, size).map_err(BedError::Io)?;
                        }
                    }

                    // 3. Output leading gap on current chromosome
                    if record.start() > 0 {
                        writeln!(buf_output, "{}\t0\t{}", chrom, record.start())
                            .map_err(BedError::Io)?;
                    }

                    current_chrom_idx = Some(chrom_idx);
                    last_end = record.end().min(chrom_size);
                }
                Some(_) => {
                    // Same chromosome - output gap if there's space
                    if record.start() > last_end {
                        writeln!(buf_output, "{}\t{}\t{}", chrom, last_end, record.start())
                            .map_err(BedError::Io)?;
                    }
                    last_end = last_end.max(record.end().min(chrom_size));
                }
            }
        }

        // Handle remaining chromosomes
        match current_chrom_idx {
            Some(last_idx) => {
                // Trailing gap for last processed chromosome
                let last_chrom = chroms[last_idx];
                let last_size = genome.chrom_size(last_chrom).unwrap();
                if last_end < last_size {
                    writeln!(buf_output, "{}\t{}\t{}", last_chrom, last_end, last_size)
                        .map_err(BedError::Io)?;
                }

                // Full chromosomes after the last one
                for i in (last_idx + 1)..chroms.len() {
                    let c = chroms[i];
                    let size = genome.chrom_size(c).unwrap();
                    if size > 0 {
                        writeln!(buf_output, "{}\t0\t{}", c, size).map_err(BedError::Io)?;
                    }
                }
            }
            None => {
                // No intervals at all - entire genome is complement
                for c in &chroms {
                    let size = genome.chrom_size(c).unwrap();
                    if size > 0 {
                        writeln!(buf_output, "{}\t0\t{}", c, size).map_err(BedError::Io)?;
                    }
                }
            }
        }

        buf_output.flush().map_err(BedError::Io)?;
        Ok(())
    }

    /// Run complement on a file with streaming output.
    ///
    /// Uses fast raw byte parsing for maximum performance.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        input: P,
        genome: &Genome,
        output: &mut W,
    ) -> Result<(), BedError> {
        let file = File::open(input)?;
        self.complement_fast(file, genome, output)
    }

    /// Fast complement using raw byte parsing.
    /// O(n) streaming with O(1) memory per chromosome.
    fn complement_fast<R: Read, W: Write>(
        &self,
        input: R,
        genome: &Genome,
        output: &mut W,
    ) -> Result<(), BedError> {
        let mut reader = BufReader::with_capacity(256 * 1024, input);
        let mut buf_output = BufWriter::with_capacity(256 * 1024, output);

        // Get genome chromosomes as ordered list
        let chroms: Vec<&String> = genome.chromosomes().collect();
        let chrom_indices: std::collections::HashMap<&[u8], usize> = chroms
            .iter()
            .enumerate()
            .map(|(i, c)| (c.as_bytes(), i))
            .collect();

        // State
        let mut current_chrom_idx: Option<usize> = None;
        let mut last_end: u64 = 0;
        let mut line_buf = String::with_capacity(1024);

        // Reusable output buffer for itoa
        let mut itoa_buf = itoa::Buffer::new();

        loop {
            line_buf.clear();
            let bytes_read = reader.read_line(&mut line_buf).map_err(BedError::Io)?;
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

            let chrom_size = genome.chrom_size(chroms[chrom_idx]).unwrap();

            match current_chrom_idx {
                None => {
                    // First interval - output full chromosomes before this one
                    for i in 0..chrom_idx {
                        let c = chroms[i];
                        let size = genome.chrom_size(c).unwrap();
                        if size > 0 {
                            Self::write_interval_fast(
                                &mut buf_output,
                                c.as_bytes(),
                                0,
                                size,
                                &mut itoa_buf,
                            )?;
                        }
                    }
                    // Output leading gap on current chromosome
                    if start > 0 {
                        Self::write_interval_fast(&mut buf_output, chrom, 0, start, &mut itoa_buf)?;
                    }
                    current_chrom_idx = Some(chrom_idx);
                    last_end = end.min(chrom_size);
                }
                Some(prev_idx) if chrom_idx != prev_idx => {
                    // Chromosome changed
                    let prev_chrom = chroms[prev_idx];
                    let prev_size = genome.chrom_size(prev_chrom).unwrap();
                    if last_end < prev_size {
                        Self::write_interval_fast(
                            &mut buf_output,
                            prev_chrom.as_bytes(),
                            last_end,
                            prev_size,
                            &mut itoa_buf,
                        )?;
                    }

                    // Output full chromosomes between prev and current
                    for i in (prev_idx + 1)..chrom_idx {
                        let c = chroms[i];
                        let size = genome.chrom_size(c).unwrap();
                        if size > 0 {
                            Self::write_interval_fast(
                                &mut buf_output,
                                c.as_bytes(),
                                0,
                                size,
                                &mut itoa_buf,
                            )?;
                        }
                    }

                    // Output leading gap on current chromosome
                    if start > 0 {
                        Self::write_interval_fast(&mut buf_output, chrom, 0, start, &mut itoa_buf)?;
                    }

                    current_chrom_idx = Some(chrom_idx);
                    last_end = end.min(chrom_size);
                }
                Some(_) => {
                    // Same chromosome - output gap if there's space
                    if start > last_end {
                        Self::write_interval_fast(
                            &mut buf_output,
                            chrom,
                            last_end,
                            start,
                            &mut itoa_buf,
                        )?;
                    }
                    last_end = last_end.max(end.min(chrom_size));
                }
            }
        }

        // Handle remaining chromosomes
        match current_chrom_idx {
            Some(last_idx) => {
                let last_chrom = chroms[last_idx];
                let last_size = genome.chrom_size(last_chrom).unwrap();
                if last_end < last_size {
                    Self::write_interval_fast(
                        &mut buf_output,
                        last_chrom.as_bytes(),
                        last_end,
                        last_size,
                        &mut itoa_buf,
                    )?;
                }

                for i in (last_idx + 1)..chroms.len() {
                    let c = chroms[i];
                    let size = genome.chrom_size(c).unwrap();
                    if size > 0 {
                        Self::write_interval_fast(
                            &mut buf_output,
                            c.as_bytes(),
                            0,
                            size,
                            &mut itoa_buf,
                        )?;
                    }
                }
            }
            None => {
                // No intervals - entire genome is complement
                for c in &chroms {
                    let size = genome.chrom_size(c).unwrap();
                    if size > 0 {
                        Self::write_interval_fast(
                            &mut buf_output,
                            c.as_bytes(),
                            0,
                            size,
                            &mut itoa_buf,
                        )?;
                    }
                }
            }
        }

        buf_output.flush().map_err(BedError::Io)?;
        Ok(())
    }

    /// Fast interval output using itoa.
    #[inline]
    fn write_interval_fast<W: Write>(
        output: &mut W,
        chrom: &[u8],
        start: u64,
        end: u64,
        itoa_buf: &mut itoa::Buffer,
    ) -> Result<(), BedError> {
        output.write_all(chrom).map_err(BedError::Io)?;
        output.write_all(b"\t").map_err(BedError::Io)?;
        output
            .write_all(itoa_buf.format(start).as_bytes())
            .map_err(BedError::Io)?;
        output.write_all(b"\t").map_err(BedError::Io)?;
        output
            .write_all(itoa_buf.format(end).as_bytes())
            .map_err(BedError::Io)?;
        output.write_all(b"\n").map_err(BedError::Io)?;
        Ok(())
    }

    /// Compute complement of intervals against genome (for testing).
    /// Returns gaps between intervals and at chromosome boundaries.
    pub fn complement(&self, intervals: &[Interval], genome: &Genome) -> Vec<Interval> {
        use std::collections::HashMap;

        // Group intervals by chromosome
        let mut by_chrom: HashMap<&str, Vec<&Interval>> = HashMap::new();
        for interval in intervals {
            by_chrom
                .entry(interval.chrom.as_str())
                .or_default()
                .push(interval);
        }

        // Sort intervals within each chromosome
        for intervals in by_chrom.values_mut() {
            intervals.sort_by_key(|i| (i.start, i.end));
        }

        let mut result = Vec::new();

        // Process each chromosome in the genome
        for chrom in genome.chromosomes() {
            let chrom_size = genome.chrom_size(chrom).unwrap();

            if let Some(intervals) = by_chrom.get(chrom.as_str()) {
                // Find gaps in this chromosome
                let gaps = self.find_gaps(chrom, intervals, chrom_size);
                result.extend(gaps);
            } else {
                // No intervals on this chromosome - entire chromosome is complement
                if chrom_size > 0 {
                    result.push(Interval::new(chrom.clone(), 0, chrom_size));
                }
            }
        }

        result
    }

    /// Find gaps in a sorted list of intervals on a single chromosome.
    fn find_gaps(&self, chrom: &str, intervals: &[&Interval], chrom_size: u64) -> Vec<Interval> {
        let mut gaps = Vec::new();
        let mut prev_end: u64 = 0;

        for interval in intervals {
            // If there's a gap before this interval, emit it
            if interval.start > prev_end {
                gaps.push(Interval::new(chrom.to_string(), prev_end, interval.start));
            }

            // Update the previous end position
            prev_end = prev_end.max(interval.end);
        }

        // Emit trailing gap if the last interval doesn't reach chrom end
        if prev_end < chrom_size {
            gaps.push(Interval::new(chrom.to_string(), prev_end, chrom_size));
        }

        gaps
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

    fn make_interval(chrom: &str, start: u64, end: u64) -> Interval {
        Interval::new(chrom.to_string(), start, end)
    }

    #[test]
    fn test_complement_basic() {
        let genome = make_genome();
        let intervals = vec![
            make_interval("chr1", 100, 200),
            make_interval("chr1", 300, 400),
        ];

        let cmd = ComplementCommand::new();
        let result = cmd.complement(&intervals, &genome);

        // Should have gaps: [0-100), [200-300), [400-1000), and entire chr2
        assert_eq!(result.len(), 4);

        // chr1 gaps
        assert_eq!(result[0], make_interval("chr1", 0, 100));
        assert_eq!(result[1], make_interval("chr1", 200, 300));
        assert_eq!(result[2], make_interval("chr1", 400, 1000));

        // chr2 - entire chromosome
        assert_eq!(result[3], make_interval("chr2", 0, 500));
    }

    #[test]
    fn test_complement_overlapping() {
        let genome = make_genome();
        let intervals = vec![
            make_interval("chr1", 100, 300),
            make_interval("chr1", 200, 400), // overlaps previous
        ];

        let cmd = ComplementCommand::new();
        let result = cmd.complement(&intervals, &genome);

        // Overlapping intervals should be merged for complement purpose
        // Gaps: [0-100), [400-1000), entire chr2
        assert_eq!(result.len(), 3);
        assert_eq!(result[0], make_interval("chr1", 0, 100));
        assert_eq!(result[1], make_interval("chr1", 400, 1000));
        assert_eq!(result[2], make_interval("chr2", 0, 500));
    }

    #[test]
    fn test_complement_full_coverage() {
        let genome = make_genome();
        let intervals = vec![make_interval("chr1", 0, 1000)];

        let cmd = ComplementCommand::new();
        let result = cmd.complement(&intervals, &genome);

        // chr1 fully covered, only chr2 in complement
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], make_interval("chr2", 0, 500));
    }

    #[test]
    fn test_complement_empty() {
        let genome = make_genome();
        let intervals: Vec<Interval> = vec![];

        let cmd = ComplementCommand::new();
        let result = cmd.complement(&intervals, &genome);

        // All chromosomes should be in complement
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_complement_at_boundaries() {
        let genome = make_genome();
        let intervals = vec![
            make_interval("chr1", 0, 100),    // starts at 0
            make_interval("chr1", 900, 1000), // ends at chrom end
        ];

        let cmd = ComplementCommand::new();
        let result = cmd.complement(&intervals, &genome);

        // Gap only in middle of chr1, plus entire chr2
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], make_interval("chr1", 100, 900));
        assert_eq!(result[1], make_interval("chr2", 0, 500));
    }

    // Tests for streaming sorted implementation
    #[test]
    fn test_streaming_sorted_basic() {
        let genome = make_genome();
        let bed_data = "chr1\t100\t200\nchr1\t300\t400\n";

        // Test unsorted path
        let cmd_unsorted = ComplementCommand::new();
        let mut output_unsorted = Vec::new();
        let reader = BedReader::new(bed_data.as_bytes());
        cmd_unsorted
            .complement_streaming(reader, &genome, &mut output_unsorted)
            .unwrap();

        // Test sorted path
        let cmd_sorted = ComplementCommand::new().with_assume_sorted(true);
        let mut output_sorted = Vec::new();
        let reader = BedReader::new(bed_data.as_bytes());
        cmd_sorted
            .complement_streaming_sorted(reader, &genome, &mut output_sorted)
            .unwrap();

        // Both should produce identical output
        assert_eq!(
            String::from_utf8(output_unsorted).unwrap(),
            String::from_utf8(output_sorted).unwrap()
        );
    }

    #[test]
    fn test_streaming_sorted_overlapping() {
        let genome = make_genome();
        let bed_data = "chr1\t100\t300\nchr1\t200\t400\n";

        let cmd_unsorted = ComplementCommand::new();
        let mut output_unsorted = Vec::new();
        let reader = BedReader::new(bed_data.as_bytes());
        cmd_unsorted
            .complement_streaming(reader, &genome, &mut output_unsorted)
            .unwrap();

        let cmd_sorted = ComplementCommand::new().with_assume_sorted(true);
        let mut output_sorted = Vec::new();
        let reader = BedReader::new(bed_data.as_bytes());
        cmd_sorted
            .complement_streaming_sorted(reader, &genome, &mut output_sorted)
            .unwrap();

        assert_eq!(
            String::from_utf8(output_unsorted).unwrap(),
            String::from_utf8(output_sorted).unwrap()
        );
    }

    #[test]
    fn test_streaming_sorted_empty_input() {
        let genome = make_genome();
        let bed_data = "";

        let cmd_unsorted = ComplementCommand::new();
        let mut output_unsorted = Vec::new();
        let reader = BedReader::new(bed_data.as_bytes());
        cmd_unsorted
            .complement_streaming(reader, &genome, &mut output_unsorted)
            .unwrap();

        let cmd_sorted = ComplementCommand::new().with_assume_sorted(true);
        let mut output_sorted = Vec::new();
        let reader = BedReader::new(bed_data.as_bytes());
        cmd_sorted
            .complement_streaming_sorted(reader, &genome, &mut output_sorted)
            .unwrap();

        assert_eq!(
            String::from_utf8(output_unsorted).unwrap(),
            String::from_utf8(output_sorted).unwrap()
        );
    }

    #[test]
    fn test_streaming_sorted_multi_chrom() {
        let genome = make_genome();
        let bed_data = "chr1\t100\t200\nchr1\t300\t400\nchr2\t50\t100\n";

        let cmd_unsorted = ComplementCommand::new();
        let mut output_unsorted = Vec::new();
        let reader = BedReader::new(bed_data.as_bytes());
        cmd_unsorted
            .complement_streaming(reader, &genome, &mut output_unsorted)
            .unwrap();

        let cmd_sorted = ComplementCommand::new().with_assume_sorted(true);
        let mut output_sorted = Vec::new();
        let reader = BedReader::new(bed_data.as_bytes());
        cmd_sorted
            .complement_streaming_sorted(reader, &genome, &mut output_sorted)
            .unwrap();

        assert_eq!(
            String::from_utf8(output_unsorted).unwrap(),
            String::from_utf8(output_sorted).unwrap()
        );
    }

    #[test]
    fn test_streaming_sorted_full_coverage() {
        let genome = make_genome();
        let bed_data = "chr1\t0\t1000\n";

        let cmd_unsorted = ComplementCommand::new();
        let mut output_unsorted = Vec::new();
        let reader = BedReader::new(bed_data.as_bytes());
        cmd_unsorted
            .complement_streaming(reader, &genome, &mut output_unsorted)
            .unwrap();

        let cmd_sorted = ComplementCommand::new().with_assume_sorted(true);
        let mut output_sorted = Vec::new();
        let reader = BedReader::new(bed_data.as_bytes());
        cmd_sorted
            .complement_streaming_sorted(reader, &genome, &mut output_sorted)
            .unwrap();

        assert_eq!(
            String::from_utf8(output_unsorted).unwrap(),
            String::from_utf8(output_sorted).unwrap()
        );
    }
}
