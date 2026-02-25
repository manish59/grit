//! Coverage command implementation - compute depth/breadth of coverage.
//!
//! Uses O(n + m) sweep-line algorithm per chromosome for optimal performance.
//! Event-based depth computation avoids per-base iteration for basic/mean modes.

use crate::bed::{read_records, BedError};
use crate::interval::BedRecord;
use crate::parallel::PARALLEL_THRESHOLD;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// Coverage command configuration.
#[derive(Debug, Clone)]
pub struct CoverageCommand {
    /// Report a histogram of coverage
    pub histogram: bool,
    /// Report depth at each position
    pub per_base: bool,
    /// Report mean depth
    pub mean: bool,
    /// Require same strand
    pub same_strand: bool,
    /// Require opposite strand
    pub opposite_strand: bool,
    /// Process in parallel by chromosome
    pub parallel: bool,
}

impl Default for CoverageCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl CoverageCommand {
    pub fn new() -> Self {
        Self {
            histogram: false,
            per_base: false,
            mean: false,
            same_strand: false,
            opposite_strand: false,
            parallel: true,
        }
    }

    /// Execute coverage command on files using O(n+m) sweep-line algorithm.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        a_path: P,
        b_path: P,
        output: &mut W,
    ) -> Result<(), BedError> {
        let a_records = read_records(&a_path)?;
        let b_records = read_records(&b_path)?;

        if a_records.is_empty() {
            return Ok(());
        }

        // Group by chromosome
        let a_by_chrom = Self::group_records_by_chrom(a_records);
        let b_by_chrom = Self::group_records_by_chrom(b_records);

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
                    self.coverage_chromosome_sweepline(a_list, b_list, &mut buf);
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
                        self.coverage_chromosome_sweepline(a_list, b_list, &mut buf);
                    }
                    buf
                })
                .collect();

            // Write results in chromosome order
            for buf in results {
                output.write_all(&buf).map_err(BedError::Io)?;
            }
        }

        // Add genome-wide histogram summary if -hist
        if self.histogram {
            self.write_genome_histogram(&a_by_chrom, &b_by_chrom, output)?;
        }

        Ok(())
    }

    /// O(n+m) sweep-line coverage for a single chromosome.
    ///
    /// Algorithm:
    /// 1. A and B are sorted by start position
    /// 2. Use two-pointer to find overlapping B intervals for each A
    /// 3. Use event-based depth computation (O(k log k) instead of O(L))
    fn coverage_chromosome_sweepline(
        &self,
        a_sorted: &[BedRecord],
        b_sorted: Option<&Vec<BedRecord>>,
        output: &mut Vec<u8>,
    ) {
        let b_sorted = match b_sorted {
            Some(b) if !b.is_empty() => b,
            _ => {
                // No B intervals - output A with zero coverage
                for a_rec in a_sorted {
                    self.write_zero_coverage(output, a_rec);
                }
                return;
            }
        };

        let b_len = b_sorted.len();
        let mut b_ptr: usize = 0;

        for a_rec in a_sorted {
            let a_start = a_rec.start();
            let a_end = a_rec.end();
            let a_len = a_end - a_start;

            if a_len == 0 {
                self.write_zero_coverage(output, a_rec);
                continue;
            }

            // Advance b_ptr: skip B intervals that end before A starts
            while b_ptr < b_len && b_sorted[b_ptr].end() <= a_start {
                b_ptr += 1;
            }

            // Collect overlapping B intervals
            let mut overlaps: Vec<(u64, u64)> = Vec::new();
            for b_rec in b_sorted.iter().skip(b_ptr) {
                let b_start = b_rec.start();
                let b_end = b_rec.end();

                if b_start >= a_end {
                    break; // No more overlaps possible
                }
                if b_end > a_start {
                    // Clip to A boundaries
                    let clip_start = b_start.max(a_start);
                    let clip_end = b_end.min(a_end);
                    if clip_end > clip_start {
                        overlaps.push((clip_start, clip_end));
                    }
                }
            }

            let num_overlaps = overlaps.len();

            if overlaps.is_empty() {
                self.write_zero_coverage(output, a_rec);
                continue;
            }

            // Use event-based computation for efficiency
            if self.per_base {
                self.write_per_base_coverage(output, a_rec, &overlaps);
            } else if self.histogram {
                self.write_histogram_coverage(output, a_rec, num_overlaps, &overlaps);
            } else if self.mean {
                self.write_mean_coverage(output, a_rec, &overlaps);
            } else {
                self.write_basic_coverage(output, a_rec, num_overlaps, &overlaps);
            }
        }
    }

    /// Write zero coverage output for an A interval with no B overlaps.
    #[inline]
    fn write_zero_coverage(&self, buf: &mut Vec<u8>, a_rec: &BedRecord) {
        use std::io::Write as IoWrite;
        let a_len = a_rec.len();

        if self.per_base {
            // Output each position with depth 0
            for pos in 1..=a_len {
                Self::write_record_fields(buf, a_rec);
                let _ = writeln!(buf, "\t{}\t0", pos);
            }
        } else if self.histogram {
            // All bases at depth 0
            Self::write_record_fields(buf, a_rec);
            let _ = writeln!(buf, "\t0\t{}\t{}\t1.0000000", a_len, a_len);
        } else if self.mean {
            Self::write_record_fields(buf, a_rec);
            let _ = writeln!(buf, "\t0.0000000");
        } else {
            Self::write_record_fields(buf, a_rec);
            let _ = writeln!(buf, "\t0\t0\t{}\t0.0000000", a_len);
        }
    }

    /// Write all fields of a BedRecord to buffer (without newline).
    #[inline]
    fn write_record_fields(buf: &mut Vec<u8>, rec: &BedRecord) {
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
    }

    /// Compute coverage using event-based approach.
    /// Returns (bases_covered, total_depth, histogram).
    fn compute_coverage_events(
        &self,
        a_start: u64,
        a_end: u64,
        overlaps: &[(u64, u64)],
    ) -> (u64, u64, HashMap<u32, u64>) {
        // Create events: (position, delta)
        // +1 at start of B overlap, -1 at end
        let mut events: Vec<(u64, i32)> = Vec::with_capacity(overlaps.len() * 2 + 2);

        // Add A boundary events
        events.push((a_start, 0));
        events.push((a_end, 0));

        for &(start, end) in overlaps {
            events.push((start, 1));
            events.push((end, -1));
        }

        // Sort by position, with -1 events before +1 at same position
        events.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        // Sweep through events
        let mut depth: i32 = 0;
        let mut prev_pos = a_start;
        let mut bases_covered: u64 = 0;
        let mut total_depth: u64 = 0;
        let mut histogram: HashMap<u32, u64> = HashMap::new();

        for (pos, delta) in events {
            if pos > prev_pos && pos <= a_end && prev_pos >= a_start {
                let span = pos - prev_pos;
                if depth > 0 {
                    bases_covered += span;
                }
                total_depth += span * (depth as u64);
                *histogram.entry(depth as u32).or_insert(0) += span;
            }
            depth += delta;
            prev_pos = pos;
        }

        (bases_covered, total_depth, histogram)
    }

    /// Write basic coverage output: A + num_overlaps + bases_covered + length + fraction
    #[inline]
    fn write_basic_coverage(
        &self,
        buf: &mut Vec<u8>,
        a_rec: &BedRecord,
        num_overlaps: usize,
        overlaps: &[(u64, u64)],
    ) {
        use std::io::Write as IoWrite;

        let a_start = a_rec.start();
        let a_end = a_rec.end();
        let a_len = a_end - a_start;

        let (bases_covered, _, _) = self.compute_coverage_events(a_start, a_end, overlaps);
        // Use f32 to match bedtools precision (bedtools uses float internally)
        let fraction: f32 = bases_covered as f32 / a_len as f32;

        Self::write_record_fields(buf, a_rec);
        let _ = writeln!(
            buf,
            "\t{}\t{}\t{}\t{:.7}",
            num_overlaps, bases_covered, a_len, fraction
        );
    }

    /// Write mean coverage output: A + mean_depth
    #[inline]
    fn write_mean_coverage(&self, buf: &mut Vec<u8>, a_rec: &BedRecord, overlaps: &[(u64, u64)]) {
        use std::io::Write as IoWrite;

        let a_start = a_rec.start();
        let a_end = a_rec.end();
        let a_len = a_end - a_start;

        let (_, total_depth, _) = self.compute_coverage_events(a_start, a_end, overlaps);
        // Use f32 to match bedtools precision (bedtools uses float internally)
        let mean_depth: f32 = total_depth as f32 / a_len as f32;

        Self::write_record_fields(buf, a_rec);
        let _ = writeln!(buf, "\t{:.7}", mean_depth);
    }

    /// Write histogram coverage output: A + depth + count + length + fraction (for each depth)
    #[inline]
    fn write_histogram_coverage(
        &self,
        buf: &mut Vec<u8>,
        a_rec: &BedRecord,
        _num_overlaps: usize,
        overlaps: &[(u64, u64)],
    ) {
        use std::io::Write as IoWrite;

        let a_start = a_rec.start();
        let a_end = a_rec.end();
        let a_len = a_end - a_start;

        let (_, _, histogram) = self.compute_coverage_events(a_start, a_end, overlaps);

        // Sort depths for deterministic output
        let mut depths: Vec<_> = histogram.into_iter().collect();
        depths.sort_by_key(|&(d, _)| d);

        for (depth, count) in depths {
            // Use f32 to match bedtools precision (bedtools uses float internally)
            let fraction: f32 = count as f32 / a_len as f32;
            Self::write_record_fields(buf, a_rec);
            let _ = writeln!(buf, "\t{}\t{}\t{}\t{:.7}", depth, count, a_len, fraction);
        }
    }

    /// Write per-base coverage output: A + 1-based-position + depth
    #[inline]
    fn write_per_base_coverage(
        &self,
        buf: &mut Vec<u8>,
        a_rec: &BedRecord,
        overlaps: &[(u64, u64)],
    ) {
        use std::io::Write as IoWrite;

        let a_start = a_rec.start();
        let a_end = a_rec.end();

        // Create events for efficient depth computation
        let mut events: Vec<(u64, i32)> = Vec::with_capacity(overlaps.len() * 2);
        for &(start, end) in overlaps {
            events.push((start, 1));
            events.push((end, -1));
        }
        events.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        // Sweep through and output per-position depth
        let mut depth: i32 = 0;
        let mut event_idx = 0;

        for pos in a_start..a_end {
            // Process all events at this position
            while event_idx < events.len() && events[event_idx].0 <= pos {
                depth += events[event_idx].1;
                event_idx += 1;
            }

            let one_based_pos = pos - a_start + 1;
            Self::write_record_fields(buf, a_rec);
            let _ = writeln!(buf, "\t{}\t{}", one_based_pos, depth);
        }
    }

    /// Write genome-wide histogram summary (all depths across all intervals).
    fn write_genome_histogram<W: Write>(
        &self,
        a_by_chrom: &HashMap<String, Vec<BedRecord>>,
        b_by_chrom: &HashMap<String, Vec<BedRecord>>,
        output: &mut W,
    ) -> Result<(), BedError> {
        let mut total_histogram: HashMap<u32, u64> = HashMap::new();
        let mut total_length: u64 = 0;

        for (chrom, a_records) in a_by_chrom {
            let b_records = b_by_chrom.get(chrom);

            for a_rec in a_records {
                let a_start = a_rec.start();
                let a_end = a_rec.end();
                let a_len = a_end - a_start;
                total_length += a_len;

                if let Some(b_list) = b_records {
                    // Find overlapping B intervals
                    let overlaps: Vec<(u64, u64)> = b_list
                        .iter()
                        .filter(|b| b.end() > a_start && b.start() < a_end)
                        .map(|b| {
                            let clip_start = b.start().max(a_start);
                            let clip_end = b.end().min(a_end);
                            (clip_start, clip_end)
                        })
                        .filter(|&(s, e)| e > s)
                        .collect();

                    if overlaps.is_empty() {
                        *total_histogram.entry(0).or_insert(0) += a_len;
                    } else {
                        let (_, _, hist) = self.compute_coverage_events(a_start, a_end, &overlaps);
                        for (depth, count) in hist {
                            *total_histogram.entry(depth).or_insert(0) += count;
                        }
                    }
                } else {
                    *total_histogram.entry(0).or_insert(0) += a_len;
                }
            }
        }

        // Output sorted histogram
        let mut depths: Vec<_> = total_histogram.into_iter().collect();
        depths.sort_by_key(|&(d, _)| d);

        for (depth, count) in depths {
            // Use f32 to match bedtools precision (bedtools uses float internally)
            let fraction: f32 = count as f32 / total_length as f32;
            writeln!(
                output,
                "all\t{}\t{}\t{}\t{:.7}",
                depth, count, total_length, fraction
            )
            .map_err(BedError::Io)?;
        }

        Ok(())
    }

    /// Group records by chromosome, returning sorted records.
    fn group_records_by_chrom(records: Vec<BedRecord>) -> HashMap<String, Vec<BedRecord>> {
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

    fn make_record(chrom: &str, start: u64, end: u64) -> BedRecord {
        BedRecord::new(chrom, start, end)
    }

    #[test]
    fn test_basic_coverage() {
        let cmd = CoverageCommand::new();
        let a = vec![make_record("chr1", 100, 200)];
        let b = vec![make_record("chr1", 100, 150), make_record("chr1", 125, 175)];

        let mut output = Vec::new();
        let a_by_chrom = CoverageCommand::group_records_by_chrom(a);
        let b_by_chrom = CoverageCommand::group_records_by_chrom(b);

        cmd.coverage_chromosome_sweepline(
            a_by_chrom.get("chr1").unwrap(),
            b_by_chrom.get("chr1"),
            &mut output,
        );

        let result = String::from_utf8(output).unwrap();
        // Should have 2 overlaps, 75 bases covered (100-175)
        assert!(result.contains("2\t75\t100"));
    }

    #[test]
    fn test_coverage_no_overlap() {
        let cmd = CoverageCommand::new();
        let a = vec![make_record("chr1", 100, 200)];
        let b = vec![make_record("chr1", 300, 400)];

        let mut output = Vec::new();
        let a_by_chrom = CoverageCommand::group_records_by_chrom(a);
        let b_by_chrom = CoverageCommand::group_records_by_chrom(b);

        cmd.coverage_chromosome_sweepline(
            a_by_chrom.get("chr1").unwrap(),
            b_by_chrom.get("chr1"),
            &mut output,
        );

        let result = String::from_utf8(output).unwrap();
        assert!(result.contains("0\t0\t100\t0.0000000"));
    }

    #[test]
    fn test_mean_coverage() {
        let mut cmd = CoverageCommand::new();
        cmd.mean = true;

        let a = vec![make_record("chr1", 100, 200)];
        let b = vec![
            make_record("chr1", 100, 200), // Full coverage
            make_record("chr1", 100, 200), // Double coverage
        ];

        let mut output = Vec::new();
        let a_by_chrom = CoverageCommand::group_records_by_chrom(a);
        let b_by_chrom = CoverageCommand::group_records_by_chrom(b);

        cmd.coverage_chromosome_sweepline(
            a_by_chrom.get("chr1").unwrap(),
            b_by_chrom.get("chr1"),
            &mut output,
        );

        let result = String::from_utf8(output).unwrap();
        // Mean depth should be 2.0
        assert!(result.contains("2.0000000"));
    }

    #[test]
    fn test_event_based_coverage() {
        let cmd = CoverageCommand::new();

        // Test overlapping intervals
        let overlaps = vec![(100, 150), (125, 175)];
        let (bases_covered, total_depth, histogram) =
            cmd.compute_coverage_events(100, 200, &overlaps);

        // 100-125: depth 1 (25 bases)
        // 125-150: depth 2 (25 bases)
        // 150-175: depth 1 (25 bases)
        // 175-200: depth 0 (25 bases)
        assert_eq!(bases_covered, 75);
        assert_eq!(total_depth, 25 * 1 + 25 * 2 + 25 * 1 + 25 * 0); // 100
        assert_eq!(histogram.get(&0), Some(&25));
        assert_eq!(histogram.get(&1), Some(&50));
        assert_eq!(histogram.get(&2), Some(&25));
    }
}
