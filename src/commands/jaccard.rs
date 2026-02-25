//! Jaccard command implementation.
//!
//! Computes Jaccard similarity coefficient between two BED files.
//! Uses true streaming merge-sweep algorithm with O(k) memory.

use crate::bed::BedError;
use crate::streaming::parsing::{parse_bed3_bytes, should_skip_line};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

/// Format a float like C's %g: 6 significant figures, trailing zeros trimmed
fn format_g(val: f64) -> String {
    if val == 0.0 {
        return "0".to_string();
    }

    let abs_val = val.abs();

    let precision = if abs_val >= 1.0 {
        let digits_before_decimal = abs_val.log10().floor() as i32 + 1;
        (6 - digits_before_decimal).max(0) as usize
    } else {
        let leading_zeros = (-abs_val.log10()).floor() as usize;
        leading_zeros + 6
    };

    let formatted = format!("{:.prec$}", val, prec = precision);
    let trimmed = formatted.trim_end_matches('0').trim_end_matches('.');

    if trimmed.is_empty() || trimmed == "-" {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

/// Jaccard command configuration.
#[derive(Debug, Clone)]
pub struct JaccardCommand {
    pub strand: bool,
    pub fraction_a: Option<f64>,
    pub fraction_b: Option<f64>,
    pub reciprocal: bool,
}

impl Default for JaccardCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl JaccardCommand {
    pub fn new() -> Self {
        Self {
            strand: false,
            fraction_a: None,
            fraction_b: None,
            reciprocal: false,
        }
    }

    /// Run jaccard analysis between two files.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        input_a: P,
        input_b: P,
        output: &mut W,
    ) -> Result<(), BedError> {
        self.jaccard_streaming(input_a.as_ref(), input_b.as_ref(), output)
    }

    /// Read the next valid BED record from a buffered reader.
    /// Returns None if EOF, Some((chrom, start, end)) otherwise.
    #[inline]
    fn read_next_record(
        reader: &mut BufReader<File>,
        line_buf: &mut String,
    ) -> Result<Option<(Vec<u8>, u64, u64)>, BedError> {
        loop {
            line_buf.clear();
            let bytes_read = reader.read_line(line_buf)?;
            if bytes_read == 0 {
                return Ok(None);
            }

            let line = line_buf.trim_end();
            let line_bytes = line.as_bytes();

            if should_skip_line(line_bytes) {
                continue;
            }

            if let Some((chrom, start, end)) = parse_bed3_bytes(line_bytes) {
                return Ok(Some((chrom.to_vec(), start, end)));
            }
        }
    }

    /// True streaming jaccard implementation with O(k) memory.
    /// Uses active set tracking for overlapping intervals within each file.
    pub fn jaccard_streaming<W: Write>(
        &self,
        a_path: &Path,
        b_path: &Path,
        output: &mut W,
    ) -> Result<(), BedError> {
        let file_a = File::open(a_path)?;
        let file_b = File::open(b_path)?;

        let mut reader_a = BufReader::with_capacity(256 * 1024, file_a);
        let mut reader_b = BufReader::with_capacity(256 * 1024, file_b);

        let mut line_buf_a = String::with_capacity(1024);
        let mut line_buf_b = String::with_capacity(1024);

        // Pending intervals (current interval being processed from each file)
        // (chrom, start, end)
        let mut pending_a = Self::read_next_record(&mut reader_a, &mut line_buf_a)?;
        let mut pending_b = Self::read_next_record(&mut reader_b, &mut line_buf_b)?;

        // Active sets: store end positions of intervals that have started but not ended
        // For O(k) memory, we use a Vec sorted by end position
        let mut active_a: Vec<u64> = Vec::with_capacity(64);
        let mut active_b: Vec<u64> = Vec::with_capacity(64);

        // Global accumulators
        let mut total_intersection: u64 = 0;
        let mut total_union: u64 = 0;
        let mut total_n_intersections: u64 = 0;

        // Current chromosome being processed
        let mut current_chrom: Vec<u8> = Vec::new();

        // Sweep state
        let mut prev_pos: u64 = 0;
        let mut in_overlap = false;

        // Main event loop
        loop {
            // Determine the next event position and type
            // Events can be: start of A, end of A, start of B, end of B

            // Find minimum end position in active sets
            let min_end_a = active_a.first().copied();
            let min_end_b = active_b.first().copied();

            // Find start positions from pending intervals (if on current chromosome)
            let start_a = pending_a
                .as_ref()
                .filter(|(c, _, _)| *c == current_chrom)
                .map(|(_, s, _)| *s);
            let start_b = pending_b
                .as_ref()
                .filter(|(c, _, _)| *c == current_chrom)
                .map(|(_, s, _)| *s);

            // Check if we need to switch chromosomes
            let need_new_chrom = active_a.is_empty()
                && active_b.is_empty()
                && start_a.is_none()
                && start_b.is_none();

            if need_new_chrom {
                // Flush any remaining overlap from previous chromosome
                if in_overlap {
                    total_n_intersections += 1;
                    in_overlap = false;
                }

                // Find next chromosome to process
                let next_chrom = match (&pending_a, &pending_b) {
                    (None, None) => break, // All done
                    (Some((c, _, _)), None) => c.clone(),
                    (None, Some((c, _, _))) => c.clone(),
                    (Some((ca, _, _)), Some((cb, _, _))) => {
                        if ca <= cb {
                            ca.clone()
                        } else {
                            cb.clone()
                        }
                    }
                };

                current_chrom = next_chrom;
                prev_pos = 0;
                continue;
            }

            // Find the next event position (minimum of all possible events)
            // Event ordering: at same position, END before START (BED half-open semantics)

            // Collect all candidate events as (position, is_end, is_a)
            // is_end=true means END event, is_end=false means START event
            // For sorting: (position, is_end order, is_a order)
            // is_end order: END=0, START=1 (END before START)

            let mut next_pos = u64::MAX;
            let mut next_is_end = false;
            let mut next_is_a = false;

            // Check end events (priority 0, 1)
            if let Some(end_a) = min_end_a {
                if end_a < next_pos || (end_a == next_pos && !next_is_end) {
                    next_pos = end_a;
                    next_is_end = true;
                    next_is_a = true;
                }
            }
            if let Some(end_b) = min_end_b {
                if end_b < next_pos || (end_b == next_pos && !next_is_end) {
                    next_pos = end_b;
                    next_is_end = true;
                    next_is_a = false;
                } else if end_b == next_pos && next_is_end && next_is_a {
                    // EndA before EndB at same position - keep EndA
                }
            }

            // Check start events (priority 2, 3)
            if let Some(start) = start_a {
                if start < next_pos {
                    next_pos = start;
                    next_is_end = false;
                    next_is_a = true;
                } else if start == next_pos && !next_is_end {
                    // At same position as another START, prefer A
                    next_is_a = true;
                }
            }
            if let Some(start) = start_b {
                if start < next_pos {
                    next_pos = start;
                    next_is_end = false;
                    next_is_a = false;
                }
                // If same position as StartA, keep StartA (it comes first)
            }

            if next_pos == u64::MAX {
                break; // No more events
            }

            // Accumulate spans before processing this event
            if next_pos > prev_pos {
                let depth_a = active_a.len();
                let depth_b = active_b.len();

                // Check if we exited an overlap region
                if in_overlap && !(depth_a > 0 && depth_b > 0) {
                    total_n_intersections += 1;
                    in_overlap = false;
                }

                let span = next_pos - prev_pos;

                if depth_a > 0 && depth_b > 0 {
                    total_intersection += span;
                }
                if depth_a > 0 || depth_b > 0 {
                    total_union += span;
                }
            }

            // Process the event
            if next_is_end {
                // End event - remove from active set
                if next_is_a {
                    // Remove the first element (minimum end)
                    if !active_a.is_empty() {
                        active_a.remove(0);
                    }
                } else if !active_b.is_empty() {
                    active_b.remove(0);
                }
            } else {
                // Start event - add to active set and read next interval
                if next_is_a {
                    if let Some((_, _, end)) = pending_a.as_ref() {
                        // Insert end position maintaining sorted order
                        let end = *end;
                        let pos = active_a.partition_point(|&e| e < end);
                        active_a.insert(pos, end);
                    }
                    // Read next A interval
                    pending_a = Self::read_next_record(&mut reader_a, &mut line_buf_a)?;
                } else {
                    if let Some((_, _, end)) = pending_b.as_ref() {
                        let end = *end;
                        let pos = active_b.partition_point(|&e| e < end);
                        active_b.insert(pos, end);
                    }
                    // Read next B interval
                    pending_b = Self::read_next_record(&mut reader_b, &mut line_buf_b)?;
                }
            }

            // Enter overlap state when both have depth > 0
            if !active_a.is_empty() && !active_b.is_empty() {
                in_overlap = true;
            }

            prev_pos = next_pos;
        }

        // Final overlap count if still in one
        if in_overlap {
            total_n_intersections += 1;
        }

        // Compute Jaccard coefficient
        let jaccard = if total_union > 0 {
            total_intersection as f64 / total_union as f64
        } else {
            0.0
        };

        let jaccard_str = format_g(jaccard);

        writeln!(output, "intersection\tunion\tjaccard\tn_intersections")?;
        writeln!(
            output,
            "{}\t{}\t{}\t{}",
            total_intersection, total_union, jaccard_str, total_n_intersections
        )?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as IoWrite;
    use tempfile::NamedTempFile;

    fn create_temp_bed(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_jaccard_basic_streaming() {
        // A: [100, 200), [150, 250), [300, 400) - note overlapping A intervals
        // B: [120, 180), [350, 450)
        let a_content = "chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n";
        let b_content = "chr1\t120\t180\nchr1\t350\t450\n";

        let a_file = create_temp_bed(a_content);
        let b_file = create_temp_bed(b_content);

        let cmd = JaccardCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = output_str.lines().collect();

        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "intersection\tunion\tjaccard\tn_intersections");

        let parts: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(parts[0], "110"); // intersection
        assert_eq!(parts[1], "300"); // union
        assert_eq!(parts[3], "2"); // n_intersections
    }

    #[test]
    fn test_jaccard_no_overlap() {
        let a_content = "chr1\t100\t200\n";
        let b_content = "chr1\t300\t400\n";

        let a_file = create_temp_bed(a_content);
        let b_file = create_temp_bed(b_content);

        let cmd = JaccardCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = output_str.lines().collect();
        let parts: Vec<&str> = lines[1].split('\t').collect();

        assert_eq!(parts[0], "0"); // intersection
        assert_eq!(parts[1], "200"); // union = 100 + 100
        assert_eq!(parts[3], "0"); // n_intersections
    }

    #[test]
    fn test_jaccard_complete_overlap() {
        let a_content = "chr1\t100\t200\n";
        let b_content = "chr1\t100\t200\n";

        let a_file = create_temp_bed(a_content);
        let b_file = create_temp_bed(b_content);

        let cmd = JaccardCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = output_str.lines().collect();
        let parts: Vec<&str> = lines[1].split('\t').collect();

        assert_eq!(parts[0], "100"); // intersection
        assert_eq!(parts[1], "100"); // union
        assert_eq!(parts[2], "1"); // jaccard = 1.0
        assert_eq!(parts[3], "1"); // n_intersections
    }

    #[test]
    fn test_jaccard_empty_file() {
        let a_content = "chr1\t100\t200\n";
        let b_content = "";

        let a_file = create_temp_bed(a_content);
        let b_file = create_temp_bed(b_content);

        let cmd = JaccardCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = output_str.lines().collect();
        let parts: Vec<&str> = lines[1].split('\t').collect();

        assert_eq!(parts[0], "0"); // intersection
        assert_eq!(parts[1], "100"); // union = just A
        assert_eq!(parts[3], "0"); // n_intersections
    }

    #[test]
    fn test_jaccard_multiple_chromosomes() {
        let a_content = "chr1\t100\t200\nchr2\t100\t200\n";
        let b_content = "chr1\t150\t250\nchr2\t150\t250\n";

        let a_file = create_temp_bed(a_content);
        let b_file = create_temp_bed(b_content);

        let cmd = JaccardCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = output_str.lines().collect();
        let parts: Vec<&str> = lines[1].split('\t').collect();

        // Each chromosome: A=[100,200), B=[150,250)
        // Intersection per chrom: [150,200) = 50bp, total = 100bp
        // Union per chrom: [100,250) = 150bp, total = 300bp
        assert_eq!(parts[0], "100"); // intersection
        assert_eq!(parts[1], "300"); // union
        assert_eq!(parts[3], "2"); // n_intersections (one per chromosome)
    }

    #[test]
    fn test_jaccard_disjoint_chromosomes() {
        let a_content = "chr1\t100\t200\n";
        let b_content = "chr2\t100\t200\n";

        let a_file = create_temp_bed(a_content);
        let b_file = create_temp_bed(b_content);

        let cmd = JaccardCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = output_str.lines().collect();
        let parts: Vec<&str> = lines[1].split('\t').collect();

        assert_eq!(parts[0], "0"); // intersection
        assert_eq!(parts[1], "200"); // union = 100 + 100
        assert_eq!(parts[3], "0"); // n_intersections
    }

    #[test]
    fn test_jaccard_nested_intervals() {
        // B is completely inside A
        let a_content = "chr1\t100\t400\n";
        let b_content = "chr1\t150\t250\n";

        let a_file = create_temp_bed(a_content);
        let b_file = create_temp_bed(b_content);

        let cmd = JaccardCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = output_str.lines().collect();
        let parts: Vec<&str> = lines[1].split('\t').collect();

        // Intersection: [150, 250) = 100bp
        // Union: [100, 400) = 300bp
        assert_eq!(parts[0], "100"); // intersection
        assert_eq!(parts[1], "300"); // union
        assert_eq!(parts[3], "1"); // n_intersections
    }

    #[test]
    fn test_jaccard_back_to_back() {
        // Intervals that touch but don't overlap (BED half-open semantics)
        let a_content = "chr1\t100\t200\n";
        let b_content = "chr1\t200\t300\n";

        let a_file = create_temp_bed(a_content);
        let b_file = create_temp_bed(b_content);

        let cmd = JaccardCommand::new();
        let mut output = Vec::new();
        cmd.run(a_file.path(), b_file.path(), &mut output).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = output_str.lines().collect();
        let parts: Vec<&str> = lines[1].split('\t').collect();

        // No overlap - A ends at 200, B starts at 200, half-open means no intersection
        assert_eq!(parts[0], "0"); // intersection
        assert_eq!(parts[1], "200"); // union = 100 + 100
        assert_eq!(parts[3], "0"); // n_intersections
    }
}
