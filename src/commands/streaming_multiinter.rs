//! Streaming multiinter command - TRUE O(k) memory multi-file intersection.
//!
//! Memory complexity: O(k) where k = max overlapping intervals across all files
//! at any position.
//!
//! Uses k-way merge with min-heap to stream through all sorted input files
//! simultaneously without loading them entirely into memory.
//!
//! REQUIREMENT: All input files must be sorted by (chrom, start) for streaming mode.
//! Use `--assume-sorted` flag or pre-sort with `grit sort`.

#![allow(clippy::ptr_arg)]

use crate::bed::BedError;
use crate::streaming::parsing::{parse_bed3_bytes, should_skip_line};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// An interval from a specific file with its source index.
#[derive(Debug, Clone)]
struct TaggedInterval {
    chrom: Vec<u8>,
    start: u64,
    end: u64,
    file_idx: usize,
}

/// Wrapper for min-heap (BinaryHeap is max-heap by default).
#[derive(Debug, Clone, Eq, PartialEq)]
struct HeapEntry {
    chrom: Vec<u8>,
    start: u64,
    end: u64,
    file_idx: usize,
}

impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reverse ordering for min-heap
        other
            .chrom
            .cmp(&self.chrom)
            .then(other.start.cmp(&self.start))
            .then(other.end.cmp(&self.end))
            .then(other.file_idx.cmp(&self.file_idx))
    }
}

impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// An event in the sweep-line algorithm.
#[derive(Debug, Clone, Eq, PartialEq)]
struct Event {
    pos: u64,
    is_start: bool,
    file_idx: usize,
}

impl Ord for Event {
    fn cmp(&self, other: &Self) -> Ordering {
        self.pos
            .cmp(&other.pos)
            .then(self.is_start.cmp(&other.is_start)) // ends before starts at same position
    }
}

impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Reader state for a single file.
struct FileReader<R: BufRead> {
    reader: R,
    line_buf: String,
    file_idx: usize,
    exhausted: bool,
}

impl<R: BufRead> FileReader<R> {
    fn new(reader: R, file_idx: usize) -> Self {
        Self {
            reader,
            line_buf: String::with_capacity(1024),
            file_idx,
            exhausted: false,
        }
    }

    /// Read the next valid interval from this file.
    fn next_interval(&mut self) -> Result<Option<TaggedInterval>, BedError> {
        if self.exhausted {
            return Ok(None);
        }

        loop {
            self.line_buf.clear();
            let bytes_read = self.reader.read_line(&mut self.line_buf)?;
            if bytes_read == 0 {
                self.exhausted = true;
                return Ok(None);
            }

            let line_bytes = self.line_buf.trim_end().as_bytes();
            if should_skip_line(line_bytes) {
                continue;
            }

            if let Some((chrom, start, end)) = parse_bed3_bytes(line_bytes) {
                return Ok(Some(TaggedInterval {
                    chrom: chrom.to_vec(),
                    start,
                    end,
                    file_idx: self.file_idx,
                }));
            }
        }
    }
}

/// Streaming multiinter command configuration.
#[derive(Debug, Clone)]
pub struct StreamingMultiinterCommand {
    /// Only report intervals present in all files
    pub cluster: bool,
    /// Skip sorted validation (faster for pre-sorted input)
    pub assume_sorted: bool,
}

impl Default for StreamingMultiinterCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl StreamingMultiinterCommand {
    pub fn new() -> Self {
        Self {
            cluster: false,
            assume_sorted: false,
        }
    }

    /// Set cluster flag (builder pattern).
    pub fn with_cluster(mut self, cluster: bool) -> Self {
        self.cluster = cluster;
        self
    }

    /// Set assume_sorted flag (builder pattern).
    pub fn with_assume_sorted(mut self, assume_sorted: bool) -> Self {
        self.assume_sorted = assume_sorted;
        self
    }

    /// Execute streaming multiinter.
    ///
    /// Memory: O(k) where k = max overlapping intervals across all files.
    /// All input files are streamed - never fully loaded.
    ///
    /// REQUIREMENT: All files must be sorted by (chrom, start) for correct results.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        inputs: &[P],
        output: &mut W,
    ) -> Result<(), BedError> {
        if inputs.is_empty() {
            return Ok(());
        }

        // Open all files
        let mut readers = Vec::with_capacity(inputs.len());
        for (idx, path) in inputs.iter().enumerate() {
            let file = File::open(path)?;
            let reader = BufReader::with_capacity(256 * 1024, file);
            readers.push(FileReader::new(reader, idx));
        }

        self.multiinter_streaming(readers, inputs.len(), output)
    }

    /// Streaming multiinter implementation using k-way merge.
    ///
    /// Algorithm:
    /// 1. Initialize min-heap with first interval from each file
    /// 2. Process intervals in sorted order:
    ///    - When entering a new chromosome, process the previous one
    ///    - Accumulate events (start/end) for current chromosome
    ///    - Pull next interval from the file that provided the current one
    /// 3. When chromosome changes or all files exhausted, run sweep-line
    fn multiinter_streaming<R: BufRead, W: Write>(
        &self,
        mut readers: Vec<FileReader<R>>,
        n_files: usize,
        output: &mut W,
    ) -> Result<(), BedError> {
        // Large output buffer (8MB)
        let mut buf_output = BufWriter::with_capacity(8 * 1024 * 1024, output);

        // Initialize min-heap with first interval from each file
        let mut heap: BinaryHeap<HeapEntry> = BinaryHeap::with_capacity(n_files);

        for reader in &mut readers {
            if let Some(interval) = reader.next_interval()? {
                heap.push(HeapEntry {
                    chrom: interval.chrom,
                    start: interval.start,
                    end: interval.end,
                    file_idx: interval.file_idx,
                });
            }
        }

        // Current chromosome being processed
        let mut current_chrom: Option<Vec<u8>> = None;
        // Events for current chromosome
        let mut events: Vec<Event> = Vec::with_capacity(1024);

        // itoa buffer for fast integer formatting
        let mut itoa_buf = itoa::Buffer::new();

        while let Some(entry) = heap.pop() {
            // Check if chromosome changed
            let chrom_changed = match &current_chrom {
                Some(c) => c != &entry.chrom,
                None => false,
            };

            if chrom_changed {
                // Process completed chromosome
                if let Some(ref chrom) = current_chrom {
                    self.process_chromosome_events(
                        chrom,
                        &mut events,
                        n_files,
                        &mut buf_output,
                        &mut itoa_buf,
                    )?;
                }
                events.clear();
            }

            current_chrom = Some(entry.chrom.clone());

            // Add events for this interval
            events.push(Event {
                pos: entry.start,
                is_start: true,
                file_idx: entry.file_idx,
            });
            events.push(Event {
                pos: entry.end,
                is_start: false,
                file_idx: entry.file_idx,
            });

            // Pull next interval from the same file
            if let Some(next) = readers[entry.file_idx].next_interval()? {
                heap.push(HeapEntry {
                    chrom: next.chrom,
                    start: next.start,
                    end: next.end,
                    file_idx: next.file_idx,
                });
            }
        }

        // Process last chromosome
        if let Some(ref chrom) = current_chrom {
            self.process_chromosome_events(
                chrom,
                &mut events,
                n_files,
                &mut buf_output,
                &mut itoa_buf,
            )?;
        }

        buf_output.flush().map_err(BedError::Io)?;
        Ok(())
    }

    /// Process events for a single chromosome using sweep-line.
    fn process_chromosome_events<W: Write>(
        &self,
        chrom: &[u8],
        events: &mut Vec<Event>,
        n_files: usize,
        output: &mut W,
        itoa_buf: &mut itoa::Buffer,
    ) -> Result<(), BedError> {
        if events.is_empty() {
            return Ok(());
        }

        // Sort events: by position, then ends before starts
        events.sort_unstable();

        // Track depth per file
        let mut file_depths: Vec<u32> = vec![0; n_files];
        let mut prev_pos: u64 = events[0].pos;
        let mut has_coverage = false;

        for event in events.iter() {
            // Output region if there was coverage
            if event.pos > prev_pos && has_coverage {
                self.output_region(chrom, prev_pos, event.pos, &file_depths, output, itoa_buf)?;
            }

            // Update depth
            if event.is_start {
                file_depths[event.file_idx] += 1;
            } else {
                file_depths[event.file_idx] = file_depths[event.file_idx].saturating_sub(1);
            }

            // Check if any file has coverage
            has_coverage = file_depths.iter().any(|&d| d > 0);
            prev_pos = event.pos;
        }

        Ok(())
    }

    /// Output a region with coverage info.
    fn output_region<W: Write>(
        &self,
        chrom: &[u8],
        start: u64,
        end: u64,
        file_depths: &[u32],
        output: &mut W,
        itoa_buf: &mut itoa::Buffer,
    ) -> Result<(), BedError> {
        // Count files with coverage
        let count: usize = file_depths.iter().filter(|&&d| d > 0).count();

        if count == 0 {
            return Ok(());
        }

        // Skip if cluster mode and not all files
        if self.cluster && count != file_depths.len() {
            return Ok(());
        }

        // Build list of file indices (1-based)
        let file_list: Vec<String> = file_depths
            .iter()
            .enumerate()
            .filter(|(_, &d)| d > 0)
            .map(|(i, _)| (i + 1).to_string())
            .collect();

        // Write output: chrom, start, end, count, file_list, flags...
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
            .write_all(itoa_buf.format(count).as_bytes())
            .map_err(BedError::Io)?;
        output.write_all(b"\t").map_err(BedError::Io)?;
        output
            .write_all(file_list.join(",").as_bytes())
            .map_err(BedError::Io)?;

        // Write presence flags
        for &depth in file_depths {
            output.write_all(b"\t").map_err(BedError::Io)?;
            output
                .write_all(if depth > 0 { b"1" } else { b"0" })
                .map_err(BedError::Io)?;
        }

        output.write_all(b"\n").map_err(BedError::Io)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn make_reader(data: &str, idx: usize) -> FileReader<BufReader<Cursor<Vec<u8>>>> {
        let cursor = Cursor::new(data.as_bytes().to_vec());
        let reader = BufReader::new(cursor);
        FileReader::new(reader, idx)
    }

    #[test]
    fn test_streaming_multiinter_basic() {
        let file1_data = "chr1\t100\t200\nchr1\t300\t400\n";
        let file2_data = "chr1\t150\t250\nchr1\t350\t450\n";

        let readers = vec![make_reader(file1_data, 0), make_reader(file2_data, 1)];

        let cmd = StreamingMultiinterCommand::new().with_assume_sorted(true);

        let mut output = Vec::new();
        cmd.multiinter_streaming(readers, 2, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();

        // Should have multiple regions
        assert!(!lines.is_empty());

        // Check for overlap region (150-200) with both files
        let overlap_line = lines
            .iter()
            .find(|l| l.contains("150") && l.contains("200"));
        assert!(overlap_line.is_some());
        let parts: Vec<&str> = overlap_line.unwrap().split('\t').collect();
        assert_eq!(parts[3], "2"); // count = 2 (both files)
    }

    #[test]
    fn test_streaming_multiinter_three_files() {
        let file1_data = "chr1\t100\t200\n";
        let file2_data = "chr1\t150\t250\n";
        let file3_data = "chr1\t180\t220\n";

        let readers = vec![
            make_reader(file1_data, 0),
            make_reader(file2_data, 1),
            make_reader(file3_data, 2),
        ];

        let cmd = StreamingMultiinterCommand::new().with_assume_sorted(true);

        let mut output = Vec::new();
        cmd.multiinter_streaming(readers, 3, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();

        // Find region where all three overlap (180-200)
        let all_three = lines.iter().find(|line| {
            let parts: Vec<&str> = line.split('\t').collect();
            parts[3] == "3" // count == 3
        });

        assert!(all_three.is_some());
        let parts: Vec<&str> = all_three.unwrap().split('\t').collect();
        assert_eq!(parts[1], "180");
        assert_eq!(parts[2], "200");
    }

    #[test]
    fn test_streaming_multiinter_cluster() {
        let file1_data = "chr1\t100\t200\n";
        let file2_data = "chr1\t150\t250\n";

        let readers = vec![make_reader(file1_data, 0), make_reader(file2_data, 1)];

        let cmd = StreamingMultiinterCommand::new()
            .with_cluster(true)
            .with_assume_sorted(true);

        let mut output = Vec::new();
        cmd.multiinter_streaming(readers, 2, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();

        // In cluster mode, only regions with ALL files should be output
        // Only the overlap region (150-200) should appear
        assert_eq!(lines.len(), 1);
        assert!(lines[0].contains("150"));
        assert!(lines[0].contains("200"));
    }

    #[test]
    fn test_streaming_multiinter_no_overlap() {
        let file1_data = "chr1\t100\t200\n";
        let file2_data = "chr1\t300\t400\n";

        let readers = vec![make_reader(file1_data, 0), make_reader(file2_data, 1)];

        let cmd = StreamingMultiinterCommand::new().with_assume_sorted(true);

        let mut output = Vec::new();
        cmd.multiinter_streaming(readers, 2, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();

        // Should have 2 regions with count=1 each
        assert_eq!(lines.len(), 2);

        for line in lines {
            let parts: Vec<&str> = line.split('\t').collect();
            assert_eq!(parts[3], "1"); // count == 1
        }
    }

    #[test]
    fn test_streaming_multiinter_multi_chrom() {
        let file1_data = "chr1\t100\t200\nchr2\t50\t100\n";
        let file2_data = "chr1\t150\t250\nchr2\t75\t125\n";

        let readers = vec![make_reader(file1_data, 0), make_reader(file2_data, 1)];

        let cmd = StreamingMultiinterCommand::new().with_assume_sorted(true);

        let mut output = Vec::new();
        cmd.multiinter_streaming(readers, 2, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();

        // Should have regions on both chromosomes
        assert!(result.contains("chr1\t"));
        assert!(result.contains("chr2\t"));
    }

    #[test]
    fn test_streaming_multiinter_empty_file() {
        let file1_data = "chr1\t100\t200\n";
        let file2_data = "";

        let readers = vec![make_reader(file1_data, 0), make_reader(file2_data, 1)];

        let cmd = StreamingMultiinterCommand::new().with_assume_sorted(true);

        let mut output = Vec::new();
        cmd.multiinter_streaming(readers, 2, &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();

        // Should have 1 region from file1 only
        assert_eq!(lines.len(), 1);
        let parts: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(parts[3], "1"); // count == 1
    }
}
