//! Ultra-fast BED file sorting using radix sort and memory-mapped I/O.
//!
//! Optimizations:
//! - Memory-mapped file I/O for zero-copy access
//! - Zero-allocation parsing (byte slices, no String allocation)
//! - LSD Radix Sort for (chrom_index, start, end) - O(n)
//! - Parallel parsing with Rayon
//! - Fast integer parsing and output
//!
//! Sort order (matches `LC_ALL=C sort -k1,1 -k2,2n -k3,3n`):
//! 1. Chromosome (lexicographic order)
//! 2. Start coordinate (ascending, numeric)
//! 3. End coordinate (ascending, numeric)
//! 4. Input order preserved for ties (stable sort)

use crate::bed::BedError;
use memchr::memchr;
use memmap2::Mmap;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

/// Buffer size for I/O operations (256KB for better throughput)
const BUF_SIZE: usize = 256 * 1024;

/// Minimum file size to use mmap (smaller files use buffered I/O)
const MMAP_THRESHOLD: usize = 64 * 1024;

/// Minimum records to trigger parallel parsing
const PARALLEL_THRESHOLD: usize = 10_000;

/// Minimum records to use radix sort (smaller uses comparison sort)
const RADIX_THRESHOLD: usize = 256;

/// Entry for sorting - designed for cache-efficient radix sort.
/// Layout optimized for sequential memory access during radix passes.
#[derive(Clone, Copy, Debug)]
#[repr(C)]
struct SortEntry {
    /// Chromosome index (from lexicographic ordering)
    chrom_index: u16,
    /// Padding for alignment
    _padding: u16,
    /// Start coordinate (u32 is sufficient for genomic coordinates)
    start: u32,
    /// End coordinate
    end: u32,
    /// Offset into the data buffer where this line starts
    line_start: u32,
    /// Length of the line (excluding newline)
    line_len: u32,
}

/// Statistics from fast sort operation.
#[derive(Debug, Default, Clone)]
pub struct FastSortStats {
    pub records_read: usize,
    pub unique_chroms: usize,
    pub used_radix_sort: bool,
    pub used_mmap: bool,
}

impl std::fmt::Display for FastSortStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Records: {}, Chroms: {}, Radix: {}, Mmap: {}",
            self.records_read,
            self.unique_chroms,
            if self.used_radix_sort { "yes" } else { "no" },
            if self.used_mmap { "yes" } else { "no" }
        )
    }
}

/// Fast sort command with optimized algorithms.
#[derive(Debug, Clone)]
pub struct FastSortCommand {
    /// Use radix sort (default: true for large files)
    pub use_radix: bool,
    /// Reverse sort order
    pub reverse: bool,
    /// Genome-based chromosome ordering (chrom bytes -> index)
    genome_order: Option<HashMap<Vec<u8>, u16>>,
}

impl Default for FastSortCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl FastSortCommand {
    pub fn new() -> Self {
        Self {
            use_radix: true,
            reverse: false,
            genome_order: None,
        }
    }

    /// Set genome-based chromosome ordering.
    /// Chromosomes will be sorted in the order they appear in the genome file.
    /// Unknown chromosomes are placed after all known chromosomes.
    pub fn with_genome(mut self, genome: &crate::genome::Genome) -> Self {
        let order: HashMap<Vec<u8>, u16> = genome
            .chromosomes()
            .enumerate()
            .map(|(i, chrom)| (chrom.as_bytes().to_vec(), i as u16))
            .collect();
        self.genome_order = Some(order);
        self
    }

    /// Run fast sort on a file.
    pub fn run<P: AsRef<Path>, W: Write>(
        &self,
        input_path: P,
        output: &mut W,
    ) -> Result<FastSortStats, BedError> {
        let path = input_path.as_ref();
        let file = File::open(path)?;
        let metadata = file.metadata()?;
        let file_size = metadata.len() as usize;

        if file_size >= MMAP_THRESHOLD {
            // Use memory-mapped I/O for large files
            let mmap = unsafe { Mmap::map(&file)? };
            self.sort_mmap(&mmap, output)
        } else {
            // Use buffered I/O for small files
            self.sort_buffered(file, output)
        }
    }

    /// Run fast sort from stdin.
    pub fn run_stdin<W: Write>(&self, output: &mut W) -> Result<FastSortStats, BedError> {
        let stdin = io::stdin();
        self.sort_buffered(stdin.lock(), output)
    }

    /// Sort using memory-mapped file (zero-copy).
    fn sort_mmap<W: Write>(&self, data: &[u8], output: &mut W) -> Result<FastSortStats, BedError> {
        let mut stats = FastSortStats {
            used_mmap: true,
            ..Default::default()
        };

        // Phase 1: Scan for line boundaries and count records
        let line_offsets = find_line_offsets(data);
        let num_records = line_offsets.len();

        if num_records == 0 {
            return Ok(stats);
        }

        // Phase 2: Build chromosome index (genome order if provided, else lexicographic)
        let chrom_index = build_chrom_index(data, &line_offsets, self.genome_order.as_ref());
        stats.unique_chroms = chrom_index.len();

        // Phase 3: Parse all records into sort entries (parallel for large files)
        let entries = if num_records >= PARALLEL_THRESHOLD {
            parse_entries_parallel(data, &line_offsets, &chrom_index)
        } else {
            parse_entries_sequential(data, &line_offsets, &chrom_index)
        };

        stats.records_read = entries.len();

        // Phase 4: Sort using LSD radix sort or comparison sort
        let sorted_entries = if self.use_radix && entries.len() >= RADIX_THRESHOLD {
            stats.used_radix_sort = true;
            radix_sort_lsd(entries)
        } else {
            comparison_sort_entries(entries)
        };

        // Phase 5: Output sorted records
        let mut writer = BufWriter::with_capacity(BUF_SIZE, output);
        if self.reverse {
            for entry in sorted_entries.iter().rev() {
                let start = entry.line_start as usize;
                let end = start + entry.line_len as usize;
                writer.write_all(&data[start..end])?;
                writer.write_all(b"\n")?;
            }
        } else {
            for entry in &sorted_entries {
                let start = entry.line_start as usize;
                let end = start + entry.line_len as usize;
                writer.write_all(&data[start..end])?;
                writer.write_all(b"\n")?;
            }
        }
        writer.flush()?;

        Ok(stats)
    }

    /// Sort using buffered I/O (for stdin or small files).
    pub fn sort_buffered<R: Read, W: Write>(
        &self,
        mut reader: R,
        output: &mut W,
    ) -> Result<FastSortStats, BedError> {
        let mut stats = FastSortStats::default();

        // Read all data into memory
        let mut data = Vec::with_capacity(BUF_SIZE);
        reader.read_to_end(&mut data)?;

        if data.is_empty() {
            return Ok(stats);
        }

        // Use the same algorithm as mmap but on owned Vec
        let line_offsets = find_line_offsets(&data);
        let num_records = line_offsets.len();

        if num_records == 0 {
            return Ok(stats);
        }

        let chrom_index = build_chrom_index(&data, &line_offsets, self.genome_order.as_ref());
        stats.unique_chroms = chrom_index.len();

        let entries = if num_records >= PARALLEL_THRESHOLD {
            parse_entries_parallel(&data, &line_offsets, &chrom_index)
        } else {
            parse_entries_sequential(&data, &line_offsets, &chrom_index)
        };

        stats.records_read = entries.len();

        let sorted_entries = if self.use_radix && entries.len() >= RADIX_THRESHOLD {
            stats.used_radix_sort = true;
            radix_sort_lsd(entries)
        } else {
            comparison_sort_entries(entries)
        };

        let mut writer = BufWriter::with_capacity(BUF_SIZE, output);
        if self.reverse {
            for entry in sorted_entries.iter().rev() {
                let start = entry.line_start as usize;
                let end = start + entry.line_len as usize;
                writer.write_all(&data[start..end])?;
                writer.write_all(b"\n")?;
            }
        } else {
            for entry in &sorted_entries {
                let start = entry.line_start as usize;
                let end = start + entry.line_len as usize;
                writer.write_all(&data[start..end])?;
                writer.write_all(b"\n")?;
            }
        }
        writer.flush()?;

        Ok(stats)
    }
}

/// Find all line start/end offsets in the data.
/// Returns Vec of (line_start, line_end) positions.
fn find_line_offsets(data: &[u8]) -> Vec<(usize, usize)> {
    let mut offsets = Vec::with_capacity(data.len() / 50); // Estimate ~50 bytes per line
    let mut pos = 0;

    while pos < data.len() {
        // Skip leading whitespace/empty lines
        while pos < data.len() && (data[pos] == b'\n' || data[pos] == b'\r') {
            pos += 1;
        }

        if pos >= data.len() {
            break;
        }

        let line_start = pos;

        // Find end of line
        if let Some(newline_offset) = memchr(b'\n', &data[pos..]) {
            let mut line_end = pos + newline_offset;
            // Handle \r\n
            if line_end > line_start && data[line_end - 1] == b'\r' {
                line_end -= 1;
            }

            // Skip headers and empty lines
            if line_end > line_start
                && data[line_start] != b'#'
                && !data[line_start..].starts_with(b"track")
                && !data[line_start..].starts_with(b"browser")
            {
                offsets.push((line_start, line_end));
            }

            pos += newline_offset + 1;
        } else {
            // Last line without newline
            let line_end = data.len();
            if line_end > line_start
                && data[line_start] != b'#'
                && !data[line_start..].starts_with(b"track")
                && !data[line_start..].starts_with(b"browser")
            {
                offsets.push((line_start, line_end));
            }
            break;
        }
    }

    offsets
}

/// Build chromosome to index mapping.
/// If genome_order is provided, use that order (unknown chroms placed at end).
/// Otherwise, sort chromosomes lexicographically (matching `sort -k1,1`).
fn build_chrom_index(
    data: &[u8],
    line_offsets: &[(usize, usize)],
    genome_order: Option<&HashMap<Vec<u8>, u16>>,
) -> HashMap<Vec<u8>, u16> {
    let mut chroms: Vec<Vec<u8>> = Vec::new();

    for &(start, end) in line_offsets {
        let line = &data[start..end];
        if let Some(tab_pos) = memchr(b'\t', line) {
            let chrom = &line[..tab_pos];
            // Use a simple linear search for small number of chromosomes
            if !chroms.iter().any(|c| c.as_slice() == chrom) {
                chroms.push(chrom.to_vec());
            }
        }
    }

    if let Some(genome_order) = genome_order {
        // Use genome file order: known chroms first, then unknown chroms lexicographically
        let max_known = genome_order
            .values()
            .max()
            .map(|v| *v as usize)
            .unwrap_or(0);

        // Separate known and unknown chromosomes
        let known: Vec<_> = chroms
            .iter()
            .filter(|c| genome_order.contains_key(c.as_slice()))
            .cloned()
            .collect();
        let mut unknown: Vec<_> = chroms
            .iter()
            .filter(|c| !genome_order.contains_key(c.as_slice()))
            .cloned()
            .collect();

        // Sort unknown chromosomes lexicographically for deterministic ordering
        unknown.sort();

        // Build the final index: known chroms use genome order, unknown chroms follow
        let mut result = HashMap::new();
        for chrom in known {
            if let Some(&idx) = genome_order.get(&chrom) {
                result.insert(chrom, idx);
            }
        }
        for (i, chrom) in unknown.into_iter().enumerate() {
            result.insert(chrom, (max_known + 1 + i) as u16);
        }
        result
    } else {
        // Sort chromosomes lexicographically (matching `sort -k1,1`)
        chroms.sort();

        chroms
            .into_iter()
            .enumerate()
            .map(|(i, c)| (c, i as u16))
            .collect()
    }
}

/// Parse BED3 fields from a line slice.
/// Returns (chrom, start, end) or None if invalid.
#[inline(always)]
fn parse_bed3(line: &[u8]) -> Option<(&[u8], u32, u32)> {
    let tab1 = memchr(b'\t', line)?;
    let chrom = &line[..tab1];

    let rest1 = &line[tab1 + 1..];
    let tab2 = memchr(b'\t', rest1)?;
    let start_bytes = &rest1[..tab2];

    let rest2 = &rest1[tab2 + 1..];
    let end_bytes = if let Some(tab3) = memchr(b'\t', rest2) {
        &rest2[..tab3]
    } else {
        rest2
    };

    let start = parse_u32_fast(start_bytes)?;
    let end = parse_u32_fast(end_bytes)?;

    Some((chrom, start, end))
}

/// Fast u32 parsing without allocation.
#[inline(always)]
fn parse_u32_fast(bytes: &[u8]) -> Option<u32> {
    if bytes.is_empty() {
        return None;
    }

    let mut result: u32 = 0;
    for &b in bytes {
        let digit = b.wrapping_sub(b'0');
        if digit > 9 {
            return None;
        }
        result = result.wrapping_mul(10).wrapping_add(digit as u32);
    }
    Some(result)
}

/// Parse entries sequentially.
fn parse_entries_sequential(
    data: &[u8],
    line_offsets: &[(usize, usize)],
    chrom_index: &HashMap<Vec<u8>, u16>,
) -> Vec<SortEntry> {
    let mut entries = Vec::with_capacity(line_offsets.len());

    for &(line_start, line_end) in line_offsets {
        let line = &data[line_start..line_end];
        if let Some((chrom, start, end)) = parse_bed3(line) {
            if let Some(&chrom_idx) = chrom_index.get(chrom) {
                entries.push(SortEntry {
                    chrom_index: chrom_idx,
                    _padding: 0,
                    start,
                    end,
                    line_start: line_start as u32,
                    line_len: (line_end - line_start) as u32,
                });
            }
        }
    }

    entries
}

/// Parse entries in parallel chunks.
fn parse_entries_parallel(
    data: &[u8],
    line_offsets: &[(usize, usize)],
    chrom_index: &HashMap<Vec<u8>, u16>,
) -> Vec<SortEntry> {
    // Determine chunk size for parallel processing
    let num_threads = rayon::current_num_threads();
    let chunk_size = (line_offsets.len() / num_threads).max(1000);

    line_offsets
        .par_chunks(chunk_size)
        .flat_map(|chunk| {
            let mut entries = Vec::with_capacity(chunk.len());
            for &(line_start, line_end) in chunk {
                let line = &data[line_start..line_end];
                if let Some((chrom, start, end)) = parse_bed3(line) {
                    if let Some(&chrom_idx) = chrom_index.get(chrom) {
                        entries.push(SortEntry {
                            chrom_index: chrom_idx,
                            _padding: 0,
                            start,
                            end,
                            line_start: line_start as u32,
                            line_len: (line_end - line_start) as u32,
                        });
                    }
                }
            }
            entries
        })
        .collect()
}

/// Comparison-based stable sort (for smaller datasets).
/// Sorts by (chrom, start, end), preserves input order for ties.
fn comparison_sort_entries(mut entries: Vec<SortEntry>) -> Vec<SortEntry> {
    entries.sort_by(|a, b| {
        a.chrom_index
            .cmp(&b.chrom_index)
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| a.end.cmp(&b.end))
            .then_with(|| a.line_start.cmp(&b.line_start))
    });
    entries
}

/// LSD Radix Sort for SortEntry.
///
/// Sorts by (chrom_index, start, end, line_start) using Least Significant Digit first.
/// This ensures stable sorting matching `LC_ALL=C sort -k1,1 -k2,2n -k3,3n`:
/// 1. Primary: chromosome (lexicographic index)
/// 2. Secondary: start coordinate
/// 3. Tertiary: end coordinate
/// 4. Quaternary: input order (line_start) for stable tie-breaking
///
/// Radix passes (8-bit radix, 256 buckets):
/// - Passes 1-4: line_start bytes 0-3 (least significant, for stability)
/// - Passes 5-8: end bytes 0-3
/// - Passes 9-12: start bytes 0-3
/// - Passes 13-14: chrom_index bytes 0-1 (most significant)
///
/// Total: 14 passes max, optimized by skipping passes where all values have same byte.
fn radix_sort_lsd(entries: Vec<SortEntry>) -> Vec<SortEntry> {
    if entries.len() < RADIX_THRESHOLD {
        return comparison_sort_entries(entries);
    }

    let n = entries.len();
    let mut src = entries;
    let mut dst = vec![
        SortEntry {
            chrom_index: 0,
            _padding: 0,
            start: 0,
            end: 0,
            line_start: 0,
            line_len: 0,
        };
        n
    ];

    // LSD radix sort: process from least significant to most significant
    // Order: line_start -> end -> start -> chrom_index

    // Pass 1-4: Sort by line_start (for deterministic ordering of identical records)
    for shift in (0u32..32).step_by(8) {
        if !radix_pass_line_start(&mut src, &mut dst, shift) {
            // All bytes were same, skip swap
            continue;
        }
        std::mem::swap(&mut src, &mut dst);
    }

    // Pass 5-8: Sort by end coordinate
    for shift in (0u32..32).step_by(8) {
        if !radix_pass_end(&mut src, &mut dst, shift) {
            continue;
        }
        std::mem::swap(&mut src, &mut dst);
    }

    // Pass 9-12: Sort by start coordinate
    for shift in (0u32..32).step_by(8) {
        if !radix_pass_start(&mut src, &mut dst, shift) {
            continue;
        }
        std::mem::swap(&mut src, &mut dst);
    }

    // Pass 13-14: Sort by chrom_index (most significant)
    for shift in (0u32..16).step_by(8) {
        if !radix_pass_chrom(&mut src, &mut dst, shift) {
            continue;
        }
        std::mem::swap(&mut src, &mut dst);
    }

    src
}

/// Single radix pass by line_start field. Returns false if all bytes are same (can skip).
#[inline]
fn radix_pass_line_start(src: &mut [SortEntry], dst: &mut [SortEntry], shift: u32) -> bool {
    let mut count = [0usize; 257];

    // Count occurrences
    for entry in src.iter() {
        let byte = ((entry.line_start >> shift) & 0xFF) as usize;
        count[byte + 1] += 1;
    }

    // Check if all values have same byte (can skip this pass)
    let mut non_zero_buckets = 0;
    for &c in &count[1..] {
        if c > 0 {
            non_zero_buckets += 1;
        }
    }
    if non_zero_buckets <= 1 {
        return false;
    }

    // Compute prefix sums
    for i in 1..257 {
        count[i] += count[i - 1];
    }

    // Distribute to destination
    for entry in src.iter() {
        let byte = ((entry.line_start >> shift) & 0xFF) as usize;
        dst[count[byte]] = *entry;
        count[byte] += 1;
    }

    true
}

/// Single radix pass by start field. Returns false if all bytes are same (can skip).
#[inline]
fn radix_pass_start(src: &mut [SortEntry], dst: &mut [SortEntry], shift: u32) -> bool {
    let mut count = [0usize; 257];

    for entry in src.iter() {
        let byte = ((entry.start >> shift) & 0xFF) as usize;
        count[byte + 1] += 1;
    }

    let mut non_zero_buckets = 0;
    for &c in &count[1..] {
        if c > 0 {
            non_zero_buckets += 1;
        }
    }
    if non_zero_buckets <= 1 {
        return false;
    }

    for i in 1..257 {
        count[i] += count[i - 1];
    }

    for entry in src.iter() {
        let byte = ((entry.start >> shift) & 0xFF) as usize;
        dst[count[byte]] = *entry;
        count[byte] += 1;
    }

    true
}

/// Single radix pass by end field. Returns false if all bytes are same (can skip).
#[inline]
fn radix_pass_end(src: &mut [SortEntry], dst: &mut [SortEntry], shift: u32) -> bool {
    let mut count = [0usize; 257];

    for entry in src.iter() {
        let byte = ((entry.end >> shift) & 0xFF) as usize;
        count[byte + 1] += 1;
    }

    let mut non_zero_buckets = 0;
    for &c in &count[1..] {
        if c > 0 {
            non_zero_buckets += 1;
        }
    }
    if non_zero_buckets <= 1 {
        return false;
    }

    for i in 1..257 {
        count[i] += count[i - 1];
    }

    for entry in src.iter() {
        let byte = ((entry.end >> shift) & 0xFF) as usize;
        dst[count[byte]] = *entry;
        count[byte] += 1;
    }

    true
}

/// Single radix pass by chrom_index field. Returns false if all bytes are same (can skip).
#[inline]
fn radix_pass_chrom(src: &mut [SortEntry], dst: &mut [SortEntry], shift: u32) -> bool {
    let mut count = [0usize; 257];

    for entry in src.iter() {
        let byte = ((entry.chrom_index >> shift) & 0xFF) as usize;
        count[byte + 1] += 1;
    }

    let mut non_zero_buckets = 0;
    for &c in &count[1..] {
        if c > 0 {
            non_zero_buckets += 1;
        }
    }
    if non_zero_buckets <= 1 {
        return false;
    }

    for i in 1..257 {
        count[i] += count[i - 1];
    }

    for entry in src.iter() {
        let byte = ((entry.chrom_index >> shift) & 0xFF) as usize;
        dst[count[byte]] = *entry;
        count[byte] += 1;
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_u32_fast() {
        assert_eq!(parse_u32_fast(b"0"), Some(0));
        assert_eq!(parse_u32_fast(b"123"), Some(123));
        assert_eq!(parse_u32_fast(b"100000000"), Some(100_000_000));
        assert_eq!(parse_u32_fast(b""), None);
        assert_eq!(parse_u32_fast(b"abc"), None);
    }

    #[test]
    fn test_parse_bed3() {
        let line = b"chr1\t100\t200";
        let (chrom, start, end) = parse_bed3(line).unwrap();
        assert_eq!(chrom, b"chr1");
        assert_eq!(start, 100);
        assert_eq!(end, 200);
    }

    #[test]
    fn test_parse_bed3_with_extra_fields() {
        let line = b"chr1\t100\t200\tname\t500\t+";
        let (chrom, start, end) = parse_bed3(line).unwrap();
        assert_eq!(chrom, b"chr1");
        assert_eq!(start, 100);
        assert_eq!(end, 200);
    }

    #[test]
    fn test_fast_sort_basic() {
        let input = b"chr2\t100\t200\nchr1\t150\t250\nchr1\t100\t200\n";
        let cmd = FastSortCommand::new();
        let mut output = Vec::new();

        let stats = cmd.sort_buffered(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 3);
        assert!(lines[0].starts_with("chr1\t100\t200"));
        assert!(lines[1].starts_with("chr1\t150\t250"));
        assert!(lines[2].starts_with("chr2\t100\t200"));
        assert_eq!(stats.records_read, 3);
    }

    #[test]
    fn test_fast_sort_reverse() {
        let input = b"chr1\t100\t200\nchr1\t200\t300\nchr2\t100\t200\n";
        let mut cmd = FastSortCommand::new();
        cmd.reverse = true;
        let mut output = Vec::new();

        cmd.sort_buffered(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 3);
        assert!(lines[0].starts_with("chr2"));
        assert!(lines[2].starts_with("chr1\t100"));
    }

    #[test]
    fn test_fast_sort_lexicographic() {
        let input = b"chr2\t100\t200\nchr10\t100\t200\nchr1\t100\t200\n";
        let cmd = FastSortCommand::new();
        let mut output = Vec::new();

        cmd.sort_buffered(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        // Lexicographic: chr1 < chr10 < chr2
        assert_eq!(lines.len(), 3);
        assert!(lines[0].starts_with("chr1\t"));
        assert!(lines[1].starts_with("chr10\t"));
        assert!(lines[2].starts_with("chr2\t"));
    }

    #[test]
    fn test_fast_sort_same_start_different_end() {
        // Test case: intervals with same chrom and start but different end
        // Sorted by end ascending (matches sort -k1,1 -k2,2n -k3,3n)
        let input = b"chr1\t100\t300\nchr1\t100\t200\nchr1\t100\t250\n";
        let cmd = FastSortCommand::new();
        let mut output = Vec::new();

        cmd.sort_buffered(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        // Same chrom and start - sorted by end ascending (200, 250, 300)
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[0], "chr1\t100\t200");
        assert_eq!(lines[1], "chr1\t100\t250");
        assert_eq!(lines[2], "chr1\t100\t300");
    }

    #[test]
    fn test_fast_sort_same_coordinates_different_extra_columns() {
        // Test case: intervals with same (chrom, start, end) but different extra columns
        // Should preserve input order (stable sort)
        let input = b"chr1\t100\t200\tgeneA\nchr1\t100\t200\tgeneB\nchr1\t100\t200\tgeneC\n";
        let cmd = FastSortCommand::new();
        let mut output = Vec::new();

        cmd.sort_buffered(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        // Same (chrom, start, end) - should maintain input order
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[0], "chr1\t100\t200\tgeneA");
        assert_eq!(lines[1], "chr1\t100\t200\tgeneB");
        assert_eq!(lines[2], "chr1\t100\t200\tgeneC");
    }

    #[test]
    fn test_fast_sort_mixed_chromosomes() {
        // Test mixed chromosome ordering with various sort keys
        // Input order: chr2:50, chr1:100(200), chr10:100, chr1:100(150), chr2:100, chr1:50
        let input = b"chr2\t50\t100\nchr1\t100\t200\nchr10\t100\t150\nchr1\t100\t150\nchr2\t100\t200\nchr1\t50\t100\n";
        let cmd = FastSortCommand::new();
        let mut output = Vec::new();

        cmd.sort_buffered(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        // Expected order (lexicographic chrom, then start, then end ascending):
        // chr1:50-100 (only one at start=50)
        // chr1:100-150 (end=150 < end=200)
        // chr1:100-200
        // chr10:100-150
        // chr2:50-100
        // chr2:100-200
        assert_eq!(lines.len(), 6);
        assert_eq!(lines[0], "chr1\t50\t100");
        assert_eq!(lines[1], "chr1\t100\t150"); // end=150 sorts before end=200
        assert_eq!(lines[2], "chr1\t100\t200");
        assert_eq!(lines[3], "chr10\t100\t150");
        assert_eq!(lines[4], "chr2\t50\t100");
        assert_eq!(lines[5], "chr2\t100\t200");
    }

    #[test]
    fn test_skip_headers() {
        let input = b"#header line\ntrack name=test\nbrowser position\nchr1\t100\t200\n";
        let cmd = FastSortCommand::new();
        let mut output = Vec::new();

        let stats = cmd.sort_buffered(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        assert_eq!(lines.len(), 1);
        assert!(lines[0].starts_with("chr1"));
        assert_eq!(stats.records_read, 1);
    }

    #[test]
    fn test_empty_input() {
        let input = b"";
        let cmd = FastSortCommand::new();
        let mut output = Vec::new();

        let stats = cmd.sort_buffered(&input[..], &mut output).unwrap();

        assert_eq!(stats.records_read, 0);
        assert!(output.is_empty());
    }

    #[test]
    fn test_radix_sort_large_dataset() {
        // Generate a dataset large enough to trigger radix sort
        let mut input = Vec::new();
        for i in (0..500).rev() {
            input.extend_from_slice(format!("chr1\t{}\t{}\n", i * 100, i * 100 + 50).as_bytes());
        }

        let cmd = FastSortCommand::new();
        let mut output = Vec::new();

        let stats = cmd.sort_buffered(&input[..], &mut output).unwrap();

        assert_eq!(stats.records_read, 500);
        assert!(stats.used_radix_sort);

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        // Verify sorted order
        assert_eq!(lines.len(), 500);
        for i in 0..500 {
            let expected = format!("chr1\t{}\t{}", i * 100, i * 100 + 50);
            assert_eq!(lines[i], expected, "Mismatch at index {}", i);
        }
    }

    #[test]
    fn test_radix_sort_entries_basic() {
        let entries = vec![
            SortEntry {
                chrom_index: 2,
                _padding: 0,
                start: 100,
                end: 200,
                line_start: 0,
                line_len: 10,
            },
            SortEntry {
                chrom_index: 1,
                _padding: 0,
                start: 200,
                end: 300,
                line_start: 10,
                line_len: 10,
            },
            SortEntry {
                chrom_index: 1,
                _padding: 0,
                start: 100,
                end: 200,
                line_start: 20,
                line_len: 10,
            },
        ];

        // For small inputs, radix sort falls back to comparison sort
        let sorted = radix_sort_lsd(entries);

        // Verify order: (chrom1, 100, 200), (chrom1, 200, 300), (chrom2, 100, 200)
        assert_eq!(sorted[0].line_start, 20); // chr1:100-200
        assert_eq!(sorted[1].line_start, 10); // chr1:200-300
        assert_eq!(sorted[2].line_start, 0); // chr2:100-200
    }

    #[test]
    fn test_deterministic_ordering() {
        // Test that identical records produce consistent output
        let input = b"chr1\t100\t200\tA\nchr1\t100\t200\tB\nchr1\t100\t200\tC\n";

        // Run sort multiple times and verify same output
        for _ in 0..10 {
            let cmd = FastSortCommand::new();
            let mut output = Vec::new();
            cmd.sort_buffered(&input[..], &mut output).unwrap();

            let result = String::from_utf8(output).unwrap();
            let lines: Vec<_> = result.lines().collect();

            // Order should always be A, B, C (input order preserved)
            assert_eq!(lines[0], "chr1\t100\t200\tA");
            assert_eq!(lines[1], "chr1\t100\t200\tB");
            assert_eq!(lines[2], "chr1\t100\t200\tC");
        }
    }

    #[test]
    fn test_end_as_tiebreaker() {
        // Verify that end is used as tiebreaker when (chrom, start) are same
        let input = b"chr1\t100\t500\nchr1\t100\t200\nchr1\t100\t300\nchr1\t100\t100\n";
        let cmd = FastSortCommand::new();
        let mut output = Vec::new();

        cmd.sort_buffered(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        // End ascending: 100, 200, 300, 500
        assert_eq!(lines.len(), 4);
        assert_eq!(lines[0], "chr1\t100\t100");
        assert_eq!(lines[1], "chr1\t100\t200");
        assert_eq!(lines[2], "chr1\t100\t300");
        assert_eq!(lines[3], "chr1\t100\t500");
    }

    #[test]
    fn test_build_chrom_index_with_genome() {
        // Create a genome order: chr2, chr1, chr10 (non-lexicographic)
        let mut genome_order = HashMap::new();
        genome_order.insert(b"chr2".to_vec(), 0u16);
        genome_order.insert(b"chr1".to_vec(), 1u16);
        genome_order.insert(b"chr10".to_vec(), 2u16);

        let data = b"chr1\t100\t200\nchr10\t50\t100\nchr2\t100\t200\n";
        let line_offsets = find_line_offsets(data);
        let chrom_index = build_chrom_index(data, &line_offsets, Some(&genome_order));

        // Verify genome ordering is used
        assert_eq!(chrom_index.get(b"chr2".as_slice()), Some(&0u16));
        assert_eq!(chrom_index.get(b"chr1".as_slice()), Some(&1u16));
        assert_eq!(chrom_index.get(b"chr10".as_slice()), Some(&2u16));
    }

    #[test]
    fn test_fast_sort_genome_order() {
        use crate::genome::Genome;

        // Create a genome with custom order: chr2, chr1, chr10
        let mut genome = Genome::new();
        genome.insert("chr2".to_string(), 1000);
        genome.insert("chr1".to_string(), 1000);
        genome.insert("chr10".to_string(), 1000);

        let input = b"chr1\t100\t200\nchr10\t50\t100\nchr2\t100\t200\n";
        let cmd = FastSortCommand::new().with_genome(&genome);
        let mut output = Vec::new();

        cmd.sort_buffered(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        // Expected order: chr2, chr1, chr10 (genome file order)
        assert_eq!(lines.len(), 3);
        assert!(lines[0].starts_with("chr2\t"));
        assert!(lines[1].starts_with("chr1\t"));
        assert!(lines[2].starts_with("chr10\t"));
    }

    #[test]
    fn test_fast_sort_genome_unknown_chrom() {
        use crate::genome::Genome;

        // Create a genome with only chr2 and chr1 (chr10 is unknown)
        let mut genome = Genome::new();
        genome.insert("chr2".to_string(), 1000);
        genome.insert("chr1".to_string(), 1000);

        let input = b"chr1\t100\t200\nchr10\t50\t100\nchr2\t100\t200\nchrX\t100\t200\n";
        let cmd = FastSortCommand::new().with_genome(&genome);
        let mut output = Vec::new();

        cmd.sort_buffered(&input[..], &mut output).unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<_> = result.lines().collect();

        // Expected order: chr2, chr1 (known in genome order), then chr10, chrX (unknown, lexicographic)
        assert_eq!(lines.len(), 4);
        assert!(lines[0].starts_with("chr2\t"));
        assert!(lines[1].starts_with("chr1\t"));
        assert!(lines[2].starts_with("chr10\t"));
        assert!(lines[3].starts_with("chrX\t"));
    }
}
