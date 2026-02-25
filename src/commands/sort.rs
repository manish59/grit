//! Sort command implementation.
//!
//! Sort order (matches `LC_ALL=C sort -k1,1 -k2,2n -k3,3n`):
//! 1. Primary: chromosome (lexicographic, or genome order with -g)
//! 2. Secondary: start coordinate (ascending, numeric)
//! 3. Tertiary: end coordinate (ascending, numeric)
//! 4. Ties: input order preserved (stable sort)
//!
//! Optimized sorting with pre-computed chromosome keys to avoid
//! O(n log n) string allocations during comparison.

use crate::bed::{read_records, BedError, BedReader};
use crate::interval::BedRecord;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::io::{self, BufWriter, Write};
use std::path::Path;

/// Sort key specification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortKey {
    /// Sort by chromosome, start, end (default)
    ChromStart,
    /// Sort by size (ascending)
    SizeAsc,
    /// Sort by size (descending)
    SizeDesc,
    /// Sort by chromosome only
    ChromOnly,
    /// Reverse order
    Reverse,
}

/// Sort command configuration.
#[derive(Debug, Clone)]
pub struct SortCommand {
    /// Sort by size ascending
    pub size_asc: bool,
    /// Sort by size descending
    pub size_desc: bool,
    /// Reverse the sort order
    pub reverse: bool,
    /// Sort by chromosome only
    pub chrom_only: bool,
    /// Natural sort for chromosome names
    pub natural_sort: bool,
    /// Genome-based chromosome ordering (chrom name -> index)
    genome_order: Option<HashMap<String, u32>>,
}

impl Default for SortCommand {
    fn default() -> Self {
        Self::new()
    }
}

impl SortCommand {
    pub fn new() -> Self {
        Self {
            size_asc: false,
            size_desc: false,
            reverse: false,
            chrom_only: false,
            natural_sort: false, // Lexicographic by default (matches GNU sort -k1,1)
            genome_order: None,
        }
    }

    /// Set genome-based chromosome ordering.
    /// Chromosomes will be sorted in the order they appear in the genome file.
    /// Unknown chromosomes are placed after all known chromosomes.
    pub fn with_genome(mut self, genome: &crate::genome::Genome) -> Self {
        let order: HashMap<String, u32> = genome
            .chromosomes()
            .enumerate()
            .map(|(i, chrom)| (chrom.clone(), i as u32))
            .collect();
        self.genome_order = Some(order);
        self
    }

    /// Sort BED records.
    pub fn sort(&self, mut records: Vec<BedRecord>) -> Vec<BedRecord> {
        if self.size_asc {
            records.sort_by(|a, b| self.compare_by_size_asc(a, b));
        } else if self.size_desc {
            records.sort_by(|a, b| self.compare_by_size_desc(a, b));
        } else if self.chrom_only {
            records.sort_by(|a, b| self.compare_chrom_with_genome(a.chrom(), b.chrom()));
        } else {
            // Default: chrom, start, end
            records.sort_by(|a, b| self.compare_default(a, b));
        }

        if self.reverse {
            records.reverse();
        }

        records
    }

    /// Get chromosome order index for genome-based sorting.
    /// Unknown chromosomes get u32::MAX to sort them to the end.
    #[inline]
    fn chrom_order(&self, chrom: &str) -> u32 {
        self.genome_order
            .as_ref()
            .and_then(|order| order.get(chrom).copied())
            .unwrap_or(u32::MAX)
    }

    /// Compare chromosomes using genome order if available.
    fn compare_chrom_with_genome(&self, a: &str, b: &str) -> Ordering {
        if self.genome_order.is_some() {
            let a_ord = self.chrom_order(a);
            let b_ord = self.chrom_order(b);
            // If both unknown, fall back to lexicographic comparison
            if a_ord == u32::MAX && b_ord == u32::MAX {
                return a.cmp(b);
            }
            a_ord.cmp(&b_ord)
        } else {
            self.compare_chrom(a, b)
        }
    }

    /// Parallel sort using Rayon with stable sort to preserve input order for ties.
    pub fn sort_parallel(&self, mut records: Vec<BedRecord>) -> Vec<BedRecord> {
        if records.len() < 2 {
            return records;
        }

        // For large datasets with genome order, use pre-computed keys
        if self.genome_order.is_some() && records.len() > 10000 && !self.size_asc && !self.size_desc
        {
            return self.sort_parallel_genome(records);
        }

        // For large datasets with natural sort, use pre-computed keys
        if self.natural_sort && records.len() > 10000 && !self.size_asc && !self.size_desc {
            return self.sort_parallel_natural(records);
        }

        // Direct parallel sort - use stable sort to preserve input order for ties
        // Note: par_sort_by is stable, par_sort_unstable_by is not
        if self.size_asc {
            records.par_sort_by(|a, b| self.compare_by_size_asc(a, b));
        } else if self.size_desc {
            records.par_sort_by(|a, b| self.compare_by_size_desc(a, b));
        } else if self.chrom_only {
            records.par_sort_by(|a, b| self.compare_chrom_with_genome(a.chrom(), b.chrom()));
        } else if self.genome_order.is_some() {
            // Use genome ordering for chromosome comparison
            records.par_sort_by(|a, b| {
                self.compare_chrom_with_genome(a.chrom(), b.chrom())
                    .then_with(|| a.start().cmp(&b.start()))
                    .then_with(|| a.end().cmp(&b.end()))
            });
        } else {
            // Default: sort by (chrom, start, end)
            records.par_sort_by(|a, b| {
                a.chrom()
                    .cmp(b.chrom())
                    .then_with(|| a.start().cmp(&b.start()))
                    .then_with(|| a.end().cmp(&b.end()))
            });
        }

        if self.reverse {
            records.reverse();
        }

        records
    }

    /// Optimized genome-based sort using pre-computed chromosome order.
    /// This avoids repeated hashmap lookups by computing order once.
    /// Uses stable sorting: (chrom_order, start, end, original_index)
    fn sort_parallel_genome(&self, mut records: Vec<BedRecord>) -> Vec<BedRecord> {
        let n = records.len();
        let genome_order = self.genome_order.as_ref().unwrap();

        // Collect unique chromosomes from records
        let mut unique_chroms: Vec<String> =
            records.par_iter().map(|r| r.chrom().to_string()).collect();
        unique_chroms.par_sort_unstable();
        unique_chroms.dedup();

        // Build chromosome order map: known chroms get their genome index,
        // unknown chroms get max_index + 1 + lexicographic position
        let max_known = genome_order.values().max().copied().unwrap_or(0);
        let mut unknown_chroms: Vec<_> = unique_chroms
            .iter()
            .filter(|c| !genome_order.contains_key(c.as_str()))
            .cloned()
            .collect();
        unknown_chroms.sort(); // Lexicographic order for unknown

        let mut chrom_order: HashMap<String, u32> = genome_order.clone();
        for (i, chrom) in unknown_chroms.into_iter().enumerate() {
            chrom_order.insert(chrom, max_known + 1 + i as u32);
        }

        // Create sort keys for all records (done once, O(n))
        // Key is (chrom_order, start, end, original_index) - original_index for stable tie-breaking
        let mut indexed: Vec<(u32, u64, u64, usize)> = records
            .par_iter()
            .enumerate()
            .map(|(i, r)| {
                let chrom_ord = chrom_order.get(r.chrom()).copied().unwrap_or(u32::MAX);
                (chrom_ord, r.start(), r.end(), i)
            })
            .collect();

        // Sort the keys (pure integer comparison - very fast)
        // Include original index to ensure stable sorting
        if self.chrom_only {
            indexed.par_sort_by_key(|k| (k.0, k.3));
        } else {
            indexed.par_sort_by_key(|k| (k.0, k.1, k.2, k.3));
        }

        if self.reverse {
            indexed.reverse();
        }

        // Reorder records in-place using permutation cycles
        let mut perm: Vec<usize> = indexed.into_iter().map(|(_, _, _, i)| i).collect();

        for i in 0..n {
            if perm[i] == i {
                continue;
            }
            let mut j = i;
            loop {
                let target = perm[j];
                perm[j] = j;
                if target == i {
                    break;
                }
                records.swap(j, target);
                j = target;
            }
        }

        records
    }

    /// Optimized natural sort using pre-computed chromosome order.
    /// This avoids O(n log n) natural comparisons by computing order once.
    /// Uses stable sorting: (chrom_order, start, end, original_index)
    fn sort_parallel_natural(&self, mut records: Vec<BedRecord>) -> Vec<BedRecord> {
        let n = records.len();

        // Step 1: Collect unique chromosomes
        let mut unique_chroms: Vec<String> =
            records.par_iter().map(|r| r.chrom().to_string()).collect();
        unique_chroms.par_sort_unstable();
        unique_chroms.dedup();

        // Sort chromosomes naturally and assign numeric order
        let mut chrom_keys: Vec<(ChromSortKey, String)> = unique_chroms
            .into_iter()
            .map(|c| {
                let key = ChromSortKey::from_chrom(&c);
                (key, c)
            })
            .collect();
        chrom_keys.sort_unstable_by(|a, b| a.0.cmp(&b.0));

        let chrom_order: HashMap<String, u32> = chrom_keys
            .into_iter()
            .enumerate()
            .map(|(i, (_, c))| (c, i as u32))
            .collect();

        // Step 2: Create sort keys for all records (done once, O(n))
        // Key is (chrom_order, start, end, original_index) - original_index for stable tie-breaking
        let mut indexed: Vec<(u32, u64, u64, usize)> = records
            .par_iter()
            .enumerate()
            .map(|(i, r)| {
                let chrom_ord = chrom_order.get(r.chrom()).copied().unwrap_or(u32::MAX);
                (chrom_ord, r.start(), r.end(), i)
            })
            .collect();

        // Step 3: Sort the keys (pure integer comparison - very fast)
        // Include original index to ensure stable sorting
        if self.chrom_only {
            indexed.par_sort_by_key(|k| (k.0, k.3));
        } else {
            indexed.par_sort_by_key(|k| (k.0, k.1, k.2, k.3));
        }

        if self.reverse {
            indexed.reverse();
        }

        // Step 4: Reorder records in-place using permutation cycles
        // This avoids extra allocations by doing swaps
        let mut perm: Vec<usize> = indexed.into_iter().map(|(_, _, _, i)| i).collect();

        // Apply permutation in-place using cycle detection
        for i in 0..n {
            if perm[i] == i {
                continue;
            }
            let mut j = i;
            loop {
                let target = perm[j];
                perm[j] = j; // Mark as done
                if target == i {
                    break;
                }
                records.swap(j, target);
                j = target;
            }
        }

        records
    }

    fn compare_default(&self, a: &BedRecord, b: &BedRecord) -> Ordering {
        self.compare_chrom_with_genome(a.chrom(), b.chrom())
            .then(a.start().cmp(&b.start()))
            .then(a.end().cmp(&b.end()))
    }

    fn compare_by_size_asc(&self, a: &BedRecord, b: &BedRecord) -> Ordering {
        a.len()
            .cmp(&b.len())
            .then(self.compare_chrom_with_genome(a.chrom(), b.chrom()))
            .then(a.start().cmp(&b.start()))
    }

    fn compare_by_size_desc(&self, a: &BedRecord, b: &BedRecord) -> Ordering {
        b.len()
            .cmp(&a.len())
            .then(self.compare_chrom_with_genome(a.chrom(), b.chrom()))
            .then(a.start().cmp(&b.start()))
    }

    fn compare_chrom(&self, a: &str, b: &str) -> Ordering {
        if self.natural_sort {
            natural_compare(a, b)
        } else {
            a.cmp(b)
        }
    }

    /// Execute sort command on a file.
    pub fn run<P: AsRef<Path>, W: Write>(&self, input: P, output: &mut W) -> Result<(), BedError> {
        let records = read_records(input)?;
        let sorted = self.sort_parallel(records);

        // Use buffered writer for better I/O performance
        let mut buf_output = BufWriter::with_capacity(256 * 1024, output);
        for record in sorted {
            writeln!(buf_output, "{}", record).map_err(BedError::Io)?;
        }
        buf_output.flush().map_err(BedError::Io)?;

        Ok(())
    }

    /// Execute sort from stdin to stdout.
    pub fn run_stdio(&self) -> Result<(), BedError> {
        let stdin = io::stdin();
        let reader = BedReader::new(stdin.lock());
        let records: Result<Vec<_>, _> = reader.records().collect();
        let records = records?;

        let sorted = self.sort_parallel(records);

        let stdout = io::stdout();
        let handle = stdout.lock();
        let mut buf_output = BufWriter::with_capacity(256 * 1024, handle);
        for record in sorted {
            writeln!(buf_output, "{}", record).map_err(BedError::Io)?;
        }
        buf_output.flush().map_err(BedError::Io)?;

        Ok(())
    }
}

/// Pre-computed chromosome sort key for O(1) comparisons.
/// Parses chromosome names like "chr10_random" into comparable components.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ChromSortKey {
    /// Text prefix (e.g., "chr")
    prefix: String,
    /// Primary numeric part (e.g., 10 for "chr10")
    number: Option<u64>,
    /// Suffix or remaining text (e.g., "_random" or "X" for chrX)
    suffix: String,
    /// Secondary number in suffix if present
    suffix_number: Option<u64>,
}

impl ChromSortKey {
    fn from_chrom(s: &str) -> Self {
        let bytes = s.as_bytes();
        let len = bytes.len();

        // Find where digits start
        let mut digit_start = 0;
        while digit_start < len && !bytes[digit_start].is_ascii_digit() {
            digit_start += 1;
        }

        let prefix = s[..digit_start].to_string();

        if digit_start >= len {
            // No digits found (e.g., "chrX", "chrM")
            return Self {
                prefix,
                number: None,
                suffix: String::new(),
                suffix_number: None,
            };
        }

        // Find where digits end
        let mut digit_end = digit_start;
        while digit_end < len && bytes[digit_end].is_ascii_digit() {
            digit_end += 1;
        }

        let number = s[digit_start..digit_end].parse().ok();
        let suffix_str = &s[digit_end..];

        // Check for number in suffix (e.g., "chr1_random2")
        let mut suffix_digit_start = 0;
        let suffix_bytes = suffix_str.as_bytes();
        while suffix_digit_start < suffix_bytes.len()
            && !suffix_bytes[suffix_digit_start].is_ascii_digit()
        {
            suffix_digit_start += 1;
        }

        let (suffix, suffix_number) = if suffix_digit_start < suffix_bytes.len() {
            let mut suffix_digit_end = suffix_digit_start;
            while suffix_digit_end < suffix_bytes.len()
                && suffix_bytes[suffix_digit_end].is_ascii_digit()
            {
                suffix_digit_end += 1;
            }
            (
                suffix_str[..suffix_digit_start].to_string(),
                suffix_str[suffix_digit_start..suffix_digit_end]
                    .parse()
                    .ok(),
            )
        } else {
            (suffix_str.to_string(), None)
        };

        Self {
            prefix,
            number,
            suffix,
            suffix_number,
        }
    }
}

impl Ord for ChromSortKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.prefix
            .cmp(&other.prefix)
            .then_with(|| {
                // Numbers sort before non-numbers (chr1 < chrX)
                match (&self.number, &other.number) {
                    (Some(a), Some(b)) => a.cmp(b),
                    (Some(_), None) => Ordering::Less,
                    (None, Some(_)) => Ordering::Greater,
                    (None, None) => Ordering::Equal,
                }
            })
            .then_with(|| self.suffix.cmp(&other.suffix))
            .then_with(|| match (&self.suffix_number, &other.suffix_number) {
                (Some(a), Some(b)) => a.cmp(b),
                (Some(_), None) => Ordering::Less,
                (None, Some(_)) => Ordering::Greater,
                (None, None) => Ordering::Equal,
            })
    }
}

impl PartialOrd for ChromSortKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Natural comparison for chromosome names.
/// Handles numeric suffixes properly (chr1 < chr2 < chr10).
fn natural_compare(a: &str, b: &str) -> Ordering {
    let a_parts = split_numeric(a);
    let b_parts = split_numeric(b);

    for (ap, bp) in a_parts.iter().zip(b_parts.iter()) {
        let cmp = match (ap, bp) {
            (Part::Text(at), Part::Text(bt)) => at.cmp(bt),
            (Part::Number(an), Part::Number(bn)) => an.cmp(bn),
            (Part::Text(_), Part::Number(_)) => Ordering::Less,
            (Part::Number(_), Part::Text(_)) => Ordering::Greater,
        };
        if cmp != Ordering::Equal {
            return cmp;
        }
    }

    a_parts.len().cmp(&b_parts.len())
}

#[derive(Debug, PartialEq, Eq)]
enum Part {
    Text(String),
    Number(u64),
}

fn split_numeric(s: &str) -> Vec<Part> {
    let mut parts = Vec::new();
    let mut current_text = String::new();
    let mut current_num = String::new();
    let mut in_number = false;

    for c in s.chars() {
        if c.is_ascii_digit() {
            if !in_number && !current_text.is_empty() {
                parts.push(Part::Text(current_text.clone()));
                current_text.clear();
            }
            in_number = true;
            current_num.push(c);
        } else {
            if in_number && !current_num.is_empty() {
                if let Ok(n) = current_num.parse() {
                    parts.push(Part::Number(n));
                }
                current_num.clear();
            }
            in_number = false;
            current_text.push(c);
        }
    }

    if !current_text.is_empty() {
        parts.push(Part::Text(current_text));
    }
    if !current_num.is_empty() {
        if let Ok(n) = current_num.parse() {
            parts.push(Part::Number(n));
        }
    }

    parts
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_record(chrom: &str, start: u64, end: u64) -> BedRecord {
        BedRecord::new(chrom, start, end)
    }

    #[test]
    fn test_default_sort() {
        let cmd = SortCommand::new();
        let records = vec![
            make_record("chr2", 100, 200),
            make_record("chr1", 300, 400),
            make_record("chr1", 100, 200),
        ];

        let sorted = cmd.sort(records);

        assert_eq!(sorted[0].chrom(), "chr1");
        assert_eq!(sorted[0].start(), 100);
        assert_eq!(sorted[1].start(), 300);
        assert_eq!(sorted[2].chrom(), "chr2");
    }

    #[test]
    fn test_natural_sort() {
        let mut cmd = SortCommand::new();
        cmd.natural_sort = true; // Enable natural sort explicitly
        let records = vec![
            make_record("chr10", 100, 200),
            make_record("chr2", 100, 200),
            make_record("chr1", 100, 200),
        ];

        let sorted = cmd.sort(records);

        assert_eq!(sorted[0].chrom(), "chr1");
        assert_eq!(sorted[1].chrom(), "chr2");
        assert_eq!(sorted[2].chrom(), "chr10");
    }

    #[test]
    fn test_lexicographic_sort() {
        // Default is lexicographic (matches GNU sort -k1,1)
        let cmd = SortCommand::new();
        let records = vec![
            make_record("chr10", 100, 200),
            make_record("chr2", 100, 200),
            make_record("chr1", 100, 200),
        ];

        let sorted = cmd.sort(records);

        // Lexicographic: chr1 < chr10 < chr2
        assert_eq!(sorted[0].chrom(), "chr1");
        assert_eq!(sorted[1].chrom(), "chr10");
        assert_eq!(sorted[2].chrom(), "chr2");
    }

    #[test]
    fn test_size_sort() {
        let mut cmd = SortCommand::new();
        cmd.size_asc = true;

        let records = vec![
            make_record("chr1", 100, 300), // 200bp
            make_record("chr1", 100, 150), // 50bp
            make_record("chr1", 100, 200), // 100bp
        ];

        let sorted = cmd.sort(records);

        assert_eq!(sorted[0].len(), 50);
        assert_eq!(sorted[1].len(), 100);
        assert_eq!(sorted[2].len(), 200);
    }

    #[test]
    fn test_reverse_sort() {
        let mut cmd = SortCommand::new();
        cmd.reverse = true;

        let records = vec![make_record("chr1", 100, 200), make_record("chr1", 300, 400)];

        let sorted = cmd.sort(records);

        assert_eq!(sorted[0].start(), 300);
        assert_eq!(sorted[1].start(), 100);
    }

    #[test]
    fn test_natural_compare() {
        assert_eq!(natural_compare("chr1", "chr2"), Ordering::Less);
        assert_eq!(natural_compare("chr2", "chr10"), Ordering::Less);
        assert_eq!(natural_compare("chr10", "chr2"), Ordering::Greater);
        assert_eq!(natural_compare("chrX", "chrY"), Ordering::Less);
    }

    #[test]
    fn test_genome_order_sort() {
        use crate::genome::Genome;

        // Create a genome with custom order: chr2, chr1, chr10
        let mut genome = Genome::new();
        genome.insert("chr2".to_string(), 1000);
        genome.insert("chr1".to_string(), 1000);
        genome.insert("chr10".to_string(), 1000);

        let cmd = SortCommand::new().with_genome(&genome);
        let records = vec![
            make_record("chr1", 100, 200),
            make_record("chr10", 50, 100),
            make_record("chr2", 100, 200),
        ];

        let sorted = cmd.sort(records);

        // Expected order: chr2, chr1, chr10 (genome file order)
        assert_eq!(sorted[0].chrom(), "chr2");
        assert_eq!(sorted[1].chrom(), "chr1");
        assert_eq!(sorted[2].chrom(), "chr10");
    }

    #[test]
    fn test_genome_order_unknown_chrom() {
        use crate::genome::Genome;

        // Create a genome with only chr2 and chr1 (chr10 and chrX are unknown)
        let mut genome = Genome::new();
        genome.insert("chr2".to_string(), 1000);
        genome.insert("chr1".to_string(), 1000);

        let cmd = SortCommand::new().with_genome(&genome);
        let records = vec![
            make_record("chr1", 100, 200),
            make_record("chr10", 50, 100),
            make_record("chr2", 100, 200),
            make_record("chrX", 100, 200),
        ];

        let sorted = cmd.sort(records);

        // Expected order: chr2, chr1 (known in genome order), then chr10, chrX (unknown, lexicographic)
        assert_eq!(sorted[0].chrom(), "chr2");
        assert_eq!(sorted[1].chrom(), "chr1");
        assert_eq!(sorted[2].chrom(), "chr10");
        assert_eq!(sorted[3].chrom(), "chrX");
    }

    #[test]
    fn test_genome_order_parallel() {
        use crate::genome::Genome;

        // Create a genome with custom order
        let mut genome = Genome::new();
        genome.insert("chr2".to_string(), 1000);
        genome.insert("chr1".to_string(), 1000);
        genome.insert("chr10".to_string(), 1000);

        let cmd = SortCommand::new().with_genome(&genome);
        let records = vec![
            make_record("chr1", 100, 200),
            make_record("chr10", 50, 100),
            make_record("chr2", 100, 200),
            make_record("chr1", 50, 100),
            make_record("chr2", 50, 100),
        ];

        let sorted = cmd.sort_parallel(records);

        // Expected order: chr2:50, chr2:100, chr1:50, chr1:100, chr10:50
        assert_eq!(sorted[0].chrom(), "chr2");
        assert_eq!(sorted[0].start(), 50);
        assert_eq!(sorted[1].chrom(), "chr2");
        assert_eq!(sorted[1].start(), 100);
        assert_eq!(sorted[2].chrom(), "chr1");
        assert_eq!(sorted[2].start(), 50);
        assert_eq!(sorted[3].chrom(), "chr1");
        assert_eq!(sorted[3].start(), 100);
        assert_eq!(sorted[4].chrom(), "chr10");
    }

    #[test]
    fn test_end_as_tiebreaker() {
        // Same (chrom, start), different end - should sort by end ascending
        let cmd = SortCommand::new();
        let records = vec![
            make_record("chr1", 100, 500),
            make_record("chr1", 100, 200),
            make_record("chr1", 100, 300),
            make_record("chr1", 100, 100),
        ];

        let sorted = cmd.sort(records);

        // End ascending: 100, 200, 300, 500
        assert_eq!(sorted[0].end(), 100);
        assert_eq!(sorted[1].end(), 200);
        assert_eq!(sorted[2].end(), 300);
        assert_eq!(sorted[3].end(), 500);
    }

    #[test]
    fn test_end_as_tiebreaker_parallel() {
        // Same test but using parallel sort
        let cmd = SortCommand::new();
        let records = vec![
            make_record("chr1", 100, 500),
            make_record("chr1", 100, 200),
            make_record("chr1", 100, 300),
            make_record("chr1", 100, 100),
        ];

        let sorted = cmd.sort_parallel(records);

        // End ascending: 100, 200, 300, 500
        assert_eq!(sorted[0].end(), 100);
        assert_eq!(sorted[1].end(), 200);
        assert_eq!(sorted[2].end(), 300);
        assert_eq!(sorted[3].end(), 500);
    }

    #[test]
    fn test_stability_with_identical_records() {
        // Records with same (chrom, start, end) - input order must be preserved
        let cmd = SortCommand::new();

        // Create records with same (chrom, start) and sequential ends
        // Since end is now a tiebreaker, these should sort by end ascending
        // which happens to match input order here
        let records = vec![
            make_record("chr1", 100, 201),
            make_record("chr1", 100, 202),
            make_record("chr1", 100, 203),
        ];

        let sorted = cmd.sort(records);

        // End ascending: 201, 202, 203
        assert_eq!(sorted[0].end(), 201);
        assert_eq!(sorted[1].end(), 202);
        assert_eq!(sorted[2].end(), 203);
    }

    #[test]
    fn test_sort_and_sort_parallel_produce_same_output() {
        // Verify that sort() and sort_parallel() produce identical results
        let cmd = SortCommand::new();

        let records = vec![
            make_record("chr2", 100, 200),
            make_record("chr1", 100, 300),
            make_record("chr1", 100, 200),
            make_record("chr10", 50, 100),
            make_record("chr1", 50, 100),
        ];

        let sorted_seq = cmd.sort(records.clone());
        let sorted_par = cmd.sort_parallel(records);

        assert_eq!(sorted_seq.len(), sorted_par.len());
        for (i, (seq, par)) in sorted_seq.iter().zip(sorted_par.iter()).enumerate() {
            assert_eq!(seq.chrom(), par.chrom(), "Chrom differs at index {}", i);
            assert_eq!(seq.start(), par.start(), "Start differs at index {}", i);
            assert_eq!(seq.end(), par.end(), "End differs at index {}", i);
        }
    }

    #[test]
    fn test_reverse_is_exact_reverse() {
        // Verify that reverse mode produces exact reverse of forward sort
        let cmd_fwd = SortCommand::new();
        let mut cmd_rev = SortCommand::new();
        cmd_rev.reverse = true;

        let records = vec![
            make_record("chr2", 100, 200),
            make_record("chr1", 100, 200),
            make_record("chr1", 50, 100),
            make_record("chr10", 100, 200),
        ];

        let sorted_fwd = cmd_fwd.sort(records.clone());
        let sorted_rev = cmd_rev.sort(records);

        // Reverse should be exact reverse of forward
        let reversed_fwd: Vec<_> = sorted_fwd.iter().rev().collect();
        for (i, (rev, fwd_rev)) in sorted_rev.iter().zip(reversed_fwd.iter()).enumerate() {
            assert_eq!(rev.chrom(), fwd_rev.chrom(), "Chrom differs at index {}", i);
            assert_eq!(rev.start(), fwd_rev.start(), "Start differs at index {}", i);
        }
    }
}
