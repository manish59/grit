//! Core interval types for genomic region representation.

use std::cmp::Ordering;
use std::fmt;

/// A genomic interval with chromosome, start, and end positions.
/// Uses 0-based, half-open coordinates (BED format).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Interval {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl Interval {
    /// Create a new interval.
    #[inline]
    pub fn new(chrom: impl Into<String>, start: u64, end: u64) -> Self {
        Self {
            chrom: chrom.into(),
            start,
            end,
        }
    }

    /// Returns the length of the interval.
    #[inline]
    pub fn len(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    /// Returns true if the interval has zero length.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }

    /// Check if this interval overlaps with another.
    #[inline]
    pub fn overlaps(&self, other: &Interval) -> bool {
        self.chrom == other.chrom && self.start < other.end && other.start < self.end
    }

    /// Check if this interval overlaps with another by at least the given fraction.
    #[inline]
    pub fn overlaps_by_fraction(&self, other: &Interval, fraction: f64) -> bool {
        if !self.overlaps(other) {
            return false;
        }
        let overlap_start = self.start.max(other.start);
        let overlap_end = self.end.min(other.end);
        let overlap_len = overlap_end - overlap_start;
        let self_len = self.len();
        if self_len == 0 {
            return false;
        }
        (overlap_len as f64 / self_len as f64) >= fraction
    }

    /// Check if intervals overlap reciprocally by at least the given fraction.
    #[inline]
    pub fn overlaps_reciprocal(&self, other: &Interval, fraction: f64) -> bool {
        self.overlaps_by_fraction(other, fraction) && other.overlaps_by_fraction(self, fraction)
    }

    /// Compute the overlap length with another interval.
    #[inline]
    pub fn overlap_length(&self, other: &Interval) -> u64 {
        if !self.overlaps(other) {
            return 0;
        }
        let overlap_start = self.start.max(other.start);
        let overlap_end = self.end.min(other.end);
        overlap_end - overlap_start
    }

    /// Compute the distance to another interval.
    /// Returns 0 if overlapping, positive distance otherwise.
    #[inline]
    pub fn distance_to(&self, other: &Interval) -> Option<u64> {
        if self.chrom != other.chrom {
            return None;
        }
        if self.overlaps(other) {
            return Some(0);
        }
        if self.end <= other.start {
            Some(other.start - self.end)
        } else {
            Some(self.start - other.end)
        }
    }

    /// Merge this interval with another, returning the union.
    #[inline]
    pub fn merge(&self, other: &Interval) -> Option<Interval> {
        if self.chrom != other.chrom {
            return None;
        }
        Some(Interval {
            chrom: self.chrom.clone(),
            start: self.start.min(other.start),
            end: self.end.max(other.end),
        })
    }

    /// Subtract another interval from this one, returning remaining pieces.
    pub fn subtract(&self, other: &Interval) -> Vec<Interval> {
        if !self.overlaps(other) {
            return vec![self.clone()];
        }

        let mut result = Vec::new();

        // Left piece
        if self.start < other.start {
            result.push(Interval {
                chrom: self.chrom.clone(),
                start: self.start,
                end: other.start,
            });
        }

        // Right piece
        if self.end > other.end {
            result.push(Interval {
                chrom: self.chrom.clone(),
                start: other.end,
                end: self.end,
            });
        }

        result
    }
}

impl fmt::Display for Interval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\t{}\t{}", self.chrom, self.start, self.end)
    }
}

impl Ord for Interval {
    fn cmp(&self, other: &Self) -> Ordering {
        self.chrom
            .cmp(&other.chrom)
            .then(self.start.cmp(&other.start))
            .then(self.end.cmp(&other.end))
    }
}

impl PartialOrd for Interval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// A BED record with all optional fields.
#[derive(Debug, Clone, PartialEq)]
pub struct BedRecord {
    pub interval: Interval,
    pub name: Option<String>,
    pub score: Option<f64>,
    pub strand: Option<Strand>,
    pub thick_start: Option<u64>,
    pub thick_end: Option<u64>,
    pub item_rgb: Option<String>,
    pub block_count: Option<u32>,
    pub block_sizes: Option<Vec<u64>>,
    pub block_starts: Option<Vec<u64>>,
    /// Additional fields beyond BED12
    pub extra_fields: Vec<String>,
}

impl BedRecord {
    /// Create a minimal BED3 record.
    pub fn new(chrom: impl Into<String>, start: u64, end: u64) -> Self {
        Self {
            interval: Interval::new(chrom, start, end),
            name: None,
            score: None,
            strand: None,
            thick_start: None,
            thick_end: None,
            item_rgb: None,
            block_count: None,
            block_sizes: None,
            block_starts: None,
            extra_fields: Vec::new(),
        }
    }

    /// Get the chromosome.
    #[inline]
    pub fn chrom(&self) -> &str {
        &self.interval.chrom
    }

    /// Get the start position.
    #[inline]
    pub fn start(&self) -> u64 {
        self.interval.start
    }

    /// Get the end position.
    #[inline]
    pub fn end(&self) -> u64 {
        self.interval.end
    }

    /// Get the interval length.
    #[inline]
    pub fn len(&self) -> u64 {
        self.interval.len()
    }

    /// Check if the interval is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.interval.is_empty()
    }
}

impl fmt::Display for BedRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.interval)?;
        if let Some(ref name) = self.name {
            write!(f, "\t{}", name)?;
            if let Some(score) = self.score {
                write!(f, "\t{}", score)?;
                if let Some(strand) = self.strand {
                    write!(f, "\t{}", strand)?;
                }
            }
        }
        Ok(())
    }
}

/// Strand orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

impl Strand {
    pub fn from_char(c: char) -> Self {
        match c {
            '+' => Strand::Plus,
            '-' => Strand::Minus,
            _ => Strand::Unknown,
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
            Strand::Unknown => write!(f, "."),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_overlap() {
        let a = Interval::new("chr1", 100, 200);
        let b = Interval::new("chr1", 150, 250);
        let c = Interval::new("chr1", 200, 300);
        let d = Interval::new("chr2", 100, 200);

        assert!(a.overlaps(&b));
        assert!(!a.overlaps(&c)); // Adjacent, not overlapping
        assert!(!a.overlaps(&d)); // Different chromosome
    }

    #[test]
    fn test_interval_overlap_fraction() {
        let a = Interval::new("chr1", 100, 200);
        let b = Interval::new("chr1", 150, 250);

        // 50 bp overlap out of 100 bp = 50%
        assert!(a.overlaps_by_fraction(&b, 0.5));
        assert!(!a.overlaps_by_fraction(&b, 0.6));
    }

    #[test]
    fn test_interval_distance() {
        let a = Interval::new("chr1", 100, 200);
        let b = Interval::new("chr1", 300, 400);

        assert_eq!(a.distance_to(&b), Some(100));
    }

    #[test]
    fn test_interval_merge() {
        let a = Interval::new("chr1", 100, 200);
        let b = Interval::new("chr1", 150, 250);

        let merged = a.merge(&b).unwrap();
        assert_eq!(merged.start, 100);
        assert_eq!(merged.end, 250);
    }

    #[test]
    fn test_interval_subtract() {
        let a = Interval::new("chr1", 100, 300);
        let b = Interval::new("chr1", 150, 200);

        let pieces = a.subtract(&b);
        assert_eq!(pieces.len(), 2);
        assert_eq!(pieces[0].start, 100);
        assert_eq!(pieces[0].end, 150);
        assert_eq!(pieces[1].start, 200);
        assert_eq!(pieces[1].end, 300);
    }

    #[test]
    fn test_interval_ordering() {
        let mut intervals = [
            Interval::new("chr2", 100, 200),
            Interval::new("chr1", 200, 300),
            Interval::new("chr1", 100, 200),
        ];
        intervals.sort();

        assert_eq!(intervals[0].chrom, "chr1");
        assert_eq!(intervals[0].start, 100);
        assert_eq!(intervals[1].start, 200);
        assert_eq!(intervals[2].chrom, "chr2");
    }
}
