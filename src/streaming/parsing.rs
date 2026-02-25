//! Zero-allocation BED parsing utilities.
//!
//! These functions provide high-performance parsing of BED records
//! without any heap allocation in the hot path.

use crate::config::normalize_end;
use memchr::memchr;

/// Fast u64 parsing - no allocation, no error formatting.
///
/// Returns None if the input is empty or contains non-digit characters.
///
/// # Performance
///
/// This is approximately 3x faster than `str::parse::<u64>()` because:
/// - No UTF-8 validation (input is already bytes)
/// - No error string formatting
/// - Inline-always for elimination of call overhead
#[inline(always)]
pub fn parse_u64_fast(bytes: &[u8]) -> Option<u64> {
    if bytes.is_empty() {
        return None;
    }
    let mut n: u64 = 0;
    for &b in bytes {
        let d = b.wrapping_sub(b'0');
        if d > 9 {
            return None;
        }
        n = n.wrapping_mul(10).wrapping_add(d as u64);
    }
    Some(n)
}

/// Parse BED3 fields using memchr - zero allocation.
///
/// Returns (chrom_bytes, start, end) or None if parsing fails.
///
/// # Performance
///
/// Uses memchr for SIMD-accelerated tab searching, avoiding
/// the overhead of splitting into a Vec.
///
/// # Bedtools Compatibility
///
/// If bedtools-compatible mode is enabled, zero-length intervals
/// (start == end) are normalized to 1bp intervals (end = start + 1).
#[inline(always)]
pub fn parse_bed3_bytes(line: &[u8]) -> Option<(&[u8], u64, u64)> {
    let tab1 = memchr(b'\t', line)?;
    let chrom = &line[..tab1];

    let rest1 = &line[tab1 + 1..];
    let tab2 = memchr(b'\t', rest1)?;
    let start = parse_u64_fast(&rest1[..tab2])?;

    let rest2 = &rest1[tab2 + 1..];
    let end_len = memchr(b'\t', rest2).unwrap_or(rest2.len());
    let end_len_trimmed = memchr(b'\n', &rest2[..end_len]).unwrap_or(end_len);
    let end = parse_u64_fast(&rest2[..end_len_trimmed])?;

    // Normalize zero-length intervals if bedtools-compatible mode is enabled
    let end = normalize_end(start, end);

    Some((chrom, start, end))
}

/// Parse BED3 fields and return the rest of line index.
///
/// Returns (chrom_bytes, start, end, rest_start_idx) where rest_start_idx
/// is the byte offset where additional columns begin.
///
/// This variant is useful when the original line needs to be preserved
/// with modified coordinates (e.g., in subtract operations).
///
/// # Bedtools Compatibility
///
/// If bedtools-compatible mode is enabled, zero-length intervals
/// (start == end) are normalized to 1bp intervals (end = start + 1).
#[inline(always)]
pub fn parse_bed3_bytes_with_rest(line: &[u8]) -> Option<(&[u8], u64, u64, usize)> {
    let tab1 = memchr(b'\t', line)?;
    let chrom = &line[..tab1];

    let rest1 = &line[tab1 + 1..];
    let tab2 = memchr(b'\t', rest1)?;
    let start = parse_u64_fast(&rest1[..tab2])?;

    let rest2 = &rest1[tab2 + 1..];
    let end_len = memchr(b'\t', rest2).unwrap_or(rest2.len());
    let end_len_trimmed = memchr(b'\n', &rest2[..end_len]).unwrap_or(end_len);
    let end = parse_u64_fast(&rest2[..end_len_trimmed])?;

    // Normalize zero-length intervals if bedtools-compatible mode is enabled
    let end = normalize_end(start, end);

    // Calculate where the rest of the line starts (after end field)
    let rest_start = tab1 + 1 + tab2 + 1 + end_len;

    Some((chrom, start, end, rest_start))
}

/// Check if a line should be skipped (empty, comment, or header).
#[inline(always)]
pub fn should_skip_line(line: &[u8]) -> bool {
    line.is_empty() || line[0] == b'#' || line.starts_with(b"track") || line.starts_with(b"browser")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_u64_fast() {
        assert_eq!(parse_u64_fast(b"12345"), Some(12345));
        assert_eq!(parse_u64_fast(b"0"), Some(0));
        assert_eq!(parse_u64_fast(b""), None);
        assert_eq!(parse_u64_fast(b"abc"), None);
        assert_eq!(parse_u64_fast(b"123abc"), None);
        assert_eq!(parse_u64_fast(b"18446744073709551615"), Some(u64::MAX));
    }

    #[test]
    fn test_parse_bed3_bytes() {
        assert_eq!(
            parse_bed3_bytes(b"chr1\t100\t200"),
            Some((&b"chr1"[..], 100, 200))
        );
        assert_eq!(
            parse_bed3_bytes(b"chr1\t100\t200\tname"),
            Some((&b"chr1"[..], 100, 200))
        );
        assert_eq!(
            parse_bed3_bytes(b"chr1\t100\t200\n"),
            Some((&b"chr1"[..], 100, 200))
        );
        assert_eq!(parse_bed3_bytes(b"chr1\t100"), None);
        assert_eq!(parse_bed3_bytes(b""), None);
    }

    #[test]
    fn test_parse_bed3_bytes_with_rest() {
        let result = parse_bed3_bytes_with_rest(b"chr1\t100\t200\tname\t50\t+");
        assert!(result.is_some());
        let (chrom, start, end, rest_start) = result.unwrap();
        assert_eq!(chrom, b"chr1");
        assert_eq!(start, 100);
        assert_eq!(end, 200);
        assert_eq!(rest_start, 12); // Position after "200"
    }

    #[test]
    fn test_should_skip_line() {
        assert!(should_skip_line(b""));
        assert!(should_skip_line(b"#comment"));
        assert!(should_skip_line(b"track name=foo"));
        assert!(should_skip_line(b"browser position chr1:1-100"));
        assert!(!should_skip_line(b"chr1\t100\t200"));
    }
}
