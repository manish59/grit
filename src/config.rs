//! Global configuration for GRIT runtime behavior.
//!
//! This module provides thread-safe global configuration that affects
//! parsing and interval semantics without adding overhead to hot loops.

use std::sync::atomic::{AtomicBool, Ordering};

/// Global flag for bedtools-compatible zero-length interval handling.
///
/// When enabled, zero-length intervals (start == end) are normalized to
/// 1bp intervals (end = start + 1) during parsing to match bedtools behavior.
///
/// This is set once at startup and read during parsing. The atomic load
/// has negligible overhead compared to the actual parsing work.
static BEDTOOLS_COMPATIBLE: AtomicBool = AtomicBool::new(false);

/// Enable bedtools-compatible mode.
///
/// When enabled, zero-length intervals (start == end) are normalized to
/// 1bp intervals during BED parsing. This matches bedtools behavior where
/// zero-length intervals still participate in overlap calculations.
///
/// # Example
///
/// ```
/// use grit_genomics::config;
///
/// // Enable at startup before any parsing
/// config::set_bedtools_compatible(true);
///
/// // Now parsing will normalize zero-length intervals
/// // chr1  100  100  ->  chr1  100  101
/// ```
#[inline]
pub fn set_bedtools_compatible(enabled: bool) {
    BEDTOOLS_COMPATIBLE.store(enabled, Ordering::Release);
}

/// Check if bedtools-compatible mode is enabled.
///
/// This function is called during interval parsing to determine whether
/// to normalize zero-length intervals.
#[inline]
pub fn is_bedtools_compatible() -> bool {
    BEDTOOLS_COMPATIBLE.load(Ordering::Acquire)
}

/// Normalize interval end position for bedtools compatibility.
///
/// If bedtools-compatible mode is enabled and start == end,
/// returns start + 1. Otherwise returns the original end value.
///
/// This should be called during parsing, not in inner loops.
#[inline]
pub fn normalize_end(start: u64, end: u64) -> u64 {
    if is_bedtools_compatible() && start == end {
        start + 1
    } else {
        end
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_strict_mode() {
        // Reset to default
        set_bedtools_compatible(false);
        assert!(!is_bedtools_compatible());
        assert_eq!(normalize_end(100, 100), 100);
    }

    #[test]
    fn test_bedtools_compatible_mode() {
        set_bedtools_compatible(true);
        assert!(is_bedtools_compatible());
        assert_eq!(normalize_end(100, 100), 101);
        assert_eq!(normalize_end(100, 200), 200); // Non-zero-length unchanged
        set_bedtools_compatible(false); // Reset
    }
}
