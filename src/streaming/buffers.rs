//! Buffer size constants for streaming operations.
//!
//! These constants control memory usage vs I/O throughput tradeoffs.
//! The default sizes balance good performance with reasonable memory usage.

/// Default output buffer size (2 MB).
/// This is large enough for efficient I/O while keeping memory low.
/// Previously was 8 MB, reduced to match bedops-level memory efficiency.
pub const DEFAULT_OUTPUT_BUFFER: usize = 2 * 1024 * 1024;

/// Low-memory output buffer size (256 KB).
/// Use this when memory is extremely constrained.
pub const LOW_MEMORY_OUTPUT_BUFFER: usize = 256 * 1024;

/// Default input buffer size (256 KB).
/// Good balance for reading sorted BED files.
pub const DEFAULT_INPUT_BUFFER: usize = 256 * 1024;

/// Low-memory input buffer size (64 KB).
/// Smaller buffer for memory-constrained environments.
pub const LOW_MEMORY_INPUT_BUFFER: usize = 64 * 1024;

/// Default line buffer capacity (1 KB).
/// Sufficient for most BED lines.
pub const DEFAULT_LINE_BUFFER: usize = 1024;

/// Default chromosome name buffer (64 bytes).
/// Handles standard chromosome names (chr1, chrX, etc).
pub const DEFAULT_CHROM_BUFFER: usize = 64;

/// Returns the appropriate output buffer size based on low_memory flag.
#[inline]
pub const fn output_buffer_size(low_memory: bool) -> usize {
    if low_memory {
        LOW_MEMORY_OUTPUT_BUFFER
    } else {
        DEFAULT_OUTPUT_BUFFER
    }
}

/// Returns the appropriate input buffer size based on low_memory flag.
#[inline]
pub const fn input_buffer_size(low_memory: bool) -> usize {
    if low_memory {
        LOW_MEMORY_INPUT_BUFFER
    } else {
        DEFAULT_INPUT_BUFFER
    }
}
