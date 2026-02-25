//! Centralized streaming utilities for GRIT.
//!
//! This module provides shared components for all streaming commands:
//! - Zero-allocation BED parsing
//! - Sort validation
//! - Efficient output formatting
//! - Active set management with automatic compaction
//!
//! All streaming commands maintain O(k) memory where k = max overlapping intervals.

pub mod active_set;
pub mod output;
pub mod parsing;
pub mod validation;

pub use active_set::{ActiveInterval, ActiveSet};
pub use output::BedWriter;
pub use parsing::{parse_bed3_bytes, parse_bed3_bytes_with_rest, parse_u64_fast, should_skip_line};
pub use validation::{
    verify_sorted, verify_sorted_reader, verify_sorted_with_genome, GenomeOrderValidator,
    SortValidator,
};
