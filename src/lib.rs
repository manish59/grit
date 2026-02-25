// Clippy allows for the whole crate
#![allow(clippy::too_many_arguments)]
#![allow(clippy::should_implement_trait)]
#![allow(clippy::type_complexity)]

//! GRIT: Genomic Range Interval Toolkit
//!
//! This library provides efficient interval operations for genomic data analysis.
//!
//! # Features
//!
//! - **Parallel processing**: Uses Rayon for multi-core parallelism
//! - **Streaming I/O**: Memory-efficient processing of large files
//! - **Drop-in compatibility**: Matches bedtools CLI arguments
//!
//! # Example
//!
//! ```rust,no_run
//! use grit_genomics::{bed, interval::Interval, commands::IntersectCommand};
//!
//! // Read BED files
//! let a = bed::read_intervals("a.bed").unwrap();
//! let b = bed::read_intervals("b.bed").unwrap();
//!
//! // Find intersections
//! let cmd = IntersectCommand::new();
//! let results = cmd.find_intersections_parallel(a, b);
//! ```

pub mod bed;
pub mod commands;
pub mod config;
pub mod genome;
pub mod index;
pub mod interval;
pub mod parallel;
pub mod streaming;

// Re-export commonly used types
pub use bed::{read_intervals, read_records, BedReader};
pub use index::IntervalIndex;
pub use interval::{BedRecord, Interval, Strand};

/// Library version
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Prelude module for convenient imports.
pub mod prelude {
    pub use crate::bed::{read_intervals, read_records, BedReader};
    pub use crate::commands::{
        ClosestCommand, CoverageCommand, IntersectCommand, MergeCommand, SortCommand,
        SubtractCommand, WindowCommand,
    };
    pub use crate::index::IntervalIndex;
    pub use crate::interval::{BedRecord, Interval, Strand};
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_basic_workflow() {
        use crate::bed::parse_intervals;
        use crate::commands::MergeCommand;

        let content = "chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n";
        let intervals = parse_intervals(content).unwrap();

        let cmd = MergeCommand::new();
        let merged = cmd.merge(intervals);

        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0].start, 100);
        assert_eq!(merged[0].end, 250);
    }

    #[test]
    fn test_intersect_workflow() {
        use crate::bed::parse_intervals;
        use crate::commands::IntersectCommand;
        use crate::index::IntervalIndex;

        let a_content = "chr1\t100\t200\nchr1\t300\t400\n";
        let b_content = "chr1\t150\t250\n";

        let a = parse_intervals(a_content).unwrap();
        let b = parse_intervals(b_content).unwrap();

        let cmd = IntersectCommand::new();
        let b_index = IntervalIndex::from_intervals(b);
        let results = cmd.find_intersections(&a, &b_index);

        assert_eq!(results.len(), 1);
    }
}
