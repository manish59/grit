//! Command implementations for bedtools-rs.

pub mod closest;
pub mod complement;
pub mod coverage;
pub mod fast_merge;
pub mod fast_sort;
pub mod generate;
pub mod genomecov;
pub mod intersect;
pub mod intersect_engine;
pub mod jaccard;
pub mod merge;
pub mod multiinter;
pub mod slop;
pub mod sort;
pub mod streaming_closest;
pub mod streaming_coverage;
pub mod streaming_genomecov;
pub mod streaming_intersect;
pub mod streaming_merge;
pub mod streaming_multiinter;
pub mod streaming_subtract;
pub mod streaming_window;
pub mod subtract;
pub mod window;

pub use crate::streaming::{
    verify_sorted, verify_sorted_reader, verify_sorted_with_genome, GenomeOrderValidator,
};
pub use closest::ClosestCommand;
pub use complement::ComplementCommand;
pub use coverage::CoverageCommand;
pub use fast_merge::{FastMergeCommand, FastMergeStats};
pub use fast_sort::{FastSortCommand, FastSortStats};
pub use generate::{
    GenerateCommand, GenerateConfig, GenerateMode, GenerateStats, SizeSpec, SortMode,
};
pub use genomecov::{GenomecovCommand, OutputMode as GenomecovOutputMode};
pub use intersect::IntersectCommand;
pub use intersect_engine::{ExecutionMode, IntersectConfig, IntersectEngine, IntersectStats};
pub use jaccard::JaccardCommand;
pub use merge::MergeCommand;
pub use multiinter::MultiinterCommand;
pub use slop::SlopCommand;
pub use sort::SortCommand;
pub use streaming_closest::{StreamingClosestCommand, StreamingClosestStats};
pub use streaming_coverage::StreamingCoverageCommand;
pub use streaming_genomecov::{StreamingGenomecovCommand, StreamingGenomecovMode};
pub use streaming_intersect::{StreamingIntersectCommand, StreamingStats};
pub use streaming_merge::{StreamingMergeCommand, StreamingMergeStats};
pub use streaming_multiinter::StreamingMultiinterCommand;
pub use streaming_subtract::{StreamingSubtractCommand, StreamingSubtractStats};
pub use streaming_window::{StreamingWindowCommand, StreamingWindowStats};
pub use subtract::SubtractCommand;
pub use window::WindowCommand;
