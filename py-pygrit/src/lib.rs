//! Python bindings for GRIT - Genomic Range Interval Toolkit.
//!
//! This module provides Python bindings for high-performance genomic interval
//! operations implemented in Rust. All core algorithms run in Rust with the
//! GIL released for maximum parallelism.

mod errors;

use errors::to_py_err;
use numpy::{PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray2};
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::path::PathBuf;

// Re-export from main crate
use grit_genomics::bed::{
    parse_intervals as rs_parse_intervals, read_intervals as rs_read_intervals, BedError,
};
use grit_genomics::commands::{
    ComplementCommand, FastSortCommand, GenerateCommand, GenerateConfig, GenerateMode,
    IntersectCommand as RsIntersectCommand, JaccardCommand, MergeCommand as RsMergeCommand,
    SizeSpec, SlopCommand, SortMode, StreamingClosestCommand, StreamingCoverageCommand,
    StreamingGenomecovCommand, StreamingGenomecovMode, StreamingIntersectCommand,
    StreamingMergeCommand, StreamingMultiinterCommand, StreamingSubtractCommand,
    StreamingWindowCommand,
};
use grit_genomics::genome::Genome;
use grit_genomics::index::IntervalIndex as RsIntervalIndex;
use grit_genomics::interval::Interval as RsInterval;

// ============================================================================
// Core Types
// ============================================================================

/// A genomic interval with chromosome, start, and end coordinates.
///
/// Coordinates are 0-based, half-open (BED format).
///
/// Example:
///     >>> iv = Interval("chr1", 100, 200)
///     >>> len(iv)
///     100
///     >>> iv.overlaps(Interval("chr1", 150, 250))
///     True
#[pyclass]
#[derive(Clone)]
pub struct Interval {
    #[pyo3(get, set)]
    pub chrom: String,
    #[pyo3(get, set)]
    pub start: u64,
    #[pyo3(get, set)]
    pub end: u64,
}

#[pymethods]
impl Interval {
    #[new]
    fn new(chrom: String, start: u64, end: u64) -> PyResult<Self> {
        if start > end {
            return Err(PyValueError::new_err(format!(
                "start ({}) must be <= end ({})",
                start, end
            )));
        }
        Ok(Self { chrom, start, end })
    }

    fn __repr__(&self) -> String {
        format!("Interval('{}', {}, {})", self.chrom, self.start, self.end)
    }

    fn __str__(&self) -> String {
        format!("{}\t{}\t{}", self.chrom, self.start, self.end)
    }

    fn __len__(&self) -> usize {
        (self.end.saturating_sub(self.start)) as usize
    }

    fn __eq__(&self, other: &Interval) -> bool {
        self.chrom == other.chrom && self.start == other.start && self.end == other.end
    }

    fn __hash__(&self) -> u64 {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        let mut hasher = DefaultHasher::new();
        self.chrom.hash(&mut hasher);
        self.start.hash(&mut hasher);
        self.end.hash(&mut hasher);
        hasher.finish()
    }

    /// Check if this interval overlaps with another.
    fn overlaps(&self, other: &Interval) -> bool {
        self.chrom == other.chrom && self.start < other.end && other.start < self.end
    }

    /// Get the overlap length with another interval.
    fn overlap_length(&self, other: &Interval) -> u64 {
        if !self.overlaps(other) {
            return 0;
        }
        let start = self.start.max(other.start);
        let end = self.end.min(other.end);
        end - start
    }

    /// Get the distance to another interval (0 if overlapping, None if different chromosomes).
    fn distance_to(&self, other: &Interval) -> Option<u64> {
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

    /// Convert to tuple (chrom, start, end).
    fn to_tuple(&self) -> (String, u64, u64) {
        (self.chrom.clone(), self.start, self.end)
    }
}

impl From<RsInterval> for Interval {
    fn from(i: RsInterval) -> Self {
        Self {
            chrom: i.chrom,
            start: i.start,
            end: i.end,
        }
    }
}

impl From<&Interval> for RsInterval {
    fn from(i: &Interval) -> Self {
        RsInterval::new(&i.chrom, i.start, i.end)
    }
}

/// A collection of genomic intervals.
///
/// Provides methods for bulk operations like merge, intersect, and sorting.
#[pyclass]
pub struct IntervalSet {
    intervals: Vec<RsInterval>,
}

#[pymethods]
impl IntervalSet {
    #[new]
    fn new() -> Self {
        Self {
            intervals: Vec::new(),
        }
    }

    /// Create an IntervalSet from a list of Interval objects.
    #[staticmethod]
    fn from_intervals(intervals: Vec<Interval>) -> Self {
        Self {
            intervals: intervals.iter().map(RsInterval::from).collect(),
        }
    }

    fn __len__(&self) -> usize {
        self.intervals.len()
    }

    fn __repr__(&self) -> String {
        format!("IntervalSet({} intervals)", self.intervals.len())
    }

    fn __getitem__(&self, idx: usize) -> PyResult<Interval> {
        self.intervals
            .get(idx)
            .map(|i| Interval::from(i.clone()))
            .ok_or_else(|| PyValueError::new_err("Index out of bounds"))
    }

    /// Add an interval.
    fn add(&mut self, interval: Interval) {
        self.intervals.push(RsInterval::from(&interval));
    }

    /// Convert to a list of Interval objects.
    fn to_list(&self) -> Vec<Interval> {
        self.intervals
            .iter()
            .map(|i| Interval::from(i.clone()))
            .collect()
    }

    /// Merge overlapping intervals.
    #[pyo3(signature = (distance = 0))]
    fn merge(&self, distance: u64) -> Self {
        let cmd = RsMergeCommand::new().with_distance(distance);
        let merged = cmd.merge(self.intervals.clone());
        Self { intervals: merged }
    }

    /// Find intersections with another IntervalSet.
    #[pyo3(signature = (other, fraction = None, reciprocal = false))]
    fn intersect(&self, other: &IntervalSet, fraction: Option<f64>, reciprocal: bool) -> Self {
        let mut cmd = RsIntersectCommand::new();
        cmd.fraction_a = fraction;
        cmd.reciprocal = reciprocal;

        let results =
            cmd.find_intersections_parallel(self.intervals.clone(), other.intervals.clone());

        let mut seen = std::collections::HashSet::new();
        let intervals: Vec<RsInterval> = results
            .into_iter()
            .filter(|r| seen.insert(r.a_index))
            .map(|r| r.a_interval)
            .collect();

        Self { intervals }
    }

    /// Find intervals with no overlap.
    fn non_overlapping(&self, other: &IntervalSet) -> Self {
        let mut cmd = RsIntersectCommand::new();
        cmd.no_overlap = true;

        let b_index = RsIntervalIndex::from_intervals(other.intervals.clone());
        let results = cmd.find_intersections(&self.intervals, &b_index);

        let intervals: Vec<RsInterval> = results.into_iter().map(|r| r.a_interval).collect();

        Self { intervals }
    }

    /// Sort intervals by chromosome and start position.
    fn sort(&mut self) {
        self.intervals.sort();
    }

    /// Convert to NumPy array (start, end only).
    fn to_numpy<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<i64>>> {
        let n = self.intervals.len();
        let mut data = Vec::with_capacity(n * 2);

        for interval in &self.intervals {
            data.push(interval.start as i64);
            data.push(interval.end as i64);
        }

        Ok(PyArray1::from_vec(py, data)
            .reshape([n, 2])
            .map_err(|e| PyValueError::new_err(format!("Failed to reshape: {}", e)))?)
    }
}

// ============================================================================
// File-Based Streaming API
// ============================================================================

/// Helper to parse BED output buffer to intervals.
fn parse_bed_output(buffer: &[u8]) -> PyResult<Vec<Interval>> {
    let content =
        std::str::from_utf8(buffer).map_err(|e| PyValueError::new_err(e.to_string()))?;

    let mut intervals = Vec::new();
    for line in content.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 3 {
            let chrom = fields[0].to_string();
            let start: u64 = fields[1]
                .parse()
                .map_err(|_| PyValueError::new_err(format!("Invalid start: {}", fields[1])))?;
            let end: u64 = fields[2]
                .parse()
                .map_err(|_| PyValueError::new_err(format!("Invalid end: {}", fields[2])))?;
            intervals.push(Interval { chrom, start, end });
        }
    }
    Ok(intervals)
}

/// Intersect two BED files using streaming algorithm.
///
/// This is the recommended way to intersect large BED files. Uses O(k) memory
/// where k = maximum number of overlapping intervals at any point.
///
/// Args:
///     a: Path to first BED file
///     b: Path to second BED file
///     output: Optional output file path. If None, returns list of intervals.
///     write_a: Include original A record in output (-wa flag)
///     write_b: Include original B record in output (-wb flag)
///     fraction: Minimum overlap fraction for A (-f flag)
///     reciprocal: Require reciprocal fraction overlap (-r flag)
///     count: Report overlap count instead of intervals (-c flag)
///     unique: Report each A interval only once (-u flag)
///     no_overlap: Report A intervals with no overlap (-v flag)
///
/// Returns:
///     List of Interval objects if output is None, otherwise None.
///
/// Example:
///     >>> results = pygrit.intersect("a.bed", "b.bed")
///     >>> pygrit.intersect("a.bed", "b.bed", output="out.bed")  # writes to file
#[pyfunction]
#[pyo3(signature = (
    a,
    b,
    output = None,
    write_a = false,
    write_b = false,
    fraction = None,
    reciprocal = false,
    count = false,
    unique = false,
    no_overlap = false
))]
pub fn intersect(
    py: Python<'_>,
    a: &str,
    b: &str,
    output: Option<&str>,
    write_a: bool,
    write_b: bool,
    fraction: Option<f64>,
    reciprocal: bool,
    count: bool,
    unique: bool,
    no_overlap: bool,
) -> PyResult<Option<Vec<Interval>>> {
    // Release GIL for heavy computation
    let result = py
        .allow_threads(|| -> Result<Vec<u8>, BedError> {
            let a_path = PathBuf::from(a);
            let b_path = PathBuf::from(b);

            let mut cmd = StreamingIntersectCommand::new();
            cmd.write_a = write_a;
            cmd.write_b = write_b;
            cmd.fraction_a = fraction;
            cmd.reciprocal = reciprocal;
            cmd.count = count;
            cmd.unique = unique;
            cmd.no_overlap = no_overlap;
            cmd.assume_sorted = true;

            let mut buffer = Vec::new();
            cmd.run(&a_path, &b_path, &mut buffer)?;
            Ok(buffer)
        })
        .map_err(to_py_err)?;

    if let Some(output_path) = output {
        std::fs::write(output_path, &result).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(None)
    } else if count {
        // Count mode returns different format - return as string instead
        Err(PyValueError::new_err(
            "count=True requires output file path",
        ))
    } else {
        let intervals = parse_bed_output(&result)?;
        Ok(Some(intervals))
    }
}

/// Merge overlapping intervals in a BED file.
///
/// Uses streaming algorithm with O(k) memory complexity.
///
/// Args:
///     input: Path to input BED file
///     output: Optional output file path. If None, returns list of intervals.
///     distance: Maximum distance between intervals to merge (default: 0)
///     strand: Merge only intervals on the same strand
///
/// Returns:
///     List of Interval objects if output is None, otherwise None.
///
/// Example:
///     >>> merged = pygrit.merge("input.bed", distance=100)
///     >>> pygrit.merge("input.bed", output="merged.bed")
#[pyfunction]
#[pyo3(signature = (input, output = None, distance = 0, strand = false))]
pub fn merge(
    py: Python<'_>,
    input: &str,
    output: Option<&str>,
    distance: u64,
    strand: bool,
) -> PyResult<Option<Vec<Interval>>> {
    let result = py
        .allow_threads(|| -> Result<Vec<u8>, BedError> {
            let input_path = PathBuf::from(input);

            let mut cmd = StreamingMergeCommand::new();
            cmd.distance = distance;
            cmd.strand_specific = strand;

            let mut buffer = Vec::new();
            cmd.run(&input_path, &mut buffer)?;
            Ok(buffer)
        })
        .map_err(to_py_err)?;

    if let Some(output_path) = output {
        std::fs::write(output_path, &result).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(None)
    } else {
        let intervals = parse_bed_output(&result)?;
        Ok(Some(intervals))
    }
}

/// Subtract B intervals from A intervals.
///
/// Args:
///     a: Path to file A
///     b: Path to file B
///     output: Optional output file path
///     remove_entire: Remove entire A interval if any overlap (-A flag)
///     fraction: Minimum overlap fraction
///     reciprocal: Require reciprocal fraction overlap
///
/// Returns:
///     List of Interval objects if output is None, otherwise None.
#[pyfunction]
#[pyo3(signature = (a, b, output = None, remove_entire = false, fraction = None, reciprocal = false))]
pub fn subtract(
    py: Python<'_>,
    a: &str,
    b: &str,
    output: Option<&str>,
    remove_entire: bool,
    fraction: Option<f64>,
    reciprocal: bool,
) -> PyResult<Option<Vec<Interval>>> {
    let result = py
        .allow_threads(|| -> Result<Vec<u8>, BedError> {
            let a_path = PathBuf::from(a);
            let b_path = PathBuf::from(b);

            let mut cmd = StreamingSubtractCommand::new();
            cmd.remove_entire = remove_entire;
            cmd.fraction = fraction;
            cmd.reciprocal = reciprocal;

            let mut buffer = Vec::new();
            cmd.run(&a_path, &b_path, &mut buffer)?;
            Ok(buffer)
        })
        .map_err(to_py_err)?;

    if let Some(output_path) = output {
        std::fs::write(output_path, &result).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(None)
    } else {
        let intervals = parse_bed_output(&result)?;
        Ok(Some(intervals))
    }
}

/// Calculate coverage of A intervals by B intervals.
///
/// Args:
///     a: Path to file A (regions)
///     b: Path to file B (reads/features)
///     output: Optional output file path
///     histogram: Report depth histogram
///     mean: Report mean depth
///
/// Returns:
///     Coverage output as string if output is None, otherwise None.
#[pyfunction]
#[pyo3(signature = (a, b, output = None, histogram = false, mean = false))]
pub fn coverage(
    py: Python<'_>,
    a: &str,
    b: &str,
    output: Option<&str>,
    histogram: bool,
    mean: bool,
) -> PyResult<Option<String>> {
    let result = py
        .allow_threads(|| -> Result<Vec<u8>, BedError> {
            let a_path = PathBuf::from(a);
            let b_path = PathBuf::from(b);

            let mut cmd = StreamingCoverageCommand::new();
            cmd.histogram = histogram;
            cmd.mean = mean;

            let mut buffer = Vec::new();
            cmd.run(a_path, b_path, &mut buffer)?;
            Ok(buffer)
        })
        .map_err(to_py_err)?;

    if let Some(output_path) = output {
        std::fs::write(output_path, &result).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(None)
    } else {
        let output_str =
            String::from_utf8(result).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Some(output_str))
    }
}

/// Find closest interval in B for each interval in A.
///
/// Args:
///     a: Path to file A
///     b: Path to file B
///     output: Optional output file path
///     ignore_overlaps: Don't report overlapping intervals
///     ignore_upstream: Ignore upstream intervals
///     ignore_downstream: Ignore downstream intervals
///
/// Returns:
///     Closest output as string if output is None, otherwise None.
#[pyfunction]
#[pyo3(signature = (a, b, output = None, ignore_overlaps = false, ignore_upstream = false, ignore_downstream = false))]
pub fn closest(
    py: Python<'_>,
    a: &str,
    b: &str,
    output: Option<&str>,
    ignore_overlaps: bool,
    ignore_upstream: bool,
    ignore_downstream: bool,
) -> PyResult<Option<String>> {
    let result = py
        .allow_threads(|| -> Result<Vec<u8>, BedError> {
            let a_path = PathBuf::from(a);
            let b_path = PathBuf::from(b);

            let mut cmd = StreamingClosestCommand::new();
            cmd.ignore_overlaps = ignore_overlaps;
            cmd.ignore_upstream = ignore_upstream;
            cmd.ignore_downstream = ignore_downstream;

            let mut buffer = Vec::new();
            cmd.run(a_path, b_path, &mut buffer)?;
            Ok(buffer)
        })
        .map_err(to_py_err)?;

    if let Some(output_path) = output {
        std::fs::write(output_path, &result).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(None)
    } else {
        let output_str =
            String::from_utf8(result).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Some(output_str))
    }
}

/// Find intervals within a window distance.
///
/// Args:
///     a: Path to file A
///     b: Path to file B
///     output: Optional output file path
///     window: Window size in base pairs (default: 1000)
///     left: Left window size (overrides window)
///     right: Right window size (overrides window)
///     count: Report count of overlaps
///     no_overlap: Report only non-overlapping
///
/// Returns:
///     Window output as string if output is None, otherwise None.
#[pyfunction]
#[pyo3(signature = (a, b, output = None, window = 1000, left = None, right = None, count = false, no_overlap = false))]
pub fn window(
    py: Python<'_>,
    a: &str,
    b: &str,
    output: Option<&str>,
    window: u64,
    left: Option<u64>,
    right: Option<u64>,
    count: bool,
    no_overlap: bool,
) -> PyResult<Option<String>> {
    let result = py
        .allow_threads(|| -> Result<Vec<u8>, BedError> {
            let a_path = PathBuf::from(a);
            let b_path = PathBuf::from(b);

            let mut cmd = StreamingWindowCommand::new();
            cmd.window = window;
            cmd.left = left;
            cmd.right = right;
            cmd.count = count;
            cmd.no_overlap = no_overlap;

            let mut buffer = Vec::new();
            cmd.run(a_path, b_path, &mut buffer)?;
            Ok(buffer)
        })
        .map_err(to_py_err)?;

    if let Some(output_path) = output {
        std::fs::write(output_path, &result).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(None)
    } else {
        let output_str =
            String::from_utf8(result).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Some(output_str))
    }
}

/// Sort a BED file by chromosome and position.
///
/// Args:
///     input: Path to input BED file
///     output: Optional output file path
///     genome: Optional genome file for chromosome ordering
///     reverse: Reverse the sort order
///
/// Returns:
///     Sorted output as string if output is None, otherwise None.
#[pyfunction]
#[pyo3(signature = (input, output = None, genome = None, reverse = false))]
pub fn sort(
    py: Python<'_>,
    input: &str,
    output: Option<&str>,
    genome: Option<&str>,
    reverse: bool,
) -> PyResult<Option<String>> {
    let result = py
        .allow_threads(|| -> Result<Vec<u8>, BedError> {
            let input_path = PathBuf::from(input);

            let mut cmd = FastSortCommand::new();
            cmd.reverse = reverse;

            let cmd = if let Some(genome_path) = genome {
                let genome_data = Genome::from_file(genome_path)?;
                cmd.with_genome(&genome_data)
            } else {
                cmd
            };

            let mut buffer = Vec::new();
            cmd.run(&input_path, &mut buffer)?;
            Ok(buffer)
        })
        .map_err(to_py_err)?;

    if let Some(output_path) = output {
        std::fs::write(output_path, &result).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(None)
    } else {
        let output_str =
            String::from_utf8(result).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Some(output_str))
    }
}

/// Extend intervals by a given number of bases.
///
/// Args:
///     input: Path to input BED file
///     genome: Path to genome file (chromosome sizes)
///     output: Optional output file path
///     both: Extend both sides by this many bases
///     left: Extend left/upstream by this many bases
///     right: Extend right/downstream by this many bases
///     strand: Use strand info (left=upstream, right=downstream)
///     pct: Interpret values as fraction of interval size
///
/// Returns:
///     Slop output as string if output is None, otherwise None.
#[pyfunction]
#[pyo3(signature = (input, genome, output = None, both = 0.0, left = None, right = None, strand = false, pct = false))]
pub fn slop(
    py: Python<'_>,
    input: &str,
    genome: &str,
    output: Option<&str>,
    both: f64,
    left: Option<f64>,
    right: Option<f64>,
    strand: bool,
    pct: bool,
) -> PyResult<Option<String>> {
    let result = py
        .allow_threads(|| -> Result<Vec<u8>, BedError> {
            let input_path = PathBuf::from(input);
            let genome_data = Genome::from_file(genome)?;

            let mut cmd = SlopCommand::new();
            cmd.both = both;
            cmd.left = left;
            cmd.right = right;
            cmd.strand = strand;
            cmd.pct = pct;

            let mut buffer = Vec::new();
            cmd.run(&input_path, &genome_data, &mut buffer)?;
            Ok(buffer)
        })
        .map_err(to_py_err)?;

    if let Some(output_path) = output {
        std::fs::write(output_path, &result).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(None)
    } else {
        let output_str =
            String::from_utf8(result).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Some(output_str))
    }
}

/// Return intervals NOT covered by the input BED file.
///
/// Args:
///     input: Path to input BED file
///     genome: Path to genome file (chromosome sizes)
///     output: Optional output file path
///
/// Returns:
///     Complement output as string if output is None, otherwise None.
#[pyfunction]
#[pyo3(signature = (input, genome, output = None))]
pub fn complement(
    py: Python<'_>,
    input: &str,
    genome: &str,
    output: Option<&str>,
) -> PyResult<Option<String>> {
    let result = py
        .allow_threads(|| -> Result<Vec<u8>, BedError> {
            let input_path = PathBuf::from(input);
            let genome_data = Genome::from_file(genome)?;

            let cmd = ComplementCommand::new().with_assume_sorted(true);

            let file = std::fs::File::open(&input_path)?;
            let reader = grit_genomics::bed::BedReader::new(file);

            let mut buffer = Vec::new();
            cmd.complement_streaming(reader, &genome_data, &mut buffer)?;
            Ok(buffer)
        })
        .map_err(to_py_err)?;

    if let Some(output_path) = output {
        std::fs::write(output_path, &result).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(None)
    } else {
        let output_str =
            String::from_utf8(result).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Some(output_str))
    }
}

/// Compute genome-wide coverage.
///
/// Args:
///     input: Path to input BED file
///     genome: Path to genome file (chromosome sizes)
///     output: Optional output file path
///     per_base: Report depth at each position (1-based)
///     bg: Report BedGraph format (non-zero only)
///     bga: Report BedGraph format (including zero coverage)
///     scale: Scale depth by factor
///
/// Returns:
///     Genomecov output as string if output is None, otherwise None.
#[pyfunction]
#[pyo3(signature = (input, genome, output = None, per_base = false, bg = false, bga = false, scale = 1.0))]
pub fn genomecov(
    py: Python<'_>,
    input: &str,
    genome: &str,
    output: Option<&str>,
    per_base: bool,
    bg: bool,
    bga: bool,
    scale: f64,
) -> PyResult<Option<String>> {
    let result = py
        .allow_threads(|| -> Result<Vec<u8>, BedError> {
            let input_path = PathBuf::from(input);
            let genome_data = Genome::from_file(genome)?;

            let mode = if per_base {
                StreamingGenomecovMode::PerBase
            } else if bg {
                StreamingGenomecovMode::BedGraph
            } else if bga {
                StreamingGenomecovMode::BedGraphAll
            } else {
                StreamingGenomecovMode::Histogram
            };

            let cmd = StreamingGenomecovCommand::new()
                .with_mode(mode)
                .with_scale(scale)
                .with_assume_sorted(true);

            let mut buffer = Vec::new();
            cmd.run(&input_path, &genome_data, &mut buffer)?;
            Ok(buffer)
        })
        .map_err(to_py_err)?;

    if let Some(output_path) = output {
        std::fs::write(output_path, &result).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(None)
    } else {
        let output_str =
            String::from_utf8(result).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Some(output_str))
    }
}

/// Calculate Jaccard similarity between two BED files.
///
/// Args:
///     a: Path to file A
///     b: Path to file B
///     output: Optional output file path
///
/// Returns:
///     Jaccard output as string if output is None, otherwise None.
///     Format: intersection\tunion\tjaccard\tn_intersections
#[pyfunction]
#[pyo3(signature = (a, b, output = None))]
pub fn jaccard(
    py: Python<'_>,
    a: &str,
    b: &str,
    output: Option<&str>,
) -> PyResult<Option<String>> {
    let result = py
        .allow_threads(|| -> Result<Vec<u8>, BedError> {
            let a_path = PathBuf::from(a);
            let b_path = PathBuf::from(b);

            let cmd = JaccardCommand::new();

            let mut buffer = Vec::new();
            cmd.run(&a_path, &b_path, &mut buffer)?;
            Ok(buffer)
        })
        .map_err(to_py_err)?;

    if let Some(output_path) = output {
        std::fs::write(output_path, &result).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(None)
    } else {
        let output_str =
            String::from_utf8(result).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Some(output_str))
    }
}

/// Identify common intervals across multiple BED files.
///
/// Args:
///     inputs: List of input BED file paths
///     output: Optional output file path
///     cluster: Only output intervals found in all files
///
/// Returns:
///     Multiinter output as string if output is None, otherwise None.
#[pyfunction]
#[pyo3(signature = (inputs, output = None, cluster = false))]
pub fn multiinter(
    py: Python<'_>,
    inputs: Vec<String>,
    output: Option<&str>,
    cluster: bool,
) -> PyResult<Option<String>> {
    if inputs.len() < 2 {
        return Err(PyValueError::new_err(
            "multiinter requires at least 2 input files",
        ));
    }

    let result = py
        .allow_threads(|| -> Result<Vec<u8>, BedError> {
            let input_paths: Vec<PathBuf> = inputs.iter().map(PathBuf::from).collect();

            let mut cmd = StreamingMultiinterCommand::new();
            cmd.cluster = cluster;
            cmd.assume_sorted = true;

            let mut buffer = Vec::new();
            cmd.run(&input_paths, &mut buffer)?;
            Ok(buffer)
        })
        .map_err(to_py_err)?;

    if let Some(output_path) = output {
        std::fs::write(output_path, &result).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(None)
    } else {
        let output_str =
            String::from_utf8(result).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Some(output_str))
    }
}

/// Generate synthetic BED datasets for benchmarking.
///
/// Args:
///     output_dir: Output directory for generated files
///     num_intervals: Number of intervals to generate
///     seed: Random seed for reproducibility
///     mode: Generation mode (balanced, skewed-a-gt-b, skewed-b-gt-a, identical, clustered)
///     sorted: Whether to sort output (yes, no, auto)
///     len_min: Minimum interval length
///     len_max: Maximum interval length
///
/// Returns:
///     Dictionary with generation statistics.
#[pyfunction]
#[pyo3(signature = (output_dir, num_intervals = 1000000, seed = 42, mode = "balanced", sorted = "auto", len_min = 50, len_max = 1000))]
pub fn generate(
    py: Python<'_>,
    output_dir: &str,
    num_intervals: u64,
    seed: u64,
    mode: &str,
    sorted: &str,
    len_min: u32,
    len_max: u32,
) -> PyResult<pyo3::Py<pyo3::types::PyDict>> {
    let gen_mode = GenerateMode::from_str(mode).ok_or_else(|| {
        PyValueError::new_err(format!(
            "Invalid mode: {}. Use: balanced, skewed-a-gt-b, skewed-b-gt-a, identical, clustered",
            mode
        ))
    })?;

    let sort_mode = SortMode::from_str(sorted).ok_or_else(|| {
        PyValueError::new_err(format!(
            "Invalid sorted value: {}. Use: yes, no, auto",
            sorted
        ))
    })?;

    let stats = py
        .allow_threads(|| -> Result<grit_genomics::commands::GenerateStats, BedError> {
            let config = GenerateConfig {
                output_dir: PathBuf::from(output_dir),
                sizes: vec![SizeSpec { count: num_intervals }],
                seed,
                mode: gen_mode,
                sorted: sort_mode,
                custom_a: None,
                custom_b: None,
                hotspot_frac: 0.05,
                hotspot_weight: 0.80,
                len_min,
                len_max,
                force: true,
            };

            let cmd = GenerateCommand::new(config);
            cmd.run()
        })
        .map_err(to_py_err)?;

    // Convert stats to Python dict
    let dict = pyo3::types::PyDict::new(py);
    dict.set_item("files_generated", stats.total_files)?;
    dict.set_item("total_intervals", stats.total_intervals)?;
    dict.set_item("elapsed_seconds", stats.elapsed_secs)?;

    Ok(dict.into())
}

// ============================================================================
// I/O Utilities
// ============================================================================

/// Read intervals from a BED file.
///
/// Args:
///     path: Path to BED file
///
/// Returns:
///     IntervalSet containing all intervals from the file.
#[pyfunction]
fn read_bed(path: &str) -> PyResult<IntervalSet> {
    let intervals = rs_read_intervals(path)
        .map_err(|e| PyIOError::new_err(format!("Failed to read BED file: {}", e)))?;
    Ok(IntervalSet { intervals })
}

/// Parse intervals from a string.
///
/// Args:
///     content: BED-formatted string content
///
/// Returns:
///     IntervalSet containing parsed intervals.
#[pyfunction]
fn parse_bed(content: &str) -> PyResult<IntervalSet> {
    let intervals = rs_parse_intervals(content)
        .map_err(|e| PyValueError::new_err(format!("Failed to parse BED content: {}", e)))?;
    Ok(IntervalSet { intervals })
}

/// Create an IntervalSet from a NumPy array.
///
/// Args:
///     chrom: Chromosome name for all intervals
///     arr: NumPy array with shape (n, 2) containing start and end coordinates
///
/// Returns:
///     IntervalSet with intervals from the array.
#[pyfunction]
fn from_numpy(_py: Python<'_>, chrom: &str, arr: PyReadonlyArray2<i64>) -> PyResult<IntervalSet> {
    let arr = arr.as_array();
    let shape = arr.shape();

    if shape.len() != 2 || shape[1] != 2 {
        return Err(PyValueError::new_err(
            "Array must have shape (n, 2) with start and end columns",
        ));
    }

    let mut intervals = Vec::with_capacity(shape[0]);
    for i in 0..shape[0] {
        let start = arr[[i, 0]] as u64;
        let end = arr[[i, 1]] as u64;
        intervals.push(RsInterval::new(chrom, start, end));
    }

    Ok(IntervalSet { intervals })
}

// ============================================================================
// Module Definition
// ============================================================================

/// GRIT: Genomic Range Interval Toolkit
///
/// High-performance genomic interval operations implemented in Rust.
///
/// Example:
///     >>> import pygrit
///     >>>
///     >>> # File-based streaming (recommended for large files)
///     >>> results = pygrit.intersect("a.bed", "b.bed")
///     >>> pygrit.merge("input.bed", output="merged.bed", distance=100)
///     >>>
///     >>> # In-memory operations
///     >>> iv = pygrit.Interval("chr1", 100, 200)
///     >>> intervals = pygrit.read_bed("file.bed")
#[pymodule]
fn pygrit(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Core types
    m.add_class::<Interval>()?;
    m.add_class::<IntervalSet>()?;

    // File-based streaming functions
    m.add_function(wrap_pyfunction!(intersect, m)?)?;
    m.add_function(wrap_pyfunction!(merge, m)?)?;
    m.add_function(wrap_pyfunction!(subtract, m)?)?;
    m.add_function(wrap_pyfunction!(coverage, m)?)?;
    m.add_function(wrap_pyfunction!(closest, m)?)?;
    m.add_function(wrap_pyfunction!(window, m)?)?;
    m.add_function(wrap_pyfunction!(sort, m)?)?;
    m.add_function(wrap_pyfunction!(slop, m)?)?;
    m.add_function(wrap_pyfunction!(complement, m)?)?;
    m.add_function(wrap_pyfunction!(genomecov, m)?)?;
    m.add_function(wrap_pyfunction!(jaccard, m)?)?;
    m.add_function(wrap_pyfunction!(multiinter, m)?)?;
    m.add_function(wrap_pyfunction!(generate, m)?)?;

    // I/O utilities
    m.add_function(wrap_pyfunction!(read_bed, m)?)?;
    m.add_function(wrap_pyfunction!(parse_bed, m)?)?;
    m.add_function(wrap_pyfunction!(from_numpy, m)?)?;

    // Version
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
