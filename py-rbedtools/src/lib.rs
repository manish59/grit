//! Python bindings for rbedtools.

use numpy::{PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray2};
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::path::PathBuf;

// Re-export from main crate
use rbedtools::bed::{read_intervals as rs_read_intervals, parse_intervals as rs_parse_intervals};
use rbedtools::commands::{
    ClosestCommand as RsClosestCommand, CoverageCommand as RsCoverageCommand,
    IntersectCommand as RsIntersectCommand, MergeCommand as RsMergeCommand,
};
use rbedtools::index::IntervalIndex as RsIntervalIndex;
use rbedtools::interval::Interval as RsInterval;

/// A genomic interval.
#[pyclass]
#[derive(Clone)]
struct Interval {
    #[pyo3(get, set)]
    chrom: String,
    #[pyo3(get, set)]
    start: u64,
    #[pyo3(get, set)]
    end: u64,
}

#[pymethods]
impl Interval {
    #[new]
    fn new(chrom: String, start: u64, end: u64) -> Self {
        Self { chrom, start, end }
    }

    fn __repr__(&self) -> String {
        format!("Interval('{}', {}, {})", self.chrom, self.start, self.end)
    }

    fn __str__(&self) -> String {
        format!("{}\t{}\t{}", self.chrom, self.start, self.end)
    }

    /// Get the length of the interval.
    fn __len__(&self) -> u64 {
        self.end.saturating_sub(self.start)
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

    /// Get the distance to another interval.
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
#[pyclass]
struct IntervalSet {
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

    /// Get the number of intervals.
    fn __len__(&self) -> usize {
        self.intervals.len()
    }

    fn __repr__(&self) -> String {
        format!("IntervalSet({} intervals)", self.intervals.len())
    }

    /// Get an interval by index.
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
        self.intervals.iter().map(|i| Interval::from(i.clone())).collect()
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

        let results = cmd.find_intersections_parallel(
            self.intervals.clone(),
            other.intervals.clone(),
        );

        // Collect unique overlapping intervals from A
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

    /// Sort intervals.
    fn sort(&mut self) {
        self.intervals.sort();
    }

    /// Convert to NumPy array (start, end only, single chromosome).
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

/// Read intervals from a BED file.
#[pyfunction]
fn read_bed(path: &str) -> PyResult<IntervalSet> {
    let intervals = rs_read_intervals(path)
        .map_err(|e| PyIOError::new_err(format!("Failed to read BED file: {}", e)))?;
    Ok(IntervalSet { intervals })
}

/// Parse intervals from a string.
#[pyfunction]
fn parse_bed(content: &str) -> PyResult<IntervalSet> {
    let intervals = rs_parse_intervals(content)
        .map_err(|e| PyValueError::new_err(format!("Failed to parse BED content: {}", e)))?;
    Ok(IntervalSet { intervals })
}

/// Create an IntervalSet from a NumPy array.
#[pyfunction]
fn from_numpy(py: Python<'_>, chrom: &str, arr: PyReadonlyArray2<i64>) -> PyResult<IntervalSet> {
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

/// Merge intervals with optional distance parameter.
#[pyfunction]
#[pyo3(signature = (intervals, distance = 0))]
fn merge(intervals: &IntervalSet, distance: u64) -> IntervalSet {
    intervals.merge(distance)
}

/// Find intersections between two interval sets.
#[pyfunction]
#[pyo3(signature = (a, b, fraction = None, reciprocal = false))]
fn intersect(a: &IntervalSet, b: &IntervalSet, fraction: Option<f64>, reciprocal: bool) -> IntervalSet {
    a.intersect(b, fraction, reciprocal)
}

/// Python module definition.
#[pymodule]
fn rbedtools(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Interval>()?;
    m.add_class::<IntervalSet>()?;
    m.add_function(wrap_pyfunction!(read_bed, m)?)?;
    m.add_function(wrap_pyfunction!(parse_bed, m)?)?;
    m.add_function(wrap_pyfunction!(from_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(merge, m)?)?;
    m.add_function(wrap_pyfunction!(intersect, m)?)?;

    // Add version
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
