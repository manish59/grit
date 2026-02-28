//! Error handling for Python bindings.
//!
//! Converts Rust errors to appropriate Python exceptions.

use grit_genomics::bed::BedError;
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::PyErr;

/// Convert BedError to appropriate Python exception.
pub fn to_py_err(e: BedError) -> PyErr {
    match e {
        BedError::Io(io_err) => PyIOError::new_err(io_err.to_string()),
        BedError::Parse { line, message } => {
            PyValueError::new_err(format!("Parse error at line {}: {}", line, message))
        }
        BedError::InvalidFormat(msg) => PyValueError::new_err(msg),
    }
}
