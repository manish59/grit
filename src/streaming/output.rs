//! Efficient output formatting for streaming operations.
//!
//! Uses itoa for integer formatting and ryu for float formatting
//! to avoid allocation in the hot path.

use crate::bed::BedError;
use std::io::{BufWriter, Write};

/// Buffer size for BedWriter (8MB default).
const DEFAULT_BUFFER_SIZE: usize = 8 * 1024 * 1024;

/// High-performance BED output writer.
///
/// Uses large buffering and zero-allocation formatting with itoa/ryu
/// for optimal throughput.
pub struct BedWriter<W: Write> {
    writer: BufWriter<W>,
    itoa_buf: itoa::Buffer,
    ryu_buf: ryu::Buffer,
}

impl<W: Write> BedWriter<W> {
    /// Create a new BedWriter with default 8MB buffer.
    pub fn new(output: W) -> Self {
        Self::with_capacity(DEFAULT_BUFFER_SIZE, output)
    }

    /// Create a new BedWriter with specified buffer size.
    pub fn with_capacity(capacity: usize, output: W) -> Self {
        Self {
            writer: BufWriter::with_capacity(capacity, output),
            itoa_buf: itoa::Buffer::new(),
            ryu_buf: ryu::Buffer::new(),
        }
    }

    /// Write a BED3 record (chrom, start, end).
    #[inline]
    pub fn write_bed3(&mut self, chrom: &[u8], start: u64, end: u64) -> Result<(), BedError> {
        self.writer.write_all(chrom).map_err(BedError::Io)?;
        self.writer.write_all(b"\t").map_err(BedError::Io)?;
        self.writer
            .write_all(self.itoa_buf.format(start).as_bytes())
            .map_err(BedError::Io)?;
        self.writer.write_all(b"\t").map_err(BedError::Io)?;
        self.writer
            .write_all(self.itoa_buf.format(end).as_bytes())
            .map_err(BedError::Io)?;
        Ok(())
    }

    /// Write a BED3 record followed by newline.
    #[inline]
    pub fn write_bed3_line(&mut self, chrom: &[u8], start: u64, end: u64) -> Result<(), BedError> {
        self.write_bed3(chrom, start, end)?;
        self.writer.write_all(b"\n").map_err(BedError::Io)?;
        Ok(())
    }

    /// Write a BED3 record with additional columns preserved from original line.
    #[inline]
    pub fn write_bed3_with_rest(
        &mut self,
        chrom: &[u8],
        start: u64,
        end: u64,
        original_line: &[u8],
        rest_start: usize,
    ) -> Result<(), BedError> {
        self.write_bed3(chrom, start, end)?;

        // Write rest of line if present
        if rest_start < original_line.len() {
            self.writer
                .write_all(&original_line[rest_start..])
                .map_err(BedError::Io)?;
        }

        self.writer.write_all(b"\n").map_err(BedError::Io)?;
        Ok(())
    }

    /// Write a full line as-is with newline.
    #[inline]
    pub fn write_line(&mut self, line: &[u8]) -> Result<(), BedError> {
        self.writer.write_all(line).map_err(BedError::Io)?;
        self.writer.write_all(b"\n").map_err(BedError::Io)?;
        Ok(())
    }

    /// Write raw bytes.
    #[inline]
    pub fn write_bytes(&mut self, bytes: &[u8]) -> Result<(), BedError> {
        self.writer.write_all(bytes).map_err(BedError::Io)?;
        Ok(())
    }

    /// Write a tab character.
    #[inline]
    pub fn write_tab(&mut self) -> Result<(), BedError> {
        self.writer.write_all(b"\t").map_err(BedError::Io)?;
        Ok(())
    }

    /// Write a newline character.
    #[inline]
    pub fn write_newline(&mut self) -> Result<(), BedError> {
        self.writer.write_all(b"\n").map_err(BedError::Io)?;
        Ok(())
    }

    /// Write an integer using itoa.
    #[inline]
    pub fn write_int<I: itoa::Integer>(&mut self, n: I) -> Result<(), BedError> {
        self.writer
            .write_all(self.itoa_buf.format(n).as_bytes())
            .map_err(BedError::Io)?;
        Ok(())
    }

    /// Write a float using ryu.
    #[inline]
    pub fn write_float(&mut self, f: f64) -> Result<(), BedError> {
        self.writer
            .write_all(self.ryu_buf.format(f).as_bytes())
            .map_err(BedError::Io)?;
        Ok(())
    }

    /// Write a float with exactly 7 decimal places (bedtools compatibility).
    /// Uses Rust's standard formatting which matches C printf rounding behavior.
    #[inline]
    pub fn write_float_7dp(&mut self, f: f64) -> Result<(), BedError> {
        write!(self.writer, "{:.7}", f).map_err(BedError::Io)
    }

    /// Write A\tB pair (two BED lines joined by tab).
    #[inline]
    pub fn write_pair(&mut self, a_line: &[u8], b_line: &[u8]) -> Result<(), BedError> {
        self.writer.write_all(a_line).map_err(BedError::Io)?;
        self.writer.write_all(b"\t").map_err(BedError::Io)?;
        self.writer.write_all(b_line).map_err(BedError::Io)?;
        self.writer.write_all(b"\n").map_err(BedError::Io)?;
        Ok(())
    }

    /// Flush the output buffer.
    pub fn flush(&mut self) -> Result<(), BedError> {
        self.writer.flush().map_err(BedError::Io)?;
        Ok(())
    }

    /// Get mutable reference to the underlying writer.
    pub fn inner_mut(&mut self) -> &mut BufWriter<W> {
        &mut self.writer
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_bed3() {
        let mut output = Vec::new();
        {
            let mut writer = BedWriter::new(&mut output);
            writer.write_bed3_line(b"chr1", 100, 200).unwrap();
            writer.flush().unwrap();
        }
        assert_eq!(output, b"chr1\t100\t200\n");
    }

    #[test]
    fn test_write_bed3_with_rest() {
        let mut output = Vec::new();
        let original = b"chr1\t100\t200\tname\t50\t+";
        {
            let mut writer = BedWriter::new(&mut output);
            writer
                .write_bed3_with_rest(b"chr1", 150, 250, original, 12)
                .unwrap();
            writer.flush().unwrap();
        }
        assert_eq!(output, b"chr1\t150\t250\tname\t50\t+\n");
    }

    #[test]
    fn test_write_pair() {
        let mut output = Vec::new();
        {
            let mut writer = BedWriter::new(&mut output);
            writer
                .write_pair(b"chr1\t100\t200", b"chr1\t150\t250")
                .unwrap();
            writer.flush().unwrap();
        }
        assert_eq!(output, b"chr1\t100\t200\tchr1\t150\t250\n");
    }

    #[test]
    fn test_write_float_7dp() {
        let mut output = Vec::new();
        {
            let mut writer = BedWriter::new(&mut output);
            writer.write_float_7dp(0.75).unwrap();
            writer.flush().unwrap();
        }
        let result = String::from_utf8(output).unwrap();
        assert_eq!(result.len(), 9); // "0.7500000"
        assert!(result.starts_with("0.75"));
    }
}
