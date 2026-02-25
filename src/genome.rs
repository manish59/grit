//! Genome file parser for chromosome sizes.
//!
//! Parses .genome files (tab-delimited: chrom\tsize)

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::bed::BedError;

/// Genome information containing chromosome sizes.
/// Preserves chromosome order from input file.
#[derive(Debug, Clone, Default)]
pub struct Genome {
    /// Map of chromosome name to size
    sizes: HashMap<String, u64>,
    /// Chromosome order (preserves input file order)
    order: Vec<String>,
}

impl Genome {
    /// Create an empty genome.
    pub fn new() -> Self {
        Self {
            sizes: HashMap::new(),
            order: Vec::new(),
        }
    }

    /// Load genome from a file.
    /// Format: tab-delimited with chrom\tsize per line
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, BedError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut sizes = HashMap::new();
        let mut order = Vec::new();

        for (line_num, line_result) in reader.lines().enumerate() {
            let line = line_result?;
            let line = line.trim();

            // Skip empty lines and comments
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 2 {
                return Err(BedError::Parse {
                    line: line_num + 1,
                    message: "Genome file requires two columns: chrom and size".to_string(),
                });
            }

            let chrom = fields[0].to_string();
            let size: u64 = fields[1].parse().map_err(|_| BedError::Parse {
                line: line_num + 1,
                message: format!("Invalid chromosome size: {}", fields[1]),
            })?;

            if !sizes.contains_key(&chrom) {
                order.push(chrom.clone());
            }
            sizes.insert(chrom, size);
        }

        Ok(Self { sizes, order })
    }

    /// Get the size of a chromosome.
    #[inline]
    pub fn chrom_size(&self, chrom: &str) -> Option<u64> {
        self.sizes.get(chrom).copied()
    }

    /// Check if a chromosome exists.
    #[inline]
    pub fn has_chrom(&self, chrom: &str) -> bool {
        self.sizes.contains_key(chrom)
    }

    /// Get all chromosome names in order.
    pub fn chromosomes(&self) -> impl Iterator<Item = &String> {
        self.order.iter()
    }

    /// Get number of chromosomes.
    pub fn len(&self) -> usize {
        self.sizes.len()
    }

    /// Check if empty.
    pub fn is_empty(&self) -> bool {
        self.sizes.is_empty()
    }

    /// Insert a chromosome size (appends to order if new).
    pub fn insert(&mut self, chrom: String, size: u64) {
        if !self.sizes.contains_key(&chrom) {
            self.order.push(chrom.clone());
        }
        self.sizes.insert(chrom, size);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_genome_from_file() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "chr1\t1000000").unwrap();
        writeln!(file, "chr2\t500000").unwrap();
        writeln!(file, "# comment line").unwrap();
        writeln!(file, "chr3\t250000").unwrap();

        let genome = Genome::from_file(file.path()).unwrap();

        assert_eq!(genome.chrom_size("chr1"), Some(1000000));
        assert_eq!(genome.chrom_size("chr2"), Some(500000));
        assert_eq!(genome.chrom_size("chr3"), Some(250000));
        assert_eq!(genome.chrom_size("chr4"), None);
        assert_eq!(genome.len(), 3);
    }

    #[test]
    fn test_genome_bounds() {
        let mut genome = Genome::new();
        genome.insert("chr1".to_string(), 1000);

        assert!(genome.has_chrom("chr1"));
        assert!(!genome.has_chrom("chr2"));
    }
}
