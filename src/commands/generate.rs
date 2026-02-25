//! Generate synthetic BED datasets for benchmarking.
//!
//! This module provides the `grit generate` command to create synthetic BED files
//! with various distributions (uniform, clustered) for benchmarking GRIT vs bedtools.
//!
//! Features:
//! - Human genome model (23 chromosomes, weighted by size)
//! - Uniform and clustered interval distributions
//! - Multiple dataset modes: balanced, skewed, identical, clustered
//! - External sorting for large files (50M+ intervals)

#![allow(clippy::manual_is_multiple_of)]
#![allow(clippy::manual_div_ceil)]
//! - Deterministic reproducibility via seed

use crate::bed::BedError;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Buffer size for I/O operations (8MB for better throughput)
const BUF_SIZE: usize = 8 * 1024 * 1024;

/// Chunk size for external sort (5M intervals)
const CHUNK_SIZE: usize = 5_000_000;

/// Generation mode for synthetic datasets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenerateMode {
    /// Same size A and B files
    Balanced,
    /// Large A file, small B file
    SkewedAGtB,
    /// Small A file, large B file
    SkewedBGtA,
    /// A = B (identical copies)
    Identical,
    /// 80% of intervals in hotspots (5% of genome)
    Clustered,
    /// Generate all modes
    All,
}

impl GenerateMode {
    /// Parse mode from string.
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "balanced" => Some(Self::Balanced),
            "skewed-a-gt-b" | "skewed_a_gt_b" => Some(Self::SkewedAGtB),
            "skewed-b-gt-a" | "skewed_b_gt_a" => Some(Self::SkewedBGtA),
            "identical" => Some(Self::Identical),
            "clustered" => Some(Self::Clustered),
            "all" => Some(Self::All),
            _ => None,
        }
    }

    /// Get directory name for this mode.
    pub fn dir_name(&self) -> &'static str {
        match self {
            Self::Balanced => "balanced",
            Self::SkewedAGtB => "skewed_A_gt_B",
            Self::SkewedBGtA => "skewed_B_gt_A",
            Self::Identical => "identical",
            Self::Clustered => "clustered",
            Self::All => "all",
        }
    }
}

/// Sorting mode for generated files.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortMode {
    /// Always sort output
    Yes,
    /// Never sort output
    No,
    /// Sort if count > 1M
    Auto,
}

impl SortMode {
    /// Parse sort mode from string.
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "yes" | "true" | "1" => Some(Self::Yes),
            "no" | "false" | "0" => Some(Self::No),
            "auto" => Some(Self::Auto),
            _ => None,
        }
    }

    /// Determine if we should sort for a given count.
    pub fn should_sort(&self, count: u64) -> bool {
        match self {
            Self::Yes => true,
            Self::No => false,
            Self::Auto => count >= 1_000_000,
        }
    }
}

/// Size specification (parses 1K, 1M, etc.).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SizeSpec {
    pub count: u64,
}

impl SizeSpec {
    /// Parse size from string (e.g., "1K", "5M", "100").
    pub fn from_str(s: &str) -> Option<Self> {
        let s = s.trim().to_uppercase();
        if s.is_empty() {
            return None;
        }

        let (num_part, multiplier) = if s.ends_with('K') {
            (&s[..s.len() - 1], 1_000u64)
        } else if s.ends_with('M') {
            (&s[..s.len() - 1], 1_000_000u64)
        } else if s.ends_with('G') {
            (&s[..s.len() - 1], 1_000_000_000u64)
        } else {
            (s.as_str(), 1u64)
        };

        num_part.parse::<u64>().ok().map(|n| Self {
            count: n * multiplier,
        })
    }

    /// Format size for display.
    pub fn display(&self) -> String {
        if self.count >= 1_000_000_000 && self.count % 1_000_000_000 == 0 {
            format!("{}G", self.count / 1_000_000_000)
        } else if self.count >= 1_000_000 && self.count % 1_000_000 == 0 {
            format!("{}M", self.count / 1_000_000)
        } else if self.count >= 1_000 && self.count % 1_000 == 0 {
            format!("{}K", self.count / 1_000)
        } else {
            self.count.to_string()
        }
    }
}

/// Configuration for the generate command.
#[derive(Debug, Clone)]
pub struct GenerateConfig {
    pub output_dir: PathBuf,
    pub sizes: Vec<SizeSpec>,
    pub seed: u64,
    pub mode: GenerateMode,
    pub sorted: SortMode,
    pub custom_a: Option<u64>,
    pub custom_b: Option<u64>,
    pub hotspot_frac: f64,
    pub hotspot_weight: f64,
    pub len_min: u32,
    pub len_max: u32,
    pub force: bool,
}

impl Default for GenerateConfig {
    fn default() -> Self {
        Self {
            output_dir: PathBuf::from("./grit_bench_data"),
            sizes: vec![
                SizeSpec { count: 1_000_000 },
                SizeSpec { count: 5_000_000 },
                SizeSpec { count: 10_000_000 },
                SizeSpec { count: 25_000_000 },
                SizeSpec { count: 50_000_000 },
            ],
            seed: 42,
            mode: GenerateMode::All,
            sorted: SortMode::Auto,
            custom_a: None,
            custom_b: None,
            hotspot_frac: 0.05,
            hotspot_weight: 0.80,
            len_min: 50,
            len_max: 1000,
            force: false,
        }
    }
}

/// Statistics from generate operation.
#[derive(Debug, Default, Clone)]
pub struct GenerateStats {
    pub total_intervals: u64,
    pub total_files: usize,
    pub elapsed_secs: f64,
}

impl std::fmt::Display for GenerateStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} intervals in {} files ({:.1}s)",
            format_count(self.total_intervals),
            self.total_files,
            self.elapsed_secs
        )
    }
}

/// Human genome model with 23 chromosomes weighted by size.
struct HumanGenome {
    /// Chromosome names and sizes
    chromosomes: Vec<(&'static str, u64)>,
    /// Cumulative sizes for weighted sampling
    cumulative: Vec<u64>,
    /// Total genome size
    total_size: u64,
}

impl HumanGenome {
    /// Create a new human genome model (hg38 approximate sizes).
    fn new() -> Self {
        // Approximate hg38 chromosome sizes
        let chromosomes: Vec<(&'static str, u64)> = vec![
            ("chr1", 248_956_422),
            ("chr2", 242_193_529),
            ("chr3", 198_295_559),
            ("chr4", 190_214_555),
            ("chr5", 181_538_259),
            ("chr6", 170_805_979),
            ("chr7", 159_345_973),
            ("chr8", 145_138_636),
            ("chr9", 138_394_717),
            ("chr10", 133_797_422),
            ("chr11", 135_086_622),
            ("chr12", 133_275_309),
            ("chr13", 114_364_328),
            ("chr14", 107_043_718),
            ("chr15", 101_991_189),
            ("chr16", 90_338_345),
            ("chr17", 83_257_441),
            ("chr18", 80_373_285),
            ("chr19", 58_617_616),
            ("chr20", 64_444_167),
            ("chr21", 46_709_983),
            ("chr22", 50_818_468),
            ("chrX", 156_040_895),
        ];

        let mut cumulative = Vec::with_capacity(chromosomes.len());
        let mut running_total = 0u64;
        for (_, size) in &chromosomes {
            running_total += size;
            cumulative.push(running_total);
        }

        Self {
            chromosomes,
            cumulative,
            total_size: running_total,
        }
    }

    /// Sample a chromosome weighted by size.
    /// Returns (chromosome index, chromosome size).
    #[inline]
    fn sample_chromosome(&self, rng: &mut SmallRng) -> (usize, u64) {
        let target = rng.gen_range(0..self.total_size);
        let idx = self.cumulative.partition_point(|&x| x <= target);
        (idx, self.chromosomes[idx].1)
    }

    /// Get chromosome name by index.
    #[inline]
    fn chrom_name(&self, idx: usize) -> &'static str {
        self.chromosomes[idx].0
    }
}

/// Compact interval for generation (avoids String allocation).
#[derive(Clone, Copy, Debug)]
#[repr(C)]
struct RawInterval {
    chrom_idx: u16,
    start: u32,
    end: u32,
}

/// Mapping from biological chromosome index to lexicographic sort index.
/// This ensures sorted output matches `LC_ALL=C sort -k1,1 -k2,2n -k3,3n`.
const CHROM_LEX_ORDER: [u16; 23] = [
    0,  // chr1 -> 0
    11, // chr2 -> 11
    15, // chr3 -> 15
    16, // chr4 -> 16
    17, // chr5 -> 17
    18, // chr6 -> 18
    19, // chr7 -> 19
    20, // chr8 -> 20
    21, // chr9 -> 21
    1,  // chr10 -> 1
    2,  // chr11 -> 2
    3,  // chr12 -> 3
    4,  // chr13 -> 4
    5,  // chr14 -> 5
    6,  // chr15 -> 6
    7,  // chr16 -> 7
    8,  // chr17 -> 8
    9,  // chr18 -> 9
    10, // chr19 -> 10
    12, // chr20 -> 12
    13, // chr21 -> 13
    14, // chr22 -> 14
    22, // chrX -> 22
];

impl RawInterval {
    /// Sort key matching `LC_ALL=C sort -k1,1 -k2,2n -k3,3n`.
    #[inline]
    fn sort_key(&self) -> (u16, u32, u32) {
        let lex_idx = CHROM_LEX_ORDER[self.chrom_idx as usize];
        (lex_idx, self.start, self.end)
    }
}

/// Hotspot for clustered distribution.
struct Hotspot {
    chrom_idx: u16,
    center: u64,
    radius: u64,
}

/// Generate command.
pub struct GenerateCommand {
    config: GenerateConfig,
    genome: HumanGenome,
}

impl GenerateCommand {
    /// Create a new generate command with the given config.
    pub fn new(config: GenerateConfig) -> Self {
        Self {
            config,
            genome: HumanGenome::new(),
        }
    }

    /// Run the generation.
    pub fn run(&self) -> Result<GenerateStats, BedError> {
        let start = Instant::now();
        let mut stats = GenerateStats::default();

        // Create output directory
        fs::create_dir_all(&self.config.output_dir)?;

        eprintln!("Output directory: {}", self.config.output_dir.display());

        // Handle custom A/B sizes
        if self.config.custom_a.is_some() || self.config.custom_b.is_some() {
            let a_count = self.config.custom_a.unwrap_or(1_000_000);
            let b_count = self.config.custom_b.unwrap_or(1_000_000);

            // Enable clustering if mode is Clustered
            let clustered = self.config.mode == GenerateMode::Clustered;

            let dir_name = format!(
                "custom_A{}_B{}",
                format_count(a_count),
                format_count(b_count)
            );
            let dir = self.config.output_dir.join(&dir_name);
            fs::create_dir_all(&dir)?;

            eprintln!(
                "Mode: custom{}",
                if clustered { " (clustered)" } else { "" }
            );
            eprintln!(
                "Size: A={}, B={}",
                format_count(a_count),
                format_count(b_count)
            );

            self.generate_pair(&dir, a_count, b_count, clustered, &mut stats)?;

            stats.elapsed_secs = start.elapsed().as_secs_f64();
            return Ok(stats);
        }

        // Process modes
        let modes = if self.config.mode == GenerateMode::All {
            vec![
                GenerateMode::Balanced,
                GenerateMode::SkewedAGtB,
                GenerateMode::SkewedBGtA,
                GenerateMode::Identical,
                GenerateMode::Clustered,
            ]
        } else {
            vec![self.config.mode]
        };

        for mode in modes {
            self.run_mode(mode, &mut stats)?;
        }

        stats.elapsed_secs = start.elapsed().as_secs_f64();
        eprintln!("\nComplete: {}", stats);

        Ok(stats)
    }

    /// Run generation for a specific mode.
    fn run_mode(&self, mode: GenerateMode, stats: &mut GenerateStats) -> Result<(), BedError> {
        eprintln!("\nMode: {}", mode.dir_name());
        let mode_dir = self.config.output_dir.join(mode.dir_name());
        fs::create_dir_all(&mode_dir)?;

        match mode {
            GenerateMode::Balanced => {
                for size in &self.config.sizes {
                    let size_dir = mode_dir.join(size.display());
                    fs::create_dir_all(&size_dir)?;

                    eprintln!("Size: {}", size.display());
                    self.generate_pair(&size_dir, size.count, size.count, false, stats)?;
                }
            }
            GenerateMode::SkewedAGtB => {
                // Large A, small B: (25M,1M), (10M,1M), (50M,2M)
                let pairs = [
                    (25_000_000, 1_000_000),
                    (10_000_000, 1_000_000),
                    (50_000_000, 2_000_000),
                ];
                for (a, b) in pairs {
                    let dir_name = format!("{}_vs_{}", format_count(a), format_count(b));
                    let size_dir = mode_dir.join(&dir_name);
                    fs::create_dir_all(&size_dir)?;

                    eprintln!("Size: {}", dir_name);
                    self.generate_pair(&size_dir, a, b, false, stats)?;
                }
            }
            GenerateMode::SkewedBGtA => {
                // Small A, large B: (1M,25M), (1M,10M), (2M,50M)
                let pairs = [
                    (1_000_000, 25_000_000),
                    (1_000_000, 10_000_000),
                    (2_000_000, 50_000_000),
                ];
                for (a, b) in pairs {
                    let dir_name = format!("{}_vs_{}", format_count(a), format_count(b));
                    let size_dir = mode_dir.join(&dir_name);
                    fs::create_dir_all(&size_dir)?;

                    eprintln!("Size: {}", dir_name);
                    self.generate_pair(&size_dir, a, b, false, stats)?;
                }
            }
            GenerateMode::Identical => {
                // Only 1M and 25M for identical mode
                for count in [1_000_000u64, 25_000_000] {
                    let size_dir = mode_dir.join(format_count(count));
                    fs::create_dir_all(&size_dir)?;

                    eprintln!("Size: {}", format_count(count));
                    self.generate_identical(&size_dir, count, stats)?;
                }
            }
            GenerateMode::Clustered => {
                // Only 10M and 50M for clustered mode
                for count in [10_000_000u64, 50_000_000] {
                    let size_dir = mode_dir.join(format_count(count));
                    fs::create_dir_all(&size_dir)?;

                    eprintln!("Size: {}", format_count(count));
                    self.generate_pair(&size_dir, count, count, true, stats)?;
                }
            }
            GenerateMode::All => unreachable!(),
        }

        Ok(())
    }

    /// Generate a pair of BED files (A.bed, B.bed).
    fn generate_pair(
        &self,
        dir: &Path,
        a_count: u64,
        b_count: u64,
        clustered: bool,
        stats: &mut GenerateStats,
    ) -> Result<(), BedError> {
        let a_path = dir.join("A.bed");
        let b_path = dir.join("B.bed");

        // Check if files exist
        if !self.config.force && a_path.exists() && b_path.exists() {
            eprintln!("  Skipping (files exist, use --force to overwrite)");
            return Ok(());
        }

        // Generate A with seed
        let mut rng_a = SmallRng::seed_from_u64(self.config.seed);
        eprint!("  Generating A.bed... ");
        let start_a = Instant::now();
        self.generate_file(&a_path, a_count, clustered, &mut rng_a)?;
        eprintln!("done ({:.1}s)", start_a.elapsed().as_secs_f64());
        eprintln!("  Saved: {}", a_path.display());

        // Generate B with different seed
        let mut rng_b = SmallRng::seed_from_u64(self.config.seed.wrapping_add(1));
        eprint!("  Generating B.bed... ");
        let start_b = Instant::now();
        self.generate_file(&b_path, b_count, clustered, &mut rng_b)?;
        eprintln!("done ({:.1}s)", start_b.elapsed().as_secs_f64());
        eprintln!("  Saved: {}", b_path.display());

        stats.total_intervals += a_count + b_count;
        stats.total_files += 2;

        Ok(())
    }

    /// Generate identical A and B files.
    fn generate_identical(
        &self,
        dir: &Path,
        count: u64,
        stats: &mut GenerateStats,
    ) -> Result<(), BedError> {
        let a_path = dir.join("A.bed");
        let b_path = dir.join("B.bed");

        // Check if files exist
        if !self.config.force && a_path.exists() && b_path.exists() {
            eprintln!("  Skipping (files exist, use --force to overwrite)");
            return Ok(());
        }

        // Generate A
        let mut rng = SmallRng::seed_from_u64(self.config.seed);
        eprint!("  Generating A.bed... ");
        let start = Instant::now();
        self.generate_file(&a_path, count, false, &mut rng)?;
        eprintln!("done ({:.1}s)", start.elapsed().as_secs_f64());
        eprintln!("  Saved: {}", a_path.display());

        // Copy A to B
        eprint!("  Copying to B.bed... ");
        fs::copy(&a_path, &b_path)?;
        eprintln!("done");
        eprintln!("  Saved: {}", b_path.display());

        stats.total_intervals += count * 2;
        stats.total_files += 2;

        Ok(())
    }

    /// Generate a single BED file.
    fn generate_file(
        &self,
        path: &Path,
        count: u64,
        clustered: bool,
        rng: &mut SmallRng,
    ) -> Result<(), BedError> {
        let should_sort = self.config.sorted.should_sort(count);

        if should_sort && count as usize > CHUNK_SIZE {
            // External sort for large files
            self.generate_with_external_sort(path, count, clustered, rng)
        } else if should_sort {
            // In-memory sort for medium files
            self.generate_with_memory_sort(path, count, clustered, rng)
        } else {
            // No sorting needed
            self.generate_unsorted(path, count, clustered, rng)
        }
    }

    /// Generate intervals (uniform or clustered distribution).
    fn generate_intervals(
        &self,
        count: u64,
        clustered: bool,
        rng: &mut SmallRng,
    ) -> Vec<RawInterval> {
        let mut intervals = Vec::with_capacity(count as usize);

        if clustered {
            // Generate hotspots
            let hotspots = self.create_hotspots(rng);

            // Calculate how many intervals go in hotspots
            let in_hotspot = (count as f64 * self.config.hotspot_weight) as u64;
            let uniform = count - in_hotspot;

            // Generate hotspot intervals
            for _ in 0..in_hotspot {
                let hotspot = &hotspots[rng.gen_range(0..hotspots.len())];
                let len = rng.gen_range(self.config.len_min..=self.config.len_max);

                // Sample position within hotspot
                let pos = rng.gen_range(
                    hotspot.center.saturating_sub(hotspot.radius)
                        ..hotspot.center.saturating_add(hotspot.radius),
                );

                let chrom_size = self.genome.chromosomes[hotspot.chrom_idx as usize].1;
                let start = pos.min(chrom_size - len as u64) as u32;

                intervals.push(RawInterval {
                    chrom_idx: hotspot.chrom_idx,
                    start,
                    end: start + len,
                });
            }

            // Generate uniform intervals
            for _ in 0..uniform {
                intervals.push(self.generate_uniform_interval(rng));
            }
        } else {
            // Pure uniform distribution
            for _ in 0..count {
                intervals.push(self.generate_uniform_interval(rng));
            }
        }

        intervals
    }

    /// Generate a single uniform interval.
    #[inline]
    fn generate_uniform_interval(&self, rng: &mut SmallRng) -> RawInterval {
        let (chrom_idx, chrom_size) = self.genome.sample_chromosome(rng);
        let len = rng.gen_range(self.config.len_min..=self.config.len_max);

        let max_start = chrom_size.saturating_sub(len as u64);
        let start = if max_start > 0 {
            rng.gen_range(0..max_start) as u32
        } else {
            0
        };

        RawInterval {
            chrom_idx: chrom_idx as u16,
            start,
            end: start + len,
        }
    }

    /// Create hotspots covering hotspot_frac of the genome.
    fn create_hotspots(&self, rng: &mut SmallRng) -> Vec<Hotspot> {
        let target_coverage = (self.genome.total_size as f64 * self.config.hotspot_frac) as u64;
        let mut hotspots = Vec::new();
        let mut total_coverage = 0u64;

        // Average hotspot size ~1Mb
        let avg_hotspot_size = 1_000_000u64;

        while total_coverage < target_coverage {
            let (chrom_idx, chrom_size) = self.genome.sample_chromosome(rng);
            let center = rng.gen_range(0..chrom_size);
            let radius = rng.gen_range(avg_hotspot_size / 2..avg_hotspot_size * 2);

            hotspots.push(Hotspot {
                chrom_idx: chrom_idx as u16,
                center,
                radius,
            });

            total_coverage += radius * 2;
        }

        hotspots
    }

    /// Generate without sorting.
    fn generate_unsorted(
        &self,
        path: &Path,
        count: u64,
        clustered: bool,
        rng: &mut SmallRng,
    ) -> Result<(), BedError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::with_capacity(BUF_SIZE, file);

        let intervals = self.generate_intervals(count, clustered, rng);
        self.write_intervals(&intervals, &mut writer)?;

        writer.flush()?;
        Ok(())
    }

    /// Generate with in-memory sorting.
    fn generate_with_memory_sort(
        &self,
        path: &Path,
        count: u64,
        clustered: bool,
        rng: &mut SmallRng,
    ) -> Result<(), BedError> {
        let mut intervals = self.generate_intervals(count, clustered, rng);

        // Sort by (chrom, start, end) matching sort -k1,1 -k2,2n -k3,3n
        intervals.par_sort_by_key(|i| i.sort_key());

        let file = File::create(path)?;
        let mut writer = BufWriter::with_capacity(BUF_SIZE, file);
        self.write_intervals(&intervals, &mut writer)?;

        writer.flush()?;
        Ok(())
    }

    /// Generate with external sorting for very large files.
    fn generate_with_external_sort(
        &self,
        path: &Path,
        count: u64,
        clustered: bool,
        rng: &mut SmallRng,
    ) -> Result<(), BedError> {
        // Create temp directory
        let temp_dir = tempfile::tempdir()?;

        // Phase 1: Generate and sort chunks
        let mut chunk_paths = Vec::new();
        let mut remaining = count;
        let mut chunk_idx = 0;

        while remaining > 0 {
            let chunk_size = remaining.min(CHUNK_SIZE as u64);
            let mut chunk = self.generate_intervals(chunk_size, clustered, rng);

            // Sort chunk by (chrom, start, end)
            chunk.par_sort_by_key(|i| i.sort_key());

            // Write chunk to temp file
            let chunk_path = temp_dir.path().join(format!("chunk_{}.bed", chunk_idx));
            let file = File::create(&chunk_path)?;
            let mut writer = BufWriter::with_capacity(BUF_SIZE, file);
            self.write_intervals(&chunk, &mut writer)?;
            writer.flush()?;

            chunk_paths.push(chunk_path);
            remaining -= chunk_size;
            chunk_idx += 1;

            eprint!(
                "\r  Sorting chunks... {}/{}",
                chunk_idx,
                (count as usize + CHUNK_SIZE - 1) / CHUNK_SIZE
            );
        }
        eprintln!();

        // Phase 2: K-way merge
        eprint!("  Merging... ");
        self.k_way_merge(&chunk_paths, path)?;
        eprintln!("done");

        Ok(())
    }

    /// K-way merge of sorted chunk files.
    fn k_way_merge(&self, chunk_paths: &[PathBuf], output: &Path) -> Result<(), BedError> {
        // Open all chunk files
        let mut readers: Vec<_> = chunk_paths
            .iter()
            .map(|p| BufReader::with_capacity(BUF_SIZE, File::open(p).unwrap()))
            .collect();

        // Min-heap: ((chrom_idx, start, end), chunk_idx, line)
        let mut heap: BinaryHeap<std::cmp::Reverse<((u16, u32, u32), usize, String)>> =
            BinaryHeap::new();

        // Initialize heap with first line from each chunk
        for (idx, reader) in readers.iter_mut().enumerate() {
            let mut line = String::new();
            if reader.read_line(&mut line)? > 0 {
                if let Some(key) = parse_sort_key(&line) {
                    heap.push(std::cmp::Reverse((key, idx, line)));
                }
            }
        }

        // Merge
        let file = File::create(output)?;
        let mut writer = BufWriter::with_capacity(BUF_SIZE, file);

        while let Some(std::cmp::Reverse((_key, chunk_idx, line))) = heap.pop() {
            writer.write_all(line.as_bytes())?;

            // Read next line from same chunk
            let mut next_line = String::new();
            if readers[chunk_idx].read_line(&mut next_line)? > 0 {
                if let Some(key) = parse_sort_key(&next_line) {
                    heap.push(std::cmp::Reverse((key, chunk_idx, next_line)));
                }
            }
        }

        writer.flush()?;
        Ok(())
    }

    /// Write intervals to a BufWriter.
    fn write_intervals<W: Write>(
        &self,
        intervals: &[RawInterval],
        writer: &mut BufWriter<W>,
    ) -> Result<(), BedError> {
        let mut buf = itoa::Buffer::new();

        for interval in intervals {
            let chrom = self.genome.chrom_name(interval.chrom_idx as usize);
            writer.write_all(chrom.as_bytes())?;
            writer.write_all(b"\t")?;
            writer.write_all(buf.format(interval.start).as_bytes())?;
            writer.write_all(b"\t")?;
            writer.write_all(buf.format(interval.end).as_bytes())?;
            writer.write_all(b"\n")?;
        }

        Ok(())
    }
}

/// Parse sort key from a BED line: (chrom_index, start, end).
fn parse_sort_key(line: &str) -> Option<(u16, u32, u32)> {
    let line = line.trim_end();
    let mut parts = line.split('\t');
    let chrom = parts.next()?;
    let start: u32 = parts.next()?.parse().ok()?;
    let end: u32 = parts.next()?.parse().ok()?;

    let chrom_idx = chrom_to_index(chrom);

    Some((chrom_idx, start, end))
}

/// Convert chromosome name to index for sorting.
/// Uses lexicographic order to match `sort -k1,1`.
fn chrom_to_index(chrom: &str) -> u16 {
    // Lexicographic order: chr1, chr10, chr11, ..., chr19, chr2, chr20, chr21, chr22, chr3, ..., chr9, chrX
    match chrom {
        "chr1" => 0,
        "chr10" => 1,
        "chr11" => 2,
        "chr12" => 3,
        "chr13" => 4,
        "chr14" => 5,
        "chr15" => 6,
        "chr16" => 7,
        "chr17" => 8,
        "chr18" => 9,
        "chr19" => 10,
        "chr2" => 11,
        "chr20" => 12,
        "chr21" => 13,
        "chr22" => 14,
        "chr3" => 15,
        "chr4" => 16,
        "chr5" => 17,
        "chr6" => 18,
        "chr7" => 19,
        "chr8" => 20,
        "chr9" => 21,
        "chrX" => 22,
        _ => 100, // Unknown chromosome
    }
}

/// Format a count for display (e.g., 1000000 -> "1M").
fn format_count(count: u64) -> String {
    if count >= 1_000_000_000 && count % 1_000_000_000 == 0 {
        format!("{}G", count / 1_000_000_000)
    } else if count >= 1_000_000 && count % 1_000_000 == 0 {
        format!("{}M", count / 1_000_000)
    } else if count >= 1_000 && count % 1_000 == 0 {
        format!("{}K", count / 1_000)
    } else {
        count.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_size_spec_parse() {
        assert_eq!(SizeSpec::from_str("1K").unwrap().count, 1_000);
        assert_eq!(SizeSpec::from_str("5M").unwrap().count, 5_000_000);
        assert_eq!(SizeSpec::from_str("1G").unwrap().count, 1_000_000_000);
        assert_eq!(SizeSpec::from_str("100").unwrap().count, 100);
        assert_eq!(SizeSpec::from_str("  10k  ").unwrap().count, 10_000);
    }

    #[test]
    fn test_size_spec_display() {
        assert_eq!(SizeSpec { count: 1_000 }.display(), "1K");
        assert_eq!(SizeSpec { count: 5_000_000 }.display(), "5M");
        assert_eq!(
            SizeSpec {
                count: 1_000_000_000
            }
            .display(),
            "1G"
        );
        assert_eq!(SizeSpec { count: 100 }.display(), "100");
    }

    #[test]
    fn test_generate_mode_parse() {
        assert_eq!(
            GenerateMode::from_str("balanced"),
            Some(GenerateMode::Balanced)
        );
        assert_eq!(
            GenerateMode::from_str("skewed-a-gt-b"),
            Some(GenerateMode::SkewedAGtB)
        );
        assert_eq!(
            GenerateMode::from_str("skewed_a_gt_b"),
            Some(GenerateMode::SkewedAGtB)
        );
        assert_eq!(GenerateMode::from_str("all"), Some(GenerateMode::All));
        assert_eq!(GenerateMode::from_str("invalid"), None);
    }

    #[test]
    fn test_sort_mode_should_sort() {
        assert!(SortMode::Yes.should_sort(100));
        assert!(!SortMode::No.should_sort(100_000_000));
        assert!(!SortMode::Auto.should_sort(999_999));
        assert!(SortMode::Auto.should_sort(1_000_000));
    }

    #[test]
    fn test_human_genome() {
        let genome = HumanGenome::new();
        assert_eq!(genome.chromosomes.len(), 23);
        assert!(genome.total_size > 3_000_000_000); // ~3Gb
    }

    #[test]
    fn test_generate_uniform_interval() {
        let config = GenerateConfig::default();
        let cmd = GenerateCommand::new(config);
        let mut rng = SmallRng::seed_from_u64(42);

        let interval = cmd.generate_uniform_interval(&mut rng);
        assert!(interval.chrom_idx < 23);
        assert!(interval.start < interval.end);
        assert!(interval.end - interval.start >= 50);
        assert!(interval.end - interval.start <= 1000);
    }

    #[test]
    fn test_deterministic_generation() {
        let config1 = GenerateConfig {
            seed: 12345,
            ..Default::default()
        };
        let cmd1 = GenerateCommand::new(config1);
        let mut rng1 = SmallRng::seed_from_u64(12345);
        let intervals1 = cmd1.generate_intervals(100, false, &mut rng1);

        let config2 = GenerateConfig {
            seed: 12345,
            ..Default::default()
        };
        let cmd2 = GenerateCommand::new(config2);
        let mut rng2 = SmallRng::seed_from_u64(12345);
        let intervals2 = cmd2.generate_intervals(100, false, &mut rng2);

        // Same seed should produce identical intervals
        assert_eq!(intervals1.len(), intervals2.len());
        for (i1, i2) in intervals1.iter().zip(intervals2.iter()) {
            assert_eq!(i1.chrom_idx, i2.chrom_idx);
            assert_eq!(i1.start, i2.start);
            assert_eq!(i1.end, i2.end);
        }
    }

    #[test]
    fn test_chrom_to_index() {
        // Lexicographic order
        assert_eq!(chrom_to_index("chr1"), 0);
        assert_eq!(chrom_to_index("chr10"), 1);
        assert_eq!(chrom_to_index("chr2"), 11);
        assert_eq!(chrom_to_index("chr22"), 14);
        assert_eq!(chrom_to_index("chrX"), 22);
    }

    #[test]
    fn test_format_count() {
        assert_eq!(format_count(1_000), "1K");
        assert_eq!(format_count(5_000_000), "5M");
        assert_eq!(format_count(1_000_000_000), "1G");
        assert_eq!(format_count(123), "123");
        assert_eq!(format_count(1_500_000), "1500K"); // Divisible by 1000, uses K suffix
        assert_eq!(format_count(1_234_567), "1234567"); // Not divisible evenly
    }
}
