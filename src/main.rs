// Clippy allows
#![allow(clippy::too_many_arguments)]

//! GRIT: Genomic Range Interval Toolkit
//!
//! Usage: grit <COMMAND> [OPTIONS]

use clap::{Parser, Subcommand};
use std::io;
use std::path::PathBuf;
use std::process;

use grit_genomics::bed::{BedError, BedReader};
use grit_genomics::commands::{
    verify_sorted, verify_sorted_reader, verify_sorted_with_genome, ClosestCommand,
    ComplementCommand, FastMergeCommand, FastSortCommand, GenomecovCommand, GenomecovOutputMode,
    IntersectCommand, JaccardCommand, MergeCommand, MultiinterCommand, SlopCommand, SortCommand,
    StreamingClosestCommand, StreamingCoverageCommand, StreamingGenomecovCommand,
    StreamingGenomecovMode, StreamingIntersectCommand, StreamingMultiinterCommand,
    StreamingSubtractCommand, StreamingWindowCommand, SubtractCommand,
};
use grit_genomics::genome::Genome;

#[derive(Parser)]
#[command(name = "grit")]
#[command(author = "Manish Kumar Bobbili")]
#[command(version)]
#[command(about = "GRIT: Genomic Range Interval Toolkit - high-performance genomic interval operations", long_about = None)]
struct Cli {
    /// Number of threads to use (default: number of CPUs)
    #[arg(long, short = 't', global = true)]
    threads: Option<usize>,

    /// Normalize zero-length intervals (start == end) to 1bp intervals
    /// to match bedtools behavior. By default, GRIT uses strict half-open
    /// interval semantics where zero-length intervals do not overlap with
    /// themselves.
    #[arg(long, global = true)]
    bedtools_compatible: bool,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Sort a BED file by chromosome and position
    Sort {
        /// Input BED file (use - for stdin)
        #[arg(short, long)]
        input: Option<PathBuf>,

        /// Genome file for chromosome ordering
        #[arg(short = 'g', long)]
        genome: Option<PathBuf>,

        /// Sort by interval size (ascending)
        #[arg(long = "sizeA")]
        size_asc: bool,

        /// Sort by interval size (descending)
        #[arg(long = "sizeD")]
        size_desc: bool,

        /// Reverse the sort order
        #[arg(short, long)]
        reverse: bool,

        /// Sort by chromosome name only
        #[arg(long = "chrThenSizeA")]
        chrom_only: bool,

        /// Legacy flag (fast mode is now default, kept for compatibility)
        #[arg(long, hide = true)]
        fast: bool,

        /// Print sorting statistics to stderr
        #[arg(long)]
        stats: bool,
    },

    /// Merge overlapping intervals
    Merge {
        /// Input BED file (use - for stdin)
        #[arg(short, long)]
        input: Option<PathBuf>,

        /// Maximum distance between intervals to merge
        #[arg(short, long, default_value = "0")]
        distance: u64,

        /// Require strand to match for merging
        #[arg(short, long)]
        strand: bool,

        /// Use in-memory mode (loads all records, handles unsorted input)
        #[arg(long)]
        in_memory: bool,

        /// Report count of merged intervals
        #[arg(short = 'c', long)]
        count: bool,

        /// Print streaming statistics to stderr
        #[arg(long)]
        stats: bool,

        /// Skip sorted validation (faster for pre-sorted input)
        #[arg(long)]
        assume_sorted: bool,

        /// Genome file for chromosome order validation
        #[arg(short = 'g', long)]
        genome: Option<PathBuf>,
    },

    /// Find overlapping intervals between two BED files
    Intersect {
        /// Input BED file A
        #[arg(short = 'a', long)]
        file_a: PathBuf,

        /// Input BED file B
        #[arg(short = 'b', long)]
        file_b: PathBuf,

        /// Write original A entry (-wa in bedtools)
        #[arg(long = "wa")]
        write_a: bool,

        /// Write original B entry (-wb in bedtools)
        #[arg(long = "wb")]
        write_b: bool,

        /// Only report unique A intervals
        #[arg(short = 'u', long)]
        unique: bool,

        /// Only report A intervals with NO overlap
        #[arg(short = 'v', long)]
        no_overlap: bool,

        /// Minimum overlap fraction for A
        #[arg(short = 'f', long)]
        fraction: Option<f64>,

        /// Require reciprocal fraction overlap
        #[arg(short = 'r', long)]
        reciprocal: bool,

        /// Report the number of overlaps
        #[arg(short = 'c', long)]
        count: bool,

        /// Use streaming mode (constant memory, requires sorted input)
        #[arg(long)]
        streaming: bool,

        /// Print streaming statistics to stderr
        #[arg(long)]
        stats: bool,

        /// Skip sorted validation (faster for pre-sorted input)
        #[arg(long)]
        assume_sorted: bool,

        /// Allow unsorted input (loads into memory and re-sorts, uses O(n) memory)
        #[arg(long)]
        allow_unsorted: bool,

        /// Genome file for chromosome order validation (streaming mode)
        #[arg(short = 'g', long)]
        genome: Option<PathBuf>,
    },

    /// Remove intervals in A that overlap with B
    Subtract {
        /// Input BED file A
        #[arg(short = 'a', long)]
        file_a: PathBuf,

        /// Input BED file B
        #[arg(short = 'b', long)]
        file_b: PathBuf,

        /// Remove entire A feature if any overlap
        #[arg(short = 'A', long)]
        remove_entire: bool,

        /// Minimum overlap fraction required
        #[arg(short = 'f', long)]
        fraction: Option<f64>,

        /// Require reciprocal fraction overlap
        #[arg(short = 'r', long)]
        reciprocal: bool,

        /// Use streaming mode (O(k) memory, requires sorted input)
        #[arg(long)]
        streaming: bool,

        /// Print streaming statistics to stderr
        #[arg(long)]
        stats: bool,

        /// Skip sorted validation (faster for pre-sorted input)
        #[arg(long)]
        assume_sorted: bool,

        /// Allow unsorted input (loads into memory and re-sorts, uses O(n) memory)
        #[arg(long)]
        allow_unsorted: bool,

        /// Genome file for chromosome order validation (streaming mode)
        #[arg(short = 'g', long)]
        genome: Option<PathBuf>,
    },

    /// Find the closest interval in B for each interval in A
    Closest {
        /// Input BED file A
        #[arg(short = 'a', long)]
        file_a: PathBuf,

        /// Input BED file B
        #[arg(short = 'b', long)]
        file_b: PathBuf,

        /// Report distance in output
        #[arg(short = 'd', long)]
        distance: bool,

        /// Report all ties
        #[arg(short = 't', long, value_parser = ["all", "first", "last"])]
        tie: Option<String>,

        /// Ignore overlapping intervals
        #[arg(long = "io")]
        ignore_overlaps: bool,

        /// Ignore upstream intervals
        #[arg(long = "iu")]
        ignore_upstream: bool,

        /// Ignore downstream intervals
        #[arg(long = "id")]
        ignore_downstream: bool,

        /// Maximum distance to report
        #[arg(short = 'D', long)]
        max_distance: Option<u64>,

        /// Use streaming mode (O(k) memory, requires sorted input)
        #[arg(long)]
        streaming: bool,

        /// Skip sorted validation (faster for pre-sorted input)
        #[arg(long)]
        assume_sorted: bool,

        /// Allow unsorted input (loads into memory and re-sorts, uses O(n) memory)
        #[arg(long)]
        allow_unsorted: bool,

        /// Genome file for chromosome order validation (streaming mode)
        #[arg(short = 'g', long)]
        genome: Option<PathBuf>,
    },

    /// Find intervals in B that are within a window of A
    Window {
        /// Input BED file A
        #[arg(short = 'a', long)]
        file_a: PathBuf,

        /// Input BED file B
        #[arg(short = 'b', long)]
        file_b: PathBuf,

        /// Window size (both sides)
        #[arg(short = 'w', long, default_value = "1000")]
        window: u64,

        /// Left window size
        #[arg(short = 'l', long)]
        left: Option<u64>,

        /// Right window size
        #[arg(short = 'r', long)]
        right: Option<u64>,

        /// Report number of overlaps
        #[arg(short = 'c', long)]
        count: bool,

        /// Only report A intervals with no matches
        #[arg(short = 'v', long)]
        no_overlap: bool,

        /// Skip sorted validation (faster for pre-sorted input)
        #[arg(long)]
        assume_sorted: bool,

        /// Genome file for chromosome order validation
        #[arg(short = 'g', long)]
        genome: Option<PathBuf>,
    },

    /// Calculate coverage of A intervals by B intervals
    Coverage {
        /// Input BED file A (regions)
        #[arg(short = 'a', long)]
        file_a: PathBuf,

        /// Input BED file B (reads/features)
        #[arg(short = 'b', long)]
        file_b: PathBuf,

        /// Report a histogram of coverage
        #[arg(long = "hist")]
        histogram: bool,

        /// Report depth at each position
        #[arg(short = 'd', long)]
        per_base: bool,

        /// Report mean depth
        #[arg(long)]
        mean: bool,

        /// Skip sorted validation (faster for pre-sorted input)
        #[arg(long)]
        assume_sorted: bool,

        /// Genome file for chromosome order validation
        #[arg(short = 'g', long)]
        genome: Option<PathBuf>,
    },

    /// Extend intervals by a given number of bases
    Slop {
        /// Input BED file
        #[arg(short, long)]
        input: PathBuf,

        /// Genome file (chrom sizes)
        #[arg(short, long)]
        genome: PathBuf,

        /// Extend both sides by this many bases (or fraction if -pct)
        #[arg(short = 'b', long)]
        both: Option<f64>,

        /// Extend left/upstream by this many bases (or fraction if -pct)
        #[arg(short = 'l', long)]
        left: Option<f64>,

        /// Extend right/downstream by this many bases (or fraction if -pct)
        #[arg(short = 'r', long)]
        right: Option<f64>,

        /// Use strand info (left=upstream, right=downstream)
        #[arg(short = 's', long)]
        strand: bool,

        /// Interpret values as fraction of interval size
        #[arg(long)]
        pct: bool,
    },

    /// Return intervals NOT covered by the input BED file
    Complement {
        /// Input BED file
        #[arg(short, long)]
        input: PathBuf,

        /// Genome file (chrom sizes)
        #[arg(short, long)]
        genome: PathBuf,

        /// Assume input is sorted in genome order (enables O(1) memory streaming)
        #[arg(long)]
        assume_sorted: bool,
    },

    /// Compute genome-wide coverage
    Genomecov {
        /// Input BED file
        #[arg(short, long)]
        input: PathBuf,

        /// Genome file (chrom sizes)
        #[arg(short, long)]
        genome: PathBuf,

        /// Report depth at each position (1-based)
        #[arg(short = 'd', long)]
        per_base: bool,

        /// Report BedGraph format (non-zero only)
        #[arg(long = "bg")]
        bedgraph: bool,

        /// Report BedGraph format (including zero coverage)
        #[arg(long = "bga")]
        bedgraph_all: bool,

        /// Scale depth by factor
        #[arg(long, default_value = "1.0")]
        scale: f64,

        /// Use streaming mode (O(k) memory, requires sorted input)
        #[arg(long)]
        streaming: bool,

        /// Skip sorted validation (faster for pre-sorted input)
        #[arg(long)]
        assume_sorted: bool,
    },

    /// Calculate Jaccard similarity between two BED files
    Jaccard {
        /// Input BED file A
        #[arg(short = 'a', long)]
        file_a: PathBuf,

        /// Input BED file B
        #[arg(short = 'b', long)]
        file_b: PathBuf,
    },

    /// Identify common intervals across multiple BED files
    Multiinter {
        /// Input BED files
        #[arg(short = 'i', long = "input", num_args = 1..)]
        inputs: Vec<PathBuf>,

        /// Only output intervals found in all files
        #[arg(long)]
        cluster: bool,

        /// Use streaming mode (O(k) memory, requires sorted input)
        #[arg(long)]
        streaming: bool,

        /// Skip sorted validation (faster for pre-sorted input)
        #[arg(long)]
        assume_sorted: bool,
    },

    /// Generate synthetic BED datasets for benchmarking
    #[command(alias = "create")]
    Generate {
        /// Output directory
        #[arg(short, long, default_value = "./grit_bench_data")]
        output: PathBuf,

        /// Sizes to generate (comma-separated, e.g., "1M,5M,10M")
        #[arg(long, default_value = "1M,5M,10M,25M,50M")]
        sizes: String,

        /// Random seed for reproducibility
        #[arg(long, default_value = "42")]
        seed: u64,

        /// Generation mode: balanced|skewed-a-gt-b|skewed-b-gt-a|identical|clustered|all
        #[arg(long, default_value = "all")]
        mode: String,

        /// Sorting behavior: yes|no|auto
        #[arg(long, default_value = "auto")]
        sorted: String,

        /// Alias for --sorted no
        #[arg(long)]
        no_sort: bool,

        /// Custom A file size (e.g., "5M")
        #[arg(long)]
        a: Option<String>,

        /// Custom B file size (e.g., "10M")
        #[arg(long)]
        b: Option<String>,

        /// Genome fraction for hotspots (clustered mode)
        #[arg(long, default_value = "0.05")]
        hotspot_frac: f64,

        /// Interval fraction in hotspots (clustered mode)
        #[arg(long, default_value = "0.80")]
        hotspot_weight: f64,

        /// Minimum interval length
        #[arg(long, default_value = "50")]
        len_min: u32,

        /// Maximum interval length
        #[arg(long, default_value = "1000")]
        len_max: u32,

        /// Overwrite existing files
        #[arg(long)]
        force: bool,
    },
}

/// Preprocess CLI arguments to support bedtools-style flags.
/// Converts -wa to --wa and -wb to --wb for compatibility.
fn preprocess_args() -> Vec<String> {
    std::env::args()
        .map(|arg| match arg.as_str() {
            "-wa" => "--wa".to_string(),
            "-wb" => "--wb".to_string(),
            _ => arg,
        })
        .collect()
}

fn main() {
    let cli = Cli::parse_from(preprocess_args());

    // Configure bedtools-compatible mode if requested
    // This must be set before any parsing occurs
    if cli.bedtools_compatible {
        grit_genomics::config::set_bedtools_compatible(true);
    }

    // Configure thread pool if --threads specified
    if let Some(n) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to initialize thread pool");
    }

    let result = match cli.command {
        Commands::Sort {
            input,
            genome,
            size_asc,
            size_desc,
            reverse,
            chrom_only,
            fast,
            stats,
        } => run_sort(
            input, genome, size_asc, size_desc, reverse, chrom_only, fast, stats,
        ),

        Commands::Merge {
            input,
            distance,
            strand,
            in_memory,
            count,
            stats,
            assume_sorted,
            genome,
        } => run_merge(
            input,
            distance,
            strand,
            in_memory,
            count,
            stats,
            assume_sorted,
            genome,
        ),

        Commands::Intersect {
            file_a,
            file_b,
            write_a,
            write_b,
            unique,
            no_overlap,
            fraction,
            reciprocal,
            count,
            streaming,
            stats,
            assume_sorted,
            allow_unsorted,
            genome,
        } => run_intersect(
            file_a,
            file_b,
            write_a,
            write_b,
            unique,
            no_overlap,
            fraction,
            reciprocal,
            count,
            streaming,
            stats,
            assume_sorted,
            allow_unsorted,
            genome,
        ),

        Commands::Subtract {
            file_a,
            file_b,
            remove_entire,
            fraction,
            reciprocal,
            streaming,
            stats,
            assume_sorted,
            allow_unsorted,
            genome,
        } => run_subtract(
            file_a,
            file_b,
            remove_entire,
            fraction,
            reciprocal,
            streaming,
            stats,
            assume_sorted,
            allow_unsorted,
            genome,
        ),

        Commands::Closest {
            file_a,
            file_b,
            distance,
            tie,
            ignore_overlaps,
            ignore_upstream,
            ignore_downstream,
            max_distance,
            streaming,
            assume_sorted,
            allow_unsorted,
            genome,
        } => run_closest(
            file_a,
            file_b,
            distance,
            tie,
            ignore_overlaps,
            ignore_upstream,
            ignore_downstream,
            max_distance,
            streaming,
            assume_sorted,
            allow_unsorted,
            genome,
        ),

        Commands::Window {
            file_a,
            file_b,
            window,
            left,
            right,
            count,
            no_overlap,
            assume_sorted,
            genome,
        } => run_window(
            file_a,
            file_b,
            window,
            left,
            right,
            count,
            no_overlap,
            assume_sorted,
            genome,
        ),

        Commands::Coverage {
            file_a,
            file_b,
            histogram,
            per_base,
            mean,
            assume_sorted,
            genome,
        } => run_coverage(
            file_a,
            file_b,
            histogram,
            per_base,
            mean,
            assume_sorted,
            genome,
        ),

        Commands::Slop {
            input,
            genome,
            both,
            left,
            right,
            strand,
            pct,
        } => run_slop(input, genome, both, left, right, strand, pct),

        Commands::Complement {
            input,
            genome,
            assume_sorted,
        } => run_complement(input, genome, assume_sorted),

        Commands::Genomecov {
            input,
            genome,
            per_base,
            bedgraph,
            bedgraph_all,
            scale,
            streaming,
            assume_sorted,
        } => run_genomecov(
            input,
            genome,
            per_base,
            bedgraph,
            bedgraph_all,
            scale,
            streaming,
            assume_sorted,
        ),

        Commands::Jaccard { file_a, file_b } => run_jaccard(file_a, file_b),

        Commands::Multiinter {
            inputs,
            cluster,
            streaming,
            assume_sorted,
        } => run_multiinter(inputs, cluster, streaming, assume_sorted),

        Commands::Generate {
            output,
            sizes,
            seed,
            mode,
            sorted,
            no_sort,
            a,
            b,
            hotspot_frac,
            hotspot_weight,
            len_min,
            len_max,
            force,
        } => run_generate(
            output,
            sizes,
            seed,
            mode,
            sorted,
            no_sort,
            a,
            b,
            hotspot_frac,
            hotspot_weight,
            len_min,
            len_max,
            force,
        ),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}

fn run_sort(
    input: Option<PathBuf>,
    genome: Option<PathBuf>,
    size_asc: bool,
    size_desc: bool,
    reverse: bool,
    chrom_only: bool,
    _fast: bool, // Legacy flag, fast mode is now default
    stats: bool,
) -> Result<(), BedError> {
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    // Load genome file if provided
    let genome = genome.map(|p| Genome::from_file(&p)).transpose()?;

    // Use fast mode by default when no special sort modes requested
    // Fast mode uses radix sort + mmap for better performance
    // Fall back to standard sort only for --sizeA, --sizeD, --chrThenSizeA
    let use_fast = !size_asc && !size_desc && !chrom_only;

    if use_fast {
        let mut cmd = FastSortCommand::new();
        cmd.reverse = reverse;

        // Apply genome ordering if provided
        if let Some(ref g) = genome {
            cmd = cmd.with_genome(g);
        }

        let result = if let Some(path) = input {
            if path.to_string_lossy() == "-" {
                cmd.run_stdin(&mut handle)?
            } else {
                cmd.run(&path, &mut handle)?
            }
        } else {
            cmd.run_stdin(&mut handle)?
        };

        if stats {
            eprintln!("Fast sort stats: {}", result);
        }

        Ok(())
    } else {
        // Use standard sort for special sort modes
        let mut cmd = SortCommand::new();
        cmd.size_asc = size_asc;
        cmd.size_desc = size_desc;
        cmd.reverse = reverse;
        cmd.chrom_only = chrom_only;

        // Apply genome ordering if provided
        if let Some(ref g) = genome {
            cmd = cmd.with_genome(g);
        }

        if let Some(path) = input {
            if path.to_string_lossy() == "-" {
                cmd.run_stdio()
            } else {
                cmd.run(path, &mut handle)
            }
        } else {
            cmd.run_stdio()
        }
    }
}

/// Helper to validate sort order, optionally using genome file for chromosome ordering.
fn validate_sorted(path: &PathBuf, genome: Option<&Genome>) -> Result<(), BedError> {
    if let Some(g) = genome {
        verify_sorted_with_genome(path, g)
    } else {
        verify_sorted(path)
    }
}

fn run_merge(
    input: Option<PathBuf>,
    distance: u64,
    strand: bool,
    in_memory: bool,
    count: bool,
    stats: bool,
    assume_sorted: bool,
    genome_path: Option<PathBuf>,
) -> Result<(), BedError> {
    // Load genome file if provided
    let genome =
        if let Some(ref gp) = genome_path {
            Some(Genome::from_file(gp).map_err(|e| {
                BedError::InvalidFormat(format!("Failed to load genome file: {}", e))
            })?)
        } else {
            None
        };
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    if in_memory {
        // Use in-memory mode - loads all records, can handle unsorted input
        let cmd = MergeCommand::new()
            .with_distance(distance)
            .with_strand(strand);

        if let Some(path) = input {
            if path.to_string_lossy() == "-" {
                let stdin = io::stdin();
                let reader = BedReader::new(stdin.lock());
                cmd.merge_streaming(reader, &mut handle)
            } else {
                cmd.run(path, &mut handle)
            }
        } else {
            let stdin = io::stdin();
            let reader = BedReader::new(stdin.lock());
            cmd.merge_streaming(reader, &mut handle)
        }
    } else if strand {
        // Strand-specific merge not yet implemented in fast path, use standard streaming
        use grit_genomics::commands::StreamingMergeCommand;
        let mut cmd = StreamingMergeCommand::new()
            .with_distance(distance)
            .with_strand(strand);
        cmd.count = count;

        let result = if let Some(path) = input {
            if path.to_string_lossy() == "-" {
                // Stdin: validate by buffering, then process
                if !assume_sorted {
                    let stdin = io::stdin();
                    let buffer = verify_sorted_reader(stdin.lock()).map_err(|e| {
                        BedError::InvalidFormat(format!(
                            "stdin is not sorted: {}\n\n\
                             Fix: Pre-sort your input before piping.\n\
                             Or use '--assume-sorted' if you know the input is sorted.",
                            e
                        ))
                    })?;
                    let cursor = std::io::Cursor::new(buffer);
                    let reader = BedReader::new(cursor);
                    cmd.run_streaming(reader, &mut handle)?
                } else {
                    cmd.run_stdin(&mut handle)?
                }
            } else {
                // File: validate before processing
                if !assume_sorted {
                    validate_sorted(&path, genome.as_ref()).map_err(|e| {
                        BedError::InvalidFormat(format!(
                            "Input is not sorted: {}\n\n\
                             Fix: Run 'grit sort -i {}{}' first.\n\
                             Or use '--assume-sorted' if you know the input is sorted.",
                            e,
                            path.display(),
                            if genome.is_some() {
                                " -g <genome.txt>"
                            } else {
                                ""
                            }
                        ))
                    })?;
                }
                cmd.run(&path, &mut handle)?
            }
        } else {
            // No path specified: read from stdin
            if !assume_sorted {
                let stdin = io::stdin();
                let buffer = verify_sorted_reader(stdin.lock()).map_err(|e| {
                    BedError::InvalidFormat(format!(
                        "stdin is not sorted: {}\n\n\
                         Fix: Pre-sort your input before piping.\n\
                         Or use '--assume-sorted' if you know the input is sorted.",
                        e
                    ))
                })?;
                let cursor = std::io::Cursor::new(buffer);
                let reader = BedReader::new(cursor);
                cmd.run_streaming(reader, &mut handle)?
            } else {
                cmd.run_stdin(&mut handle)?
            }
        };

        if stats {
            eprintln!("Streaming merge stats: {}", result);
        }

        Ok(())
    } else {
        // Use fast streaming mode (default) - O(1) memory, zero-allocation parsing
        let mut cmd = FastMergeCommand::new().with_distance(distance);
        cmd.count = count;

        let result = if let Some(path) = input {
            if path.to_string_lossy() == "-" {
                // Stdin: validate by buffering, then process
                if !assume_sorted {
                    let stdin = io::stdin();
                    let buffer = verify_sorted_reader(stdin.lock()).map_err(|e| {
                        BedError::InvalidFormat(format!(
                            "stdin is not sorted: {}\n\n\
                             Fix: Pre-sort your input before piping.\n\
                             Or use '--assume-sorted' if you know the input is sorted.",
                            e
                        ))
                    })?;
                    let cursor = std::io::Cursor::new(buffer);
                    cmd.run_reader(cursor, &mut handle)?
                } else {
                    cmd.run_stdin(&mut handle)?
                }
            } else {
                // File: validate before processing
                if !assume_sorted {
                    validate_sorted(&path, genome.as_ref()).map_err(|e| {
                        BedError::InvalidFormat(format!(
                            "Input is not sorted: {}\n\n\
                             Fix: Run 'grit sort -i {}{}' first.\n\
                             Or use '--assume-sorted' if you know the input is sorted.",
                            e,
                            path.display(),
                            if genome.is_some() {
                                " -g <genome.txt>"
                            } else {
                                ""
                            }
                        ))
                    })?;
                }
                cmd.run(&path, &mut handle)?
            }
        } else {
            // No path specified: read from stdin
            if !assume_sorted {
                let stdin = io::stdin();
                let buffer = verify_sorted_reader(stdin.lock()).map_err(|e| {
                    BedError::InvalidFormat(format!(
                        "stdin is not sorted: {}\n\n\
                         Fix: Pre-sort your input before piping.\n\
                         Or use '--assume-sorted' if you know the input is sorted.",
                        e
                    ))
                })?;
                let cursor = std::io::Cursor::new(buffer);
                cmd.run_reader(cursor, &mut handle)?
            } else {
                cmd.run_stdin(&mut handle)?
            }
        };

        if stats {
            eprintln!("Fast merge stats: {}", result);
        }

        Ok(())
    }
}

fn run_intersect(
    file_a: PathBuf,
    file_b: PathBuf,
    write_a: bool,
    write_b: bool,
    unique: bool,
    no_overlap: bool,
    fraction: Option<f64>,
    reciprocal: bool,
    count: bool,
    streaming: bool,
    stats: bool,
    assume_sorted: bool,
    allow_unsorted: bool,
    genome_path: Option<PathBuf>,
) -> Result<(), BedError> {
    // Load genome file if provided
    let genome =
        if let Some(ref gp) = genome_path {
            Some(Genome::from_file(gp).map_err(|e| {
                BedError::InvalidFormat(format!("Failed to load genome file: {}", e))
            })?)
        } else {
            None
        };

    let stdout = io::stdout();
    let mut handle = stdout.lock();
    let genome_flag = if genome.is_some() {
        " -g <genome.txt>"
    } else {
        ""
    };

    if streaming {
        // Use streaming mode - constant memory, requires sorted input
        // Only validate sorted order if --assume-sorted is not set
        if !assume_sorted {
            validate_sorted(&file_a, genome.as_ref()).map_err(|e| {
                BedError::InvalidFormat(format!(
                    "File A is not sorted: {}\n\n\
                     Fix: Run 'grit sort -i {}{}' > sorted_a.bed first.\n\
                     Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).",
                    e,
                    file_a.display(),
                    genome_flag
                ))
            })?;
            validate_sorted(&file_b, genome.as_ref()).map_err(|e| {
                BedError::InvalidFormat(format!(
                    "File B is not sorted: {}\n\n\
                     Fix: Run 'grit sort -i {}{}' > sorted_b.bed first.\n\
                     Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).",
                    e,
                    file_b.display(),
                    genome_flag
                ))
            })?;
        }

        let mut cmd = StreamingIntersectCommand::new();
        cmd.write_a = write_a;
        cmd.write_b = write_b;
        cmd.unique = unique;
        cmd.no_overlap = no_overlap;
        cmd.fraction_a = fraction;
        cmd.reciprocal = reciprocal;
        cmd.count = count;
        // Always skip inline validation in streaming mode - we either validated above or user assumes sorted
        cmd.assume_sorted = true;

        let result = cmd.run(&file_a, &file_b, &mut handle)?;

        if stats {
            eprintln!("Streaming intersect stats: {}", result);
        }

        Ok(())
    } else {
        // Non-streaming mode: validate sorted input unless --allow-unsorted
        if !allow_unsorted {
            validate_sorted(&file_a, genome.as_ref()).map_err(|e| {
                BedError::InvalidFormat(format!(
                    "File A is not sorted: {}\n\n\
                     Fix: Run 'grit sort -i {}{}' > sorted_a.bed first.\n\
                     Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).\n\
                     Or use '--streaming' for O(k) memory with pre-sorted input.",
                    e,
                    file_a.display(),
                    genome_flag
                ))
            })?;
            validate_sorted(&file_b, genome.as_ref()).map_err(|e| {
                BedError::InvalidFormat(format!(
                    "File B is not sorted: {}\n\n\
                     Fix: Run 'grit sort -i {}{}' > sorted_b.bed first.\n\
                     Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).\n\
                     Or use '--streaming' for O(k) memory with pre-sorted input.",
                    e,
                    file_b.display(),
                    genome_flag
                ))
            })?;
        }

        // Use standard parallel mode
        let mut cmd = IntersectCommand::new();
        cmd.write_a = write_a;
        cmd.write_b = write_b;
        cmd.unique = unique;
        cmd.no_overlap = no_overlap;
        cmd.fraction_a = fraction;
        cmd.reciprocal = reciprocal;
        cmd.count = count;

        cmd.run(file_a, file_b, &mut handle)
    }
}

fn run_subtract(
    file_a: PathBuf,
    file_b: PathBuf,
    remove_entire: bool,
    fraction: Option<f64>,
    reciprocal: bool,
    streaming: bool,
    stats: bool,
    assume_sorted: bool,
    allow_unsorted: bool,
    genome_path: Option<PathBuf>,
) -> Result<(), BedError> {
    // Load genome file if provided
    let genome =
        if let Some(ref gp) = genome_path {
            Some(Genome::from_file(gp).map_err(|e| {
                BedError::InvalidFormat(format!("Failed to load genome file: {}", e))
            })?)
        } else {
            None
        };

    let stdout = io::stdout();
    let mut handle = stdout.lock();
    let genome_flag = if genome.is_some() {
        " -g <genome.txt>"
    } else {
        ""
    };

    if streaming {
        // Use streaming mode - O(k) memory, requires sorted input
        // Validate that both input files are sorted (unless --assume-sorted)
        if !assume_sorted {
            validate_sorted(&file_a, genome.as_ref()).map_err(|e| {
                BedError::InvalidFormat(format!(
                    "File A is not sorted: {}\n\n\
                     Fix: Run 'grit sort -i {}{}' > sorted_a.bed first.\n\
                     Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).",
                    e,
                    file_a.display(),
                    genome_flag
                ))
            })?;
            validate_sorted(&file_b, genome.as_ref()).map_err(|e| {
                BedError::InvalidFormat(format!(
                    "File B is not sorted: {}\n\n\
                     Fix: Run 'grit sort -i {}{}' > sorted_b.bed first.\n\
                     Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).",
                    e,
                    file_b.display(),
                    genome_flag
                ))
            })?;
        }

        let mut cmd = StreamingSubtractCommand::new();
        cmd.remove_entire = remove_entire;
        cmd.fraction = fraction;
        cmd.reciprocal = reciprocal;

        let result = cmd.run(&file_a, &file_b, &mut handle)?;

        if stats {
            eprintln!("Streaming subtract stats: {}", result);
        }

        Ok(())
    } else {
        // Non-streaming mode: validate sorted input unless --allow-unsorted
        if !allow_unsorted {
            validate_sorted(&file_a, genome.as_ref()).map_err(|e| {
                BedError::InvalidFormat(format!(
                    "File A is not sorted: {}\n\n\
                     Fix: Run 'grit sort -i {}{}' > sorted_a.bed first.\n\
                     Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).\n\
                     Or use '--streaming' for O(k) memory with pre-sorted input.",
                    e,
                    file_a.display(),
                    genome_flag
                ))
            })?;
            validate_sorted(&file_b, genome.as_ref()).map_err(|e| {
                BedError::InvalidFormat(format!(
                    "File B is not sorted: {}\n\n\
                     Fix: Run 'grit sort -i {}{}' > sorted_b.bed first.\n\
                     Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).\n\
                     Or use '--streaming' for O(k) memory with pre-sorted input.",
                    e,
                    file_b.display(),
                    genome_flag
                ))
            })?;
        }

        // Use standard mode
        let mut cmd = SubtractCommand::new();
        cmd.remove_entire = remove_entire;
        cmd.fraction = fraction;
        cmd.reciprocal = reciprocal;

        cmd.run(file_a, file_b, &mut handle)
    }
}

fn run_closest(
    file_a: PathBuf,
    file_b: PathBuf,
    _distance: bool,
    tie: Option<String>,
    ignore_overlaps: bool,
    ignore_upstream: bool,
    ignore_downstream: bool,
    _max_distance: Option<u64>,
    streaming: bool,
    assume_sorted: bool,
    allow_unsorted: bool,
    genome_path: Option<PathBuf>,
) -> Result<(), BedError> {
    // Load genome file if provided
    let genome =
        if let Some(ref gp) = genome_path {
            Some(Genome::from_file(gp).map_err(|e| {
                BedError::InvalidFormat(format!("Failed to load genome file: {}", e))
            })?)
        } else {
            None
        };

    let stdout = io::stdout();
    let mut handle = stdout.lock();
    let genome_flag = if genome.is_some() {
        " -g <genome.txt>"
    } else {
        ""
    };

    if streaming {
        // Validate that both input files are sorted (unless --assume-sorted)
        if !assume_sorted {
            validate_sorted(&file_a, genome.as_ref()).map_err(|e| {
                BedError::InvalidFormat(format!(
                    "File A is not sorted: {}\n\n\
                     Fix: Run 'grit sort -i {}{}' > sorted_a.bed first.\n\
                     Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).",
                    e,
                    file_a.display(),
                    genome_flag
                ))
            })?;
            validate_sorted(&file_b, genome.as_ref()).map_err(|e| {
                BedError::InvalidFormat(format!(
                    "File B is not sorted: {}\n\n\
                     Fix: Run 'grit sort -i {}{}' > sorted_b.bed first.\n\
                     Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).",
                    e,
                    file_b.display(),
                    genome_flag
                ))
            })?;
        }

        // Use streaming implementation (O(k) memory)
        let mut cmd = StreamingClosestCommand::new();
        cmd.ignore_overlaps = ignore_overlaps;
        cmd.ignore_upstream = ignore_upstream;
        cmd.ignore_downstream = ignore_downstream;
        cmd.report_all_ties = tie.as_ref().is_none_or(|t| t == "all");

        cmd.run(file_a, file_b, &mut handle)?;
        Ok(())
    } else {
        // Non-streaming mode: validate sorted input unless --allow-unsorted
        if !allow_unsorted {
            validate_sorted(&file_a, genome.as_ref()).map_err(|e| {
                BedError::InvalidFormat(format!(
                    "File A is not sorted: {}\n\n\
                     Fix: Run 'grit sort -i {}{}' > sorted_a.bed first.\n\
                     Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).\n\
                     Or use '--streaming' for O(k) memory with pre-sorted input.",
                    e,
                    file_a.display(),
                    genome_flag
                ))
            })?;
            validate_sorted(&file_b, genome.as_ref()).map_err(|e| {
                BedError::InvalidFormat(format!(
                    "File B is not sorted: {}\n\n\
                     Fix: Run 'grit sort -i {}{}' > sorted_b.bed first.\n\
                     Or use '--allow-unsorted' to load and re-sort in memory (uses O(n) memory).\n\
                     Or use '--streaming' for O(k) memory with pre-sorted input.",
                    e,
                    file_b.display(),
                    genome_flag
                ))
            })?;
        }

        use grit_genomics::commands::closest::TieHandling;

        let mut cmd = ClosestCommand::new();
        cmd.report_distance = _distance;
        cmd.ignore_overlaps = ignore_overlaps;
        cmd.ignore_upstream = ignore_upstream;
        cmd.ignore_downstream = ignore_downstream;
        cmd.max_distance = _max_distance;

        if let Some(t) = tie {
            cmd.tie_handling = match t.as_str() {
                "all" => TieHandling::All,
                "first" => TieHandling::First,
                "last" => TieHandling::Last,
                _ => TieHandling::All,
            };
        }

        cmd.run(file_a, file_b, &mut handle)
    }
}

fn run_window(
    file_a: PathBuf,
    file_b: PathBuf,
    window: u64,
    left: Option<u64>,
    right: Option<u64>,
    count: bool,
    no_overlap: bool,
    assume_sorted: bool,
    genome_path: Option<PathBuf>,
) -> Result<(), BedError> {
    // Load genome file if provided
    let genome =
        if let Some(ref gp) = genome_path {
            Some(Genome::from_file(gp).map_err(|e| {
                BedError::InvalidFormat(format!("Failed to load genome file: {}", e))
            })?)
        } else {
            None
        };
    let genome_flag = if genome.is_some() {
        " -g <genome.txt>"
    } else {
        ""
    };

    // Use streaming implementation for better performance
    // Validate that both input files are sorted (unless --assume-sorted)
    if !assume_sorted {
        validate_sorted(&file_a, genome.as_ref()).map_err(|e| {
            BedError::InvalidFormat(format!(
                "File A is not sorted: {}\n\nFix: Run 'grit sort -i {}{}' first.",
                e,
                file_a.display(),
                genome_flag
            ))
        })?;
        validate_sorted(&file_b, genome.as_ref()).map_err(|e| {
            BedError::InvalidFormat(format!(
                "File B is not sorted: {}\n\nFix: Run 'grit sort -i {}{}' first.",
                e,
                file_b.display(),
                genome_flag
            ))
        })?;
    }

    let mut cmd = StreamingWindowCommand::new();
    cmd.window = window;
    cmd.left = left;
    cmd.right = right;
    cmd.count = count;
    cmd.no_overlap = no_overlap;

    let stdout = io::stdout();
    let mut handle = stdout.lock();

    cmd.run(file_a, file_b, &mut handle)?;
    Ok(())
}

fn run_coverage(
    file_a: PathBuf,
    file_b: PathBuf,
    histogram: bool,
    per_base: bool,
    mean: bool,
    assume_sorted: bool,
    genome_path: Option<PathBuf>,
) -> Result<(), BedError> {
    // Load genome file if provided
    let genome =
        if let Some(ref gp) = genome_path {
            Some(Genome::from_file(gp).map_err(|e| {
                BedError::InvalidFormat(format!("Failed to load genome file: {}", e))
            })?)
        } else {
            None
        };
    let genome_flag = if genome.is_some() {
        " -g <genome.txt>"
    } else {
        ""
    };

    // Validate that both input files are sorted (unless --assume-sorted)
    if !assume_sorted {
        validate_sorted(&file_a, genome.as_ref()).map_err(|e| {
            BedError::InvalidFormat(format!(
                "File A is not sorted: {}\n\nFix: Run 'grit sort -i {}{}' first.",
                e,
                file_a.display(),
                genome_flag
            ))
        })?;
        validate_sorted(&file_b, genome.as_ref()).map_err(|e| {
            BedError::InvalidFormat(format!(
                "File B is not sorted: {}\n\nFix: Run 'grit sort -i {}{}' first.",
                e,
                file_b.display(),
                genome_flag
            ))
        })?;
    }

    // Use streaming mode by default for memory efficiency
    // Memory: O(B × 8 bytes) instead of O((A + B) × ~400 bytes)
    let mut cmd = StreamingCoverageCommand::new();
    cmd.histogram = histogram;
    cmd.per_base = per_base;
    cmd.mean = mean;

    let stdout = io::stdout();
    let mut handle = stdout.lock();

    cmd.run(file_a, file_b, &mut handle)
}

fn run_slop(
    input: PathBuf,
    genome_file: PathBuf,
    both: Option<f64>,
    left: Option<f64>,
    right: Option<f64>,
    strand: bool,
    pct: bool,
) -> Result<(), BedError> {
    let genome = Genome::from_file(&genome_file)?;

    let mut cmd = SlopCommand::new();
    cmd.both = both.unwrap_or(0.0);
    cmd.left = left;
    cmd.right = right;
    cmd.strand = strand;
    cmd.pct = pct;

    let stdout = io::stdout();
    let mut handle = stdout.lock();

    cmd.run(input, &genome, &mut handle)
}

fn run_complement(
    input: PathBuf,
    genome_file: PathBuf,
    assume_sorted: bool,
) -> Result<(), BedError> {
    let genome = Genome::from_file(&genome_file)?;
    let cmd = ComplementCommand::new().with_assume_sorted(assume_sorted);

    let stdout = io::stdout();
    let mut handle = stdout.lock();

    cmd.run(input, &genome, &mut handle)
}

fn run_genomecov(
    input: PathBuf,
    genome_file: PathBuf,
    per_base: bool,
    bedgraph: bool,
    bedgraph_all: bool,
    scale: f64,
    streaming: bool,
    assume_sorted: bool,
) -> Result<(), BedError> {
    let genome = Genome::from_file(&genome_file)?;

    let stdout = io::stdout();
    let mut handle = stdout.lock();

    if streaming || assume_sorted {
        // Use streaming implementation with O(k) memory
        let mode = if per_base {
            StreamingGenomecovMode::PerBase
        } else if bedgraph_all {
            StreamingGenomecovMode::BedGraphAll
        } else if bedgraph {
            StreamingGenomecovMode::BedGraph
        } else {
            StreamingGenomecovMode::Histogram
        };

        let cmd = StreamingGenomecovCommand::new()
            .with_mode(mode)
            .with_scale(scale)
            .with_assume_sorted(assume_sorted);

        cmd.run(input, &genome, &mut handle)
    } else {
        // Use original implementation (loads all intervals into memory)
        let mut cmd = GenomecovCommand::new();
        cmd.scale = scale;

        if per_base {
            cmd.mode = GenomecovOutputMode::PerBase;
        } else if bedgraph_all {
            cmd.mode = GenomecovOutputMode::BedGraphAll;
        } else if bedgraph {
            cmd.mode = GenomecovOutputMode::BedGraph;
        }
        // else default to Histogram

        cmd.run(input, &genome, &mut handle)
    }
}

fn run_jaccard(file_a: PathBuf, file_b: PathBuf) -> Result<(), BedError> {
    let cmd = JaccardCommand::new();

    let stdout = io::stdout();
    let mut handle = stdout.lock();

    cmd.run(file_a, file_b, &mut handle)
}

fn run_multiinter(
    inputs: Vec<PathBuf>,
    cluster: bool,
    streaming: bool,
    assume_sorted: bool,
) -> Result<(), BedError> {
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    if streaming || assume_sorted {
        // Use streaming implementation with O(k) memory and k-way merge
        let cmd = StreamingMultiinterCommand::new()
            .with_cluster(cluster)
            .with_assume_sorted(assume_sorted);

        cmd.run(&inputs, &mut handle)
    } else {
        // Use original implementation (loads all intervals into memory)
        let mut cmd = MultiinterCommand::new();
        cmd.cluster = cluster;

        cmd.run(&inputs, &mut handle)
    }
}

fn run_generate(
    output: PathBuf,
    sizes: String,
    seed: u64,
    mode: String,
    sorted: String,
    no_sort: bool,
    a: Option<String>,
    b: Option<String>,
    hotspot_frac: f64,
    hotspot_weight: f64,
    len_min: u32,
    len_max: u32,
    force: bool,
) -> Result<(), BedError> {
    use grit_genomics::commands::generate::{
        GenerateCommand, GenerateConfig, GenerateMode, SizeSpec, SortMode,
    };

    // Parse mode
    let mode = GenerateMode::from_str(&mode).ok_or_else(|| {
        BedError::InvalidFormat(format!(
            "Invalid mode '{}'. Use: balanced, skewed-a-gt-b, skewed-b-gt-a, identical, clustered, all",
            mode
        ))
    })?;

    // Parse sizes
    let sizes: Result<Vec<SizeSpec>, _> = sizes
        .split(',')
        .map(|s| {
            SizeSpec::from_str(s.trim()).ok_or_else(|| {
                BedError::InvalidFormat(format!(
                    "Invalid size '{}'. Use formats like 1K, 5M, 100",
                    s
                ))
            })
        })
        .collect();
    let sizes = sizes?;

    // Parse sort mode
    let sorted = if no_sort {
        SortMode::No
    } else {
        SortMode::from_str(&sorted).ok_or_else(|| {
            BedError::InvalidFormat(format!(
                "Invalid sorted value '{}'. Use: yes, no, auto",
                sorted
            ))
        })?
    };

    // Parse custom A/B sizes
    let custom_a = a
        .map(|s| {
            SizeSpec::from_str(&s)
                .map(|spec| spec.count)
                .ok_or_else(|| BedError::InvalidFormat(format!("Invalid --a size '{}'", s)))
        })
        .transpose()?;
    let custom_b = b
        .map(|s| {
            SizeSpec::from_str(&s)
                .map(|spec| spec.count)
                .ok_or_else(|| BedError::InvalidFormat(format!("Invalid --b size '{}'", s)))
        })
        .transpose()?;

    let config = GenerateConfig {
        output_dir: output,
        sizes,
        seed,
        mode,
        sorted,
        custom_a,
        custom_b,
        hotspot_frac,
        hotspot_weight,
        len_min,
        len_max,
        force,
    };

    let cmd = GenerateCommand::new(config);
    cmd.run()?;

    Ok(())
}
