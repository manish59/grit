"""Type stubs for pygrit - Python bindings for GRIT."""

from typing import overload
import numpy as np
import numpy.typing as npt

__version__: str
"""Package version string."""


class Interval:
    """A genomic interval with chromosome, start, and end coordinates.

    Coordinates are 0-based, half-open (BED format).

    Attributes:
        chrom: Chromosome name (e.g., "chr1", "chrX").
        start: Start position (0-based, inclusive).
        end: End position (exclusive).

    Example:
        >>> iv = Interval("chr1", 100, 200)
        >>> len(iv)
        100
        >>> iv.overlaps(Interval("chr1", 150, 250))
        True
    """

    chrom: str
    """Chromosome name."""

    start: int
    """Start position (0-based, inclusive)."""

    end: int
    """End position (exclusive)."""

    def __init__(self, chrom: str, start: int, end: int) -> None:
        """Create a new interval.

        Args:
            chrom: Chromosome name.
            start: Start position (0-based, inclusive).
            end: End position (exclusive).

        Raises:
            ValueError: If start > end.
        """
        ...

    def __len__(self) -> int:
        """Return the length of the interval (end - start)."""
        ...

    def __repr__(self) -> str:
        """Return Python representation: Interval('chr1', 100, 200)."""
        ...

    def __str__(self) -> str:
        """Return BED format string: chr1\\t100\\t200."""
        ...

    def __eq__(self, other: object) -> bool:
        """Check equality with another interval."""
        ...

    def __hash__(self) -> int:
        """Return hash for use in sets and dicts."""
        ...

    def overlaps(self, other: "Interval") -> bool:
        """Check if this interval overlaps another.

        Args:
            other: The interval to check against.

        Returns:
            True if intervals overlap, False otherwise.
        """
        ...

    def overlap_length(self, other: "Interval") -> int:
        """Calculate the number of overlapping bases.

        Args:
            other: The interval to check against.

        Returns:
            Number of overlapping bases (0 if no overlap).
        """
        ...

    def distance_to(self, other: "Interval") -> int | None:
        """Calculate the distance to another interval.

        Args:
            other: The interval to measure distance to.

        Returns:
            - 0 if intervals overlap
            - Positive integer for the gap between non-overlapping intervals
            - None if intervals are on different chromosomes
        """
        ...

    def to_tuple(self) -> tuple[str, int, int]:
        """Convert to a tuple (chrom, start, end)."""
        ...


class IntervalSet:
    """A collection of genomic intervals with bulk operations.

    Example:
        >>> intervals = IntervalSet.from_intervals([
        ...     Interval("chr1", 100, 200),
        ...     Interval("chr1", 150, 250),
        ... ])
        >>> merged = intervals.merge()
        >>> len(merged)
        1
    """

    def __init__(self) -> None:
        """Create an empty IntervalSet."""
        ...

    @staticmethod
    def from_intervals(intervals: list[Interval]) -> "IntervalSet":
        """Create an IntervalSet from a list of Interval objects.

        Args:
            intervals: List of Interval objects.

        Returns:
            A new IntervalSet containing the intervals.
        """
        ...

    def __len__(self) -> int:
        """Return the number of intervals."""
        ...

    def __getitem__(self, index: int) -> Interval:
        """Access interval by index.

        Args:
            index: Index (supports negative indexing).

        Returns:
            The interval at the given index.

        Raises:
            IndexError: If index is out of range.
        """
        ...

    def add(self, interval: Interval) -> None:
        """Add an interval to the set.

        Args:
            interval: The interval to add.
        """
        ...

    def to_list(self) -> list[Interval]:
        """Convert to a list of Interval objects."""
        ...

    def merge(self, distance: int = 0) -> "IntervalSet":
        """Merge overlapping or nearby intervals.

        Args:
            distance: Maximum gap to bridge when merging. Default is 0.

        Returns:
            A new IntervalSet with merged intervals.
        """
        ...

    def intersect(
        self,
        other: "IntervalSet",
        fraction: float | None = None,
        reciprocal: bool = False,
    ) -> "IntervalSet":
        """Find intervals that overlap with another IntervalSet.

        Args:
            other: IntervalSet to intersect with.
            fraction: Minimum overlap fraction (0.0-1.0).
            reciprocal: Require reciprocal overlap fraction.

        Returns:
            Intervals from self that overlap with other.
        """
        ...

    def non_overlapping(self, other: "IntervalSet") -> "IntervalSet":
        """Find intervals that don't overlap with another IntervalSet.

        Args:
            other: IntervalSet to check against.

        Returns:
            Intervals from self that don't overlap with other.
        """
        ...

    def sort(self) -> None:
        """Sort intervals in place by chromosome and start position."""
        ...

    def to_numpy(self) -> npt.NDArray[np.int64]:
        """Convert to a NumPy array.

        Returns:
            Array with shape (n, 2) containing [start, end] columns.
        """
        ...


# File-based streaming functions

@overload
def intersect(
    a: str,
    b: str,
    *,
    output: None = None,
    write_a: bool = False,
    write_b: bool = False,
    fraction: float | None = None,
    reciprocal: bool = False,
    count: bool = False,
    unique: bool = False,
    no_overlap: bool = False,
) -> list[Interval]:
    ...


@overload
def intersect(
    a: str,
    b: str,
    *,
    output: str,
    write_a: bool = False,
    write_b: bool = False,
    fraction: float | None = None,
    reciprocal: bool = False,
    count: bool = False,
    unique: bool = False,
    no_overlap: bool = False,
) -> None:
    ...


def intersect(
    a: str,
    b: str,
    *,
    output: str | None = None,
    write_a: bool = False,
    write_b: bool = False,
    fraction: float | None = None,
    reciprocal: bool = False,
    count: bool = False,
    unique: bool = False,
    no_overlap: bool = False,
) -> list[Interval] | None:
    """Find overlapping intervals between two BED files.

    Args:
        a: Path to file A.
        b: Path to file B.
        output: Output file path. If None, returns list.
        write_a: Include original A record in output.
        write_b: Include original B record in output.
        fraction: Minimum overlap as fraction of A (0.0-1.0).
        reciprocal: Require fraction overlap in both A and B.
        count: Report count of overlaps per A interval.
        unique: Report each A interval only once.
        no_overlap: Report A intervals with no B overlap.

    Returns:
        List of overlapping Interval objects if output is None,
        otherwise None (results written to file).

    Raises:
        IOError: File not found or I/O error.
        RuntimeError: Processing error (e.g., unsorted input).
        ValueError: Invalid parameter values.
    """
    ...


@overload
def merge(
    input: str,
    *,
    output: None = None,
    distance: int = 0,
    strand: bool = False,
) -> list[Interval]:
    ...


@overload
def merge(
    input: str,
    *,
    output: str,
    distance: int = 0,
    strand: bool = False,
) -> None:
    ...


def merge(
    input: str,
    *,
    output: str | None = None,
    distance: int = 0,
    strand: bool = False,
) -> list[Interval] | None:
    """Merge overlapping or nearby intervals.

    Args:
        input: Path to input BED file.
        output: Output file path. If None, returns list.
        distance: Maximum gap to bridge when merging.
        strand: Only merge intervals on same strand.

    Returns:
        List of merged Interval objects if output is None,
        otherwise None (results written to file).
    """
    ...


@overload
def subtract(
    a: str,
    b: str,
    *,
    output: None = None,
    remove_entire: bool = False,
    fraction: float | None = None,
    reciprocal: bool = False,
) -> list[Interval]:
    ...


@overload
def subtract(
    a: str,
    b: str,
    *,
    output: str,
    remove_entire: bool = False,
    fraction: float | None = None,
    reciprocal: bool = False,
) -> None:
    ...


def subtract(
    a: str,
    b: str,
    *,
    output: str | None = None,
    remove_entire: bool = False,
    fraction: float | None = None,
    reciprocal: bool = False,
) -> list[Interval] | None:
    """Subtract B intervals from A intervals.

    Args:
        a: Path to file A.
        b: Path to file B.
        output: Output file path. If None, returns list.
        remove_entire: Remove entire A interval if any overlap.
        fraction: Minimum overlap fraction to subtract.
        reciprocal: Require reciprocal fraction.

    Returns:
        List of resulting Interval objects if output is None,
        otherwise None (results written to file).
    """
    ...


@overload
def coverage(
    a: str,
    b: str,
    *,
    output: None = None,
    histogram: bool = False,
    mean: bool = False,
) -> str:
    ...


@overload
def coverage(
    a: str,
    b: str,
    *,
    output: str,
    histogram: bool = False,
    mean: bool = False,
) -> None:
    ...


def coverage(
    a: str,
    b: str,
    *,
    output: str | None = None,
    histogram: bool = False,
    mean: bool = False,
) -> str | None:
    """Calculate coverage of A regions by B features.

    Args:
        a: Path to regions file.
        b: Path to reads/features file.
        output: Output file path. If None, returns string.
        histogram: Report depth histogram.
        mean: Report mean depth per region.

    Returns:
        Coverage output as string if output is None,
        otherwise None (results written to file).
    """
    ...


@overload
def closest(
    a: str,
    b: str,
    *,
    output: None = None,
    ignore_overlaps: bool = False,
    ignore_upstream: bool = False,
    ignore_downstream: bool = False,
) -> str:
    ...


@overload
def closest(
    a: str,
    b: str,
    *,
    output: str,
    ignore_overlaps: bool = False,
    ignore_upstream: bool = False,
    ignore_downstream: bool = False,
) -> None:
    ...


def closest(
    a: str,
    b: str,
    *,
    output: str | None = None,
    ignore_overlaps: bool = False,
    ignore_upstream: bool = False,
    ignore_downstream: bool = False,
) -> str | None:
    """Find closest B interval for each A interval.

    Args:
        a: Path to file A.
        b: Path to file B.
        output: Output file path. If None, returns string.
        ignore_overlaps: Skip overlapping intervals.
        ignore_upstream: Only look downstream (3').
        ignore_downstream: Only look upstream (5').

    Returns:
        Closest output as string if output is None,
        otherwise None (results written to file).
    """
    ...


@overload
def window(
    a: str,
    b: str,
    *,
    output: None = None,
    window: int = 1000,
    left: int | None = None,
    right: int | None = None,
    count: bool = False,
    no_overlap: bool = False,
) -> str:
    ...


@overload
def window(
    a: str,
    b: str,
    *,
    output: str,
    window: int = 1000,
    left: int | None = None,
    right: int | None = None,
    count: bool = False,
    no_overlap: bool = False,
) -> None:
    ...


def window(
    a: str,
    b: str,
    *,
    output: str | None = None,
    window: int = 1000,
    left: int | None = None,
    right: int | None = None,
    count: bool = False,
    no_overlap: bool = False,
) -> str | None:
    """Find B intervals within window distance of A intervals.

    Args:
        a: Path to file A.
        b: Path to file B.
        output: Output file path. If None, returns string.
        window: Window size in base pairs.
        left: Left window (overrides window).
        right: Right window (overrides window).
        count: Report count of B in window.
        no_overlap: Only report non-overlapping.

    Returns:
        Window output as string if output is None,
        otherwise None (results written to file).
    """
    ...


@overload
def sort(
    input: str,
    *,
    output: None = None,
    genome: str | None = None,
    reverse: bool = False,
) -> str:
    ...


@overload
def sort(
    input: str,
    *,
    output: str,
    genome: str | None = None,
    reverse: bool = False,
) -> None:
    ...


def sort(
    input: str,
    *,
    output: str | None = None,
    genome: str | None = None,
    reverse: bool = False,
) -> str | None:
    """Sort a BED file by chromosome and start position.

    Args:
        input: Path to input BED file.
        output: Output file path. If None, returns string.
        genome: Genome file for chromosome ordering.
        reverse: Sort in reverse order.

    Returns:
        Sorted output as string if output is None,
        otherwise None (results written to file).
    """
    ...


@overload
def slop(
    input: str,
    genome: str,
    *,
    output: None = None,
    both: float | None = None,
    left: float | None = None,
    right: float | None = None,
    pct: bool = False,
) -> str:
    ...


@overload
def slop(
    input: str,
    genome: str,
    *,
    output: str,
    both: float | None = None,
    left: float | None = None,
    right: float | None = None,
    pct: bool = False,
) -> None:
    ...


def slop(
    input: str,
    genome: str,
    *,
    output: str | None = None,
    both: float | None = None,
    left: float | None = None,
    right: float | None = None,
    pct: bool = False,
) -> str | None:
    """Extend intervals by a specified amount on each side.

    Args:
        input: Path to input BED file.
        genome: Path to genome file (chromosome sizes).
        output: Output file path. If None, returns string.
        both: Extend both sides by this amount.
        left: Extend left (5') side.
        right: Extend right (3') side.
        pct: Interpret values as percentage of interval length.

    Returns:
        Extended intervals as string if output is None,
        otherwise None (results written to file).
    """
    ...


@overload
def complement(
    input: str,
    genome: str,
    *,
    output: None = None,
) -> str:
    ...


@overload
def complement(
    input: str,
    genome: str,
    *,
    output: str,
) -> None:
    ...


def complement(
    input: str,
    genome: str,
    *,
    output: str | None = None,
) -> str | None:
    """Calculate the complement of intervals (gaps between intervals).

    Args:
        input: Path to input BED file.
        genome: Path to genome file (chromosome sizes).
        output: Output file path. If None, returns string.

    Returns:
        Complement intervals as string if output is None,
        otherwise None (results written to file).
    """
    ...


@overload
def genomecov(
    input: str,
    genome: str,
    *,
    output: None = None,
    bg: bool = False,
    bga: bool = False,
    scale: float = 1.0,
) -> str:
    ...


@overload
def genomecov(
    input: str,
    genome: str,
    *,
    output: str,
    bg: bool = False,
    bga: bool = False,
    scale: float = 1.0,
) -> None:
    ...


def genomecov(
    input: str,
    genome: str,
    *,
    output: str | None = None,
    bg: bool = False,
    bga: bool = False,
    scale: float = 1.0,
) -> str | None:
    """Calculate genome-wide coverage.

    Args:
        input: Path to input BED file.
        genome: Path to genome file (chromosome sizes).
        output: Output file path. If None, returns string.
        bg: Output BedGraph format (non-zero regions only).
        bga: Output BedGraph format (all regions, including zero).
        scale: Scale coverage by this factor.

    Returns:
        Coverage output as string if output is None,
        otherwise None (results written to file).
    """
    ...


@overload
def jaccard(
    a: str,
    b: str,
    *,
    output: None = None,
) -> str:
    ...


@overload
def jaccard(
    a: str,
    b: str,
    *,
    output: str,
) -> None:
    ...


def jaccard(
    a: str,
    b: str,
    *,
    output: str | None = None,
) -> str | None:
    """Calculate Jaccard similarity between two BED files.

    Args:
        a: Path to file A.
        b: Path to file B.
        output: Output file path. If None, returns string.

    Returns:
        Jaccard statistics as string if output is None,
        otherwise None (results written to file).
    """
    ...


@overload
def multiinter(
    inputs: list[str],
    *,
    output: None = None,
    cluster: bool = False,
) -> str:
    ...


@overload
def multiinter(
    inputs: list[str],
    *,
    output: str,
    cluster: bool = False,
) -> None:
    ...


def multiinter(
    inputs: list[str],
    *,
    output: str | None = None,
    cluster: bool = False,
) -> str | None:
    """Find intervals that overlap across multiple BED files.

    Args:
        inputs: List of paths to input BED files (minimum 2).
        output: Output file path. If None, returns string.
        cluster: Cluster overlapping intervals.

    Returns:
        Multiinter output as string if output is None,
        otherwise None (results written to file).

    Raises:
        ValueError: If fewer than 2 input files provided.
    """
    ...


def generate(
    output_dir: str,
    num_intervals: int = 10000,
    num_chroms: int = 5,
    chrom_size: int = 100000000,
    len_min: int = 50,
    len_max: int = 5000,
    mode: str = "uniform",
    num_files: int = 2,
    sorted: bool = True,
    seed: int | None = None,
) -> dict[str, int | float]:
    """Generate synthetic BED files for testing and benchmarking.

    Args:
        output_dir: Output directory path.
        num_intervals: Number of intervals per file.
        num_chroms: Number of chromosomes.
        chrom_size: Size of each chromosome.
        len_min: Minimum interval length.
        len_max: Maximum interval length.
        mode: Distribution mode: "uniform", "balanced", "clustered".
        num_files: Number of files to generate.
        sorted: Generate sorted output.
        seed: Random seed for reproducibility.

    Returns:
        Dict with generation statistics:
            - total_files: Number of files generated
            - total_intervals: Total intervals across all files
            - elapsed_secs: Time taken in seconds

    Raises:
        ValueError: If invalid mode specified.
    """
    ...


# I/O functions

def read_bed(path: str) -> IntervalSet:
    """Read intervals from a BED file.

    Args:
        path: Path to the BED file.

    Returns:
        IntervalSet containing all intervals from the file.

    Raises:
        IOError: File not found or I/O error.
    """
    ...


def parse_bed(content: str) -> IntervalSet:
    """Parse intervals from a BED-formatted string.

    Args:
        content: BED-formatted string content.

    Returns:
        IntervalSet containing parsed intervals.
    """
    ...


def from_numpy(chrom: str, arr: npt.NDArray[np.int64]) -> IntervalSet:
    """Create an IntervalSet from a NumPy array.

    Args:
        chrom: Chromosome name for all intervals.
        arr: Array with shape (n, 2) containing [start, end] pairs.

    Returns:
        IntervalSet with n intervals, all on the specified chromosome.

    Raises:
        ValueError: If array shape is not (n, 2).
    """
    ...
