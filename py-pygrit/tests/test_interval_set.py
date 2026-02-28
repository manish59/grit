"""Unit tests for IntervalSet class."""

import pytest
import numpy as np
import pygrit
from pygrit import Interval, IntervalSet


class TestIntervalSetCreation:
    """Tests for IntervalSet construction."""

    def test_create_empty_interval_set(self):
        """Test creating empty IntervalSet."""
        iset = IntervalSet()
        assert len(iset) == 0

    def test_from_intervals(self):
        """Test creating IntervalSet from list of intervals."""
        intervals = [
            Interval("chr1", 100, 200),
            Interval("chr1", 300, 400),
            Interval("chr2", 100, 200),
        ]
        iset = IntervalSet.from_intervals(intervals)
        assert len(iset) == 3

    def test_add_interval(self):
        """Test adding intervals to set."""
        iset = IntervalSet()
        iset.add(Interval("chr1", 100, 200))
        iset.add(Interval("chr1", 300, 400))
        assert len(iset) == 2

    def test_repr(self):
        """Test string representation."""
        intervals = [Interval("chr1", 100, 200), Interval("chr1", 300, 400)]
        iset = IntervalSet.from_intervals(intervals)
        assert repr(iset) == "IntervalSet(2 intervals)"


class TestIntervalSetAccess:
    """Tests for IntervalSet element access."""

    def test_getitem(self):
        """Test accessing intervals by index."""
        intervals = [
            Interval("chr1", 100, 200),
            Interval("chr1", 300, 400),
        ]
        iset = IntervalSet.from_intervals(intervals)
        assert iset[0] == Interval("chr1", 100, 200)
        assert iset[1] == Interval("chr1", 300, 400)

    def test_getitem_out_of_bounds(self):
        """Test index out of bounds."""
        iset = IntervalSet.from_intervals([Interval("chr1", 100, 200)])
        with pytest.raises(ValueError, match="out of bounds"):
            _ = iset[10]

    def test_to_list(self):
        """Test converting to list."""
        intervals = [
            Interval("chr1", 100, 200),
            Interval("chr1", 300, 400),
        ]
        iset = IntervalSet.from_intervals(intervals)
        result = iset.to_list()
        assert len(result) == 2
        assert result[0] == Interval("chr1", 100, 200)


class TestIntervalSetMerge:
    """Tests for IntervalSet merge operation."""

    def test_merge_overlapping(self):
        """Test merging overlapping intervals."""
        intervals = [
            Interval("chr1", 100, 200),
            Interval("chr1", 150, 250),
            Interval("chr1", 200, 300),
        ]
        iset = IntervalSet.from_intervals(intervals)
        merged = iset.merge()
        assert len(merged) == 1
        result = merged.to_list()
        assert result[0].start == 100
        assert result[0].end == 300

    def test_merge_non_overlapping(self):
        """Test merging non-overlapping intervals."""
        intervals = [
            Interval("chr1", 100, 200),
            Interval("chr1", 300, 400),
        ]
        iset = IntervalSet.from_intervals(intervals)
        merged = iset.merge()
        assert len(merged) == 2

    def test_merge_with_distance(self):
        """Test merging with distance parameter."""
        intervals = [
            Interval("chr1", 100, 200),
            Interval("chr1", 250, 350),  # 50bp gap
        ]
        iset = IntervalSet.from_intervals(intervals)

        # Without distance - should not merge
        merged = iset.merge(distance=0)
        assert len(merged) == 2

        # With distance - should merge
        merged = iset.merge(distance=100)
        assert len(merged) == 1
        result = merged.to_list()
        assert result[0].start == 100
        assert result[0].end == 350

    def test_merge_multiple_chroms(self):
        """Test merging intervals on different chromosomes."""
        intervals = [
            Interval("chr1", 100, 200),
            Interval("chr1", 150, 250),
            Interval("chr2", 100, 200),
            Interval("chr2", 150, 250),
        ]
        iset = IntervalSet.from_intervals(intervals)
        merged = iset.merge()
        assert len(merged) == 2  # One merged interval per chromosome


class TestIntervalSetIntersect:
    """Tests for IntervalSet intersect operation."""

    def test_intersect_basic(self):
        """Test basic intersection."""
        set_a = IntervalSet.from_intervals([
            Interval("chr1", 100, 200),
            Interval("chr1", 300, 400),
        ])
        set_b = IntervalSet.from_intervals([
            Interval("chr1", 150, 250),  # Overlaps first
            Interval("chr1", 500, 600),  # No overlap
        ])
        result = set_a.intersect(set_b)
        assert len(result) == 1  # Only first interval overlaps

    def test_intersect_no_overlap(self):
        """Test intersection with no overlaps."""
        set_a = IntervalSet.from_intervals([
            Interval("chr1", 100, 200),
        ])
        set_b = IntervalSet.from_intervals([
            Interval("chr1", 300, 400),
        ])
        result = set_a.intersect(set_b)
        assert len(result) == 0

    def test_intersect_with_fraction(self):
        """Test intersection with minimum fraction."""
        set_a = IntervalSet.from_intervals([
            Interval("chr1", 100, 200),  # 100bp
        ])
        set_b = IntervalSet.from_intervals([
            Interval("chr1", 150, 160),  # 10bp overlap = 10%
        ])
        # With 50% fraction requirement, should not match
        result = set_a.intersect(set_b, fraction=0.5)
        assert len(result) == 0


class TestIntervalSetNonOverlapping:
    """Tests for IntervalSet non_overlapping operation."""

    def test_non_overlapping(self):
        """Test finding non-overlapping intervals."""
        set_a = IntervalSet.from_intervals([
            Interval("chr1", 100, 200),
            Interval("chr1", 300, 400),
            Interval("chr1", 500, 600),
        ])
        set_b = IntervalSet.from_intervals([
            Interval("chr1", 150, 250),  # Overlaps first
        ])
        result = set_a.non_overlapping(set_b)
        # Should return 300-400 and 500-600
        assert len(result) == 2


class TestIntervalSetSort:
    """Tests for IntervalSet sort operation."""

    def test_sort(self):
        """Test sorting intervals."""
        iset = IntervalSet()
        iset.add(Interval("chr1", 300, 400))
        iset.add(Interval("chr1", 100, 200))
        iset.add(Interval("chr2", 50, 100))
        iset.sort()
        result = iset.to_list()
        assert result[0].start == 100
        assert result[1].start == 300
        assert result[2].chrom == "chr2"


class TestIntervalSetNumpy:
    """Tests for IntervalSet NumPy integration."""

    def test_to_numpy(self):
        """Test conversion to NumPy array."""
        intervals = [
            Interval("chr1", 100, 200),
            Interval("chr1", 300, 400),
        ]
        iset = IntervalSet.from_intervals(intervals)
        arr = iset.to_numpy()
        assert arr.shape == (2, 2)
        assert arr[0, 0] == 100
        assert arr[0, 1] == 200
        assert arr[1, 0] == 300
        assert arr[1, 1] == 400

    def test_to_numpy_empty(self):
        """Test conversion of empty set to NumPy."""
        iset = IntervalSet()
        arr = iset.to_numpy()
        assert arr.shape == (0, 2)
