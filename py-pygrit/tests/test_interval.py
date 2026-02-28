"""Unit tests for Interval class."""

import pytest
import pygrit
from pygrit import Interval


class TestIntervalCreation:
    """Tests for Interval construction and validation."""

    def test_create_interval(self):
        """Test basic interval creation."""
        iv = Interval("chr1", 100, 200)
        assert iv.chrom == "chr1"
        assert iv.start == 100
        assert iv.end == 200

    def test_interval_length(self):
        """Test interval length calculation."""
        iv = Interval("chr1", 100, 200)
        assert len(iv) == 100

    def test_zero_length_interval(self):
        """Test zero-length interval (point)."""
        iv = Interval("chr1", 100, 100)
        assert len(iv) == 0

    def test_invalid_interval_start_greater_than_end(self):
        """Test that start > end raises error."""
        with pytest.raises(ValueError, match="start.*must be <= end"):
            Interval("chr1", 200, 100)

    def test_interval_repr(self):
        """Test interval string representation."""
        iv = Interval("chr1", 100, 200)
        assert repr(iv) == "Interval('chr1', 100, 200)"

    def test_interval_str(self):
        """Test interval BED format string."""
        iv = Interval("chr1", 100, 200)
        assert str(iv) == "chr1\t100\t200"


class TestIntervalComparison:
    """Tests for Interval comparison and hashing."""

    def test_interval_equality(self):
        """Test interval equality."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr1", 100, 200)
        assert iv1 == iv2

    def test_interval_inequality_chrom(self):
        """Test inequality with different chromosome."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr2", 100, 200)
        assert iv1 != iv2

    def test_interval_inequality_start(self):
        """Test inequality with different start."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr1", 150, 200)
        assert iv1 != iv2

    def test_interval_inequality_end(self):
        """Test inequality with different end."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr1", 100, 250)
        assert iv1 != iv2

    def test_interval_hash(self):
        """Test interval hashing for use in sets/dicts."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr1", 100, 200)
        iv3 = Interval("chr1", 100, 300)

        s = {iv1, iv2, iv3}
        assert len(s) == 2  # iv1 and iv2 should be the same

    def test_interval_in_dict(self):
        """Test using interval as dictionary key."""
        iv = Interval("chr1", 100, 200)
        d = {iv: "test"}
        assert d[Interval("chr1", 100, 200)] == "test"


class TestIntervalOverlap:
    """Tests for Interval overlap calculations."""

    def test_overlaps_true(self):
        """Test overlapping intervals."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr1", 150, 250)
        assert iv1.overlaps(iv2)
        assert iv2.overlaps(iv1)

    def test_overlaps_false_adjacent(self):
        """Test adjacent (non-overlapping) intervals."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr1", 200, 300)
        assert not iv1.overlaps(iv2)
        assert not iv2.overlaps(iv1)

    def test_overlaps_false_separate(self):
        """Test non-overlapping intervals."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr1", 300, 400)
        assert not iv1.overlaps(iv2)

    def test_overlaps_false_different_chrom(self):
        """Test intervals on different chromosomes."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr2", 100, 200)
        assert not iv1.overlaps(iv2)

    def test_overlaps_contained(self):
        """Test one interval contained within another."""
        iv1 = Interval("chr1", 100, 300)
        iv2 = Interval("chr1", 150, 250)
        assert iv1.overlaps(iv2)
        assert iv2.overlaps(iv1)

    def test_overlap_length(self):
        """Test overlap length calculation."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr1", 150, 250)
        assert iv1.overlap_length(iv2) == 50
        assert iv2.overlap_length(iv1) == 50

    def test_overlap_length_no_overlap(self):
        """Test overlap length when no overlap."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr1", 300, 400)
        assert iv1.overlap_length(iv2) == 0

    def test_overlap_length_contained(self):
        """Test overlap length when one contained in other."""
        iv1 = Interval("chr1", 100, 300)
        iv2 = Interval("chr1", 150, 250)
        assert iv1.overlap_length(iv2) == 100
        assert iv2.overlap_length(iv1) == 100


class TestIntervalDistance:
    """Tests for Interval distance calculations."""

    def test_distance_overlapping(self):
        """Test distance is 0 for overlapping intervals."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr1", 150, 250)
        assert iv1.distance_to(iv2) == 0

    def test_distance_upstream(self):
        """Test distance to upstream interval."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr1", 300, 400)
        assert iv1.distance_to(iv2) == 100  # 300 - 200

    def test_distance_downstream(self):
        """Test distance to downstream interval."""
        iv1 = Interval("chr1", 300, 400)
        iv2 = Interval("chr1", 100, 200)
        assert iv1.distance_to(iv2) == 100  # 300 - 200

    def test_distance_adjacent(self):
        """Test distance for adjacent intervals."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr1", 200, 300)
        assert iv1.distance_to(iv2) == 0  # Adjacent = 0 distance

    def test_distance_different_chrom(self):
        """Test distance returns None for different chromosomes."""
        iv1 = Interval("chr1", 100, 200)
        iv2 = Interval("chr2", 100, 200)
        assert iv1.distance_to(iv2) is None


class TestIntervalConversion:
    """Tests for Interval conversion methods."""

    def test_to_tuple(self):
        """Test conversion to tuple."""
        iv = Interval("chr1", 100, 200)
        assert iv.to_tuple() == ("chr1", 100, 200)
