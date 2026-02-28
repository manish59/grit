"""Unit tests for file-based streaming functions."""

import pytest
import pygrit
from pygrit import Interval


class TestIntersect:
    """Tests for pygrit.intersect function."""

    def test_intersect_basic(self, sample_bed_a, sample_bed_b):
        """Test basic intersection."""
        results = pygrit.intersect(str(sample_bed_a), str(sample_bed_b))
        assert isinstance(results, list)
        assert len(results) > 0
        assert all(isinstance(iv, Interval) for iv in results)

    def test_intersect_to_file(self, sample_bed_a, sample_bed_b, temp_dir):
        """Test intersection with output to file."""
        output = temp_dir / "out.bed"
        result = pygrit.intersect(str(sample_bed_a), str(sample_bed_b), output=str(output))
        assert result is None  # Returns None when writing to file
        assert output.exists()
        content = output.read_text()
        assert len(content) > 0

    def test_intersect_write_a(self, sample_bed_a, sample_bed_b, temp_dir):
        """Test intersection with -wa flag."""
        output = temp_dir / "out.bed"
        pygrit.intersect(str(sample_bed_a), str(sample_bed_b), output=str(output), write_a=True)
        content = output.read_text()
        # Should contain original A coordinates
        assert "100\t200" in content or "300\t400" in content

    def test_intersect_write_b(self, sample_bed_a, sample_bed_b, temp_dir):
        """Test intersection with -wb flag."""
        output = temp_dir / "out.bed"
        pygrit.intersect(str(sample_bed_a), str(sample_bed_b), output=str(output), write_b=True)
        content = output.read_text()
        # Should contain B coordinates
        assert "150\t250" in content or "550\t650" in content

    def test_intersect_no_overlap(self, sample_bed_a, sample_bed_b):
        """Test intersection -v flag (no overlap)."""
        results = pygrit.intersect(str(sample_bed_a), str(sample_bed_b), no_overlap=True)
        # Should return A intervals with no overlap
        assert isinstance(results, list)

    def test_intersect_unique(self, sample_bed_a, sample_bed_b):
        """Test intersection -u flag (unique)."""
        results = pygrit.intersect(str(sample_bed_a), str(sample_bed_b), unique=True)
        assert isinstance(results, list)

    def test_intersect_fraction(self, temp_dir):
        """Test intersection with fraction requirement."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t100\t200\n")
        b_file.write_text("chr1\t150\t160\n")  # 10bp overlap = 10%

        # Without fraction - should match
        results = pygrit.intersect(str(a_file), str(b_file))
        assert len(results) >= 1

        # With 50% fraction - should not match
        results = pygrit.intersect(str(a_file), str(b_file), fraction=0.5)
        assert len(results) == 0

    def test_intersect_reciprocal(self, temp_dir):
        """Test intersection with reciprocal fraction."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t100\t200\n")  # 100bp
        b_file.write_text("chr1\t50\t250\n")   # 200bp, overlap=100bp

        # A overlap fraction = 100%, B overlap fraction = 50%
        # With reciprocal 60%, should not match (B doesn't meet)
        results = pygrit.intersect(str(a_file), str(b_file), fraction=0.6, reciprocal=True)
        assert len(results) == 0

    def test_intersect_empty_file(self, sample_bed_a, empty_bed):
        """Test intersection with empty file."""
        results = pygrit.intersect(str(sample_bed_a), str(empty_bed))
        assert results == []

    def test_intersect_file_not_found(self, temp_dir):
        """Test intersection with missing file."""
        with pytest.raises(IOError):
            pygrit.intersect(str(temp_dir / "nonexistent.bed"), str(temp_dir / "also_missing.bed"))


class TestMerge:
    """Tests for pygrit.merge function."""

    def test_merge_basic(self, overlapping_bed):
        """Test basic merge."""
        results = pygrit.merge(str(overlapping_bed))
        assert isinstance(results, list)
        # 100-200, 150-250, 200-300 should merge into one
        # 500-600 stays separate
        assert len(results) == 2

    def test_merge_to_file(self, overlapping_bed, temp_dir):
        """Test merge with output to file."""
        output = temp_dir / "merged.bed"
        result = pygrit.merge(str(overlapping_bed), output=str(output))
        assert result is None
        assert output.exists()

    def test_merge_with_distance(self, temp_dir):
        """Test merge with distance parameter."""
        bed_file = temp_dir / "gaps.bed"
        bed_file.write_text("chr1\t100\t200\nchr1\t250\t350\n")  # 50bp gap

        # Without distance
        results = pygrit.merge(str(bed_file), distance=0)
        assert len(results) == 2

        # With distance=100 (covers 50bp gap)
        results = pygrit.merge(str(bed_file), distance=100)
        assert len(results) == 1
        assert results[0].start == 100
        assert results[0].end == 350

    def test_merge_strand(self, stranded_bed):
        """Test merge with strand-specific mode."""
        # Without strand - should merge all overlapping
        results = pygrit.merge(str(stranded_bed), strand=False)
        all_merged_count = len(results)

        # With strand - should keep strands separate
        results = pygrit.merge(str(stranded_bed), strand=True)
        strand_merged_count = len(results)

        # Strand-specific should produce more intervals
        assert strand_merged_count >= all_merged_count

    def test_merge_empty_file(self, empty_bed):
        """Test merge of empty file."""
        results = pygrit.merge(str(empty_bed))
        assert results == []


class TestSubtract:
    """Tests for pygrit.subtract function."""

    def test_subtract_basic(self, sample_bed_a, sample_bed_b):
        """Test basic subtraction."""
        results = pygrit.subtract(str(sample_bed_a), str(sample_bed_b))
        assert isinstance(results, list)

    def test_subtract_to_file(self, sample_bed_a, sample_bed_b, temp_dir):
        """Test subtraction with output to file."""
        output = temp_dir / "subtracted.bed"
        result = pygrit.subtract(str(sample_bed_a), str(sample_bed_b), output=str(output))
        assert result is None
        assert output.exists()

    def test_subtract_complete_overlap(self, temp_dir):
        """Test subtraction where A is completely covered by B."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t100\t200\n")
        b_file.write_text("chr1\t50\t250\n")

        results = pygrit.subtract(str(a_file), str(b_file))
        assert len(results) == 0

    def test_subtract_partial_overlap(self, temp_dir):
        """Test subtraction with partial overlap."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t100\t300\n")
        b_file.write_text("chr1\t150\t200\n")

        results = pygrit.subtract(str(a_file), str(b_file))
        # Should produce two intervals: 100-150 and 200-300
        assert len(results) == 2

    def test_subtract_remove_entire(self, temp_dir):
        """Test subtraction with -A flag (remove entire)."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t100\t200\nchr1\t300\t400\n")
        b_file.write_text("chr1\t120\t140\n")  # Partially overlaps first

        results = pygrit.subtract(str(a_file), str(b_file), remove_entire=True)
        # First interval should be entirely removed
        assert len(results) == 1
        assert results[0].start == 300


class TestCoverage:
    """Tests for pygrit.coverage function."""

    def test_coverage_basic(self, sample_bed_a, sample_bed_b):
        """Test basic coverage."""
        result = pygrit.coverage(str(sample_bed_a), str(sample_bed_b))
        assert isinstance(result, str)
        assert len(result) > 0

    def test_coverage_to_file(self, sample_bed_a, sample_bed_b, temp_dir):
        """Test coverage with output to file."""
        output = temp_dir / "coverage.bed"
        result = pygrit.coverage(str(sample_bed_a), str(sample_bed_b), output=str(output))
        assert result is None
        assert output.exists()

    def test_coverage_histogram(self, sample_bed_a, sample_bed_b):
        """Test coverage with histogram output."""
        result = pygrit.coverage(str(sample_bed_a), str(sample_bed_b), histogram=True)
        assert isinstance(result, str)

    def test_coverage_mean(self, sample_bed_a, sample_bed_b):
        """Test coverage with mean depth."""
        result = pygrit.coverage(str(sample_bed_a), str(sample_bed_b), mean=True)
        assert isinstance(result, str)


class TestClosest:
    """Tests for pygrit.closest function."""

    def test_closest_basic(self, sample_bed_a, sample_bed_b):
        """Test basic closest."""
        result = pygrit.closest(str(sample_bed_a), str(sample_bed_b))
        assert isinstance(result, str)
        assert len(result) > 0

    def test_closest_to_file(self, sample_bed_a, sample_bed_b, temp_dir):
        """Test closest with output to file."""
        output = temp_dir / "closest.bed"
        result = pygrit.closest(str(sample_bed_a), str(sample_bed_b), output=str(output))
        assert result is None
        assert output.exists()

    def test_closest_ignore_overlaps(self, temp_dir):
        """Test closest with ignore_overlaps."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t100\t200\n")
        b_file.write_text("chr1\t150\t250\nchr1\t300\t400\n")

        # Without ignore_overlaps - should find overlapping
        result = pygrit.closest(str(a_file), str(b_file))
        assert "150\t250" in result

        # With ignore_overlaps - should find 300-400
        result = pygrit.closest(str(a_file), str(b_file), ignore_overlaps=True)
        assert "300\t400" in result

    def test_closest_ignore_upstream(self, temp_dir):
        """Test closest with ignore_upstream."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t200\t300\n")
        b_file.write_text("chr1\t100\t150\nchr1\t350\t400\n")

        result = pygrit.closest(str(a_file), str(b_file), ignore_upstream=True)
        assert "350\t400" in result

    def test_closest_ignore_downstream(self, temp_dir):
        """Test closest with ignore_downstream."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t200\t300\n")
        b_file.write_text("chr1\t100\t150\nchr1\t350\t400\n")

        result = pygrit.closest(str(a_file), str(b_file), ignore_downstream=True)
        assert "100\t150" in result


class TestWindow:
    """Tests for pygrit.window function."""

    def test_window_basic(self, sample_bed_a, sample_bed_b):
        """Test basic window."""
        result = pygrit.window(str(sample_bed_a), str(sample_bed_b))
        assert isinstance(result, str)

    def test_window_to_file(self, sample_bed_a, sample_bed_b, temp_dir):
        """Test window with output to file."""
        output = temp_dir / "window.bed"
        result = pygrit.window(str(sample_bed_a), str(sample_bed_b), output=str(output))
        assert result is None
        assert output.exists()

    def test_window_size(self, temp_dir):
        """Test window with different sizes."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t100\t200\n")
        b_file.write_text("chr1\t500\t600\n")  # 300bp away

        # With small window - should not find
        result = pygrit.window(str(a_file), str(b_file), window=100)
        # Result should be empty or not contain the B interval

        # With large window - should find
        result = pygrit.window(str(a_file), str(b_file), window=500)
        assert "500\t600" in result

    def test_window_asymmetric(self, temp_dir):
        """Test window with left/right parameters."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t500\t600\n")
        b_file.write_text("chr1\t100\t200\nchr1\t800\t900\n")

        # Only look left
        result = pygrit.window(str(a_file), str(b_file), left=500, right=0)
        assert "100\t200" in result


class TestReadBed:
    """Tests for pygrit.read_bed function."""

    def test_read_bed_basic(self, sample_bed_a):
        """Test reading BED file."""
        result = pygrit.read_bed(str(sample_bed_a))
        assert len(result) == 5

    def test_read_bed_empty(self, empty_bed):
        """Test reading empty BED file."""
        result = pygrit.read_bed(str(empty_bed))
        assert len(result) == 0


class TestParseBed:
    """Tests for pygrit.parse_bed function."""

    def test_parse_bed_string(self):
        """Test parsing BED content from string."""
        content = "chr1\t100\t200\nchr1\t300\t400\n"
        result = pygrit.parse_bed(content)
        assert len(result) == 2
        assert result[0].chrom == "chr1"
        assert result[0].start == 100
        assert result[0].end == 200


class TestFromNumpy:
    """Tests for pygrit.from_numpy function."""

    def test_from_numpy_basic(self):
        """Test creating intervals from NumPy array."""
        import numpy as np
        arr = np.array([[100, 200], [300, 400]], dtype=np.int64)
        result = pygrit.from_numpy("chr1", arr)
        assert len(result) == 2
        assert result[0].chrom == "chr1"
        assert result[0].start == 100
        assert result[0].end == 200
