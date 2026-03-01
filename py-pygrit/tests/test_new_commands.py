"""Unit tests for new pygrit commands: sort, slop, complement, genomecov, jaccard, multiinter, generate."""

import os
import tempfile
import pytest
import pygrit


@pytest.fixture
def genome_file(temp_dir):
    """Create a genome file with chromosome sizes."""
    genome_path = temp_dir / "genome.txt"
    genome_path.write_text("chr1\t1000\nchr2\t800\nchr3\t600\n")
    return genome_path


@pytest.fixture
def unsorted_bed(temp_dir):
    """Create an unsorted BED file."""
    bed_path = temp_dir / "unsorted.bed"
    bed_path.write_text(
        "chr2\t100\t200\n"
        "chr1\t300\t400\n"
        "chr1\t100\t200\n"
        "chr3\t50\t150\n"
        "chr2\t200\t300\n"
    )
    return bed_path


class TestSort:
    """Tests for pygrit.sort function."""

    def test_sort_basic(self, unsorted_bed):
        """Test basic sorting."""
        result = pygrit.sort(str(unsorted_bed))
        lines = result.strip().split("\n")

        # Should be sorted by chrom, then start
        assert len(lines) == 5
        assert lines[0].startswith("chr1\t100\t200")
        assert lines[1].startswith("chr1\t300\t400")
        assert lines[2].startswith("chr2\t100\t200")

    def test_sort_to_file(self, unsorted_bed, temp_dir):
        """Test sorting with output to file."""
        output = temp_dir / "sorted.bed"
        result = pygrit.sort(str(unsorted_bed), output=str(output))

        assert result is None
        assert output.exists()
        content = output.read_text()
        assert "chr1\t100\t200" in content

    def test_sort_with_genome(self, unsorted_bed, genome_file):
        """Test sorting with genome file for chromosome ordering."""
        result = pygrit.sort(str(unsorted_bed), genome=str(genome_file))
        lines = result.strip().split("\n")

        # Should be sorted according to genome order
        assert len(lines) == 5
        # chr1 should come first based on genome order
        assert lines[0].startswith("chr1")

    def test_sort_reverse(self, unsorted_bed):
        """Test reverse sorting."""
        result = pygrit.sort(str(unsorted_bed), reverse=True)
        lines = result.strip().split("\n")

        # Should be reverse sorted
        assert len(lines) == 5
        # Last chromosome (chr3) should come first
        assert lines[0].startswith("chr3")


class TestSlop:
    """Tests for pygrit.slop function."""

    def test_slop_both(self, sample_bed_a, genome_file):
        """Test extending both sides."""
        result = pygrit.slop(str(sample_bed_a), str(genome_file), both=50.0)
        assert result is not None
        # Check that intervals were extended
        lines = result.strip().split("\n")
        assert len(lines) > 0

    def test_slop_left_right(self, temp_dir, genome_file):
        """Test extending left and right differently."""
        bed_file = temp_dir / "test.bed"
        bed_file.write_text("chr1\t100\t200\n")

        result = pygrit.slop(str(bed_file), str(genome_file), left=10.0, right=20.0)
        lines = result.strip().split("\n")

        # Should be extended: start=90, end=220
        assert "90\t220" in lines[0]

    def test_slop_to_file(self, sample_bed_a, genome_file, temp_dir):
        """Test slop with output to file."""
        output = temp_dir / "slopped.bed"
        result = pygrit.slop(str(sample_bed_a), str(genome_file), both=10.0, output=str(output))

        assert result is None
        assert output.exists()

    def test_slop_boundary(self, temp_dir, genome_file):
        """Test that slop respects chromosome boundaries."""
        bed_file = temp_dir / "test.bed"
        bed_file.write_text("chr1\t10\t50\n")  # Close to start

        result = pygrit.slop(str(bed_file), str(genome_file), both=20.0)
        lines = result.strip().split("\n")

        # Start should be clamped to 0
        assert lines[0].startswith("chr1\t0\t")

    def test_slop_percentage(self, temp_dir, genome_file):
        """Test slop with percentage mode."""
        bed_file = temp_dir / "test.bed"
        bed_file.write_text("chr1\t100\t200\n")  # 100bp interval

        result = pygrit.slop(str(bed_file), str(genome_file), both=0.5, pct=True)
        lines = result.strip().split("\n")

        # 50% of 100bp = 50bp extension each side
        # New: start=50, end=250
        assert "50\t250" in lines[0]


class TestComplement:
    """Tests for pygrit.complement function."""

    def test_complement_basic(self, temp_dir, genome_file):
        """Test basic complement."""
        bed_file = temp_dir / "test.bed"
        bed_file.write_text("chr1\t100\t200\nchr1\t300\t400\n")

        result = pygrit.complement(str(bed_file), str(genome_file))

        # Should have gaps: 0-100, 200-300, 400-1000
        assert "chr1\t0\t100" in result
        assert "chr1\t200\t300" in result
        assert "chr1\t400\t1000" in result

    def test_complement_to_file(self, temp_dir, genome_file):
        """Test complement with output to file."""
        bed_file = temp_dir / "test.bed"
        bed_file.write_text("chr1\t100\t200\n")
        output = temp_dir / "complement.bed"

        result = pygrit.complement(str(bed_file), str(genome_file), output=str(output))

        assert result is None
        assert output.exists()

    def test_complement_full_coverage(self, temp_dir, genome_file):
        """Test complement when chromosome is fully covered."""
        bed_file = temp_dir / "test.bed"
        bed_file.write_text("chr1\t0\t1000\n")

        result = pygrit.complement(str(bed_file), str(genome_file))

        # chr1 is fully covered, should not appear in complement
        # but chr2 and chr3 should appear
        assert "chr2\t0\t800" in result
        assert "chr3\t0\t600" in result


class TestGenomecov:
    """Tests for pygrit.genomecov function."""

    def test_genomecov_histogram(self, temp_dir, genome_file):
        """Test histogram output mode."""
        bed_file = temp_dir / "test.bed"
        bed_file.write_text("chr1\t100\t200\nchr1\t150\t250\n")

        result = pygrit.genomecov(str(bed_file), str(genome_file))

        # Should contain histogram data
        assert result is not None
        assert len(result) > 0

    def test_genomecov_bedgraph(self, temp_dir, genome_file):
        """Test BedGraph output mode."""
        bed_file = temp_dir / "test.bed"
        bed_file.write_text("chr1\t100\t200\n")

        result = pygrit.genomecov(str(bed_file), str(genome_file), bg=True)

        # Should be in BedGraph format (non-zero only)
        assert "chr1" in result

    def test_genomecov_to_file(self, temp_dir, genome_file):
        """Test genomecov with output to file."""
        bed_file = temp_dir / "test.bed"
        bed_file.write_text("chr1\t100\t200\n")
        output = temp_dir / "coverage.txt"

        result = pygrit.genomecov(str(bed_file), str(genome_file), output=str(output))

        assert result is None
        assert output.exists()


class TestJaccard:
    """Tests for pygrit.jaccard function."""

    def test_jaccard_identical(self, temp_dir):
        """Test Jaccard similarity of identical files."""
        bed_file = temp_dir / "test.bed"
        bed_file.write_text("chr1\t100\t200\nchr1\t300\t400\n")

        result = pygrit.jaccard(str(bed_file), str(bed_file))

        # Identical files should have Jaccard = 1.0
        assert "1" in result  # Jaccard index should be 1

    def test_jaccard_no_overlap(self, temp_dir):
        """Test Jaccard similarity with no overlap."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t100\t200\n")
        b_file.write_text("chr1\t300\t400\n")

        result = pygrit.jaccard(str(a_file), str(b_file))

        # No overlap should have Jaccard = 0
        lines = result.strip().split("\n")
        assert len(lines) >= 1

    def test_jaccard_to_file(self, temp_dir):
        """Test Jaccard with output to file."""
        bed_file = temp_dir / "test.bed"
        bed_file.write_text("chr1\t100\t200\n")
        output = temp_dir / "jaccard.txt"

        result = pygrit.jaccard(str(bed_file), str(bed_file), output=str(output))

        assert result is None
        assert output.exists()


class TestMultiinter:
    """Tests for pygrit.multiinter function."""

    def test_multiinter_basic(self, temp_dir):
        """Test basic multiinter with two files."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t100\t200\nchr1\t300\t400\n")
        b_file.write_text("chr1\t150\t250\nchr1\t350\t450\n")

        result = pygrit.multiinter([str(a_file), str(b_file)])

        assert result is not None
        assert "chr1" in result

    def test_multiinter_three_files(self, temp_dir):
        """Test multiinter with three files."""
        files = []
        for i, name in enumerate(["a", "b", "c"]):
            f = temp_dir / f"{name}.bed"
            f.write_text(f"chr1\t{100 + i*50}\t{200 + i*50}\n")
            files.append(str(f))

        result = pygrit.multiinter(files)

        assert result is not None

    def test_multiinter_cluster(self, temp_dir):
        """Test multiinter with cluster mode."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t100\t200\n")
        b_file.write_text("chr1\t100\t200\n")  # Same interval

        result = pygrit.multiinter([str(a_file), str(b_file)], cluster=True)

        assert result is not None

    def test_multiinter_to_file(self, temp_dir):
        """Test multiinter with output to file."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        a_file.write_text("chr1\t100\t200\n")
        b_file.write_text("chr1\t150\t250\n")
        output = temp_dir / "multiinter.bed"

        result = pygrit.multiinter([str(a_file), str(b_file)], output=str(output))

        assert result is None
        assert output.exists()

    def test_multiinter_requires_two_files(self):
        """Test that multiinter requires at least 2 files."""
        with pytest.raises(ValueError):
            pygrit.multiinter(["single_file.bed"])


class TestGenerate:
    """Tests for pygrit.generate function."""

    def test_generate_basic(self, temp_dir):
        """Test basic generation."""
        output_dir = temp_dir / "generated"

        stats = pygrit.generate(
            str(output_dir),
            num_intervals=1000,
            seed=42,
        )

        assert "files_generated" in stats
        assert "total_intervals" in stats
        assert "elapsed_seconds" in stats
        assert stats["total_intervals"] > 0

    def test_generate_modes(self, temp_dir):
        """Test different generation modes."""
        for mode in ["balanced", "clustered"]:
            output_dir = temp_dir / f"generated_{mode}"

            stats = pygrit.generate(
                str(output_dir),
                num_intervals=500,
                mode=mode,
                seed=42,
            )

            assert stats["total_intervals"] > 0

    def test_generate_reproducible(self, temp_dir):
        """Test that generation is reproducible with same seed."""
        dir1 = temp_dir / "gen1"
        dir2 = temp_dir / "gen2"

        stats1 = pygrit.generate(str(dir1), num_intervals=100, seed=12345)
        stats2 = pygrit.generate(str(dir2), num_intervals=100, seed=12345)

        assert stats1["total_intervals"] == stats2["total_intervals"]

    def test_generate_custom_lengths(self, temp_dir):
        """Test generation with custom interval lengths."""
        output_dir = temp_dir / "generated"

        stats = pygrit.generate(
            str(output_dir),
            num_intervals=500,
            len_min=100,
            len_max=500,
            seed=42,
        )

        assert stats["total_intervals"] > 0

    def test_generate_invalid_mode(self, temp_dir):
        """Test that invalid mode raises error."""
        with pytest.raises(ValueError):
            pygrit.generate(str(temp_dir), mode="invalid_mode")


class TestNewCommandsIntegration:
    """Integration tests for new commands."""

    def test_sort_then_merge(self, unsorted_bed, temp_dir):
        """Test sorting then merging."""
        # Sort
        sorted_output = temp_dir / "sorted.bed"
        pygrit.sort(str(unsorted_bed), output=str(sorted_output))

        # Merge
        merged = pygrit.merge(str(sorted_output))

        assert len(merged) > 0

    def test_slop_then_merge(self, temp_dir, genome_file):
        """Test slop then merge workflow."""
        bed_file = temp_dir / "test.bed"
        bed_file.write_text("chr1\t100\t110\nchr1\t120\t130\n")

        # Slop to make them overlap
        slopped = temp_dir / "slopped.bed"
        pygrit.slop(str(bed_file), str(genome_file), both=10.0, output=str(slopped))

        # Then merge
        merged = pygrit.merge(str(slopped))

        # Should be merged into one interval
        assert len(merged) == 1
