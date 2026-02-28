"""Integration tests comparing pygrit output with grit CLI."""

import subprocess
import pytest
import pygrit
from pathlib import Path


def run_grit_cli(args: list[str]) -> str:
    """Run grit CLI command and return output."""
    result = subprocess.run(
        ["grit"] + args,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"grit failed: {result.stderr}")
    return result.stdout


def grit_available() -> bool:
    """Check if grit CLI is available."""
    try:
        subprocess.run(["grit", "--version"], capture_output=True)
        return True
    except FileNotFoundError:
        return False


# Skip all tests if grit CLI not available
pytestmark = pytest.mark.skipif(
    not grit_available(),
    reason="grit CLI not available"
)


class TestIntersectParity:
    """Test pygrit.intersect produces same results as grit intersect."""

    def test_intersect_parity(self, sample_bed_a, sample_bed_b, temp_dir):
        """Compare intersect output between Python and CLI."""
        # Run with Python
        py_output = temp_dir / "py_intersect.bed"
        pygrit.intersect(
            str(sample_bed_a),
            str(sample_bed_b),
            output=str(py_output)
        )

        # Run with CLI
        cli_output = run_grit_cli([
            "intersect",
            "--streaming",
            "--assume-sorted",
            "-a", str(sample_bed_a),
            "-b", str(sample_bed_b),
        ])

        # Compare
        py_content = py_output.read_text()
        assert py_content == cli_output, f"Output mismatch:\nPython:\n{py_content}\nCLI:\n{cli_output}"

    def test_intersect_wa_wb_parity(self, sample_bed_a, sample_bed_b, temp_dir):
        """Compare intersect -wa -wb output."""
        py_output = temp_dir / "py_intersect.bed"
        pygrit.intersect(
            str(sample_bed_a),
            str(sample_bed_b),
            output=str(py_output),
            write_a=True,
            write_b=True,
        )

        cli_output = run_grit_cli([
            "intersect",
            "--streaming",
            "--assume-sorted",
            "-a", str(sample_bed_a),
            "-b", str(sample_bed_b),
            "--wa", "--wb",
        ])

        py_content = py_output.read_text()
        assert py_content == cli_output

    def test_intersect_unique_parity(self, sample_bed_a, sample_bed_b, temp_dir):
        """Compare intersect -u output."""
        py_output = temp_dir / "py_intersect.bed"
        pygrit.intersect(
            str(sample_bed_a),
            str(sample_bed_b),
            output=str(py_output),
            unique=True,
        )

        cli_output = run_grit_cli([
            "intersect",
            "--streaming",
            "--assume-sorted",
            "-a", str(sample_bed_a),
            "-b", str(sample_bed_b),
            "-u",
        ])

        py_content = py_output.read_text()
        assert py_content == cli_output

    def test_intersect_no_overlap_parity(self, sample_bed_a, sample_bed_b, temp_dir):
        """Compare intersect -v output."""
        py_output = temp_dir / "py_intersect.bed"
        pygrit.intersect(
            str(sample_bed_a),
            str(sample_bed_b),
            output=str(py_output),
            no_overlap=True,
        )

        cli_output = run_grit_cli([
            "intersect",
            "--streaming",
            "--assume-sorted",
            "-a", str(sample_bed_a),
            "-b", str(sample_bed_b),
            "-v",
        ])

        py_content = py_output.read_text()
        assert py_content == cli_output


class TestMergeParity:
    """Test pygrit.merge produces same results as grit merge."""

    def test_merge_parity(self, overlapping_bed, temp_dir):
        """Compare merge output."""
        py_output = temp_dir / "py_merge.bed"
        pygrit.merge(str(overlapping_bed), output=str(py_output))

        cli_output = run_grit_cli([
            "merge",
            "--assume-sorted",
            "-i", str(overlapping_bed),
        ])

        py_content = py_output.read_text()
        assert py_content == cli_output

    def test_merge_distance_parity(self, temp_dir):
        """Compare merge with distance."""
        bed_file = temp_dir / "gaps.bed"
        bed_file.write_text("chr1\t100\t200\nchr1\t250\t350\n")

        py_output = temp_dir / "py_merge.bed"
        pygrit.merge(str(bed_file), output=str(py_output), distance=100)

        cli_output = run_grit_cli([
            "merge",
            "--assume-sorted",
            "-i", str(bed_file),
            "-d", "100",
        ])

        py_content = py_output.read_text()
        assert py_content == cli_output


class TestSubtractParity:
    """Test pygrit.subtract produces same results as grit subtract."""

    def test_subtract_parity(self, sample_bed_a, sample_bed_b, temp_dir):
        """Compare subtract output."""
        py_output = temp_dir / "py_subtract.bed"
        pygrit.subtract(
            str(sample_bed_a),
            str(sample_bed_b),
            output=str(py_output)
        )

        cli_output = run_grit_cli([
            "subtract",
            "--streaming",
            "--assume-sorted",
            "-a", str(sample_bed_a),
            "-b", str(sample_bed_b),
        ])

        py_content = py_output.read_text()
        assert py_content == cli_output

    def test_subtract_remove_entire_parity(self, sample_bed_a, sample_bed_b, temp_dir):
        """Compare subtract -A output."""
        py_output = temp_dir / "py_subtract.bed"
        pygrit.subtract(
            str(sample_bed_a),
            str(sample_bed_b),
            output=str(py_output),
            remove_entire=True,
        )

        cli_output = run_grit_cli([
            "subtract",
            "--streaming",
            "--assume-sorted",
            "-a", str(sample_bed_a),
            "-b", str(sample_bed_b),
            "-A",
        ])

        py_content = py_output.read_text()
        assert py_content == cli_output


class TestClosestParity:
    """Test pygrit.closest produces same results as grit closest."""

    def test_closest_parity(self, sample_bed_a, sample_bed_b, temp_dir):
        """Compare closest output."""
        py_output = temp_dir / "py_closest.bed"
        pygrit.closest(
            str(sample_bed_a),
            str(sample_bed_b),
            output=str(py_output)
        )

        cli_output = run_grit_cli([
            "closest",
            "--streaming",
            "--assume-sorted",
            "-a", str(sample_bed_a),
            "-b", str(sample_bed_b),
        ])

        py_content = py_output.read_text()
        assert py_content == cli_output


class TestCoverageParity:
    """Test pygrit.coverage produces same results as grit coverage."""

    def test_coverage_parity(self, sample_bed_a, sample_bed_b, temp_dir):
        """Compare coverage output."""
        py_output = temp_dir / "py_coverage.bed"
        pygrit.coverage(
            str(sample_bed_a),
            str(sample_bed_b),
            output=str(py_output)
        )

        cli_output = run_grit_cli([
            "coverage",
            "--assume-sorted",
            "-a", str(sample_bed_a),
            "-b", str(sample_bed_b),
        ])

        py_content = py_output.read_text()
        assert py_content == cli_output


class TestLargeFileParity:
    """Test parity with larger files."""

    def test_intersect_large_parity(self, large_bed_a, large_bed_b, temp_dir):
        """Compare intersect on larger files."""
        py_output = temp_dir / "py_intersect.bed"
        pygrit.intersect(
            str(large_bed_a),
            str(large_bed_b),
            output=str(py_output)
        )

        cli_output = run_grit_cli([
            "intersect",
            "--streaming",
            "--assume-sorted",
            "-a", str(large_bed_a),
            "-b", str(large_bed_b),
        ])

        py_content = py_output.read_text()
        assert py_content == cli_output

    def test_merge_large_parity(self, large_bed_a, temp_dir):
        """Compare merge on larger file."""
        py_output = temp_dir / "py_merge.bed"
        pygrit.merge(str(large_bed_a), output=str(py_output))

        cli_output = run_grit_cli([
            "merge",
            "--assume-sorted",
            "-i", str(large_bed_a),
        ])

        py_content = py_output.read_text()
        assert py_content == cli_output
