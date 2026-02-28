"""Pytest configuration and fixtures for pygrit tests."""

import os
import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_bed_a(temp_dir):
    """Create a sample BED file A with sorted intervals."""
    path = temp_dir / "a.bed"
    content = """chr1\t100\t200
chr1\t300\t400
chr1\t500\t600
chr2\t100\t200
chr2\t300\t400
"""
    path.write_text(content)
    return path


@pytest.fixture
def sample_bed_b(temp_dir):
    """Create a sample BED file B with sorted intervals."""
    path = temp_dir / "b.bed"
    content = """chr1\t150\t250
chr1\t550\t650
chr2\t50\t150
chr2\t350\t450
"""
    path.write_text(content)
    return path


@pytest.fixture
def overlapping_bed(temp_dir):
    """Create a BED file with overlapping intervals for merge testing."""
    path = temp_dir / "overlapping.bed"
    content = """chr1\t100\t200
chr1\t150\t250
chr1\t200\t300
chr1\t500\t600
"""
    path.write_text(content)
    return path


@pytest.fixture
def stranded_bed(temp_dir):
    """Create a BED file with strand information."""
    path = temp_dir / "stranded.bed"
    content = """chr1\t100\t200\tgene1\t0\t+
chr1\t150\t250\tgene2\t0\t+
chr1\t180\t280\tgene3\t0\t-
chr1\t220\t320\tgene4\t0\t-
"""
    path.write_text(content)
    return path


@pytest.fixture
def empty_bed(temp_dir):
    """Create an empty BED file."""
    path = temp_dir / "empty.bed"
    path.write_text("")
    return path


@pytest.fixture
def large_bed_a(temp_dir):
    """Create a larger BED file for performance testing."""
    path = temp_dir / "large_a.bed"
    lines = []
    for chrom in range(1, 3):
        for i in range(1000):
            start = i * 1000
            end = start + 500
            lines.append(f"chr{chrom}\t{start}\t{end}")
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def large_bed_b(temp_dir):
    """Create a larger BED file B for performance testing."""
    path = temp_dir / "large_b.bed"
    lines = []
    for chrom in range(1, 3):
        for i in range(500):
            start = i * 2000 + 250
            end = start + 500
            lines.append(f"chr{chrom}\t{start}\t{end}")
    path.write_text("\n".join(lines) + "\n")
    return path
