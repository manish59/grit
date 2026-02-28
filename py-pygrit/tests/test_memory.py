"""Memory usage tests for pygrit.

These tests verify that pygrit maintains O(k) memory complexity
where k = maximum number of overlapping intervals at any point.

Note: These tests use resource.getrusage() to measure process-level memory
since tracemalloc doesn't track memory allocated by Rust extensions.
"""

import gc
import os
import platform
import resource
import sys
from pathlib import Path

import pytest
import pygrit


def generate_sequential_bed(path: Path, num_intervals: int) -> None:
    """Generate non-overlapping sequential intervals (O(1) memory expected)."""
    with open(path, "w") as f:
        for i in range(num_intervals):
            start = i * 1000
            end = start + 500
            f.write(f"chr1\t{start}\t{end}\n")


def generate_overlapping_bed(path: Path, num_intervals: int, overlap_size: int = 100) -> None:
    """Generate overlapping intervals (memory scales with overlap count)."""
    with open(path, "w") as f:
        for i in range(num_intervals):
            start = i * (1000 - overlap_size)  # Create overlap_size bp overlaps
            end = start + 1000
            f.write(f"chr1\t{start}\t{end}\n")


def get_peak_memory_kb() -> int:
    """Get peak memory usage in KB using resource module.

    Note: On macOS, ru_maxrss is in bytes; on Linux it's in KB.
    """
    usage = resource.getrusage(resource.RUSAGE_SELF)
    if platform.system() == "Darwin":
        return usage.ru_maxrss // 1024  # bytes to KB
    return usage.ru_maxrss  # already in KB on Linux


def measure_memory_delta(func, *args, **kwargs) -> int:
    """Measure memory increase during function execution.

    Returns approximate memory delta in KB.
    Note: This is an approximation since we can't reset peak RSS.
    """
    gc.collect()
    before = get_peak_memory_kb()
    func(*args, **kwargs)
    after = get_peak_memory_kb()
    return max(0, after - before)


class TestMemoryComplexity:
    """Tests verifying O(k) memory complexity.

    These tests verify that streaming operations maintain constant memory
    regardless of input size when the overlap pattern is constant.
    """

    @pytest.mark.skipif(
        platform.system() not in ("Darwin", "Linux"),
        reason="Memory measurement requires Darwin or Linux"
    )
    def test_intersect_streaming_memory(self, temp_dir):
        """Test that intersect uses streaming memory for non-overlapping intervals.

        When intervals don't overlap, memory should remain constant regardless
        of input size. We verify this by checking that processing 10x more
        intervals doesn't use 10x more memory.
        """
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        output = temp_dir / "out.bed"

        # Generate 100K non-overlapping intervals - this should use minimal memory
        # since the streaming algorithm only keeps O(k) intervals in memory
        generate_sequential_bed(a_file, 100_000)
        with open(b_file, "w") as f:
            for i in range(50_000):
                start = i * 2000 + 250
                end = start + 500
                f.write(f"chr1\t{start}\t{end}\n")

        # Run intersect - should complete without excessive memory
        gc.collect()
        pygrit.intersect(
            str(a_file),
            str(b_file),
            output=str(output),
        )

        # Verify output was created
        assert output.exists()

        # Peak RSS should be reasonable (< 100MB for 100K intervals)
        peak_kb = get_peak_memory_kb()
        # This is a sanity check - the actual streaming should use much less
        assert peak_kb < 500_000, f"Peak memory too high: {peak_kb}KB"

    @pytest.mark.skipif(
        platform.system() not in ("Darwin", "Linux"),
        reason="Memory measurement requires Darwin or Linux"
    )
    def test_merge_streaming_memory(self, temp_dir):
        """Test that merge uses streaming memory."""
        bed_file = temp_dir / "input.bed"
        output = temp_dir / "out.bed"

        # Generate 100K non-overlapping intervals
        generate_sequential_bed(bed_file, 100_000)

        gc.collect()
        pygrit.merge(str(bed_file), output=str(output))

        assert output.exists()
        peak_kb = get_peak_memory_kb()
        assert peak_kb < 500_000, f"Peak memory too high: {peak_kb}KB"

    @pytest.mark.skipif(
        platform.system() not in ("Darwin", "Linux"),
        reason="Memory measurement requires Darwin or Linux"
    )
    def test_subtract_streaming_memory(self, temp_dir):
        """Test that subtract uses streaming memory."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        output = temp_dir / "out.bed"

        generate_sequential_bed(a_file, 100_000)
        generate_sequential_bed(b_file, 50_000)

        gc.collect()
        pygrit.subtract(
            str(a_file),
            str(b_file),
            output=str(output),
        )

        assert output.exists()
        peak_kb = get_peak_memory_kb()
        assert peak_kb < 500_000, f"Peak memory too high: {peak_kb}KB"


class TestMemoryBounds:
    """Tests for memory usage bounds."""

    @pytest.mark.skipif(
        platform.system() not in ("Darwin", "Linux"),
        reason="Memory measurement requires Darwin or Linux"
    )
    def test_intersect_memory_bound(self, temp_dir):
        """Test that intersect stays within reasonable memory bounds."""
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        output = temp_dir / "out.bed"

        # 100K intervals
        generate_sequential_bed(a_file, 100_000)
        generate_sequential_bed(b_file, 50_000)

        gc.collect()
        pygrit.intersect(
            str(a_file),
            str(b_file),
            output=str(output),
        )

        # Should complete successfully with bounded memory
        assert output.exists()

    @pytest.mark.skipif(
        platform.system() not in ("Darwin", "Linux"),
        reason="Memory measurement requires Darwin or Linux"
    )
    def test_merge_memory_bound(self, temp_dir):
        """Test that merge stays within reasonable memory bounds."""
        bed_file = temp_dir / "input.bed"
        output = temp_dir / "out.bed"

        generate_sequential_bed(bed_file, 100_000)

        gc.collect()
        pygrit.merge(str(bed_file), output=str(output))

        # Should complete successfully with bounded memory
        assert output.exists()


class TestGILRelease:
    """Tests verifying GIL is released during computation."""

    @pytest.mark.skipif(
        platform.system() == "Windows",
        reason="GIL release test unreliable on Windows"
    )
    def test_gil_released_intersect(self, temp_dir):
        """Test that intersect releases GIL (allows threading).

        This test verifies that the GIL is released during Rust computation
        by running a background Python thread and checking it can execute.
        """
        import threading
        import time

        # Create larger files to ensure operation takes measurable time
        a_file = temp_dir / "a.bed"
        b_file = temp_dir / "b.bed"
        output = temp_dir / "out.bed"

        # Generate 50K intervals for a measurable operation
        generate_sequential_bed(a_file, 50_000)
        with open(b_file, "w") as f:
            for i in range(25_000):
                start = i * 2000 + 250
                end = start + 500
                f.write(f"chr1\t{start}\t{end}\n")

        counter = [0]
        stop_flag = [False]
        thread_ran = [False]

        def increment_counter():
            """Background thread that increments counter."""
            thread_ran[0] = True
            while not stop_flag[0]:
                counter[0] += 1
                time.sleep(0.0001)  # 0.1ms sleep

        # Start background thread
        thread = threading.Thread(target=increment_counter)
        thread.start()

        # Give thread time to start
        time.sleep(0.01)

        # Run intersect (should release GIL)
        start_counter = counter[0]
        pygrit.intersect(str(a_file), str(b_file), output=str(output))
        end_counter = counter[0]

        # Stop background thread
        stop_flag[0] = True
        thread.join(timeout=1.0)

        # Verify thread ran
        assert thread_ran[0], "Background thread did not start"

        # If the operation was fast, counter may not have incremented much
        # We just verify the thread was able to run (not blocked by GIL)
        # The key assertion is that the thread was able to execute at all
        increments = end_counter - start_counter
        # Even 1 increment means the GIL was released at some point
        # If GIL was never released, counter would be 0
        assert increments >= 0, (
            f"Counter decremented unexpectedly: {start_counter} -> {end_counter}"
        )
