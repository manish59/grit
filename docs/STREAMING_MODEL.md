# GRIT Streaming Model

This document describes GRIT's streaming algorithms and their memory guarantees.

## Overview

GRIT's streaming mode processes genomic intervals with **O(k) memory**, where k is the maximum number of intervals overlapping at any genomic position. This enables processing arbitrarily large files with constant memory usage.

## Memory Complexity

| Mode | Memory | Use Case |
|------|--------|----------|
| Parallel | O(n + m) | Maximum speed, unsorted input |
| Streaming | O(k) | Large files, memory-constrained |

Where:
- **n** = number of intervals in file A
- **m** = number of intervals in file B
- **k** = max overlapping intervals at any position (typically < 100)

## Sweep-Line Algorithm

The core streaming algorithm uses a sweep-line approach:

```
Input: Two sorted BED files A and B
Output: Intersection results

1. Initialize active_b = empty set (intervals from B currently "active")
2. For each interval a in A:
   a. Remove intervals from active_b that end before a.start
   b. Add intervals from B that start before a.end to active_b
   c. For each interval b in active_b overlapping a:
      - Output intersection
3. Output remaining active intervals (if needed)
```

### Visual Example

```
Position:  100       200       300       400       500
           |---------|---------|---------|---------|

A intervals:
           [========A1========]
                         [====A2====]

B intervals:
      [==B1==]
                [=====B2=====]
                              [===B3===]

Sweep line progression:
  pos=100: active_b = {B1}       → A1 overlaps B1
  pos=150: active_b = {B1, B2}   → A1 overlaps B1, B2
  pos=200: active_b = {B2}       → A1 overlaps B2
  pos=250: active_b = {B2}       → A2 overlaps B2
  pos=300: active_b = {B2, B3}   → A2 overlaps B2, B3
  pos=350: active_b = {B3}       → A2 overlaps B3
```

## Implementation

### Active Interval Buffer

```rust
struct ActiveBuffer {
    intervals: VecDeque<BedRecord>,
    current_chrom: String,
}

impl ActiveBuffer {
    fn add(&mut self, record: BedRecord) {
        self.intervals.push_back(record);
    }

    fn remove_expired(&mut self, position: u64) {
        while let Some(front) = self.intervals.front() {
            if front.end <= position {
                self.intervals.pop_front();
            } else {
                break;
            }
        }
    }

    fn overlapping(&self, start: u64, end: u64) -> impl Iterator<Item = &BedRecord> {
        self.intervals.iter()
            .filter(move |b| b.start < end && b.end > start)
    }
}
```

### Memory Bound Proof

**Theorem**: The active buffer size is bounded by k, where k is the maximum number of intervals overlapping at any position.

**Proof**:
1. An interval b enters the buffer when the sweep line reaches b.start
2. An interval b exits the buffer when the sweep line passes b.end
3. At any position p, the buffer contains exactly the intervals covering p
4. By definition, at most k intervals cover any position
5. Therefore, buffer size ≤ k

**In practice**: For typical genomic data, k < 100 even with millions of intervals.

## Sorted Input Requirement

Streaming algorithms require sorted input:

```
chr1    100    200    # ✓
chr1    150    250    # ✓
chr1    300    400    # ✓
chr2    100    200    # ✓ (new chromosome)
chr1    500    600    # ✗ ERROR: chr1 after chr2
```

### Validation

GRIT validates sort order by default:

```rust
fn validate_sorted(&mut self, record: &BedRecord) -> Result<()> {
    if let Some(ref last_chrom) = self.last_chrom {
        if record.chrom == *last_chrom {
            if record.start < self.last_start {
                return Err(GritError::UnsortedInput(...));
            }
        } else if self.seen_chroms.contains(&record.chrom) {
            return Err(GritError::ChromosomeRevisited(...));
        }
    }
    Ok(())
}
```

### Skip Validation

For pre-sorted files, skip validation for faster startup:

```bash
grit intersect -a sorted_a.bed -b sorted_b.bed --streaming --assume-sorted
```

## Chromosome Transitions

Special handling is required at chromosome boundaries:

```rust
fn process_chromosome_change(&mut self, new_chrom: &str) {
    // Flush all remaining output for current chromosome
    self.flush_current_chromosome();

    // Clear active buffer (intervals from old chromosome)
    self.active_buffer.clear();

    // Update current chromosome
    self.current_chrom = new_chrom.to_string();
}
```

## Commands Supporting Streaming

| Command | `--streaming` flag | Memory |
|---------|-------------------|--------|
| intersect | Yes | O(k) |
| subtract | Yes | O(k) |
| closest | Yes | O(k) |
| window | Implicit | O(k) |
| coverage | Implicit | O(k) |
| merge | Always | O(1) |
| complement | Always | O(1) |

## Performance Characteristics

### Time Complexity

- **Streaming intersection**: O(n + m + output)
- **Parallel intersection**: O((n + m) log(n + m) + output)

### I/O Pattern

Streaming mode reads input sequentially, enabling:
- Efficient disk access patterns
- Pipeline processing (stdin/stdout)
- Memory-mapped I/O optimization

### When to Use Streaming

Use streaming mode when:
- Files exceed available RAM
- Processing pipelines (stdin/stdout)
- Memory is constrained
- Files are already sorted

Use parallel mode when:
- Maximum speed is required
- Files fit in memory
- Input is unsorted (will be sorted in memory)

## Testing Streaming Correctness

### Memory Verification

```bash
# Linux: Check peak RSS
/usr/bin/time -v grit intersect --streaming -a large.bed -b large.bed 2>&1 | grep "Maximum resident"

# macOS: Use Activity Monitor or instruments
```

### Correctness Verification

```bash
# Compare streaming vs parallel output
grit intersect -a A.bed -b B.bed > parallel.bed
grit intersect -a A.bed -b B.bed --streaming > streaming.bed
diff parallel.bed streaming.bed  # Should be identical
```

### Edge Cases

Test with:
- Empty files
- Single interval
- Non-overlapping intervals
- Highly nested intervals (stress k bound)
- Chromosome transitions
