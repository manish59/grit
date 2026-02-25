//! Active set management for streaming operations.
//!
//! The ActiveSet maintains a collection of intervals that could potentially
//! overlap the current position, with automatic compaction to bound memory.

/// Compaction threshold - trigger when head_idx exceeds this value.
const COMPACTION_THRESHOLD: usize = 4096;

/// Active interval - stores only coordinates (8 bytes total).
///
/// The chromosome is tracked separately to avoid per-interval allocation.
#[derive(Debug, Clone, Copy)]
pub struct ActiveInterval {
    pub start: u32,
    pub end: u32,
}

impl ActiveInterval {
    /// Create a new active interval.
    #[inline]
    pub fn new(start: u64, end: u64) -> Self {
        Self {
            start: start as u32,
            end: end as u32,
        }
    }
}

/// Active set with automatic compaction.
///
/// Uses Vec + head_idx pattern for better cache locality than VecDeque.
/// Elements before head_idx are logically removed but not deallocated
/// until compaction is triggered.
///
/// # Memory Complexity
///
/// O(k) where k = max number of overlapping intervals at any position.
/// Periodic compaction ensures memory doesn't grow unbounded.
#[derive(Debug)]
pub struct ActiveSet<T> {
    /// Storage for active elements.
    data: Vec<T>,
    /// Index of the first logically active element.
    head_idx: usize,
    /// Maximum observed active size (for statistics).
    max_active: usize,
}

impl<T> Default for ActiveSet<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> ActiveSet<T> {
    /// Create a new empty active set.
    pub fn new() -> Self {
        Self::with_capacity(1024)
    }

    /// Create a new active set with specified initial capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            data: Vec::with_capacity(capacity),
            head_idx: 0,
            max_active: 0,
        }
    }

    /// Add an element to the active set.
    #[inline]
    pub fn push(&mut self, value: T) {
        self.data.push(value);
        self.update_max_active();
    }

    /// Get the number of logically active elements.
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len() - self.head_idx
    }

    /// Check if the active set is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.head_idx >= self.data.len()
    }

    /// Get a slice of all active elements.
    #[inline]
    pub fn as_slice(&self) -> &[T] {
        &self.data[self.head_idx..]
    }

    /// Get a mutable slice of all active elements.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data[self.head_idx..]
    }

    /// Advance the head index, logically removing the first element.
    #[inline]
    pub fn advance_head(&mut self) {
        if self.head_idx < self.data.len() {
            self.head_idx += 1;
        }
    }

    /// Advance head while condition is true for the front element.
    ///
    /// Returns the number of elements removed.
    #[inline]
    pub fn advance_while<F>(&mut self, mut condition: F) -> usize
    where
        F: FnMut(&T) -> bool,
    {
        let start_idx = self.head_idx;
        while self.head_idx < self.data.len() && condition(&self.data[self.head_idx]) {
            self.head_idx += 1;
        }
        self.head_idx - start_idx
    }

    /// Compact the internal storage if needed.
    ///
    /// This is called automatically but can be invoked manually.
    pub fn compact_if_needed(&mut self) {
        if self.head_idx > COMPACTION_THRESHOLD && self.head_idx * 2 > self.data.len() {
            self.data.drain(0..self.head_idx);
            self.head_idx = 0;
        }
    }

    /// Clear all elements and reset state.
    pub fn clear(&mut self) {
        self.data.clear();
        self.head_idx = 0;
    }

    /// Get the maximum active size observed (for statistics).
    pub fn max_active(&self) -> usize {
        self.max_active
    }

    /// Update max_active tracking.
    #[inline]
    fn update_max_active(&mut self) {
        let current = self.len();
        if current > self.max_active {
            self.max_active = current;
        }
    }

    /// Iterator over active elements.
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.data[self.head_idx..].iter()
    }

    /// Get reference to element at logical index.
    #[inline]
    pub fn get(&self, index: usize) -> Option<&T> {
        self.data.get(self.head_idx + index)
    }

    /// Get mutable reference to element at logical index.
    #[inline]
    pub fn get_mut(&mut self, index: usize) -> Option<&mut T> {
        self.data.get_mut(self.head_idx + index)
    }

    /// Get reference to first active element.
    #[inline]
    pub fn front(&self) -> Option<&T> {
        self.data.get(self.head_idx)
    }
}

/// Active set specialized for intervals with expiration by end position.
impl ActiveSet<ActiveInterval> {
    /// Remove expired intervals (those ending at or before the given position).
    ///
    /// Returns the number of intervals removed.
    #[inline]
    pub fn expire_before(&mut self, position: u64) -> usize {
        let count = self.advance_while(|b| (b.end as u64) <= position);
        self.compact_if_needed();
        count
    }

    /// Count intervals that overlap the given range [start, end).
    #[inline]
    pub fn count_overlapping(&self, start: u64, end: u64) -> usize {
        self.iter()
            .filter(|b| (b.end as u64) > start && (b.start as u64) < end)
            .count()
    }

    /// Iterate over intervals that overlap the given range [start, end).
    pub fn iter_overlapping(&self, start: u64, end: u64) -> impl Iterator<Item = &ActiveInterval> {
        self.iter()
            .filter(move |b| (b.end as u64) > start && (b.start as u64) < end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_active_interval_size() {
        assert_eq!(std::mem::size_of::<ActiveInterval>(), 8);
    }

    #[test]
    fn test_active_set_basic() {
        let mut set: ActiveSet<u32> = ActiveSet::new();
        set.push(1);
        set.push(2);
        set.push(3);

        assert_eq!(set.len(), 3);
        assert_eq!(set.as_slice(), &[1, 2, 3]);
    }

    #[test]
    fn test_active_set_advance() {
        let mut set: ActiveSet<u32> = ActiveSet::new();
        set.push(1);
        set.push(2);
        set.push(3);

        set.advance_head();
        assert_eq!(set.len(), 2);
        assert_eq!(set.as_slice(), &[2, 3]);
    }

    #[test]
    fn test_active_set_advance_while() {
        let mut set: ActiveSet<u32> = ActiveSet::new();
        for i in 0..10 {
            set.push(i);
        }

        let removed = set.advance_while(|&x| x < 5);
        assert_eq!(removed, 5);
        assert_eq!(set.len(), 5);
        assert_eq!(set.front(), Some(&5));
    }

    #[test]
    fn test_active_set_clear() {
        let mut set: ActiveSet<u32> = ActiveSet::new();
        set.push(1);
        set.push(2);

        set.clear();
        assert!(set.is_empty());
        assert_eq!(set.len(), 0);
    }

    #[test]
    fn test_interval_expire_before() {
        let mut set: ActiveSet<ActiveInterval> = ActiveSet::new();
        set.push(ActiveInterval::new(100, 200));
        set.push(ActiveInterval::new(150, 250));
        set.push(ActiveInterval::new(200, 300));

        let expired = set.expire_before(200);
        assert_eq!(expired, 1);
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn test_interval_count_overlapping() {
        let mut set: ActiveSet<ActiveInterval> = ActiveSet::new();
        set.push(ActiveInterval::new(100, 200));
        set.push(ActiveInterval::new(150, 250));
        set.push(ActiveInterval::new(300, 400));

        assert_eq!(set.count_overlapping(175, 225), 2);
        assert_eq!(set.count_overlapping(350, 450), 1);
        assert_eq!(set.count_overlapping(500, 600), 0);
    }
}
