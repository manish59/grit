# Contributing to GRIT

Thank you for your interest in contributing to GRIT (Genomic Range Interval Toolkit)!

## Table of Contents

- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Code Style](#code-style)
- [Testing Requirements](#testing-requirements)
- [Streaming Correctness](#streaming-correctness)
- [Performance Guidelines](#performance-guidelines)
- [Float Precision Policy](#float-precision-policy)
- [Pull Request Process](#pull-request-process)
- [Commit Conventions](#commit-conventions)

## Getting Started

### Prerequisites

- **Rust 1.82+** (stable) - Check with `rustc --version`
- **bedtools** - For benchmark comparisons and SHA256 parity testing
- **Git** - For version control

### Development Setup

```bash
# Clone the repository
git clone https://github.com/manish59/grit.git
cd grit

# Build in debug mode (fast compilation)
cargo build

# Build in release mode (optimized, for benchmarking)
cargo build --release

# Run all tests
cargo test --all

# Run linting
cargo clippy --all-targets -- -D warnings

# Check formatting
cargo fmt --all -- --check
```

### IDE Setup

**VS Code:**
- Install `rust-analyzer` extension
- Enable format on save

**CLion/IntelliJ:**
- Install Rust plugin
- Enable rustfmt integration

## Development Workflow

### Branch Naming

Use descriptive branch names with prefixes:

| Prefix | Use Case | Example |
|--------|----------|---------|
| `feature/` | New features | `feature/add-flank-command` |
| `fix/` | Bug fixes | `fix/streaming-memory-leak` |
| `perf/` | Performance improvements | `perf/optimize-intersection` |
| `docs/` | Documentation | `docs/streaming-guide` |
| `refactor/` | Code refactoring | `refactor/extract-parser` |
| `ci/` | CI/tooling changes | `ci/add-coverage` |

### Development Cycle

1. **Create a branch** from `main`
2. **Make changes** following code style guidelines
3. **Write tests** for new functionality
4. **Run checks locally**:
   ```bash
   cargo fmt --all
   cargo clippy --all-targets -- -D warnings
   cargo test --all
   ```
5. **Verify bedtools parity** (if applicable)
6. **Submit PR** with clear description

## Code Style

### Rust Guidelines

- Follow [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)
- Use `rustfmt` defaults (no custom configuration)
- Fix all `clippy` warnings
- Prefer explicit types over inference for public APIs
- Document all public items with `///` comments

### Naming Conventions

```rust
// Types: PascalCase
struct IntervalTree { ... }

// Functions/methods: snake_case
fn merge_intervals(intervals: &[Interval]) -> Vec<Interval> { ... }

// Constants: SCREAMING_SNAKE_CASE
const DEFAULT_WINDOW_SIZE: u64 = 1000;

// Modules: snake_case
mod streaming_merge;
```

### Error Handling

- Use `thiserror` for library errors
- Provide context in error messages
- Don't panic in library code; return `Result`

```rust
#[derive(Error, Debug)]
pub enum GritError {
    #[error("file not sorted: {0} comes after {1} on {2}")]
    UnsortedInput(u64, u64, String),
}
```

## Testing Requirements

### Test Checklist

- [ ] Unit tests for new functions
- [ ] Integration tests for new commands
- [ ] Edge cases (empty files, single interval, large files)
- [ ] SHA256 parity with bedtools (where applicable)
- [ ] Streaming mode correctness
- [ ] Memory bounds verification

### Running Tests

```bash
# All tests
cargo test --all

# Specific module
cargo test streaming

# With output
cargo test -- --nocapture

# Release mode tests (for timing-sensitive tests)
cargo test --release
```

### SHA256 Parity Testing

For commands that should match bedtools output:

```bash
# Generate test data
./benchmarks/bench.sh data 100K 50K

# Run parity check
diff <(bedtools intersect -a A.bed -b B.bed | sort) \
     <(grit intersect -a A.bed -b B.bed | sort)

# Or use SHA256
sha256sum <(bedtools intersect -a A.bed -b B.bed | sort)
sha256sum <(grit intersect -a A.bed -b B.bed | sort)
```

## Streaming Correctness

GRIT's streaming algorithms use **O(k) memory** where k is the maximum number of overlapping intervals at any position.

### Requirements for Streaming Code

1. **Memory Bound**: Never buffer entire files; maintain O(k) memory
2. **Single Pass**: Process input sequentially without random access
3. **Sorted Input**: Require and validate sorted input
4. **Deterministic Output**: Same input always produces same output

### Validation Checklist

- [ ] Memory usage stays constant regardless of file size
- [ ] Works with stdin pipes (`cat file | grit merge -i -`)
- [ ] Validates sort order by default
- [ ] `--assume-sorted` skips validation correctly
- [ ] Handles chromosome transitions correctly

### Testing Streaming Mode

```bash
# Verify memory usage is constant
/usr/bin/time -v grit intersect --streaming -a large.bed -b large.bed

# Test with pipes
cat sorted.bed | grit merge -i - --assume-sorted

# Test sort validation
grit merge -i unsorted.bed  # Should fail with helpful error
```

## Performance Guidelines

### General Principles

1. **Zero-allocation hot paths**: Avoid heap allocations in inner loops
2. **Cache-friendly access**: Process data sequentially when possible
3. **Minimize copies**: Use references and borrows
4. **Profile before optimizing**: Use `cargo flamegraph` or `perf`

### Benchmarking Changes

```bash
# Generate benchmark data
./benchmarks/bench.sh data 10M 5M

# Run specific benchmark
./benchmarks/bench.sh run 10M 5M intersect

# Compare with bedtools
time bedtools intersect -a A.bed -b B.bed > /dev/null
time grit intersect -a A.bed -b B.bed > /dev/null
```

### Document Complexity

For new algorithms, document:
- Time complexity: O(?)
- Space complexity: O(?)
- I/O pattern: sequential/random

## Float Precision Policy

GRIT uses specific float formatting for bedtools compatibility:

### Rules

1. **Use `ryu` for float formatting**: Fast and deterministic
2. **Match bedtools output precision**: Typically 6 decimal places for fractions
3. **Avoid floating-point comparisons**: Use integer arithmetic where possible
4. **Document precision requirements**: Note when exact match is required

### Implementation

```rust
use ryu::Buffer;

fn format_fraction(value: f64) -> String {
    let mut buffer = Buffer::new();
    buffer.format(value).to_string()
}
```

## Pull Request Process

### Before Submitting

1. **Rebase on main**: `git fetch origin && git rebase origin/main`
2. **Run all checks**:
   ```bash
   cargo fmt --all
   cargo clippy --all-targets -- -D warnings
   cargo test --all
   ```
3. **Update documentation** if adding new features
4. **Update CHANGELOG.md** for user-facing changes

### PR Description

Include:
- Clear description of changes
- Type of change (bug fix, feature, etc.)
- Testing performed
- Performance impact (if applicable)
- Link to related issues

### Review Process

1. CI must pass
2. At least one maintainer review
3. All discussions resolved
4. Squash merge preferred for clean history

## Commit Conventions

Use [Conventional Commits](https://www.conventionalcommits.org/):

### Format

```
<type>: <description>

[optional body]

[optional footer]
```

### Types

| Type | Description |
|------|-------------|
| `feat` | New feature |
| `fix` | Bug fix |
| `perf` | Performance improvement |
| `docs` | Documentation only |
| `style` | Formatting, no code change |
| `refactor` | Code change that neither fixes bug nor adds feature |
| `test` | Adding or updating tests |
| `chore` | Maintenance tasks |
| `ci` | CI/CD changes |

### Examples

```
feat: add streaming mode for window command

- Implement O(k) memory streaming algorithm
- Add --assume-sorted flag for pre-sorted inputs
- Verify SHA256 parity with bedtools

Closes #42
```

```
fix: correct off-by-one error in complement

The complement command was excluding the last base of each
chromosome. Fixed by using inclusive end coordinates.

Fixes #55
```

```
perf: optimize intersection with SIMD parsing

- Use memchr for field detection
- Reduce allocations in hot loop
- 1.3x speedup on 10M interval benchmark
```

## Questions?

- **General questions**: Open a [Discussion](https://github.com/manish59/grit/discussions)
- **Bug reports**: Use [Issue templates](https://github.com/manish59/grit/issues/new/choose)
- **Feature ideas**: Use [Feature request template](https://github.com/manish59/grit/issues/new/choose)

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
