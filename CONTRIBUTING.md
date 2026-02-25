# Contributing to GRIT

Thank you for your interest in contributing to GRIT (Genomic Range Interval Toolkit)!

## Getting Started

### Development Setup

```bash
# Clone the repository
git clone https://github.com/manish59/grit.git
cd grit

# Build in release mode
cargo build --release

# Run tests
cargo test

# Run clippy
cargo clippy
```

### Prerequisites

- Rust 1.70+ (stable)
- For benchmarking: `bedtools` installed

## How to Contribute

### Reporting Bugs

1. Check [existing issues](https://github.com/manish59/grit/issues) to avoid duplicates
2. Use the bug report template
3. Include:
   - GRIT version (`grit --version`)
   - Operating system
   - Minimal reproduction steps
   - Input files (if small) or description of data characteristics

### Suggesting Features

1. Check if the feature aligns with GRIT's goals (fast, memory-efficient genomic operations)
2. Use the feature request template
3. Explain the use case and expected behavior

### Pull Requests

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Make your changes following the code style guidelines below
4. Add tests for new functionality
5. Ensure all checks pass:
   ```bash
   cargo test
   cargo clippy
   cargo fmt --check
   ```
6. Submit a PR with a clear description

## Code Style

### General Guidelines

- Follow Rust idioms and conventions
- Run `cargo fmt` before committing
- Run `cargo clippy` and fix all warnings
- Add documentation for public APIs
- Write tests for new features

### Performance Considerations

GRIT prioritizes performance. When contributing:

- Prefer zero-allocation algorithms where possible
- Use streaming/O(k) memory patterns for large file operations
- Benchmark changes against bedtools for comparison
- Document time and space complexity

### Testing Requirements

- Unit tests for new functions
- Integration tests for new commands
- SHA256 parity tests with bedtools where applicable
- Run the full test suite: `cargo test`

## Commit Messages

Use conventional commits:

- `feat:` new features
- `fix:` bug fixes
- `perf:` performance improvements
- `docs:` documentation changes
- `chore:` maintenance tasks
- `refactor:` code refactoring
- `test:` test additions/changes

Example:
```
feat: add streaming mode for window command

- Implement O(k) memory streaming algorithm
- Add --assume-sorted flag for pre-sorted inputs
- Verify SHA256 parity with bedtools
```

## Running Benchmarks

```bash
# Generate test data and run benchmarks
cd benchmarks
./bench.sh data 10M 5M
./bench.sh run 10M 5M all
```

## Questions?

Open an issue for discussion or questions about contributing.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
