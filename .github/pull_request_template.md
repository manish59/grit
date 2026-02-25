## Description

<!-- Brief description of the changes in this PR -->

## Type of Change

- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds functionality)
- [ ] Performance improvement
- [ ] Breaking change (fix or feature that would cause existing functionality to change)
- [ ] Documentation update
- [ ] Refactoring (no functional changes)
- [ ] CI/tooling improvement

## Testing

- [ ] Added unit tests for new functionality
- [ ] Added integration tests (if applicable)
- [ ] Verified SHA256 parity with bedtools (if applicable)
- [ ] Tested with large files (if performance-related)
- [ ] Tested streaming mode (if applicable)

## Checklist

- [ ] Code follows project style guidelines
- [ ] `cargo fmt --all` produces no changes
- [ ] `cargo clippy --all-targets -- -D warnings` produces no warnings
- [ ] `cargo test --all` passes
- [ ] Documentation updated (if applicable)
- [ ] CHANGELOG.md updated (for user-facing changes)

## Performance Impact

<!-- If this PR affects performance, describe: -->
<!-- - Benchmarks before/after -->
<!-- - Memory usage changes -->
<!-- - Time complexity changes -->

N/A or describe impact...

## Streaming Correctness

<!-- If this PR affects streaming algorithms: -->
<!-- - Verified O(k) memory bound maintained -->
<!-- - Tested with sorted input validation -->
<!-- - Checked edge cases (empty files, single interval, etc.) -->

N/A or describe verification...

## Related Issues

<!-- Link related issues: Closes #XXX, Fixes #XXX, Related to #XXX -->
