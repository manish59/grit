---
layout: default
title: Release Process
nav_exclude: true
---

# GRIT Release Process

This document describes the release process for GRIT.

## Version Numbering

GRIT follows [Semantic Versioning](https://semver.org/):

- **MAJOR.MINOR.PATCH** (e.g., `1.2.3`)
- **MAJOR**: Breaking API/CLI changes
- **MINOR**: New features, backward compatible
- **PATCH**: Bug fixes, backward compatible

Pre-release versions: `1.0.0-alpha.1`, `1.0.0-beta.1`, `1.0.0-rc.1`

## Release Checklist

### 1. Pre-Release Verification

```bash
# Ensure all tests pass
cargo test --all

# Run clippy with no warnings
cargo clippy --all-targets -- -D warnings

# Check formatting
cargo fmt --all -- --check

# Run benchmarks (optional)
./benchmarks/bench.sh run 1M 500K all

# Verify SHA256 parity with bedtools
./benchmarks/bench.sh truth 1M 500K all
```

### 2. Update Version

Update version in `Cargo.toml`:

```toml
[package]
version = "X.Y.Z"
```

### 3. Update CHANGELOG

Add release notes to `CHANGELOG.md`:

```markdown
## [X.Y.Z] - YYYY-MM-DD

### Added
- New feature description

### Changed
- Changed behavior description

### Fixed
- Bug fix description

### Performance
- Performance improvement description
```

### 4. Update Documentation

- Verify README is accurate
- Update benchmark numbers if changed
- Check all documentation links

### 5. Create Release Commit

```bash
git add Cargo.toml CHANGELOG.md
git commit -m "chore: release vX.Y.Z"
```

### 6. Create and Push Tag

```bash
git tag -a vX.Y.Z -m "Release vX.Y.Z"
git push origin main
git push origin vX.Y.Z
```

### 7. Automated Release

Pushing a tag triggers the release workflow:

1. **Build binaries** for:
   - Linux x86_64
   - macOS x86_64
   - macOS aarch64

2. **Create GitHub Release** with:
   - Compiled binaries
   - SHA256 checksums
   - Auto-generated release notes

3. **Update Homebrew tap** (if configured)

4. **Build Conda package** (if configured)

### 8. Publish to crates.io

```bash
# Dry run first
cargo publish --dry-run

# Publish
cargo publish
```

### 9. Post-Release

- Announce on relevant channels
- Update project homepage
- Monitor for issues

## Automated Workflows

### Release Workflow (`.github/workflows/release.yml`)

Triggered by: `git push origin vX.Y.Z`

Jobs:
1. `build`: Compile binaries for all platforms
2. `release`: Create GitHub release with artifacts
3. `update-homebrew`: Update Homebrew formula
4. `build-conda`: Build Conda packages

### Homebrew Update

The release workflow automatically updates the Homebrew tap:

1. Calculates SHA256 of source tarball
2. Updates `Formula/grit.rb` in `homebrew-grit` repo
3. Commits and pushes changes

**Requirements**:
- `HOMEBREW_TAP_TOKEN` secret with repo access

### crates.io Publishing

Manual step to ensure quality control:

```bash
cargo publish
```

**Requirements**:
- `CARGO_REGISTRY_TOKEN` for automated publishing
- Verified email on crates.io

## Hotfix Releases

For urgent bug fixes:

```bash
# Create hotfix branch from tag
git checkout -b hotfix/vX.Y.Z vX.Y.Z

# Apply fix
# ... make changes ...

# Bump patch version
# Update Cargo.toml and CHANGELOG.md

# Commit and tag
git commit -m "fix: critical bug description"
git tag -a vX.Y.Z+1 -m "Hotfix release"

# Push
git push origin hotfix/vX.Y.Z
git push origin vX.Y.Z+1

# Merge back to main
git checkout main
git merge hotfix/vX.Y.Z
git push origin main
```

## Release Notes Template

```markdown
## What's New in vX.Y.Z

### Highlights

Brief summary of major changes.

### New Features

- **Feature name**: Description (#PR)

### Bug Fixes

- Fixed issue with X (#PR)

### Performance

- Improved Y by Z% (#PR)

### Breaking Changes

- Changed behavior of X (migration guide below)

### Migration Guide

If upgrading from vA.B.C:

1. Step one
2. Step two

### Full Changelog

https://github.com/manish59/grit/compare/vA.B.C...vX.Y.Z
```

## Rollback Procedure

If a release has critical issues:

1. **Yank from crates.io** (if published):
   ```bash
   cargo yank --vers X.Y.Z
   ```

2. **Mark release as pre-release** on GitHub

3. **Communicate** via:
   - GitHub release notes
   - Issue for tracking

4. **Prepare hotfix** following hotfix process

## Security Releases

For security vulnerabilities:

1. **Do not disclose** details publicly
2. **Prepare fix** in private
3. **Coordinate disclosure** with reporter
4. **Release fix** with minimal details
5. **Disclose details** after users have time to update
6. **Credit reporter** (with permission)

See [SECURITY.md](../SECURITY.md) for vulnerability reporting.
