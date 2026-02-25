# GRIT Packaging

This directory contains packaging recipes for distributing GRIT.

## Homebrew (macOS/Linux)

### Option 1: Create Your Own Tap (Recommended)

1. **Create a new repository** named `homebrew-grit`:
   ```bash
   gh repo create manish59/homebrew-grit --public
   git clone https://github.com/manish59/homebrew-grit.git
   cd homebrew-grit
   ```

2. **Copy the formula**:
   ```bash
   mkdir -p Formula
   cp /path/to/grit/packaging/homebrew/grit.rb Formula/
   git add Formula/grit.rb
   git commit -m "Add grit formula v0.1.0"
   git push
   ```

3. **Users can then install with**:
   ```bash
   brew install manish59/grit/grit
   ```

### Option 2: Submit to Homebrew Core

For wider distribution, submit a PR to [homebrew-core](https://github.com/Homebrew/homebrew-core):
1. Fork homebrew-core
2. Add the formula to `Formula/g/grit.rb`
3. Submit a PR following their [contribution guidelines](https://docs.brew.sh/How-To-Open-a-Homebrew-Pull-Request)

---

## Bioconda (Conda)

### Submitting to Bioconda

1. **Fork bioconda-recipes**:
   ```bash
   gh repo fork bioconda/bioconda-recipes
   git clone https://github.com/YOUR_USERNAME/bioconda-recipes.git
   cd bioconda-recipes
   ```

2. **Create a new recipe**:
   ```bash
   mkdir -p recipes/grit-genomics
   cp /path/to/grit/packaging/bioconda/meta.yaml recipes/grit-genomics/
   cp /path/to/grit/packaging/bioconda/build.sh recipes/grit-genomics/
   ```

3. **Test locally** (optional but recommended):
   ```bash
   # Install bioconda-utils
   conda install -c conda-forge bioconda-utils

   # Test the recipe
   bioconda-utils lint --packages grit-genomics
   bioconda-utils build --packages grit-genomics
   ```

4. **Submit PR**:
   ```bash
   git checkout -b add-grit-genomics
   git add recipes/grit-genomics/
   git commit -m "Add grit-genomics v0.1.0"
   git push origin add-grit-genomics
   ```

   Then open a PR at https://github.com/bioconda/bioconda-recipes/pulls

5. **After merge, users install with**:
   ```bash
   conda install -c bioconda grit-genomics
   ```

---

## Installation Methods Summary

| Method | Command |
|--------|---------|
| **crates.io** | `cargo install grit-genomics` |
| **Homebrew** | `brew install manish59/grit/grit` |
| **Bioconda** | `conda install -c bioconda grit-genomics` |
| **From source** | `git clone ... && cargo build --release` |

All methods install the `grit` command-line tool.
