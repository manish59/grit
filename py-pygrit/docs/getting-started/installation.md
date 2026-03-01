# Installation

## Requirements

- Python 3.9 or later
- NumPy 1.20 or later

## Install from PyPI

```bash
pip install grit-genomics
```

The package is imported as `pygrit`:

```python
import pygrit
print(pygrit.__version__)
```

## Install from Source

Building from source requires [Rust](https://rustup.rs/) and [maturin](https://github.com/PyO3/maturin).

### 1. Install Rust

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
```

### 2. Clone the Repository

```bash
git clone https://github.com/manish59/grit.git
cd grit/py-pygrit
```

### 3. Create a Virtual Environment

```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

### 4. Install maturin and Build

```bash
pip install maturin
maturin develop --release
```

## Verify Installation

```python
import pygrit

# Check version
print(pygrit.__version__)

# Test basic functionality
iv = pygrit.Interval("chr1", 100, 200)
print(iv)  # chr1	100	200
```

## Platform Support

| Platform | Architecture | Pre-built Wheel | From Source |
|----------|--------------|-----------------|-------------|
| macOS | ARM64 (Apple Silicon) | Yes | Yes |
| macOS | x86_64 | No | Yes |
| Linux | x86_64 | No | Yes |
| Linux | ARM64 | No | Yes |
| Windows | x86_64 | No | Untested |

For platforms without pre-built wheels, pip will attempt to build from source. This requires a Rust toolchain to be installed.

## Troubleshooting

### ImportError: No module named 'pygrit'

Ensure you have installed grit-genomics:

```bash
pip install grit-genomics
```

### Build Errors (Source Install)

Make sure you have the latest Rust toolchain:

```bash
rustup update
```

### NumPy Version Issues

```bash
pip install --upgrade numpy
pip install --force-reinstall grit-genomics
```
