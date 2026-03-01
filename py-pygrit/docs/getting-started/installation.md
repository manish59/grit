# Installation

## Requirements

- Python 3.9 or later
- NumPy 1.20 or later

## Install from PyPI

The recommended way to install pygrit is via pip:

```bash
pip install pygrit
```

## Install with Optional Dependencies

### Development Tools

For running tests and type checking:

```bash
pip install pygrit[dev]
```

### Documentation Tools

For building documentation locally:

```bash
pip install pygrit[docs]
```

## Install from Source

To install from source, you need Rust and maturin:

### 1. Install Rust

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
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
print(pygrit.__version__)
```

## Platform Support

| Platform | Architecture | Status |
|----------|--------------|--------|
| Linux | x86_64 | Supported |
| Linux | aarch64 | Supported |
| macOS | x86_64 | Supported |
| macOS | arm64 (Apple Silicon) | Supported |
| Windows | x86_64 | Experimental |

## Troubleshooting

### ImportError: No module named 'pygrit'

Ensure you have installed pygrit in your active Python environment:

```bash
pip install pygrit
```

### Compilation Errors on Source Install

Make sure you have the latest Rust toolchain:

```bash
rustup update
```

### NumPy Version Conflicts

If you encounter NumPy compatibility issues:

```bash
pip install --upgrade numpy
pip install --force-reinstall pygrit
```
