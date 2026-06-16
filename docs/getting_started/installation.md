# Installation

## Requirements

- Python ≥ 3.8
- No C/C++ compiler required
- No conda required

## Install from PyPI

```bash
pip install chematic
```

That's it. `chematic` ships as a pre-built binary wheel for:

| Platform | x86_64 | aarch64 (Apple M-series / AWS Graviton) |
|---|---|---|
| Linux | yes | yes |
| macOS | yes | yes |
| Windows | yes | — |

## Verify the installation

```python
import chematic
print(chematic.__version__)

mol = chematic.from_smiles("c1ccccc1")
print(mol.formula)   # C6H6
```

## Optional dependencies

`chematic` itself has no required runtime dependencies beyond Python and numpy.

For the cookbook examples you may also want:

```bash
pip install pandas scikit-learn matplotlib
```

## Build from source

If a pre-built wheel is not available for your platform:

```bash
pip install maturin
git clone https://github.com/kent-tokyo/chematic
cd chematic/crates/chematic-py
maturin develop --release
```

Requires Rust ≥ 1.75 (`rustup install stable`).
