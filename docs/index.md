# chematic

**Pure-Rust cheminformatics for Python** — no C/C++ dependencies, WASM-native, 190+ descriptors.

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions)
[![PyPI](https://img.shields.io/pypi/v/chematic)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic)](https://crates.io/crates/chematic)

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
print(mol.mw, mol.logp, mol.tpsa)   # 180.16  1.31  63.6
print(mol.lipinski_passes)           # True
print(mol.admet())                   # BBB, Caco-2, hERG, CYP3A4 in one call
```

## Why chematic?

| | chematic | RDKit |
|---|---|---|
| Install | `pip install chematic` | conda or complex build |
| C/C++ deps | **Zero** | Required |
| WASM | **Yes** (~1.1 MB gzip) | No (30–50 MB) |
| pKa prediction | **Built-in** | External tool |
| ADMET profile | **Built-in** | External tool |
| Pure Python wheel | **Yes** | No |

## Features

- **SMILES / SDF / MOL / InChI** parsing and writing
- **190+ molecular descriptors** (MW, LogP, TPSA, QED, SA Score, pKa, ADMET …)
- **Fingerprints**: ECFP4/6, FCFP4/6, MACCS, AtomPair, Torsion, Layered
- **SMARTS** substructure search with full recursive SMARTS support
- **Reactions**: SMIRKS application, reaction SMARTS matching, MDL RXN I/O
- **Standardisation**: salt removal, charge neutralisation, tautomer normalisation
- **Scaffolds**: Murcko, generic Murcko, BRICS fragmentation
- **3D**: distance geometry, UFF/MMFF94 minimisation, SASA, conformer ensembles
- **Visualisation**: 2D SVG with CPK colours, atom highlighting, grid depiction
- **Bulk / parallel**: Rayon-powered batch descriptors, fingerprint matrices, similarity search
- **LSH index**: approximate nearest-neighbour search for large libraries
- **`mol.describe()`**: natural-language property summary for LLM / MCP agents
- **`mol.diff(other)`**: element-level and descriptor-level structural diff

## Installation

```bash
pip install chematic
```

Requires Python ≥ 3.8. No conda, no RDKit, no C compiler needed.

## Quick links

- [Cookbook](cookbook.md) — 20 copy-paste-ready tasks
- [Use cases](use-cases/) — AI agent workflows, notebooks, browser apps, Rust servers, batch analysis
- [Benchmark](benchmark.md) — performance vs RDKit, descriptor accuracy
- [RDKit migration guide](rdkit_cheatsheet.md) — side-by-side API comparison
- [API Reference](api/chematic.md) — full Python API
- [GitHub](https://github.com/kent-tokyo/chematic)
- [crates.io](https://crates.io/crates/chematic)
- [PyPI](https://pypi.org/project/chematic/)
