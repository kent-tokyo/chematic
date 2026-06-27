# chematic vs RDKit — Detailed Comparison

This page gives a direct comparison between chematic and RDKit for teams evaluating which library to use.

---

## TL;DR

| | chematic | RDKit |
|---|---|---|
| **Install** | `pip install chematic` | conda / cmake required |
| **Browser / WASM** | Yes — 504 KB | No |
| **C++ dependency** | None (default) | Required |
| **Batch fingerprint speed** | 3.6 µs/mol | 20–50 µs/mol |
| **AI agent integration** | MCP server built-in | None |
| **Ecosystem maturity** | Growing (2024–) | Established (2006–) |

---

## 1. Infrastructure

### Installation

| | chematic | RDKit |
|---|---|---|
| Python | `pip install chematic` | `conda install -c conda-forge rdkit` |
| C++ compiler | Not required | Required (Boost, CMake) |
| Docker image delta | ~4 MB | ~200 MB+ |
| GitHub Actions | `pip install chematic` | Separate conda setup step |
| Cloudflare Workers | Yes | No |
| AWS Lambda | Yes | Difficult (binary size) |
| Embedded / no-std | Partial | No |

### WASM deployment

| Library | Bundle size | Build toolchain |
|---|---|---|
| **chematic** | **504 KB gzip** | `wasm-pack build` only |
| RDKit.js | ~30 MB | Emscripten SDK + cmake |
| Indigo WASM | ~40 MB | Emscripten SDK + cmake |

chematic compiles to `wasm32-unknown-unknown` natively — no Emscripten, no cmake, no clang.

---

## 2. Performance

All measurements: Python 3.12, Apple M-series, chematic v0.4.22, RDKit 2026.03.3.

### Import time (cold process)

| | chematic | RDKit |
|---|---|---|
| `import` only | ~35 ms | ~400 ms |
| `import` + first parse | ~38 ms | ~430 ms |
| **Speedup** | **~11×** | baseline |

Measured by spawning a fresh subprocess per sample (module cache excluded).

Reproduce:
```bash
python scripts/bench_startup.py --runs 5 --rdkit
```

### SMILES parse throughput

| N molecules | chematic | RDKit | Speedup |
|---|---|---|---|
| 1,000 | ~1 ms | ~10 ms | ~10× |
| 5,000 | ~5 ms | ~50 ms | ~10× |
| 10,000 | ~10 ms | ~100 ms | ~10× |

Per-molecule: **~1 µs/mol** (chematic) vs ~10 µs/mol (RDKit).

Reproduce:
```bash
python scripts/bench_smiles_parse.py --n 5000 --rdkit
```

### ECFP4 fingerprint generation (batch)

| N molecules | chematic (`bulk.ecfp4`) | RDKit (Python loop) | Speedup |
|---|---|---|---|
| 100 | 0.36 ms | 2 ms | 5× |
| 1,000 | 3.6 ms | 20 ms | 5× |
| 10,000 | 36 ms | ~500 ms | **~14×** |

chematic uses Rayon for parallel batch processing. Speedup grows with batch size.

Reproduce:
```bash
python scripts/bench_smiles_parse.py --n 10000 --rdkit
```

### Where RDKit is faster or better

| Task | RDKit advantage |
|---|---|
| Publication-quality 3D structures | ETKDGv3 with ML torsion corrections (chematic: good for screening) |
| Exotic molecule handling | 20 years of edge-case fixes |
| Large SDF file streaming | Optimized C++ reader |

---

## 3. Feature coverage

### chematic has, RDKit does not

| Feature | chematic | RDKit |
|---|---|---|
| Native WASM (no Emscripten) | Yes | No |
| MCP server (AI agent API) | 15 tools | None |
| pKa prediction (built-in) | 15 SMARTS rules | External tool required |
| ADMET profile (built-in) | BBB / Caco-2 / hERG / CYP3A4 | External tool required |
| MAP4 fingerprint | Yes (Minervini 2020) | No (external package) |
| UFF force field for metals | Yes (Zn, Fe, Cu, …) | No |
| IUPAC name generation (offline) | 25+ compound classes | No |
| Retrosynthesis (template-based) | 60 retro-SMIRKS built-in | External tool required |
| `pip install` anywhere | Yes | No (conda/cmake) |

### RDKit has, chematic does not (or is weaker)

| Feature | RDKit advantage |
|---|---|
| Publication-quality 3D conformers | ETKDGv3 with ML torsion corrections; chematic uses chair/envelope + MMFF94 (good for virtual screening) |
| File format support | 100+ formats (chematic: ~20) |
| Validated production docking | Years of benchmarking |
| Community plugins | Large ecosystem |
| Exact InChI (default) | C lib bundled by default |

---

## 4. API comparison (Python)

Most common operations map directly:

| Operation | chematic | RDKit |
|---|---|---|
| Parse SMILES | `chematic.from_smiles(s)` | `Chem.MolFromSmiles(s)` |
| Molecular weight | `mol.mw` | `Descriptors.MolWt(mol)` |
| LogP | `mol.logp` | `Descriptors.MolLogP(mol)` |
| TPSA | `mol.tpsa` | `Descriptors.TPSA(mol)` |
| ECFP4 | `chematic.bulk.ecfp4(smiles)` | `AllChem.GetMorganFingerprintAsBitVect(mol, 2)` |
| Substructure | `mol.has_substructure("[OH]")` | `mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]"))` |
| Batch descriptors | `chematic.descriptors_df(smiles)` | `PandasTools` + manual loop |
| Drug-likeness | `mol.lipinski_passes` | `Descriptors.rdMolDescriptors.CalcNumHBD(mol) <= 5 and …` |
| Canonical SMILES | `mol.smiles` | `Chem.MolToSmiles(mol)` |

---

## 5. When to choose

**Choose chematic if:**

- You want chemistry in the browser (WASM, 504 KB, no server)
- You need a pure Rust stack with no C++ toolchain
- You deploy to Lambda, Cloudflare Workers, or other constrained environments
- You build AI agents and want native MCP tool integration
- You need fast batch processing (ECFP4: 5–14× faster)
- You want `pip install` to just work — anywhere

**Choose RDKit if:**

- You need maximum ecosystem compatibility and 20+ years of production validation
- You need publication-quality 3D structures with ML-assisted torsion corrections (ETKDGv3)
- You rely on community plugins written against the RDKit Python API
- You need bit-exact standard InChI without an extra feature flag
- You work with exotic file formats or unusual molecule types

---

## 6. Migration quick-reference

If you have existing RDKit code, these are the most common substitutions:

```python
# RDKit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem

mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
mw  = Descriptors.MolWt(mol)
fp  = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
```

```python
# chematic
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
mw  = mol.mw
fp  = chematic.bulk.ecfp4(["CC(=O)Oc1ccccc1C(=O)O"])[0]  # numpy uint8 array
```

For large-scale migration, chematic's Python API is designed to be familiar to
RDKit users while adding batch-first and browser-native capabilities.

---

*Benchmark data: Apple M-series, Python 3.12, chematic v0.4.22, RDKit 2026.03.3.*  
*Reproduce all benchmarks: see [benchmark details](benchmark.md).*
