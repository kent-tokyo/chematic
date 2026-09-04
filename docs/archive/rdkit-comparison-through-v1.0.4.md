# chematic vs RDKit — Detailed Comparison

This page gives a direct comparison between chematic and RDKit for teams evaluating which library to use.

The performance figures are dated snapshots, not current-release guarantees.
Use [`benchmark.md`](../benchmark.md) for the latest reproducible run and its
exact environment.

See also: [`rdkit-migration.md`](../rdkit-migration.md) for a feature-by-feature Supported/Partial/Not-supported breakdown.

---

## TL;DR

| | chematic | RDKit |
|---|---|---|
| **Install** | `pip install chematic` | `pip install rdkit` (official prebuilt wheels) or conda |
| **Browser / WASM** | Yes — 3.30 MB raw / 1.21 MB gzip | Not applicable (Python/C++ library); RDKit.js is a separate community project — 6.91 MB raw |
| **C++ dependency** | None (default) | Required |
| **Batch fingerprint speed** | ~78 µs/mol (2–3× faster, diverse corpus) | ~160–235 µs/mol |
| **AI agent integration** | MCP server built-in | None |
| **Ecosystem maturity** | Growing (2024–) | Established (2006–) |

**Descriptor accuracy caveat**: `Mol.descriptors()` returns 194 values, but bulk (4,999-molecule) RDKit-agreement testing currently covers a subset — the physicochemical core (MW/HBA/HBD/TPSA/LogP/MR/Fsp3/ring & stereocenter counts/etc.) is a verified 100% (or near-100%) match. A larger set of descriptors that are *named* to match RDKit 1:1 (Kappa/HallKierAlpha/BertzCT/BalabanJ/BCUT2D/VSA descriptor families/MQN/SA Score) were found to diverge substantially once measured at corpus scale. Treat these as chematic-specific implementations unless a public validation table explicitly reports agreement.

---

## 1. Infrastructure

### Installation

| | chematic | RDKit |
|---|---|---|
| Python | `pip install chematic` | `pip install rdkit` (official prebuilt wheels) or `conda install -c conda-forge rdkit` |
| C++ compiler | Not required (pure Rust, even building from source) | Not required when installing the official prebuilt wheel; required (Boost, CMake) when building RDKit itself from source |
| Docker image delta | ~4 MB (approximate; not independently re-measured this cycle) | ~200 MB+ (approximate; not independently re-measured this cycle) |
| GitHub Actions | `pip install chematic` | Separate conda setup step |
| Cloudflare Workers | Yes | No |
| AWS Lambda | Yes | Difficult (binary size) |
| Embedded / no-std | Partial | No |

### WASM deployment

| Library | Bundle size | Build toolchain |
|---|---|---|
| **chematic** | **3.30 MB raw / 1.21 MB gzip** | `wasm-pack build` only |
| RDKit.js (`@rdkit/rdkit`) | 6.91 MB raw (`RDKit_minimal.wasm`, gzip not independently measured) | Emscripten SDK + cmake (upstream; not needed by consumers of the published npm package) |
| Indigo (Ketcher-oriented build, `indigo-ketcher`) | 11.24 MB raw (main `.wasm`; a "norender" variant is 7.17 MB) | Emscripten SDK + cmake (upstream) |

The chematic size was measured 2026-09-04 from the v1.0.2 release candidate using
`wasm-pack 0.13.1` + `wasm-opt 130 -O3`. The pinned historical comparators were obtained from
`@rdkit/rdkit@2025.3.4-1.0.0`'s file listing on unpkg.com and
`indigo-ketcher@1.45.1`'s file listing on jsDelivr. chematic's raw binary is
currently about 2.1× smaller than RDKit.js's and about 3.8× smaller than Indigo's Ketcher-oriented
build. See the exact chematic artifact digest in the
[2026-09-04 record](../../benchmarks/2026-09-04-wasm-size.md). chematic compiles to
`wasm32-unknown-unknown` natively — no Emscripten, no cmake, no clang.

---

## 2. Performance

All measurements: Python 3.13.6, Apple M4, chematic v0.4.29, RDKit 2026.03.3 (2026-07-17;
see [`benchmarks/2026-07-17.md`](../../benchmarks/2026-07-17.md) for full methodology). Import
time and SMILES parse throughput below were not remeasured this cycle.

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

Small repeated fixture (same 20-molecule set `scripts/benchmark_vs_rdkit.py` has always used):

| N molecules | chematic (`bulk.ecfp4`) | RDKit (Python loop) | Speedup |
|---|---|---|---|
| 100 | 1.4 ms | 8 ms | 6.1× |
| 1,000 | 10 ms | 83 ms | 8.0× |
| 10,000 | 126 ms | 839 ms | **6.7×** |

On a large, structurally diverse 5,000-molecule ChEMBL corpus the margin is narrower:
**~78 µs/mol** (chematic) vs ~160–235 µs/mol (RDKit), ~2–3×. chematic uses Rayon for
parallel batch processing across all CPU cores. Neither number reproduces the
previously-reported 3.6 µs/mol / 5–14× figures, even on the identical small-fixture
script — see [`benchmarks/2026-07-17.md`](../../benchmarks/2026-07-17.md) for what was
measured and what remains an open question.

Reproduce:
```bash
python scripts/benchmark_vs_rdkit.py --rdkit
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
| MCP server (AI agent API) | 20 tools (stdio only) | None |
| pKa prediction (built-in, rule-based screening — not for clinical use) | 23 SMARTS rules | External tool required |
| ADMET profile (built-in, rule-based screening — not for clinical use) | BBB / Caco-2 / hERG / CYP3A4 | External tool required |
| MAP4 fingerprint | Yes (Minervini 2020) | No (external package) |
| UFF force field for metals | Yes (Zn, Fe, Cu, …) | No |
| IUPAC name generation (offline) | 25+ compound classes | No |
| Retrosynthesis (template-based) | 60 retro-SMIRKS built-in | External tool required |
| Pure Rust, zero C/C++ toolchain even building from source | Yes | No (RDKit itself is a C++ library; official wheels avoid a *local* build, but the library's own build still requires cmake/Boost) |

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

- You want chemistry in the browser (WASM, 1.21 MB gzip, no server)
- You need a pure Rust stack with no C++ toolchain
- You deploy to Lambda, Cloudflare Workers, or other constrained environments
- You build AI agents and want native MCP tool integration
- You need fast batch processing (ECFP4: 2–3× faster, Rayon-parallel)
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

*Benchmark data: Apple M4, Python 3.13.6, chematic v0.18.0, RDKit 2026.03.3 (2026-07-17 snapshot;
see [`benchmarks/2026-07-17.md`](../../benchmarks/2026-07-17.md) for the underlying methodology —
descriptor calculation paths unchanged since, not re-measured this cycle). The chematic bundle-size
figure above was measured 2026-09-04 from the v1.0.2 candidate; comparator sizes are pinned historical
figures (see the WASM deployment table).*
*Reproduce all benchmarks: see [benchmark details](../benchmark.md).*
