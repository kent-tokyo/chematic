# Benchmark

Measured environment: Python 3.12, Apple M-series, chematic v0.4.23, RDKit 2026.03.3.

---

## Summary

| Metric | chematic | RDKit |
|--------|----------|-------|
| Import time | **~35 ms** | ~400 ms (11×) |
| SMILES parse — 5,000 mol | **~5 ms** | ~50 ms (10×) |
| ECFP4 batch — 10,000 mol | **36 ms** | ~500 ms (14×) |
| Descriptor accuracy vs RDKit | **19 metrics 100%** — MW/HBA/HBD/TPSA/LogP/ARC/RotB/Spiro/Bridge/… (4,999-mol) | baseline |
| Install | `pip install chematic` | conda or cmake |
| C/C++ dependencies | **Zero** | Required |
| WASM binary size | **504 KB** | ~30 MB |

---

## 1. Startup Time (import latency)

Cold-process import time measured by spawning a fresh Python subprocess per sample
(5 samples, median reported).  No module-cache warm-up.

| | chematic | RDKit |
|---|---|---|
| `import` only | **~35 ms** | ~400 ms |
| `import` + first parse | **~38 ms** | ~430 ms |
| **Speedup** | **~11×** | baseline |

**chematic**

```bash
python scripts/bench_startup.py --runs 5
```

**Why chematic is faster**: chematic is a single PyO3 extension module with no
transitive Python dependencies. RDKit initialises multiple C++ modules, reads
SMARTS definition files, and triggers Boost data-structure setup on first import.

---

## 2. SMILES Parse Throughput

Timed on the built-in 20-molecule diverse corpus repeated to 5,000 total parses.
Warm-up pass excluded.

| N | chematic | RDKit | Speedup |
|---|---|---|---|
| 1,000 | ~1 ms | ~10 ms | ~10× |
| 5,000 | **~5 ms** | **~50 ms** | **~10×** |
| 10,000 | ~10 ms | ~100 ms | ~10× |

Per-molecule: **~1 µs/mol** (chematic) vs ~10 µs/mol (RDKit).

**chematic**

```python
import chematic
mols = [chematic.from_smiles(s) for s in smiles_list]
# or batch:
mols = chematic.from_smiles_list(smiles_list)
```

**RDKit**

```python
from rdkit import Chem
mols = [Chem.MolFromSmiles(s) for s in smiles_list]
```

### How to reproduce

```bash
python scripts/bench_smiles_parse.py --n 5000 --rdkit
```

---

## 3. Speed — ECFP4 Fingerprint Generation (batch)

Rayon parallelism across all CPU cores; speedup grows with batch size.

| Molecules (N) | chematic (`bulk.ecfp4`) | RDKit (Python loop) | Speedup |
|---------------|------------------------|---------------------|---------|
| 100 | 0.36 ms | 2 ms | 5× |
| 1,000 | 3.6 ms | 20 ms | 5× |
| 10,000 | 36 ms | ~500 ms | **~14×** |

Per-molecule: **3.6 µs/mol** (chematic) vs 20–50 µs/mol (RDKit).

**chematic**

```python
import chematic
fps = chematic.bulk.ecfp4(smiles_list)  # (N, 2048) uint8 numpy array
```

**RDKit**

```python
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(
           Chem.MolFromSmiles(s), 2, 2048)
       for s in smiles_list]
```

### How to reproduce

```bash
python scripts/benchmark_vs_rdkit.py --rdkit
```

---

## 4. Descriptor Accuracy vs RDKit

Tested on a 4,999-molecule ChEMBL-derived SMILES corpus (`scripts/bench5k.py`). See [Validation](validation.md) for full per-metric breakdown.

| Descriptor | Agreement | Tolerance |
|-----------|-----------|-----------|
| Molecular weight | **100%** | exact |
| Heavy atom count | **100%** | exact |
| H-bond donors (HBD) | **100%** | exact |
| H-bond acceptors (HBA) | **100%** | exact |
| TPSA | **100%** | ±0.1 Å² |
| LogP (Crippen) | **100%** | exact* |
| MR (molar refractivity) | **100%** | ±0.01 |
| Fsp3 | **100%** | ±0.001 |
| Aromatic ring count | **100%** | exact |
| Aliphatic ring count | **100%** | exact |
| Saturated ring count | **100%** | exact |
| Rotatable bonds | **100%** | exact |
| Num heteroatoms | **100%** | exact |
| Num spiro atoms | **100%** | exact |
| Num bridgehead atoms | **100%** | exact |
| Num amide bonds | **100%** | exact |
| Aromatic/aliphatic heterocycles | **100%** | exact |
| Num stereocenters | **99.98%** | exact |
| [nH] SMARTS match | **100%** | precision/recall |

19 of 19 metrics reach ≥99.98% on the 4,999-molecule ChEMBL corpus (RDKit 2026.03.3).
Stereocenters: 4998/4999 vs `CalcNumAtomStereoCenters` (legacy RDKit CIP); the 1 discrepancy
is a highly-symmetric polyester where chematic and RDKit's new CIP agree on 4 stereocenters
while the legacy function returns 2 (under-count in legacy).

### How to reproduce

```bash
# Requires RDKit and the 5k SMILES file
python scripts/bench5k.py path/to/SMILES.csv --detail
```

---

## 5. Installation & Deployment

| | chematic | RDKit |
|---|----------|-------|
| Python | `pip install chematic` | `conda install -c conda-forge rdkit` |
| C/C++ compiler | Not required | Required (Boost) |
| Docker image size delta | ~4 MB | ~200 MB+ |
| GitHub Actions | Single pip line | Separate conda setup step |
| JavaScript / WASM | `npm install @kent-tokyo/chematic` (504 KB) | No official package |
| Browser deployment | Yes | No |

---

## 6. Feature Comparison

| Feature | chematic | RDKit |
|---------|----------|-------|
| pKa prediction | Built-in (15 SMARTS rules) | External tool required |
| ADMET profile (BBB, Caco-2, hERG, CYP3A4) | Built-in | External tool required |
| MCP server (AI agent integration) | 15 tools | Not available |
| LSH approximate nearest-neighbour index | Built-in | Not available |
| IUPAC name generation | Built-in (offline) | Not available |
| Browser / WASM deployment | Yes (504 KB) | No |
| ECFP4 batch speed | 5–14× faster | Baseline |
| SMARTS atom map `:N` | Yes | Yes |
| Retrosynthesis (template-based) | 60 retro-SMIRKS built-in | External tool |
| File formats | 20+ | 100+ |
| 3D conformer quality | Good (ETKDG rules) | Better (ML-assisted) |
| Community & publications | Growing | Established (20+ years) |

---

## 7. Batch Descriptor Computation

`chematic.bulk.descriptors` returns 55+ descriptors per molecule including ADMET and pKa — all in parallel.
`chematic.bulk.descriptors_array` returns selected columns as numpy arrays (~25% faster for column-oriented access).

| N | chematic (`bulk.descriptors`) | Descriptors per call |
|---|-------------------------------|---------------------|
| 100 | ~10 ms | 55+ (incl. pKa, ADMET) |
| 1,000 | ~50 ms | 55+ (incl. pKa, ADMET) |

```python
import chematic
import pandas as pd

# list-of-dicts (general purpose)
df = pd.DataFrame(chematic.bulk.descriptors(smiles_list))

# columnar numpy arrays (faster for specific columns)
result = chematic.bulk.descriptors_array(smiles_list, ["mw", "logp", "tpsa", "hba"])
df = pd.DataFrame(result)   # float64 / bool arrays, no per-molecule dict overhead
```

### Compound screening

```python
# One call bundles lipinski / veber / pains / brenk / qed / sa_score:
results = chematic.screen(smiles_list, profile="druglike")
passing = [r for r in results if r["overall_pass"]]
```

### Large SDF files (streaming)

```python
# iter_sdf() streams one record at a time — no full-file load:
for rec in chematic.iter_sdf("large.sdf"):
    print(rec.smiles, rec.get("Activity"))

# batch pipeline:
for batch in chematic.iter_sdf_batched("large.sdf", batch_size=1000):
    descs = chematic.bulk.descriptors([r.smiles for r in batch])
```

RDKit's `rdkit.Chem.Descriptors.CalcMolDescriptors` covers ~200 descriptors but does not include pKa or ADMET.
