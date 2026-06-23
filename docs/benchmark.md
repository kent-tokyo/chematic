# Benchmark

Measured environment: Python 3.12, Apple M-series, chematic v0.4.17, RDKit 2024.09.

---

## Summary

| Metric | chematic | RDKit |
|--------|----------|-------|
| ECFP4 batch — 10,000 mol | **36 ms** | ~500 ms |
| Descriptor accuracy vs RDKit | **100%** on 5,000-mol corpus | baseline |
| Install | `pip install chematic` | conda or cmake |
| C/C++ dependencies | **Zero** | Required |
| WASM binary size | **~550 KB** | ~30 MB |

---

## 1. Speed — ECFP4 Fingerprint Generation (batch)

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

## 2. Descriptor Accuracy vs RDKit

Tested on a 5,000-molecule ChEMBL-like SMILES corpus (`scripts/bench5k.py`).

| Descriptor | Agreement | Tolerance |
|-----------|-----------|-----------|
| Molecular weight | 100% | exact |
| Heavy atom count | 100% | exact |
| H-bond donors (HBD) | 100% | exact |
| H-bond acceptors (HBA) | 100% | exact |
| TPSA | 100% | ±0.1 Å² |
| LogP (Crippen) | 100% | ±0.3 |
| Aromatic ring count | 100% | exact |

All metrics reached 100% agreement with `rdkit.Chem.Descriptors` on the 5,000-molecule corpus as of v0.4.14.

### How to reproduce

```bash
# Requires RDKit and the 5k SMILES file
python scripts/bench5k.py path/to/SMILES.csv --detail
```

---

## 3. Installation & Deployment

| | chematic | RDKit |
|---|----------|-------|
| Python | `pip install chematic` | `conda install -c conda-forge rdkit` |
| C/C++ compiler | Not required | Required (Boost) |
| Docker image size delta | ~4 MB | ~200 MB+ |
| GitHub Actions | Single pip line | Separate conda setup step |
| JavaScript / WASM | `npm install @kent-tokyo/chematic` (~550 KB) | No official package |
| Browser deployment | Yes | No |

---

## 4. Feature Comparison

| Feature | chematic | RDKit |
|---------|----------|-------|
| pKa prediction | Built-in (15 SMARTS rules) | External tool required |
| ADMET profile (BBB, Caco-2, hERG, CYP3A4) | Built-in | External tool required |
| MCP server (AI agent integration) | 15 tools | Not available |
| LSH approximate nearest-neighbour index | Built-in | Not available |
| IUPAC name generation | Built-in (offline) | Not available |
| Browser / WASM deployment | Yes (~550 KB) | No |
| ECFP4 batch speed | 5–14× faster | Baseline |
| SMARTS atom map `:N` | Yes | Yes |
| Retrosynthesis (template-based) | 60 retro-SMIRKS built-in | External tool |
| File formats | 20+ | 100+ |
| 3D conformer quality | Good (ETKDG rules) | Better (ML-assisted) |
| Community & publications | Growing | Established (20+ years) |

---

## 5. Batch Descriptor Computation

`chematic.bulk.descriptors` returns 55+ descriptors per molecule including ADMET and pKa — all in parallel.

| N | chematic (`bulk.descriptors`) | Descriptors per call |
|---|-------------------------------|---------------------|
| 100 | ~10 ms | 55+ (incl. pKa, ADMET) |
| 1,000 | ~50 ms | 55+ (incl. pKa, ADMET) |

```python
import chematic
import pandas as pd

df = pd.DataFrame(chematic.bulk.descriptors(smiles_list))
# One call returns: mw, logp, tpsa, hbd, hba, qed, sa_score,
#                   pka_acid, pka_base, bbb_score, caco2, herg_risk, ...
```

RDKit's `rdkit.Chem.Descriptors.CalcMolDescriptors` covers ~200 descriptors but does not include pKa or ADMET.
