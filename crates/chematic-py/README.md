# chematic

Pure-Rust cheminformatics library for Python — SMILES parsing, 190+ descriptor values (71 functions), fingerprints, pKa prediction, ADMET profiling, and template-based retrosynthesis.

## Installation

```bash
pip install chematic
```

## Quick Start

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin

print(mol.mw)      # 180.16
print(mol.logp)    # 1.31
print(mol.tpsa)    # 63.6
print(mol.qed)     # 0.55

# New descriptors
print(mol.vabc)              # van der Waals volume (no 3D needed)
print(mol.schultz_mti)       # Schultz MTI
print(mol.gutman_mti)        # Gutman MTI*
print(mol.gravitational_index)  # gravitational index

# pKa prediction
print(mol.pka())   # {"most_acidic": 3.49, "most_basic": None}

# ADMET profile
print(mol.admet())
# {"bbb": False, "bbb_score": ..., "caco2": ..., "herg_risk": ..., "cyp3a4_risk": ...}

# Fingerprints (bytes, 2048-bit ECFP4)
fp = mol.ecfp4()

# Tanimoto similarity
mol2 = chematic.from_smiles("c1ccccc1")
sim = chematic.tanimoto(mol.ecfp4(), mol2.ecfp4())

# Natural-language property summary (for LLM / MCP agents)
print(mol.describe())

# Structural diff between two molecules
ibuprofen = chematic.from_smiles("CC(C)Cc1ccc(CC(C)C(=O)O)cc1")
d = mol.diff(ibuprofen)  # {"summary": "...", "delta_mw": 66.1, "delta_logp": 2.75, ...}

# SVG / PDF / EPS depiction
svg = mol.to_svg()
pdf_bytes = mol.to_pdf()   # bytes; requires pdf feature
eps_str   = mol.to_eps()   # PostScript string

# ChemicalJSON (Avogadro 2 / MolSSI)
cjson_str = mol.to_cjson(coords=[])   # coords: list of (x,y,z) tuples, optional
mol2, coords = chematic.from_cjson(cjson_str)

# Template-based retrosynthesis (60 retro-SMIRKS templates)
mol3 = chematic.from_smiles("CC(=O)Nc1ccccc1")  # acetanilide
results = mol3.retro_disconnect(max_results=5)
for r in results:
    print(r["template"], "→", r["precursors"])
# amide_secondary → ['CC(=O)O', 'Nc1ccccc1']

# Filter by reaction class
amides = mol3.retro_disconnect(reaction_class="AmideBond")

# Bulk substructure match against a pre-parsed Mol list (returns indices)
mols = [chematic.from_smiles(s) for s in ["CCO", "c1ccccc1O", "CC(=O)O"]]
hits = chematic.bulk.substructure_match("[OH]", mols)  # → [0, 1, 2]

# All descriptors as a dict (for Pandas)
import pandas as pd
smiles = ["CCO", "c1ccccc1", "CC(=O)O"]
df = pd.DataFrame([chematic.from_smiles(s).descriptors() for s in smiles])
```

## Features

- **Zero C/C++ dependencies** — pure Rust, no RDKit or OpenBabel required
- **SMILES / MOL / SDF / ChemicalJSON** parsing and writing
- **190+ descriptor values** (71 functions; MQN returns 42 values, BCUT2D / autocorr2d / geary / moran return multi-value arrays): MW, LogP (±0.01, 96.5% of 4,999-mol ChEMBL subset), TPSA (±0.1 Å², 98.1%), QED, Fsp3, SA Score, HBD (100% vs RDKit, incl. S-H), `vabc`, `schultz_mti`, `gutman_mti`, `gravitational_index`
- **14 fingerprint algorithms**: ECFP2/4/6, FCFP4/6, MACCS, AtomPair, Torsion, …
- **pKa prediction** (15 SMARTS rules — unique to chematic)
- **ADMET profile**: BBB, Caco-2, hERG, CYP3A4
- **Template-based retrosynthesis**: `mol.retro_disconnect()` — 60 retro-SMIRKS templates, SA Score ranked
- **SMARTS substructure search** — `chematic.smarts_match()` and `bulk.substructure_match(smarts, mols)` (pre-parsed Mol list, returns indices)
- **SVG / PDF / EPS depiction**: `mol.to_svg()`, `mol.to_pdf()`, `mol.to_eps()`
- **ChemicalJSON**: `mol.to_cjson(coords=[])`, `chematic.from_cjson(s)` — Avogadro 2 / MolSSI compatible

## RDKit compatibility

`chematic.rdkit_compat` provides a lightweight RDKit-compatible subset for environments where RDKit is unavailable (WASM, serverless, conda-free CI):

```python
from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import Descriptors, rdMolDescriptors, DataStructs

mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")

# Descriptors
Descriptors.MolWt(mol)          # 180.16
rdMolDescriptors.CalcTPSA(mol)  # 63.6

# Fingerprint (ExplicitBitVect) with bitInfo
bitInfo = {}
fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048, bitInfo=bitInfo)
fp.GetNumBits()                         # 2048
bitInfo                                 # {bit: ((atom_idx, radius), ...)}
DataStructs.TanimotoSimilarity(fp, fp)  # 1.0
DataStructs.BulkTanimotoSimilarity(fp, [fp])  # [1.0]

import numpy as np
arr = DataStructs.ConvertToNumpyArray(fp)  # (2048,) int8 for sklearn / PyTorch

# Atom / Bond traversal
for atom in mol.GetAtoms():
    atom.GetSymbol(), atom.GetAtomicNum(), atom.IsInRing()
for bond in mol.GetBonds():
    bond.GetBondType(), bond.GetBondTypeAsDouble(), bond.IsInRing()

# Ring information
ri = mol.GetRingInfo()
ri.NumRings()       # 1
ri.AtomRings()      # tuple of tuples of atom indices
ri.NumAtomRings(0)  # rings containing atom 0

# SDF I/O with SD properties
with Chem.SDWriter("out.sdf") as w:
    mol.SetProp("ID", "aspirin")
    w.write(mol)
for m in Chem.SDMolSupplier("out.sdf"):
    print(m.GetProp("ID"))
```

Unsupported options raise `NotImplementedError` or `TypeError` — they are never silently ignored.

### Compatibility matrix

| Area | Status | Notes |
|------|--------|-------|
| SMILES I/O | ✅ Supported | `MolFromSmiles` (aromaticity perceived when `sanitize=True`) / `MolToSmiles` |
| SDF I/O | ✅ Supported | `SDMolSupplier` / `SDWriter` + SD properties |
| Mol properties | ✅ Supported | `Get`/`Set`/`Has`/`Clear`Prop, typed setters, `GetPropsAsDict` |
| Mol / Atom / Bond | ✅ Supported | read-only traversal (`GetAtoms`/`GetBonds`/`GetAtomWithIdx`/…) |
| RingInfo | ✅ Supported | SSSR-based; `NumRings`/`AtomRings`/`BondRings`/`NumAtomRings`/`NumBondRings` |
| Substructure | 🟡 Partial | SMARTS via chematic; match **order** may differ from RDKit (use set comparison) |
| Descriptors | ✅ Supported | MW/HBA/HBD exact, TPSA ±1.0, LogP ±0.5 vs RDKit (differential-tested) |
| Morgan fingerprint | 🟡 Partial | `nBits` folding + `bitInfo` shape-/origin-consistent, **not RDKit bit-identical** (FNV-1a vs MurmurHash) |
| DataStructs | ✅ Supported | `TanimotoSimilarity`/`DiceSimilarity`/`BulkTanimotoSimilarity`/`ConvertToNumpyArray` |
| RWMol / editing | ❌ Unsupported | read-only layer |
| `useFeatures`, `useBondTypes=False` | 🔊 Fails loudly | raise `NotImplementedError` instead of silently ignoring |

A live differential suite (`tests/test_rdkit_diff.py`, auto-skipped when RDKit is
absent) compares chematic against RDKit across descriptors, ring counts, SMARTS
match counts, SDF round-trips, and Morgan self-similarity, writing an explainable
diff to `validation/results/rdkit_diff.jsonl`.

`chematic.rdkit_compat` is **not a full RDKit clone** — it is a lightweight
RDKit-compatible subset for common 2D cheminformatics workflows. See the full
[RDKit compatibility guide](https://github.com/kent-tokyo/chematic/blob/main/docs/rdkit_compat.md)
(compatibility matrix, differential-validation results, known divergences, and runnable examples).

## License

MIT OR Apache-2.0
