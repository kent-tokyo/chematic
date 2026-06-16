# chematic

Pure-Rust cheminformatics library for Python — SMILES parsing, 70+ molecular descriptors, fingerprints, pKa prediction, and ADMET profiling.

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

# All descriptors as a dict (for Pandas)
import pandas as pd
smiles = ["CCO", "c1ccccc1", "CC(=O)O"]
df = pd.DataFrame([chematic.from_smiles(s).descriptors() for s in smiles])
```

## Features

- **Zero C/C++ dependencies** — pure Rust, no RDKit or OpenBabel required
- **SMILES / MOL / SDF** parsing and writing
- **70+ descriptors**: MW, LogP, TPSA, QED, Fsp3, SA Score, and more
- **14 fingerprint algorithms**: ECFP2/4/6, FCFP4/6, MACCS, AtomPair, Torsion, …
- **pKa prediction** (15 SMARTS rules — unique to chematic)
- **ADMET profile**: BBB, Caco-2, hERG, CYP3A4
- **SMARTS substructure search**
- **SVG depiction** for Jupyter notebooks

## License

MIT OR Apache-2.0
