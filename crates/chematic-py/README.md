# chematic

Pure-Rust cheminformatics library for Python — SMILES parsing, 70+ molecular descriptors, fingerprints, pKa prediction, ADMET profiling, and template-based retrosynthesis.

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
- **70+ descriptors**: MW, LogP (±0.3 vs RDKit), TPSA (±1.0 Å²), QED, Fsp3, SA Score, HBD (100% vs RDKit, incl. S-H), `vabc`, `schultz_mti`, `gutman_mti`, `gravitational_index`
- **14 fingerprint algorithms**: ECFP2/4/6, FCFP4/6, MACCS, AtomPair, Torsion, …
- **pKa prediction** (15 SMARTS rules — unique to chematic)
- **ADMET profile**: BBB, Caco-2, hERG, CYP3A4
- **Template-based retrosynthesis**: `mol.retro_disconnect()` — 60 retro-SMIRKS templates, SA Score ranked
- **SMARTS substructure search** — `chematic.smarts_match()` and `bulk.substructure_match(smarts, mols)` (pre-parsed Mol list, returns indices)
- **SVG / PDF / EPS depiction**: `mol.to_svg()`, `mol.to_pdf()`, `mol.to_eps()`
- **ChemicalJSON**: `mol.to_cjson(coords=[])`, `chematic.from_cjson(s)` — Avogadro 2 / MolSSI compatible

## License

MIT OR Apache-2.0
