# chematic cookbook — 20 common tasks

Each task is copy-paste-ready. Examples use aspirin (`CC(=O)Oc1ccccc1C(=O)O`) as the reference molecule.

```python
import chematic
```

---

## 1. Get basic properties from a SMILES string

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin

print(mol.mw)              # 180.16
print(mol.logp)            # 1.31
print(mol.tpsa)            # 63.6
print(mol.hbd, mol.hba)   # 1 4
print(mol.qed)             # 0.55
print(mol.formula)         # C9H8O4
```

## 2. Check drug-likeness filters

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

print(mol.lipinski_passes)  # True
print(mol.veber_passes)     # True
print(mol.ghose_passes)     # True
print(mol.pains_passes)     # True
print(mol.brenk_passes)     # True
```

## 3. Predict pKa

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

pka = mol.pka()
print(pka["most_acidic"])   # 3.49  (carboxylic acid)
print(pka["most_basic"])    # None  (no basic site)
```

## 4. Get an ADMET profile

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

profile = mol.admet()
# {
#   "bbb": False,
#   "bbb_score": -1.3,
#   "caco2": -4.9,         # log Papp (Caco-2)
#   "herg_risk": 0.12,     # 0–1 risk score
#   "cyp3a4_risk": 0.21,
# }
```

## 5. Compute fingerprint similarity

```python
aspirin   = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
ibuprofen = chematic.from_smiles("CC(C)Cc1ccc(CC(C)C(=O)O)cc1")

sim = chematic.tanimoto(aspirin.ecfp4(), ibuprofen.ecfp4())
print(f"Tanimoto (ECFP4): {sim:.3f}")  # 0.148
```

## 6. Read an SDF file into a DataFrame

```python
import pandas as pd

records = list(chematic.iter_sdf("library.sdf"))
df = pd.DataFrame({
    "name":   [r.name for r in records],
    "smiles": [r.mol.smiles for r in records],
    "mw":     [r.mol.mw for r in records],
    "logp":   [r.mol.logp for r in records],
    "qed":    [r.mol.qed for r in records],
})
print(df.head())
```

## 7. Compute descriptors for many SMILES in parallel

```python
import pandas as pd

smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "c1cccnc1", "CCCCCCCC"]

# Automatically parallelized across CPU cores
df = pd.DataFrame(chematic.bulk.descriptors(smiles_list))
print(df[["mw", "logp", "tpsa", "qed"]].round(2))
```

## 8. Get an ECFP4 matrix as numpy (for ML)

```python
import numpy as np

smiles_list = ["CCO", "c1ccccc1", "CC(=O)O"]

# shape: (N, 2048), dtype: uint8
X = chematic.bulk.ecfp4(smiles_list)

# Drop directly into scikit-learn
from sklearn.ensemble import RandomForestClassifier
# clf.fit(X, y)
```

## 9. Compute a Tanimoto similarity matrix

```python
smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "c1cccnc1"]

# shape: (N, N), dtype: float32
matrix = chematic.bulk.tanimoto(smiles_list, smiles_list)
print(matrix.round(2))
```

## 10. Virtual screening (similarity search)

```python
library_smiles = [...]  # tens of thousands of molecules is fine

query = "CC(=O)Oc1ccccc1C(=O)O"
scores = chematic.bulk.tanimoto_search(query, library_smiles)  # (N,) float32

# Top 10 hits
top_indices = scores.argsort()[::-1][:10]
for i in top_indices:
    print(f"{library_smiles[i]}: {scores[i]:.3f}")
```

## 11. Fast nearest-neighbour search with LSH

```python
# Designed for large libraries (hundreds of thousands of molecules)
idx = chematic.SimilarityIndex.from_smiles(library_smiles)

hits = idx.search("CC(=O)Oc1ccccc1C(=O)O", threshold=0.5, k=20)
for mol_idx, score in hits:
    print(f"{library_smiles[mol_idx]}: {score:.3f}")
```

## 12. Substructure search with SMARTS

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

# Boolean match
if chematic.smarts_match("[CX3](=O)[OX2H1]", mol):
    print("carboxylic acid found")

# Get matched atom indices
matches = chematic.smarts_find("[CX3](=O)[OX2H1]", mol)
print(matches)  # [[7, 8, 9], ...]
```

## 13. Standardise a molecule

```python
# Remove salts/solvents, neutralise charges, canonicalise tautomers
mol = chematic.from_smiles("[Na+].[O-]c1ccccc1")
clean = mol.standardize()
print(clean.smiles)  # Oc1ccccc1  (phenol)

# Individual operations
mol.largest_fragment()   # keep the largest connected fragment
mol.neutralize()         # neutralise formal charges
mol.remove_isotopes()    # strip isotope labels
mol.remove_stereo()      # remove stereochemistry
```

## 14. Extract the Murcko scaffold

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

scaffold = mol.scaffold()
print(scaffold.smiles)         # c1ccc(CC(=O)O)cc1  (Murcko)

generic = mol.generic_scaffold()
print(generic.smiles)          # generic form (all atoms → C, all bonds → single)
```

## 15. BRICS fragmentation

```python
mol = chematic.from_smiles("CC(=O)Nc1ccc(O)cc1")  # paracetamol

frags = mol.brics_fragments()
for f in frags:
    print(f.smiles)
```

## 16. Enumerate and canonicalise tautomers

```python
mol = chematic.from_smiles("OC1=CC=CC=N1")  # 2-pyridinol

canonical = mol.canonical_tautomer()
print(canonical.smiles)       # O=C1CC=CC=N1

all_tautomers = mol.enumerate_tautomers()
print(len(all_tautomers))     # typically 2–5
```

## 17. Find the Maximum Common Substructure (MCS)

```python
mol1 = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")   # aspirin
mol2 = chematic.from_smiles("CC(=O)Nc1ccc(O)cc1")       # paracetamol

mcs = chematic.find_mcs([mol1, mol2])
if mcs:
    print(mcs.smiles)   # common substructure
```

## 18. Apply a reaction with SMIRKS

```python
# Deprotonate a phenol
phenol = chematic.from_smiles("c1ccccc1O")
products = chematic.run_smirks("[OH:1]>>[O-:1]", [phenol])

for product_set in products:
    for p in product_set:
        print(p.smiles)  # [O-]c1ccccc1
```

## 19. Render SVG in Jupyter

```python
from IPython.display import SVG, display

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

# Single molecule
display(SVG(mol.svg()))

# Highlight atoms matched by a SMARTS query
matches = chematic.smarts_find("[CX3](=O)[OX2H1]", mol)
atoms = [i for match in matches for i in match]
display(SVG(mol.svg_highlighted(atoms, color="#FF6B6B")))

# Grid view
mols = [chematic.from_smiles(s) for s in ["CCO", "c1ccccc1", "CC(=O)O", "CCCC"]]
display(SVG(chematic.depict_grid(mols, cols=2)))
```

## 20. Export all descriptors + filter results to a DataFrame

```python
import pandas as pd

smiles_list = [
    "CC(=O)Oc1ccccc1C(=O)O",         # aspirin
    "CC(C)Cc1ccc(CC(C)C(=O)O)cc1",   # ibuprofen
    "c1ccc2ccccc2c1",                 # naphthalene
    "CCCCCCCCCCCCCCCCCC(=O)O",        # stearic acid
]

df = pd.DataFrame(chematic.bulk.descriptors(smiles_list))
# includes lipinski_passes, pains_passes, brenk_passes, etc.
print(df[["mw", "logp", "tpsa", "qed", "lipinski_passes", "pains_passes"]].round(2))
```

---

## Bonus: InChI / InChIKey

```python
mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")

print(mol.inchi)     # InChI=1S/C9H8O4/...
print(mol.inchikey)  # BSYNRYMUTXBXSQ-UHFFFAOYSA-N

# InChI → Mol
mol2 = chematic.from_inchi(mol.inchi)
assert mol2.smiles == mol.smiles
```
