# Quick Start

## Your first molecule

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
print(mol.mw)       # 180.16
print(mol.formula)  # C9H8O4
print(mol.logp)     # 1.31
print(mol.tpsa)     # 63.6
```

## Drug-likeness in one line

```python
print(mol.lipinski_passes)   # True
print(mol.pains_passes)      # True
```

## pKa and ADMET

```python
print(mol.pka())    # {"most_acidic": 3.49, "most_basic": None}
print(mol.admet())  # {"bbb": False, "caco2": -4.9, "herg_risk": 0.12, ...}
```

## Fingerprint similarity

```python
aspirin   = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
ibuprofen = chematic.from_smiles("CC(C)Cc1ccc(CC(C)C(=O)O)cc1")

sim = chematic.tanimoto(aspirin.ecfp4(), ibuprofen.ecfp4())
print(f"{sim:.3f}")   # 0.148
```

## Substructure search

```python
if chematic.smarts_match("[CX3](=O)[OX2H1]", mol):
    print("has carboxylic acid")

matches = chematic.smarts_find("[CX3](=O)[OX2H1]", mol)
print(matches)   # [[7, 8, 9]]
```

## Parallel descriptors for many molecules

```python
import pandas as pd

smiles = ["CCO", "c1ccccc1", "CC(=O)O", "CCCC"]
df = pd.DataFrame(chematic.bulk.descriptors(smiles))
print(df[["mw", "logp", "tpsa", "qed"]])
```

## SVG in Jupyter

```python
from IPython.display import SVG, display
display(SVG(mol.svg()))
```

## Next steps

- [Cookbook](../cookbook.md) — 20 ready-to-run recipes
- [RDKit migration guide](../rdkit_cheatsheet.md) — coming from RDKit?
- [API Reference](../api/chematic.md) — full API docs
