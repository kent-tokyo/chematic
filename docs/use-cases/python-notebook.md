# Explore molecules in a Python or Jupyter workflow

## Problem

You want to inspect a compound library, compute descriptors, build fingerprint
features, or visualize structures from a Python notebook.

## Solution

Published wheels install with `pip` on supported platforms and do not require a
local C/C++ chemistry toolchain. Use checked or streaming readers when source
row identity and rejected records matter. Descriptor and alert results are
computed properties, not experimental efficacy or safety conclusions.

## Related APIs

- `chematic.descriptors_df(smiles_list)` — one-liner pandas DataFrame of 72 descriptors
- `chematic.bulk.ecfp4(smiles_list)` — `(N, 2048)` uint8 for sklearn / PyTorch
- `chematic.SimilarityIndex` — LSH approximate nearest-neighbour search
- `mol.svg()` / `chematic.depict_grid(mols, cols=3)` — SVG rendering in Jupyter
- [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/kent-tokyo/chematic/blob/main/notebooks/quickstart.ipynb)

## Environment

```bash
pip install chematic pandas numpy scikit-learn matplotlib
```

## 1. Load a compound library from SDF

```python
import chematic
import pandas as pd

records = list(chematic.iter_sdf("library.sdf"))
df = pd.DataFrame({
    "name":  [r.name for r in records],
    "smiles": [r.mol.smiles for r in records],
    "activity": [float(r.get("IC50_nM") or 0) for r in records],
})
print(f"Loaded {len(df)} molecules")
```

## 2. Compute 190+ descriptors in parallel

```python
desc_df = pd.DataFrame(chematic.bulk.descriptors(df["smiles"].tolist()))
df = pd.concat([df, desc_df], axis=1)

# Quick filter
drug_like = df[df["lipinski_passes"] & df["pains_passes"]]
print(f"{len(drug_like)} / {len(df)} pass Lipinski + PAINS")
```

## 3. Generate ECFP4 fingerprints for ML

```python
import numpy as np

X = chematic.bulk.ecfp4(df["smiles"].tolist())   # (N, 2048) uint8
y = (df["activity"] < 100).astype(int)           # binary: IC50 < 100 nM

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

clf = RandomForestClassifier(n_estimators=100, random_state=42)
scores = cross_val_score(clf, X, y, cv=5)
print(f"CV AUC: {scores.mean():.3f} ± {scores.std():.3f}")
```

## 4. Tanimoto similarity clustering

```python
import matplotlib.pyplot as plt

sim_matrix = chematic.bulk.tanimoto(
    df["smiles"].tolist(), df["smiles"].tolist()
)  # (N, N) float32

# Hierarchical clustering
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

dist = 1.0 - sim_matrix
linkage_mat = linkage(squareform(dist), method="ward")
dendrogram(linkage_mat, no_labels=True)
plt.title("Compound cluster (Tanimoto distance)")
plt.show()
```

## 5. Nearest-neighbor search in a large library

```python
idx = chematic.SimilarityIndex.from_smiles(df["smiles"].tolist())

hits = idx.search("CC(=O)Nc1ccc(O)cc1", threshold=0.4, k=10)
for mol_idx, score in hits:
    print(f"{df['name'].iloc[mol_idx]}: Tanimoto {score:.3f}")
```

## 6. Visualize in Jupyter

```python
from IPython.display import SVG, display

mol = chematic.from_smiles("CC(=O)Nc1ccc(O)cc1")
display(SVG(mol.svg()))

# Highlight SMARTS match
matches = chematic.smarts_find("[NH]C(=O)", mol)
atoms = [i for m in matches for i in m]
display(SVG(mol.svg_highlighted(atoms, color="#FF6B6B")))

# Grid of top hits
top_mols = [chematic.from_smiles(df["smiles"].iloc[i]) for i, _ in hits[:6]]
display(SVG(chematic.depict_grid(top_mols, cols=3)))
```
