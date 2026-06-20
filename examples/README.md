# chematic — Python Examples

Runnable Python scripts showing chematic in common ML/cheminformatics workflows.

## Setup

```bash
pip install chematic scikit-learn pandas numpy
```

## Scripts

### `qsar_sklearn.py` — QSAR with scikit-learn

Builds a binary classifier (BBB penetration toy dataset) using ECFP4 fingerprints and RandomForest.
Demonstrates that `bulk.ecfp4()` → `(N, 2048) uint8 ndarray` plugs directly into sklearn with zero conversion.

```
python examples/qsar_sklearn.py
```

Key patterns:

```python
import chematic
from sklearn.ensemble import RandomForestClassifier

fps = chematic.bulk.ecfp4(smiles_list)   # (N, 2048) uint8 ndarray — sklearn-ready
clf = RandomForestClassifier().fit(fps, labels)

maccs = chematic.bulk.maccs(smiles_list)  # (N, 166) uint8 ndarray
```

### `descriptors_pandas.py` — Descriptors → pandas → ML feature matrix

Computes 55+ molecular descriptors in one call, converts to a pandas DataFrame,
applies Lipinski filtering, ranks by QED, and feeds descriptors into sklearn PCA.

```
python examples/descriptors_pandas.py
```

Key patterns:

```python
import pandas as pd
import chematic

descs = chematic.bulk.descriptors(smiles_list)  # list[dict]
df = pd.DataFrame(descs)                        # one-liner → DataFrame

# Lipinski filter
ro5 = df[(df["mw"] <= 500) & (df["logp"] <= 5) & (df["hbd"] <= 5) & (df["hba"] <= 10)]

# Pairwise Tanimoto matrix
fps = [chematic.from_smiles(s).ecfp4() for s in smiles_list]
sim = chematic.tanimoto_matrix(fps, fps)  # (N, N) list[list[float]]
```

## What chematic returns (no conversion needed)

| API | Return type | sklearn-compatible? |
|-----|-------------|---------------------|
| `bulk.ecfp4(smiles)` | `(N, 2048)` `uint8` ndarray | Yes |
| `bulk.maccs(smiles)` | `(N, 166)` `uint8` ndarray | Yes |
| `bulk.descriptors(smiles)` | `list[dict]` → `pd.DataFrame` | Yes (after DataFrame) |
| `mol.ecfp4()` | `bytes` (256 bytes = 2048 bits) | Via `np.frombuffer` |
| `tanimoto_matrix(fps_a, fps_b)` | `list[list[float]]` | Via `np.array` |
