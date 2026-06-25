# Use case: Large-scale batch analysis

## Problem

You have a SMILES list or SDF file with thousands of compounds and need to compute descriptors, generate fingerprints, or filter by drug-likeness criteria — fast, without spinning up a conda environment or waiting minutes for RDKit to process serially.

## Solution

chematic's `bulk.*` API uses Rayon internally to parallelise work across all CPU cores. A 10 k-molecule SDF → CSV pipeline runs in under a second on modern hardware. Install with `pip install chematic`; no C++ toolchain needed.

## Output / What you get

```
$ python batch.py library.sdf
Computed 10 000 molecules, 72 descriptors each → descriptors.csv
3 847 / 10 000 pass Lipinski + PAINS
Done in 0.84 s
```

## SDF → filtered CSV

```python
import chematic
import csv

# Stream through SDF without loading all molecules into memory at once
with open("filtered.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["name", "smiles", "mw", "logp", "tpsa", "qed"])
    writer.writeheader()

    for rec in chematic.iter_sdf("library.sdf"):
        mol = rec.mol
        if not mol.lipinski_passes:
            continue
        if not mol.pains_passes:
            continue
        writer.writerow({
            "name":  rec.name,
            "smiles": mol.smiles,
            "mw":    round(mol.mw, 2),
            "logp":  round(mol.logp, 2),
            "tpsa":  round(mol.tpsa, 1),
            "qed":   round(mol.qed, 3),
        })

print("Done")
```

## Parallel descriptors for a SMILES list

```python
import chematic
import pandas as pd

smiles = open("smiles.txt").read().splitlines()

# Parallelised across all CPU cores; returns list[dict]
df = pd.DataFrame(chematic.bulk.descriptors(smiles))
df.to_csv("descriptors.csv", index=False)
print(f"Computed {len(df)} molecules, {df.columns.size} descriptors each")
```

## Parallel fingerprint matrix for clustering

```python
import chematic
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

smiles = open("smiles.txt").read().splitlines()

# (N, 2048) uint8 in one call, parallelised
X = chematic.bulk.ecfp4(smiles)

# Tanimoto matrix (N, N) float32
T = chematic.bulk.tanimoto(smiles, smiles)

# Cluster at 0.6 similarity cutoff
Z = linkage(squareform(1.0 - T), method="ward")
labels = fcluster(Z, t=0.4, criterion="distance")
print(f"{labels.max()} clusters found")
```

## Parallel standardisation

```python
import chematic

raw_smiles = open("raw.smi").read().splitlines()

# Removes salts, neutralises charges, canonicalises tautomers — all in parallel
clean_mols = chematic.bulk.standardize(raw_smiles)

for mol in clean_mols:
    if mol is not None:
        print(mol.smiles)
```

## Parallel 3D generation

```python
import chematic

smiles = ["CCO", "c1ccccc1", "CC(=O)O", "c1cccnc1", "CCCC"]

# Returns list[(Mol, coords)] where coords is a flat list of (x, y, z) triples
results = chematic.bulk.generate_3d(smiles, method="etkdg")

for mol, coords in results:
    if mol is not None:
        print(mol.smiles, len(coords) // 3, "atoms")
```

## Substructure filter over a large library

```python
import chematic

smiles = open("smiles.txt").read().splitlines()

# Returns bool list — True for matches; parallelised
has_carboxylic = chematic.bulk.substructure_search("[CX3](=O)[OX2H1]", smiles)

hits = [s for s, match in zip(smiles, has_carboxylic) if match]
print(f"{len(hits)} molecules contain a carboxylic acid")
```

## Performance reference

| Task | 10 k molecules | 100 k molecules |
|------|---------------|-----------------|
| ECFP4 fingerprints | ~36 ms | ~360 ms |
| 70-descriptor batch | ~80 ms | ~800 ms |
| Tanimoto matrix (N×N) | ~200 ms | — |
| Standardisation | ~120 ms | ~1.2 s |

Measured on 8-core Apple M2. All tasks scale linearly with core count.

## Related APIs

- `chematic.bulk.descriptors(smiles)` — 72-descriptor batch, returns `list[dict]`
- `chematic.bulk.ecfp4(smiles)` — `(N, 2048)` uint8 fingerprint matrix
- `chematic.bulk.tanimoto(queries, library)` — `(M, N)` float32 similarity matrix
- `chematic.bulk.standardize(smiles)` — parallel salt stripping + tautomer canon
- `chematic.bulk.substructure_search(smarts, smiles)` — parallel VF2
