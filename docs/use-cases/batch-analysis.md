# Batch descriptors and filtering with explicit row outcomes

## Problem

You have a SMILES list or SDF file and need descriptors, fingerprints, or
substructure filters without losing the identity of rejected records.

## Solution

Use streaming SDF APIs when names, properties, source order, and parse failures
must be preserved. Use `bulk.*` for an in-memory SMILES collection. The bulk
paths use parallel Rust implementations where the operation supports them;
actual throughput depends on the operation, host, thread count, input, and
failure policy. Current bulk descriptor and array APIs skip invalid SMILES, so
validate first when source-row accounting matters.

## SDF → filtered CSV

```python
import chematic
import csv

# Stream without loading the whole file; keep one output row per source record.
with open("filtered.csv", "w", newline="") as f:
    fields = [
        "source_index", "status", "error", "name", "smiles",
        "mw", "logp", "tpsa", "qed",
    ]
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()

    for source_index, rec, error in chematic.iter_sdf_checked("library.sdf"):
        if error is not None:
            writer.writerow({
                "source_index": source_index,
                "status": "invalid",
                "error": error,
            })
            continue

        mol = rec.mol
        if not mol.lipinski_passes or not mol.pains_passes:
            writer.writerow({
                "source_index": source_index,
                "status": "filtered",
                "name": rec.name,
                "smiles": mol.smiles,
            })
            continue

        writer.writerow({
            "source_index": source_index,
            "status": "accepted",
            "name":  rec.name,
            "smiles": mol.smiles,
            "mw":    round(mol.mw, 2),
            "logp":  round(mol.logp, 2),
            "tpsa":  round(mol.tpsa, 1),
            "qed":   round(mol.qed, 3),
        })

print("Done")
```

## Parallel descriptors after input validation

```python
import chematic
import pandas as pd

source_smiles = open("smiles.txt").read().splitlines()
valid = []
rejected = []

for source_index, smiles in enumerate(source_smiles):
    try:
        chematic.from_smiles(smiles)
    except chematic.ChematicInputError as error:
        rejected.append({
            "source_index": source_index,
            "smiles": smiles,
            "error": str(error),
        })
    else:
        valid.append((source_index, smiles))

# The bulk call now receives only validated rows, preserving positional mapping.
df = pd.DataFrame(chematic.bulk.descriptors([smiles for _, smiles in valid]))
df.insert(0, "source_index", [source_index for source_index, _ in valid])
df.to_csv("descriptors.csv", index=False)
pd.DataFrame(rejected).to_csv("descriptor_errors.csv", index=False)
```

## Excel workbook → XLSX report

When a compound table is maintained in Excel, use ElixCee for workbook I/O and
chematic for molecular parsing and descriptors. This keeps XLSX/XLSM concerns
out of the chemistry core while avoiding an openpyxl-specific integration.

```bash
pip install "chematic[excel]"
python examples/xlsx_batch_analysis.py compounds.xlsx results.xlsx
```

The example streams the first worksheet by default, auto-detects a `smiles`,
`canonical_smiles`, or `structure` column, and writes typed values to XLSX.
Invalid and blank SMILES are retained with their original worksheet row and an
error. Use `--sheet`, `--smiles-column`, or `--name-column` when the workbook
uses different names. Excel I/O is provided by
[ElixCee](https://github.com/kent-tokyo/elixcee); it is an optional dependency.

## Parallel fingerprint matrix for clustering

The remaining bulk examples assume every input has already passed the same
validation step. Otherwise, preserve a source-index sidecar before calling a
bulk API that skips invalid SMILES.

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

# Uses chematic's documented bounded standardization profile
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

## Measure the operation you will deploy

Do not estimate production throughput from an unrelated descriptor or from a
different package version. The [benchmark method](../benchmark.md) separates
parsing, first use, prepared reuse, serialization, and memory, while the
[record index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks)
retains dated raw results. Re-run the same operation on your corpus and keep
invalid, unsupported, and refused rows in the denominator.

## Related APIs

- `chematic.bulk.descriptors(smiles)` — 72-descriptor batch, returns `list[dict]`
- `chematic.bulk.ecfp4(smiles)` — `(N, 2048)` uint8 fingerprint matrix
- `chematic.bulk.tanimoto(queries, library)` — `(M, N)` float32 similarity matrix
- `chematic.bulk.standardize(smiles)` — parallel salt stripping + tautomer canon
- `chematic.bulk.substructure_search(smarts, smiles)` — parallel VF2
