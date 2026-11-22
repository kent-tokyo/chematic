# chematic

Pure-Rust cheminformatics library for Python — SMILES parsing, descriptors,
fingerprints, substructure search, formats, reactions, and bounded 3D tools.

See [PyPI](https://pypi.org/project/chematic/) for available versions and
the [validation report](https://kent-tokyo.github.io/chematic/reference/validation/)
for package results.

## Installation

```bash
pip install chematic
```

## Quick Start

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
print(mol.mw, mol.tpsa)
print(chematic.tanimoto(mol.ecfp4(), chematic.from_smiles("CCO").ecfp4()))
svg = mol.to_svg()
```

See the [Python notebook guide](https://kent-tokyo.github.io/chematic/use-cases/python-notebook/)
for descriptors and batch analysis and the
[API reference](https://kent-tokyo.github.io/chematic/api/chematic/) for molecule
methods, indexed SDF, pickle and optional workflows. See
[Cleanup migration](../../docs/reference/cleanup-migration.md) for error behavior
in the unreleased source.

## Molecular GNN input

`Mol.ogb_graph_features()` returns the OGB categorical graph schema as NumPy
`int64` arrays. It works with molecules from SMILES, SDF records, and MOL2.
The keys are `node_feat` `(N, 9)`, `edge_index` `(2, 2E)`, `edge_feat`
`(2E, 3)`, and `num_nodes`. The edges and their features have the same order,
with one edge in each direction per bond. PyTorch and PyG are optional:

```python
import chematic

for source_index, record, error in chematic.iter_sdf_checked("ligands.sdf"):
    if error is not None:
        print(f"Rejected SDF record {source_index}: {error}")
        continue
    assert record is not None
    graph = record.ogb_graph_features()
    # Attach record.name and record.properties() as labels or IDs as needed.
    # With PyTorch Geometric installed:
    # import torch
    # from torch_geometric.data import Data
    # data = Data(
    #     x=torch.from_numpy(graph["node_feat"]),
    #     edge_index=torch.from_numpy(graph["edge_index"]),
    #     edge_attr=torch.from_numpy(graph["edge_feat"]),
    #     num_nodes=graph["num_nodes"],
    # )
```

The default removes RDKit-removable hydrogen graph atoms, following OGB's
default SMILES conversion. Pass `remove_hs=False` to retain the atoms already
present in the molecule; this option does not add hydrogens. Unsupported
chemistry raises `ValueError`. This export contains categorical graph features,
not coordinates or SD properties. `iter_sdf_checked` reports a rejected record
at its source index. Its default reader limits are 1 GiB and 100,000 records;
pass positive `max_input_bytes` and `max_records` overrides for larger files.
Use `SdfRecord.ogb_graph_features()` for V2000 or V3000 SDF input: it retains
explicitly unknown double-bond stereo (V2000 code 3 or V3000 `CFG=2`). The generic
`Mol.ogb_graph_features()` has no record-level stereo diagnostic to read.
For a 3D SDF, chirality inferred by RDKit from coordinates can differ from
the chirality retained in `record.mol` from the CTAB wedge representation.

## Features

- Read and write common molecule formats, including indexed SDF. See the
  [format matrix](../../docs/reference/format-capabilities.md).
- Compute descriptors and fingerprints, run SMARTS and reaction workflows,
  and render SVG, PDF or EPS depending on build features.
- Use opt-in 3D embedding and prepared MMFF94/UFF force fields. The
  [validation report](../../docs/reference/validation.md) describes
  their tested limits.

## RDKit compatibility

The operation/profile boundaries and measured oracle lanes are listed in the
[Compatibility Contract dashboard](../../docs/reference/compatibility-dashboard.md).
For persisted fingerprints, saved indices, canonical identity, stereo, and
browser/Worker migration decisions, read the
[RDKit migration guide](../../docs/reference/rdkit-migration.md) before changing a
production workflow. Native ECFP and RDKit-compatible Morgan are separate
profiles; rebuild stored fingerprints and indices under the chosen profile.

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

# Explain RDKit atom-pair / topological-torsion codes
code = mol.rdkit_atom_code(0)
chematic.explain_rdkit_atom_code(code)  # (symbol, branches, pi_electrons)
path_code = mol.rdkit_topological_torsion_path_code([0, 1, 2, 3])
chematic.explain_rdkit_topological_torsion_code(path_code)

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

# Large seekable SDF: index once, parse records only when requested
from chematic import SdfDataset
dataset = SdfDataset.open("out.sdf")
print(len(dataset), dataset.metadata(0).title)
print(dataset[-1].smiles)
```

Unsupported options raise `NotImplementedError` or `TypeError` — they are never silently ignored.

### Compatibility matrix

| Area | Status | Notes |
|------|--------|-------|
| SMILES I/O | ✅ Supported | `MolFromSmiles` (aromaticity perceived when `sanitize=True`) / `MolToSmiles` |
| SDF I/O | ✅ Supported | `SDMolSupplier` / `SDWriter` + SD properties |
| Mol properties | ✅ Supported | `GetProp(..., default=...)`, `Set`/`Has`/`Clear`Prop, typed getters/setters, `GetPropsAsDict` |
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
[RDKit migration guide](https://github.com/kent-tokyo/chematic/blob/main/docs/reference/rdkit-migration.md)
(compatibility matrix, differential-validation results, known divergences, and runnable examples).

## License

MIT OR Apache-2.0
