# `chematic.rdkit_compat` — RDKit compatibility

`chematic.rdkit_compat` is a **lightweight RDKit-compatible subset**, not a full RDKit
clone. It lets common 2D cheminformatics scripts run with minimal changes in environments
where RDKit is unavailable (WASM, serverless, pure-Rust services, conda-free Python).

**Unsupported options fail loudly** — they raise `NotImplementedError` / `TypeError` /
`ValueError` instead of being silently ignored.

```python
from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import Descriptors, rdMolDescriptors, DataStructs

mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
Descriptors.MolWt(mol)                       # 180.16
fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
DataStructs.TanimotoSimilarity(fp, fp)       # 1.0
```

---

## Compatibility matrix

| Area | Status | Notes |
|------|--------|-------|
| SMILES I/O | ✅ Supported | `MolFromSmiles` (aromaticity perceived when `sanitize=True`), `MolToSmiles`, `MolFromMolBlock`/`MolToMolBlock` |
| SDF I/O | ✅ Supported | `SDMolSupplier` (iterate, `len`, `[i]`, context manager), `SDWriter` + SD properties |
| SMILES file I/O | ✅ Supported | `SmilesMolSupplier` (title-line property columns, `[i]`), `SmilesWriter` (`SetProps`) |
| Mol properties | ✅ Supported | `Get`/`Set`/`Has`/`Clear`Prop, `SetInt/Double/BoolProp`, `GetPropsAsDict`, `GetPropNames` |
| Mol / Atom / Bond | ✅ Supported | read-only traversal: `GetAtoms`/`GetBonds`/`GetAtomWithIdx`/`GetBondWithIdx`; atom/bond getters; `BondType` |
| RingInfo | ✅ Supported | `GetRingInfo()` → `NumRings`/`AtomRings`/`BondRings`/`NumAtomRings`/`NumBondRings` (SSSR-based) |
| Descriptors | ✅ Supported | MW/HBA/HBD **exact**, TPSA ±1.0, LogP ±0.5 vs RDKit (differential-tested) |
| Aromaticity | ✅ Supported | aromatic atom/bond counts match RDKit on 99.0% / 98.3% of a 5k corpus; ring-junction carbonyls differ |
| Substructure (SMARTS) | 🟡 Partial | match **sets** agree 96.9% on a 5k corpus; ring-size `[rN]` queries differ (SSSR vs RDKit ring perception). Match **order** may differ — compare as sets |
| Morgan fingerprint | 🟡 Partial | `radius`, `nBits` (modulo folding), `bitInfo` shape-/origin-consistent — **not RDKit bit-identical** (FNV-1a vs MurmurHash) |
| DataStructs | ✅ Supported | `TanimotoSimilarity`/`DiceSimilarity`/`BulkTanimotoSimilarity`/`ConvertToNumpyArray` |
| Canonical SMILES | 🟡 Partial | 99.62% semantic round-trip vs RDKit (5k corpus); exocyclic C=N E/Z stereo not always emitted |
| RWMol / structure editing | ❌ Unsupported | read-only layer |
| `useFeatures=True`, `useBondTypes=False`, `nBits<=0`, `bitInfo+useChirality` | 🔊 Fails loudly | raise instead of silently ignoring |

---

## Differential validation

A live differential suite compares chematic against real RDKit and writes explainable
JSONL diffs. Tests auto-skip when RDKit is absent (`pytest.importorskip`), so CI without
RDKit is unaffected.

| Suite | Script | Artifact | Headline (corpus) |
|-------|--------|----------|-------------------|
| descriptors / ring / SMARTS-count / SDF / Morgan | `crates/chematic-py/tests/test_rdkit_diff.py` | `validation/results/rdkit_diff.jsonl` | MW/HBA/HBD exact; ring count exact |
| canonical SMILES round-trip + idempotency | `scripts/canonical_diff.py` | `validation/results/canonical_diff.jsonl` | 100% RDKit-parseable, 99.62% round-trip, 98.42% idempotent (5k) |
| SMARTS match sets + aromaticity counts | `scripts/rdkit_compat_diff.py` | `validation/results/rdkit_compat_diff.jsonl` | 96.9% SMARTS match-set, 99.0%/98.3% aromatic atom/bond counts (5k) |

Run the standalone differentials (RDKit required):

```bash
python scripts/canonical_diff.py --limit 2000
python scripts/rdkit_compat_diff.py --limit 2000
```

### Known divergence classes

- **Ring-size SMARTS (`[r5]`, `[r6]`, …)** — chematic SSSR ring membership differs from
  RDKit's for some fused systems (the dominant SMARTS divergence, ~93% of mismatches).
- **Aromatic ring-junction carbonyls** — atoms like `c(=O)` in fused aromatics are
  modeled differently, shifting a few `c` / `C=O` matches and aromatic counts.
- **Morgan bit positions** — FNV-1a (chematic) vs MurmurHash (RDKit); similarity ranking
  is consistent, individual bit indices are not comparable across libraries.
- **Exocyclic C=N E/Z** — canonical SMILES does not always emit the directional stereo.

---

## Examples

Runnable, self-contained task scripts live in [`examples/rdkit_compat/`](../examples/rdkit_compat/):

- `sdf_to_csv.py` — SDF → descriptor table (CSV)
- `substructure_filter.py` — filter a SMILES list by a SMARTS query
- `fingerprint_similarity.py` — Morgan fingerprint + Tanimoto ranking
- `numpy_features.py` — fingerprint matrix for scikit-learn / PyTorch

A 10-pattern RDKit→chematic migration guide is in
[`examples/rdkit_compat_migration.py`](../examples/rdkit_compat_migration.py).
