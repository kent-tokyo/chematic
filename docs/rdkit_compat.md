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

Each class below was reproduced and root-caused (not assumed). The scope decision states
why it is or isn't pursued.

| Divergence | Root cause | Scope |
|------------|-----------|-------|
| Ring-size SMARTS (`[r5]`, `[r6]`, …) — ~93% of SMARTS mismatches | chematic uses the **minimal SSSR**; RDKit uses a richer **symmetrized-SSSR** ring model. The gap is **non-monotonic** (purine: chematic 12 > RDKit 10; morphinan: chematic 12 < RDKit 17) — genuinely different ring models, not a tweak. `augmented_ring_set` recovers *smaller* XOR rings, the opposite of what bridged-cage cases need. | **Won't fix** — niche query; aligning means reverse-engineering RDKit's ring model on the core that gives 100% aromatic-ring-count + exact descriptors. Compare SMARTS results as **sets**, and avoid `[rN]` for cross-library parity. |
| Aromatic ring-junction carbonyls (`c(=O)` in fused aromatics) | chematic and RDKit model these atoms differently, shifting a few `c` / `C=O` matches and aromatic atom/bond counts. | Documented; small tail (~1–3%). |
| Morgan bit positions | FNV-1a (chematic) vs MurmurHash (RDKit). | By design — similarity **ranking** is consistent; individual bit indices are not comparable across libraries. |
| Exocyclic C=N E/Z in canonical SMILES (~0.4% round-trip) | The **parser** drops `/`,`\` directional bonds that flank an aromatic ring atom during aromatization (`crates/chematic-smiles/src/parser.rs`), *before* the canonical writer runs — so the geometry is already gone by write time. | **Deferred** — a parser + aromaticity change (broad blast radius), not a writer fix. |
| Canonical idempotency on large fused polycyclics (~1.6%) | **Aromaticity-perception round-trip inconsistency** — *not* Morgan-rank tie-breaking (the failing molecules have all-distinct ranks). A molecule and the re-parse of its own canonical SMILES can disagree on which bonds are aromatic (e.g. 16 vs 17 on a fluorene/carbazole-type linkage); because Morgan ranks weight aromatic vs single bonds differently (`bond_order_value`), this shifts the canonical atom order and the emitted string. The molecule is preserved (**InChI invariant**); only the representation differs. | **Deferred** — the fix is in the aromaticity/parser core (which delivers the 100% aromatic-ring-count + exact descriptors), gated on full descriptor regression; not a canonicalizer patch. Pairs with the exocyclic-C=N parser work above. |

**The `/`,`\` direction choice itself is deterministic and idempotent on stable skeletons**
(verified on a 12-molecule E/Z corpus — see `crates/chematic-smiles/src/canonical.rs`
tests and `tests/test_canonical_diff.py`). Canonical-SMILES E/Z is therefore reliable for the
overwhelming majority of molecules; the residue is confined to the two named structural
classes above.

**The canonical Morgan ranking itself is order-invariant for distinct-rank molecules** — the
residual idempotency failures are upstream in aromaticity perception, not in the ranking or the
writer. Canonical SMILES is idempotent for the ~98.4% of a 5k corpus whose aromaticity round-trips
consistently (guarded by `fused_aromatic_canonical_is_idempotent` and
`test_aromaticity_roundtrip_consistent`).

---

## Examples

Runnable, self-contained task scripts live in [`examples/rdkit_compat/`](../examples/rdkit_compat/):

- `sdf_to_csv.py` — SDF → descriptor table (CSV)
- `substructure_filter.py` — filter a SMILES list by a SMARTS query
- `fingerprint_similarity.py` — Morgan fingerprint + Tanimoto ranking
- `numpy_features.py` — fingerprint matrix for scikit-learn / PyTorch

A 10-pattern RDKit→chematic migration guide is in
[`examples/rdkit_compat_migration.py`](../examples/rdkit_compat_migration.py).
