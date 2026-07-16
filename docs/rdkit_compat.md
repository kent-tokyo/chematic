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
| Substructure (SMARTS) | 🟡 Partial | match **sets** agree 96.9% on a 5k corpus (frozen benchmark — SMARTS-R0 found this figure predates a since-landed SSSR fix, current figure is somewhat higher; not yet regenerated, see below); `[rN]` predicate-semantics bug found recoverable, not a ring-model gap (see "SMARTS-R0" below). Match **order** may differ — compare as sets |
| Morgan fingerprint | 🟡 Partial | `radius`, `nBits` (modulo folding), `bitInfo`, `useFeatures=True` (FCFP) — shape-/origin-consistent, **not RDKit bit-identical** (FNV-1a vs MurmurHash) |
| DataStructs | ✅ Supported | `TanimotoSimilarity`/`DiceSimilarity`/`BulkTanimotoSimilarity`/`ConvertToNumpyArray` |
| Canonical SMILES | 🟡 Partial | 99.62% semantic round-trip vs RDKit (5k corpus); exocyclic C=N E/Z stereo not always emitted |
| RWMol / structure editing | 🟡 Partial | `AddAtom`/`AddBond`/`RemoveAtom`/`RemoveBond`/`GetMol` supported; no mid-edit atom/bond iteration (call `GetMol()` first) |
| `useBondTypes=False`, `nBits<=0`, `bitInfo+useChirality`, `useFeatures+useChirality` | 🔊 Fails loudly | raise instead of silently ignoring |

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
| Ring-size SMARTS `[rN]` (`[r5]`, `[r6]`, …) — ~94% of SMARTS mismatches | **Corrected 2026-07-16 (SMARTS-R0), supersedes the row below this table previously stated.** NOT a ring-model difference — a scoped predicate-semantics bug. `crates/chematic-smarts/src/parser.rs` maps `[rN]` to the *same* `AtomPrimitive::RingSize(n)` as `[kN]`, and `match_vf2.rs`'s `RingSize` arm matches "any SSSR ring of size N contains this atom." RDKit's real `[rN]` means "this atom's *smallest* ring is exactly size N" — a materially different predicate, confirmed empirically distinct from `[kN]` (RDKit `[k6]` includes a fusion atom on a purine that `[r6]` excludes for the same atom). Chematic's own SSSR is not the driver: ring count and size multiset agree between chematic and RDKit on 98.8% of the previously-mismatched molecules. See "SMARTS-R0" below for the full diagnosis. | **Recoverable without a ring-model change** — measured, not estimated: min-ring-size semantics on chematic's existing SSSR takes agreement from 70.95% to 100.00% on the full 5,000-molecule corpus (two independent measurements). Needs `[rN]` split from `[kN]`'s primitive (they are NOT the same predicate in RDKit) — not attempted this round, diagnosis only. |
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

### SMARTS-R0 — ring-semantics mismatch diagnosis (2026-07-16, diagnosis only)

Reclassification of `scripts/rdkit_compat_diff.py`'s SMARTS match-set mismatches, per the
same 16-pattern, per-`(molecule, pattern)`-pair, set-comparison methodology the frozen
`96.9%` figure already uses. **No production code changed** — this section reports
recoverability and cost for a future fix, it does not implement one.

**Bucket breakdown**, from the committed `validation/results/rdkit_compat_diff.jsonl`
(2,510 mismatched pairs / 2,006 unique molecules):

| Bucket | Pairs | Share |
|---|---|---|
| `[r5]`/`[r6]` ring-size queries | 2,370 | 94.4% |
| Aromatic ring-junction carbonyls (`c`, `C=O`, `[#6]=[#6]` tail) | 140 | 5.6%, already documented above, same ~47-molecule root cause across all three patterns |
| `[R]`/`[R0]`/`[R1]` ring-membership-count | 0 in this corpus | not tested at all by the current 16-pattern list |

By RDKit ring topology of the mismatched molecules (fused/bridged/spiro, computed from
`GetRingInfo().AtomRings()`, not guessed from SMILES): fused systems dominate (1,793 of
2,006), which is exactly where an atom sits in rings of two different sizes — the
precondition for the semantics bug below to matter at all. Pure bridged (16) and pure
spiro (4) are rare in this corpus.

**Root cause, confirmed against real RDKit 2026.03.3, not assumed**: this is *not* an
SSSR ring-decomposition disagreement. Ring count and ring-size multiset agree between
chematic and RDKit on 98.8% of the mismatched molecules (1,982/2,006) — chematic's own
SSSR choice is not the driver. The actual bug: `crates/chematic-smarts/src/parser.rs`
(~line 847, ~line 858) parses both `[rN]` and `[kN]` into the identical
`AtomPrimitive::RingSize(n)`, and `crates/chematic-smarts/src/match_vf2.rs` (~line 315)
evaluates it as "this atom belongs to *any* SSSR ring of size N." RDKit's actual `[rN]`
semantics is "this atom's *smallest* ring is exactly size N" — genuinely different from
`[kN]`'s any-ring semantics, confirmed with a discriminating test: on a purine
(`Nc1nc(N)c2ncn(...)c2n1`), RDKit's `[k6]` matches the ring-fusion atom that RDKit's
`[r6]` explicitly excludes, for the identical atom, identical molecule. Existing unit
tests (`parser.rs`, `test_parse_k_exact_ring_size`, `test_k_vs_r_equivalent_for_cyclopentane`)
assert `[k5]==[r5]`/`[k6]==[r6]` — they encode this bug (the second test happens to pass
today only because cyclopentane has one ring, where min-ring-size and any-ring-size
coincide).

**Recoverability, measured twice independently on the full 5,000-molecule corpus, not
estimated**: computing each atom's *minimum* SSSR ring size (using chematic's existing,
unchanged SSSR — no ring-model or `find_sssr` change) instead of "any ring of size N"
takes `[rN]` agreement from **70.95% (3,547/4,999 or 3,548/5,000, both measurements
agree) to 100.00% (4,999/4,999 and 5,000/5,000)**. Zero regressions, zero remaining
exceptions, either measurement.

**Performance cost**: chematic-smarts already computes SSSR once per matching context
and passes it by reference (`match_vf2.rs`'s `rings: &RingSet`) — a min-ring-size
predicate aggregates over the *same* already-computed ring set the current any-ring
check already iterates; no new ring-finding, no new SSSR computation. Not expected to
be measurably slower than today's `[rN]` evaluation.

**What this recoverability figure does *not* cover**: `[R1]`/`[R2]` ring-*count*
predicates (not in the current 16-pattern test corpus at all) show a smaller, genuinely
different residual — 17/1,000 (1.7%) on a spot sample, traced to bridged/cage systems
where RDKit's and chematic's SSSR *bases* disagree in cardinality (e.g. an
adamantane-type cage: RDKit finds 5 SSSR rings, chematic finds 4 — a real basis-choice
disagreement, the same phenomenon this doc's morphinan example already noted). This
bucket is **not** fixed by the `[rN]` min-ring-size change and would need an actual
richer, RDKit-compatible symmetrized ring model, scoped to SMARTS ring-count predicates
only, per this doc's original proposal — kept as a separate, smaller, harder problem,
not conflated with the `[rN]` fix above.

**Bottom line**: the "Won't fix — genuinely different ring models" verdict this doc
previously gave `[rN]` was not supported by the evidence once actually measured. The
dominant SMARTS mismatch driver (94.4% of pairs) is a small, scoped predicate-semantics
bug recoverable without touching the ring model at all. The frozen 96.9% benchmark also
predates a since-landed SSSR fix (`698ba3f`, Horton's algorithm) — a fresh partial
sample measured 97.50%, higher than the frozen figure; the full corpus wasn't
regenerated this round (a multi-minute run, kept out of scope for a diagnosis-only
pass) — regenerating it is recommended as part of any future fix work, alongside
correcting the `[k5]==[r5]`/`[k6]==[r6]` unit tests that currently encode the bug.

---

## Examples

Runnable, self-contained task scripts live in [`examples/rdkit_compat/`](../examples/rdkit_compat/):

- `sdf_to_csv.py` — SDF → descriptor table (CSV)
- `substructure_filter.py` — filter a SMILES list by a SMARTS query
- `fingerprint_similarity.py` — Morgan fingerprint + Tanimoto ranking
- `numpy_features.py` — fingerprint matrix for scikit-learn / PyTorch

A 10-pattern RDKit→chematic migration guide is in
[`examples/rdkit_compat_migration.py`](../examples/rdkit_compat_migration.py).
