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
| MOL/SDF 2D wedge/hash stereo | ✅ Supported | `MolFromMolBlock`/`SDMolSupplier` (V2000 and V3000) perceive tetrahedral local parity (`Atom.chirality`) from wedge/hash bonds automatically, matching RDKit's own unconditional `assignChiralTypesFromBondDirs` — CIP-independent, never gated behind `sanitize`. Malformed/contradictory wedges never guess: see `chematic.from_mol_block_with_diagnostics`/`from_mol_v3000_with_diagnostics` for the typed rejection reasons. |
| MOL/SDF 2D double-bond E/Z | ✅ Supported (P1-S2) | `MolFromMolBlock`/`SDMolSupplier` derive SMILES `/`/`\` direction from 2D geometry automatically, mirroring RDKit's `setDoubleBondNeighborDirections` — CIP-independent, never gated on the legacy CIP-label engine. Fail-closed on: terminal alkenes, symmetric substituents, missing/non-finite/collinear coordinates, zero-length double bonds, MDL "either" stereo (V2000 code 3 / V3000 `CFG=2`), cumulenes/allenes, Kekulized aromatic-ring bonds (never stereogenic — see below), and any shared-carrier conflict (issue #149's joint-resolution problem is explicitly out of scope; a shared bond between two independent double bonds is only used when both agree, verified empirically for the ordinary conjugated-diene case). Broad-corpus validation (4,999 molecules): 622 RDKit-resolved double bonds, 276 bond-level semantic agreements, 346 abstentions, **0 semantic inversions, 0 false positives**. |
| SMILES file I/O | ✅ Supported | `SmilesMolSupplier` (title-line property columns, `[i]`), `SmilesWriter` (`SetProps`) |
| Mol properties | ✅ Supported | `Get`/`Set`/`Has`/`Clear`Prop, `SetInt/Double/BoolProp`, `GetPropsAsDict`, `GetPropNames` |
| Mol / Atom / Bond | ✅ Supported | read-only traversal: `GetAtoms`/`GetBonds`/`GetAtomWithIdx`/`GetBondWithIdx`; atom/bond getters; `BondType` |
| RingInfo | ✅ Supported | `GetRingInfo()` → `NumRings`/`AtomRings`/`BondRings`/`NumAtomRings`/`NumBondRings` (SSSR-based) |
| Descriptors | ✅ Supported | MW/HBA/HBD **exact**, TPSA ±1.0, LogP ±0.5 vs RDKit (differential-tested) |
| Aromaticity | ✅ Supported | aromatic atom/bond counts match RDKit on 99.44% / 98.82% of a 5k corpus; one bridgehead-N fused-ring scaffold over-aromatizes under `apply_aromaticity("rdkit_like")`, see "SMARTS-A0" below |
| Substructure (SMARTS) | ✅ Supported | match **sets** agree **99.93%** on a 5k corpus (up from 96.9% — SMARTS-R1 fixed a `[rN]`/`[kN]` predicate-semantics bug that was 94.4% of all mismatches, see "SMARTS-R0"/"SMARTS-R1" below); residual is 100% one bridgehead-N ring-fusion aromaticity over-extension scaffold, see "SMARTS-A0" below (previously mislabeled "carbonyl" — it isn't). Match **order** may differ — compare as sets. **Opt-in** `find_matches_rdkit_parity` (`chematic_smarts::rdkit_parity_match`, default matcher unchanged) additionally fixes the `[R1]`/`[R2]` ring-count residual on bridged/cage systems — 155,650/155,651 (99.9994%) on a 155,651-cell corpus, 0 regressions, see "SMARTS-R2" below |
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
| SMARTS match sets + aromaticity counts | `scripts/rdkit_compat_diff.py` | `validation/results/rdkit_compat_diff.jsonl` | 99.93% SMARTS match-set (SMARTS-R1, was 96.9%), 99.44%/98.82% aromatic atom/bond counts (5k) |
| MOL/SDF 2D wedge/hash reader integration | `crates/chematic-mol/examples/stereo2d_reader_integration_fixture_dump.rs` + `scripts/stereo2d_reader_diagnosis.py` | `validation/results/stereo2d_reader_fixture_dump.jsonl` + `_diagnosis_summary.json` | 12/13 fixtures agree on both raw local parity and accurate-CIP R/S (13th is a documented, characterized divergence — see `docs/stereo2d_reader_integration_rfc.md`) |
| MOL/SDF 2D E/Z direction reader integration (P1-S2) | `crates/chematic-mol/examples/stereo2d_ez_reader_fixture_dump.rs` + `scripts/stereo2d_ez_reader_diagnosis.py` (fixture-level, InChI `/b`-layer semantic comparison, negative-control-verified) + `scripts/stereo2d_ez_corpus_diagnosis.py` (broad-corpus, 4,999 molecules) | `validation/results/stereo2d_ez_reader_fixture_dump.jsonl` + `_diagnosis_summary.json`, `stereo2d_ez_corpus_diagnosis_summary.json` | 14/14 fixtures semantically agree; corpus: 622 RDKit-resolved bonds, 276 agreements, 346 abstentions, 0 inversions, 0 false positives |

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
| ~~Ring-size SMARTS `[rN]`~~ — **Fixed 2026-07-16 (SMARTS-R1)**, removed from this table | Was a predicate-semantics bug (`[rN]` wrongly aliased to `[kN]`'s any-ring semantics; RDKit's real `[rN]` is smallest-ring semantics), not a ring-model gap — see "SMARTS-R0"/"SMARTS-R1" below for the full diagnosis and fix. | **Fixed** — `[rN]` now has its own `AtomPrimitive::MinRingSize` primitive, `[kN]` unchanged. Full-corpus SMARTS match-set agreement 96.9% → 99.93%. |
| `[R1]`/`[R2]` ring-count on bridged/cage systems (~1.7% of such molecules) | Genuine SSSR-*basis*-cardinality disagreement (not a predicate bug) — e.g. an adamantane-type cage where RDKit's SSSR algorithm picks 5 basis rings and chematic's picks 4, both valid bases of the same cycle space. Distinct from the `[rN]` bug above; NOT fixed by SMARTS-R1. | **Deferred (SMARTS-R2)** — would need a richer, RDKit-compatible symmetrized ring model scoped to SMARTS ring-count predicates only; not started. |
| Bridgehead-N fused-ring aromaticity over-extension (one tricyclic scaffold, 28 substituent variants) | `apply_aromaticity("rdkit_like")` propagates aromaticity from a genuinely aromatic ring, across a bridgehead-N ring fusion, into an adjacent ring that shouldn't qualify — see "SMARTS-A0" below for the full diagnosis. **Not carbonyl-related** — an earlier version of this table mislabeled it as such; the earlier "ring-junction carbonyl" description never had a confirmed minimal reproducer. | Documented (diagnosis only, SMARTS-A0); entire remaining SMARTS residual (56/80,000 comparisons, 5k corpus) — the `[rN]` fix left nothing else. Fix requires a scope decision — see SMARTS-A0's options below. |
| Morgan bit positions | FNV-1a (chematic) vs MurmurHash (RDKit). | By design — similarity **ranking** is consistent; individual bit indices are not comparable across libraries. |
| ~~MOL/SDF double-bond cis/trans (MDL stereo code 3) not perceived~~ — **Implemented (P1-S2)**, removed from this table | Tetrahedral wedge/hash perception (P1-S1a) and E/Z direction-writing (P1-S2) are RDKit's own two structurally separate pipelines (`assignChiralTypesFromBondDirs` vs `detectBondStereochemistry`) — see `docs/stereo2d_reader_integration_rfc.md` §5b. Both are now implemented via `chematic_perception::stereo2d_ez_direction`, wired into the V2000/V3000 readers (SDF inherits via the V2000 core). | **Fixed** — see the compatibility matrix and differential validation table above. Known residual: a branch point (alkene end with 2 substituents) adjacent to another double bond, and a shared-carrier bond between two independently-stereogenic double bonds that disagree, both fail closed (abstain) rather than guess — this is Issue #149's joint-carrier-resolution problem, explicitly out of scope. Kekulized aromatic-ring bonds are excluded entirely (never stereogenic). |
| MDL bond-stereo code 4 / V3000 `CFG=2` ("either"/unspecified) round-trips as a plain, undirected bond | Code 4 is a third, distinct MDL state (RDKit: `Bond::BondDir::UNKNOWN`) from a definite wedge/hash; treating it as "no defined direction" avoids fabricating a stereocenter the file explicitly declines to specify, but nothing currently preserves "was marked either" through a write. | **Documented, accepted** — never surfaces as a wrong wedge, only as a dropped "either" annotation. |
| V3000 atom-line `CFG` (parity 1/2/3) not decoded | RDKit's own primary `MolFromMolBlock` path is bond-CFG/wedge-based (`assignChiralTypesFromBondDirs`); its atom-parity path (`assignChiralTypesFromMolParity`) had no found callers in a pinned-commit source audit. | **Out of scope** — bond `CFG` (wedge direction) is the implemented mechanism. |
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
| Aromatic ring-junction tail (`c`, `[#6]=[#6]`) — root cause diagnosed in "SMARTS-A0" below, not carbonyl-related despite the name used here at the time | 140 | 5.6%, same 28-molecule root cause across both patterns |
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

### SMARTS-R1 — `[rN]`/`[kN]` predicate split, implemented (2026-07-16)

Implements SMARTS-R0's finding: `[rN]` now has its own `AtomPrimitive::MinRingSize`
primitive (`crates/chematic-smarts/src/query.rs`), computing "the smallest SSSR ring
containing this atom" — separate from `[kN]`'s existing `AtomPrimitive::RingSize`
("any SSSR ring of size N"), which is untouched. `crates/chematic-smarts/src/parser.rs`
routes `[rN]` to `MinRingSize` and `[kN]` to `RingSize` (previously both routed to
`RingSize`). No change to ring perception itself (`find_sssr`) or to the descriptor/
aromaticity ring model, per this doc's original constraint.

**Matching**: `crates/chematic-smarts/src/match_vf2.rs`'s `EvalCtx` gained a lazily-computed,
per-atom min-ring-size cache (`RefCell<Option<Vec<Option<u8>>>>`), populated on first
`[rN]` query and reused for the rest of that pattern's match (VF2 backtracking
re-evaluates atom predicates many times per search) — never recomputed per atom, per
backtrack step, or for patterns that don't use `[rN]` at all. No new ring-finding: the
cache aggregates over the same `RingSet` already computed once per molecule and shared
across every pattern (unchanged from before this milestone).

**Tests**: 4 new (`crates/chematic-smarts/src/parser.rs`) — a fused 5/6-ring bicyclic
where `[k5]`/`[k6]` both match the fusion atom but `[r6]` does not (`[r5]` does); the
exact purine discriminating case from SMARTS-R0's diagnosis, pinned against RDKit's
real atom sets; `[rN]` false on every size for a fully acyclic molecule;
`[R]`/`[R1]`/`[R2]` unaffected by the split. 2 existing tests that asserted
`[k5]==[r5]`/`[k6]==[r6]` (encoding the bug) corrected to assert they're distinct
primitives with distinct semantics.

**Full-corpus verification** (regenerated `validation/results/rdkit_compat_diff.jsonl`,
same 5,000-molecule corpus, same RDKit version):

| Metric | Before (frozen, stale) | After (SMARTS-R1) |
|---|---|---|
| SMARTS match-set agreement | 96.9% | **99.93%** (79,944/80,000) |
| `[rN]` agreement, isolated check | 70.95% | **100.00%** (4,999/4,999) |
| `[kN]` agreement, regression check | — | 99.98% (4,998/4,999, 1 pre-existing bridged-cage mismatch confirmed present on `main` before this change too — not a regression) |
| Aromatic atom count | 99.0% | 99.44% |
| Aromatic bond count | 98.3% | 98.82% |

Zero regressions: the single remaining `[kN]` mismatch (a bridged bicyclic where
RDKit's and chematic's SSSR bases disagree in cardinality) was independently confirmed
present and identical on `main` before this change (stashed the fix, rebuilt, reran).
The entire remaining 56-row SMARTS residual (of 80,000 comparisons) is the `c`/`[#6]=[#6]`
aromatic ring-junction tail diagnosed in full in "SMARTS-A0" below — no `[rN]`
contribution left at all.

**Explicitly not done**: `[R1]`/`[R2]` on bridged/cage systems (SMARTS-R2, a genuine
SSSR-basis-cardinality gap, not a predicate bug — kept separate per SMARTS-R0's own
partition, above).

### SMARTS-R2 — opt-in RDKit-parity `[RN]` ring-count model (2026-07-25, implemented, opt-in)

Implements SMARTS-R0's deferred finding: an opt-in-only fix for the `[R1]`/`[R2]`
ring-*count* residual, entirely inside `crates/chematic-smarts/` (`rdkit_ring_model.rs` +
`rdkit_parity_match.rs`). **The default matcher (`find_matches`/`match_vf2.rs`) is
untouched — zero diff** (`git diff crates/chematic-smarts/src/match_vf2.rs` is empty);
the new mode is a deliberately-duplicated VF2 matcher (`find_matches_rdkit_parity`)
reached only by explicitly opting in, never by the existing entry points.

**Root cause, confirmed by reading RDKit source at the pinned commit
`8afba32ec539dcb2369bc84549d802aca3f7eb39`**: `[R]`/`[RN]`/`[r]`/`[rN]`/`[k]`/`[kN]`/
`[x]`/`[xN]`/ring-bond `@`/`!@` are all backed by one shared `RingInfo` object
(`Code/GraphMol/QueryOps.h:283-370` — `queryIsAtomInNRings`→`numAtomRings`,
`queryAtomMinRingSize`→`minAtomRingSize`, `queryAtomIsInRingOfSize`→`isAtomInRingOfSize`,
`queryAtomRingBondCount`/`queryIsBondInRing`→`numBondRings`; `AtomRingQuery`,
`QueryOps.h:752-789`, backs bare `[R]`/`[RN]`). RDKit's default sanitization does not
populate that `RingInfo` from a minimal SSSR alone: `MolOps::symmetrizeSSSR`
(`Code/GraphMol/FindRings.cpp:996-1093`) adds "extra" rings back in whenever a candidate
(found among RDKit's own SSSR search's rejected duplicate D2 candidates,
`findSSSRforDupCands`, same file line ~283) is the same size as some basis ring, shares
≥1 bond with it, and does not drop any bond that basis ring is the *sole* provider of
(the `bondCounts`/`replacesAllUniqueBonds` logic at `FindRings.cpp:1046-1087`). This can
make RDKit's per-atom ring **count** larger than a minimal-basis SSSR's — the "genuine
SSSR-basis-cardinality disagreement" SMARTS-R0 already named.

**A graph-theory fact narrows the fix to exactly one primitive.** An edge lies on *some*
basis cycle of a graph if and only if it lies on *any* cycle at all (a cycle's edge set is
a member of the GF(2) cycle space, so if an edge appeared in zero basis cycles it could
not appear in any combination of them either) — independent of which valid basis is
chosen. Consequently `[R]`/`[R0]` (ring membership, boolean), `[x]`/`[xN]` (ring-bond
count) and ring-bond `@`/`!@` are **provably invariant** to which SSSR basis backs them;
only `[RN]` (N ≥ 1, exact ring *count*) can move. `[rN]`/`[kN]` could in principle move
too, but both already measure ~100%/99.98% on plain SSSR alone (SMARTS-R1, above) — the
new opt-in matcher deliberately leaves them wired to plain SSSR, unchanged.

**Design.** `rdkit_ring_model.rs` builds an auxiliary, atom-indexed ring-*count* table
from chematic's own already-computed SSSR (never modifies `find_sssr`): it enumerates
simple cycles up to the basis's largest ring size via a depth-bounded DFS (depth-bounded
by ring size, not by molecule size or ring count — chosen over a GF(2) basis-subset
search specifically because a fixed subset-size cap silently misses cage topologies like
cubane, where the 6th face is the XOR of *all 5* basis rings), then re-applies RDKit's own
same-size/shares-a-bond/doesn't-drop-a-unique-bond acceptance rule to each candidate. This
is a *different candidate-generation mechanism* reaching for the *same acceptance rule* —
chematic's `find_sssr` doesn't expose RDKit's own rejected-duplicate-candidate list, so
this module can't literally replay RDKit's search, only its filter.

**Measured, full corpus (5,000-molecule `~/Downloads/SMILES.csv` + 21 hand-built
fused/bridged/spiro/cage structures, 30 patterns, 155,651 molecule×pattern cells — see
provenance in the PR body)**:

| Metric | Value |
|---|---|
| Total cells compared | 155,651 |
| Default matcher agreement with RDKit | 155,474/155,651 (99.8863%) |
| Opt-in RDKit-parity matcher agreement with RDKit | 155,650/155,651 (99.9994%) |
| `[RN]`-family cells only | 25,105 (176 fixed, 0 regressed, 24,929 already agreed) |
| Regressions (opt-in mode disagrees where default agreed) | **0** |
| Residual (both disagree with RDKit) | 1 cell — pre-existing `[k5]` bridged-bicyclic SSSR-basis mismatch, unrelated to `[RN]`, already named by SMARTS-R1 above; opt-in mode deliberately leaves `[kN]` on plain SSSR so this is unchanged, not newly introduced |

Real cage examples the new model recovers (RDKit ground truth via
`rdkit==2026.03.3` live oracle, PubChem canonical SMILES): adamantane (cycle_rank 3 →
RDKit 4 rings), bicyclo[2.2.2]octane (cycle_rank 2 → RDKit 3 rings), cubane (cycle_rank 5
→ RDKit 6 rings), **dodecahedrane** (cycle_rank 11 → RDKit 12 rings, all 20 atoms
uniformly in exactly 3 faces — used as this track's adversarial highly-symmetric
termination test; confirmed to terminate and match RDKit exactly within the default
2,000,000-candidate budget).

**Bridgehead-N (SMARTS-A0) bucket, separately measured, not conflated with the ring-count
fix above.** This crate's matching entry points (`find_matches`/`find_matches_rdkit_parity`)
never call any aromaticity re-perception themselves — they match whatever flags the input
molecule already carries. On the SMARTS-A0 bare-core reproducer
(`C1=Cc2ccccc2C2=NCCCN12`), direct parsing already carries RDKit-correct flags, so the bug
does not fire on this pipeline (0/21 hand-corpus cells). It fires only if a caller
explicitly re-perceives with `AromaticityAlgorithm::RdkitLike` first (e.g.
`rdkit_compat.MolFromSmiles(sanitize=True)`) — reproduced directly in
`rdkit_parity_match.rs`'s test suite. A third, independent engine,
`apply_aromaticity_rdkit_parity_experimental` (this opt-in matcher's separate
`use_rdkit_parity_aromaticity` flag, default off), does **not** reproduce the bug on this
reproducer — checked directly, not assumed, since the two aromaticity engines are
unrelated code paths.

**No molecule-specific allowlist** anywhere in `rdkit_ring_model.rs`/`rdkit_parity_match.rs`
— the acceptance rule and DFS bound are structural, not keyed on any SMILES/molecule
identity.

**API.** `chematic_smarts::{find_matches_rdkit_parity, has_match_rdkit_parity_bounded,
RdkitParityConfig, RdkitParityError, RdkitRingModelBudget}`. `RdkitParityError` is a typed,
non-silent failure: `RingModelBudgetExceeded` (the candidate search hit its cap — see
`RdkitRingModelBudget`) or `Aromaticity` (only when `use_rdkit_parity_aromaticity=true`
and re-perception failed). Never a silent partial match or a silent fallback to plain
SSSR under the RDKit-parity name.

### SMARTS-A0 — bridgehead-N ring-fusion over-extension diagnosis (2026-07-16, diagnosis only)

Full diagnosis of the 56-row residual SMARTS-R1 left behind (`c` and `[#6]=[#6]`, 28
rows each). **No production code changed.** Corrects an earlier label in this doc: the
residual was previously called "ring-junction carbonyl"; that name was never backed by a
confirmed reproducer. It isn't carbonyl-related at all.

**Reduces to one scaffold.** All 56 rows collapse to exactly 28 unique molecules, all
sharing one tricyclic core — a 6-6-6 linearly fused ring system (benzo ring + a
bridgehead-N-fused six-membered ring + a third ring closing back through that same
bridgehead N), varying only in an aryl/alkyl substituent hung off an exocyclic `C=C`.
Minimal bare-core reproducer (14 atoms, no substituent needed):

```
C1=Cc2ccccc2C2=NCCCN12
```

**Not a parser bug.** `chematic.from_smiles()` raw parsing matches RDKit's aromatic-atom
set exactly on all 28/28 molecules (0 mismatches). The defect is entirely in the
re-perception step `mol.apply_aromaticity("rdkit_like")` — what `rdkit_compat.py`'s
`MolFromSmiles(sanitize=True)` calls. Confirmed with a positive control: skip that call
and the SMARTS residual disappears.

**Mechanism, in `crates/chematic-perception/src/aromaticity.rs`.** The module's own
Pass-1/Pass-2 design (doc comment at the top of the file) already names the two
ingredients that combine here: a "bridgehead N" special case ("fused-ring N atoms whose
entire valence is satisfied by single σ-bonds, like indolizine's junction nitrogen") and
Pass 2's `aromatic_context` propagation ("confirmed-aromatic atoms contribute 1π
unconditionally, allowing fused rings to be recognised bottom-up"). On the reproducer:
Pass 1 correctly marks the benzo ring (atoms `2..7`) aromatic. Pass 2 then re-evaluates
the adjacent ring (atoms `0,1,2,7,8,13` — sharing the `2–7` edge with the benzo ring)
using that context: atoms `2,7` contribute 1π unconditionally per the propagation rule,
the `0=1` exocyclic-to-benzo double bond contributes 2π, and the bridgehead N at `13`
(single-bonded only, entering the third ring) is treated as aromaticity-compatible by
the bridgehead-N rule — enough for the whole adjacent ring to pass Hückel's count and
get swept into "aromatic," even though RDKit's own perception does not extend
aromaticity past the benzo ring on this scaffold.

**Isolated to two necessary-and-sufficient conditions**, by direct substitution on the
bare-core reproducer (all four combinations tested; only both together reproduce the bug):

| Variant | SMILES | Over-aromatizes? |
|---|---|---|
| Bridgehead N + exocyclic C=C (actual bug) | `C1=Cc2ccccc2C2=NCCCN12` | **Yes** — 4 extra atoms |
| All-carbon 3rd ring (N removed) | `C1=Cc2ccccc2C2=CCCC12` | No |
| N present but not at the bridgehead | `C1=Cc2ccccc2C2=CCNC12` | No |
| Bridgehead N kept, exocyclic C=C saturated | `C1Cc2ccccc2C2=NCCCN12` | No |

Three carbonyl-based mechanisms from the original taxonomy (benzene + exocyclic C=O,
fused-bicyclic + exocyclic C=O, heteroaromatic-junction C=O) were also tested directly
and **do not** reproduce any over-aromatization — ruling out the "exocyclic C=O strips
aromaticity" category entirely for this residual.

**Representation-stable**, per the acceptance criterion: the fully-Kekulized SMILES of
the same reproducer molecule (`C1=CC2=CC=CC=C2C2=NCCCN12`, no lowercase atoms at all)
produces the *identical* wrong result after `apply_aromaticity("rdkit_like")` as the
mixed aromatic/Kekulé spelling above — same 4 extra atoms. The diagnosis does not depend
on how the input was spelled.

**Fully quantified across all 28/28 molecules, 0 unclassified** (`scripts/smarts_a0_junction_diagnosis.py`,
data in `validation/smarts_a0_junction_diagnosis.jsonl`):

- Every molecule gets exactly 4 spurious aromatic atoms: 3 carbon + 1 nitrogen — 100% uniform, no variation.
- `c` query: chematic over-counts by exactly **+3** on all 28 (the 3 spurious carbons; the
  spurious nitrogen doesn't count towards `c`, which filters by element).
- `[#6]=[#6]` query: chematic under-counts by exactly **−1** on all 28 (the ring-fusion
  `C=C` bond itself gets re-typed `Double → Aromatic` once both endpoints are wrongly
  marked aromatic, so it no longer satisfies the literal `=` primitive).
- 28 + 28 = 56/56 rows fully accounted for by this single, uniform mechanism.

**Relation to already-known defects.** This is the same code region and general
mechanism family as two existing `#[ignore]`d regression tests in `aromaticity.rs`
(`test_azulene_kekulized_aromatic`, `test_purine_aromatic`, both tagged "fix belongs in
the `aromatic_context`-removal PR") — but a **different failure direction**. Those two
are false-*negatives*: a genuinely aromatic ring fails to get recognized because Pass 1
scores it non-aromatic and Pass 2 never revisits it. This is a false-*positive*: a
non-aromatic ring gets swept in because Pass 2's propagation rule doesn't check whether
the *adjacent* ring's own electron count, independent of the borrowed context, actually
supports delocalization. Same `aromatic_context` mechanism, opposite symptom — not
previously documented. Referenced design doc `greedy-hopping-crescent.md` (named in
those tests' comments) is not present in the repository; whatever plan it held did not
land as a file.

**What a fix would require — not decided here.** The defect lives in the shared
aromaticity-perception engine (`apply_aromaticity`/`assign_aromaticity_ex`), not in
`chematic-smarts`. Every consumer of that engine is affected in principle (descriptors,
canonical SMILES, fingerprints, InChI), even though this particular scaffold currently
only shows up as a SMARTS mismatch in the 5k corpus. Three options, not mutually
exclusive with the existing SMARTS-R2 deferral:

1. **Fix `apply_aromaticity`'s Pass-2 propagation** (require the adjacent ring's own
   non-context electron count to be Hückel-compatible before inheriting `aromatic_context`,
   or equivalent) — fixes SMARTS, descriptors, canonical SMILES, and fingerprints at
   once, but is a core-perception change: broad blast radius, needs a full descriptor/
   canonical/fingerprint regression pass, and likely needs to be resolved together with
   the azulene/purine `#[ignore]`d cases rather than in isolation, since they share the
   same propagation rule.
2. **SMARTS-only aromatic view** (mirroring the `[rN]`/`[kN]` scoped-predicate pattern
   from SMARTS-R1) — would keep `c`/`a` matching independent of `apply_aromaticity`'s
   own atom flags for this case. Architecturally heavier than SMARTS-R1's fix: it would
   make a molecule's SMARTS aromaticity diverge from its own `atom.aromatic` flags,
   which risks inconsistency in recursive SMARTS and aromatic-bond queries that read
   those flags directly. Not designed here, only flagged as an option.
3. **Accept 99.93% and defer** — the standing decision already in place for the
   `aromatic_context` family (azulene/purine).

No implementation in this round, per SMARTS-A0's scope.

**Follow-up decided**: Option 1 (fix the shared `apply_aromaticity` engine) was chosen
over Option 2/3 above. See `docs/aromaticity_a1_rfc.md` for the full design, the
component-level solver this points toward, and "Aromaticity-A1-0" (diagnosis/corpus,
no behavior change, landed) as the first step.

---

## Examples

Runnable, self-contained task scripts live in [`examples/rdkit_compat/`](../examples/rdkit_compat/):

- `sdf_to_csv.py` — SDF → descriptor table (CSV)
- `substructure_filter.py` — filter a SMILES list by a SMARTS query
- `fingerprint_similarity.py` — Morgan fingerprint + Tanimoto ranking
- `numpy_features.py` — fingerprint matrix for scikit-learn / PyTorch

A 10-pattern RDKit→chematic migration guide is in
[`examples/rdkit_compat_migration.py`](../examples/rdkit_compat_migration.py).
