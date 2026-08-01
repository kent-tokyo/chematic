# RFC: Canonical SMILES residual diagnosis (diagnosis only, no fix in this round)

Status: **Diagnosis complete. No production behavior change.**

> **Update (Wave 2A remeasurement, main@1bc1b63):** re-ran
> `scripts/canonical_residual_diagnosis.py` unmodified (no script changes) on
> the same 5,000-mol corpus, a week and several dozen merged PRs after the
> commit this RFC was originally written against. Fresh Check-2
> permutation-invariance failures: **91/5,000** (down from the pre-PR-#148
> baseline of 348 documented below — expected, since PR #148 landed in
> between; not re-litigated here). Breakdown:
> - **18 detected via random relabeling** (Root Cause 1's own mechanism,
>   below) — cross-checked byte-for-byte against
>   `EZ_SHARED_CANDIDATE_BOND_RESIDUALS` in
>   `crates/chematic-smiles/src/canonical.rs`: **exact match, 18/18, 0 new,
>   0 resolved.** Issue #149's residual is stable, not drifting in either
>   direction. `ez_carrier_shared_candidate_bond_residuals_never_corrupt`
>   still passes on current `main`.
> - **73 detected via idempotence only** (Root Cause 3, aromaticity
>   round-trip — a separate, deferred track this remeasurement does not
>   touch): up from the 67 reported below, and now 16/73 carry an E/Z
>   marker (vs 1/67 previously) — a compositional shift flagged here for
>   whoever owns the aromaticity track, not investigated further (out of
>   scope for an E/Z-carrier-focused remeasurement).
> - **0 detected via both probes** (down from 11) — consistent with Root
>   Cause 1 having no failures left beyond the pinned 18.
> - **91/91 semantically identical** (RDKit structural re-comparison) — 0
>   corruption, unchanged from the original finding.
>
> No production code changed for this remeasurement — only
> `validation/results/canonical_residual_diagnosis.jsonl` and
> `..._summary.json` were regenerated (same reproduce command as below).
> **Conclusion for a future Wave 2B (issue #149) fix**: the target is
> unchanged and exactly scoped — the same 18 SMILES, still abstaining for
> the same reason, still 0% corrupted. A joint/simultaneous carrier resolver
> per the issue's own spec remains the right next step; nothing here changes
> that scope.

> **Update (C1a):** Root Cause 1 (the dominant, E/Z-marker-carrier cluster
> described below) has a safe, verified **partial** fix — branch
> `fix/canonical-ez-carrier-normalization`, PR #148. 264/282 of the
> `has_ez_marker` diagnosis subset now converge to one canonical string; 18
> remain a deliberate, documented residual (two independently stereogenic
> double bonds sharing one candidate carrier bond — see tracking issue
> "canonical E/Z: jointly resolve shared carrier bonds across coupled stereo
> systems"). Zero semantic corruption in either case. Root Causes 2 and 3
> below are still untouched/deferred. See
> `validation/results/canonical_ez_carrier_c1a_reconciliation.json` for the
> full count reconciliation (including why this doc's/task history's "281"
> and the committed diagnosis artifact's "282" are not a regression against
> each other).

> **Update (C2), correction — Root Cause 2 is a verification error, not a
> chematic bug:** re-investigation (branch
> `fix/canonical-bridged-ring-ordering`) found that the sole reproduction
> for Root Cause 2 — `chematic.from_smiles("C1CC2CCC1CC2")` vs.
> `chematic.from_smiles("C1CCC2CC1CC2")`, both labeled "bicyclo[2.2.2]octane"
> below — was never checked against an independent structural oracle, unlike
> every corpus-derived finding in this RFC. RDKit `MolToInchi` shows the two
> inputs are in fact **different constitutional isomers** of C8H14:
> `C1CC2CCC1CC2` is bicyclo[2.2.2]octane (three 2-atom bridges,
> `InChI=1S/C8H14/c1-2-8-5-3-7(1)4-6-8/h7-8H,1-6H2`), while `C1CCC2CC1CC2`
> is a different bridged bicyclooctane with bridges of 1/2/3 atoms
> (`InChI=1S/C8H14/c1-2-7-4-5-8(3-1)6-7/h7-8H,1-6H2` — a different `/c`
> connectivity layer, confirmed via RDKit's own divergent canonical SMILES
> for the two inputs, `C1CC2CCC1CC2` vs. `C1CC2CCC(C1)C2`). Chematic giving
> the two inputs different canonical strings was **correct** — it was never
> a permutation-invariance bug, because the two inputs never encoded the
> same molecule. This also explains, independently, why the RFC's own
> corpus scan (below) already reported **0/232** real bridged-ring corpus
> failures attributable to this specific mechanism: it does not exist. A
> targeted follow-up probe (RDKit `RenumberAtoms`-generated same-molecule
> respellings, the same method the corpus check uses, across 22 diverse
> bridged/bicyclic/spiro/cage/fused systems — including stereocenters and a
> heteroatom bridgehead) found **zero** convergence failures. The previously
> "known gap" test (`bridged_bicyclic_canonical_gap_documentation` in
> `crates/chematic-smiles/tests/canonical_robustness.rs`) has been replaced
> with `bridged_fused_spiro_permutation_invariance`, an assertion-based
> regression test using 8 InChI-verified same-molecule pairs. No change to
> `crates/chematic-smiles/src/canonical.rs` was made or needed. See PR
> (branch `fix/canonical-bridged-ring-ordering`) for the full evidence
> trail. Root Cause 3 below remains untouched/deferred (parallel
> aromaticity-diagnosis track).

Branch: `diag/canonical-smiles-residual` (forked from `main`@`659baca221f71f135ce0e1780e71245d8770f132`).

## Scope declaration

- **Files created:** this doc (`docs/canonical_smiles_residual_rfc.md`),
  `scripts/canonical_residual_diagnosis.py`,
  `validation/canonical_residual_fixtures.jsonl`,
  `validation/results/canonical_residual_diagnosis.jsonl`,
  `validation/results/canonical_residual_diagnosis_summary.json`.
- **Files touched:** none under `crates/*/src/**` (no production code). No existing
  script/doc edited in place (`canonical_diff.py` was re-run, not modified, to
  reconcile a stale headline number — see "Reconciling with existing docs" below).
- **Explicitly out of scope / not touched:** `feat/io-mrv`, `feat/io-tdt`,
  `feat/io-smiles-supplier-writer`, `fix/smiles-bracket-implicit-h`,
  `diag/stereo-reader-integration-boundary`, `feat/stereo2d-local-parity`, or any
  other agent's branch/PR. No fix to `crates/chematic-smiles/src/canonical.rs`
  (root cause is pinpointed below to specific line ranges, but not patched).
- **Deliverables:** this RFC; a re-runnable, seeded (reproducible)
  diagnosis script; a JSON summary + JSONL per-mismatch detail under
  `validation/results/`; a frozen fixture file for hand-constructed and
  representative corpus-sampled repros.
- **Done when:** all four checks below are measured and reported separately
  (never pooled), every Check-1/Check-4 residual lands in a named bucket or
  "unclassified" with an accounted-for total, and every claim of "same
  molecule" is backed by RDKit structural introspection, not string
  comparison. All of that is satisfied as of this document.

## Prior work (read first, built on, not duplicated)

Substantial canonical-SMILES work already exists in this repo:

- `scripts/canonical_diff.py` — round-trip semantic equivalence (`RDKit(chematic
  canonical) == RDKit(original)`) + chematic-internal idempotency, on a
  5,000-mol corpus. Committed result: `validation/results/canonical_diff.jsonl`.
- `scripts/canonical_structural_correctness.py` — per-molecule check that N
  RDKit-`doRandom` respellings of the same molecule, fed through chematic,
  all land on the same RDKit-canonical molecule as the original (a semantic-
  correctness check, not a chematic-self-consistency check).
- `crates/chematic-smiles/tests/canonical_robustness.rs` — hand-picked
  regression corpus (RDKit issue #8775-inspired), including one **documented,
  non-panicking, still-open gap**: `bridged_bicyclic_canonical_gap_documentation`
  (two spellings of bicyclo[2.2.2]octane give different chematic canonical
  output).
- `crates/chematic-py/tests/test_canonical_diff.py` — records that the
  formerly-known "exocyclic C=N E/Z direction dropped near an aromatic ring"
  divergence is **fixed** (bond-direction side-channel stash).
- `docs/rdkit_compat.md` — "Known divergence classes" table, citing **99.62%
  round-trip** and **98.42% idempotent** on the 5k corpus, and root-causing
  the idempotency residual to an aromaticity-perception round-trip
  inconsistency on large fused polycyclics (deferred to the aromaticity
  track).

This diagnosis reuses the same 5,000-molecule ChEMBL-derived corpus
(`~/Downloads/SMILES.csv`, the project's standard corpus) and RDKit 2026.03.3
(verified: `.venv/bin/python -c "import rdkit; print(rdkit.__version__)"` →
`2026.03.3`), and adds the four checks requested for this round, each
measured and reported **separately** — never pooled into one "agreement %".

## Reconciling with existing docs: the 99.62%/98.42% figures are stale

Re-running the existing, unmodified `scripts/canonical_diff.py` against the
current commit (`659baca`, same as `main`) gives:

```
round-trip MATCH:   5000/5000 (100.00%)   [docs cite 99.62%]
idempotent:         4922/5000 (98.44%)    [docs cite 98.42%, within noise]
```

Round-trip parity moved from 99.62% to **100.00%** because the exocyclic-C=N
E/Z fix referenced above (`bond_directions` stash) landed after the 99.62%
figure was measured — `test_canonical_diff.py`'s own comment already says
this ("Formerly-known divergence, now fixed"), but `docs/rdkit_compat.md`'s
table was never resynced. Idempotency is unchanged within measurement noise
(78 vs the previously-reported failures — same order of magnitude, same
mechanism, see Check 3 below). **Recommendation, not performed here:** a
docs-sync pass to `docs/rdkit_compat.md`'s divergence table — out of scope
for a diagnosis-only round, and touching a shared doc risks colliding with
a concurrently-running agent's own edits to the same file.

## Method

`scripts/canonical_residual_diagnosis.py`, full corpus (5,000 mol), `--k 8`
permutation variants per molecule, `--seed 0` (byte-reproducible — verified
by running twice and diffing the JSON output). Runtime: ~19s.

Four checks, kept separate:

1. **RDKit exact canonical string parity** — `chematic.from_smiles(s).smiles
   == RDKit.MolToSmiles(RDKit.MolFromSmiles(s))`, byte-for-byte.
2. **Permutation invariance** — K reproducibly-relabeled RDKit spellings of
   the *same parsed molecule* (`Chem.RenumberAtoms` with a shuffle seeded by
   `random.Random(seed)`, deliberately **not** RDKit's own `doRandom=True`,
   which draws from RDKit's unseeded global RNG and would make every run —
   and this diagnosis's own committed JSON — non-reproducible), **plus
   chematic's own canonical output fed back through chematic** (see the note
   below), fed through chematic; chematic's own output must be identical
   across all spellings. A failure here is always a real chematic bug
   (chematic-internal self-consistency; RDKit is only used to generate
   alternate valid spellings of one molecule). **This is a lower bound at
   the tested K=8 plus one idempotence probe** — passing all of them is
   evidence against instability being found in those samples, not a proof
   of invariance under all possible relabelings.

   **Check 2 and Check 3 are not independent — Check 2 subsumes Check 3.**
   `canonical(s)` is itself a valid spelling of the same molecule `s`
   encodes, so true permutation invariance (identical output for *every*
   valid spelling) logically implies idempotence; contrapositive, **every
   idempotence failure is automatically a permutation-invariance failure**.
   The two checks are not measuring different things; they are two
   different *probes* for the same underlying property, and they have
   different blind spots: chematic's own canonical DFS traversal is one
   very specific atom ordering that a uniform-random relabeling essentially
   never reconstructs, so a mechanism that only manifests when fed
   chematic's *own* canonical spelling back in (Root cause 3 below) is
   invisible to K random relabelings alone and needs the idempotence probe
   to surface at all. This diagnosis therefore folds `canonical(canonical(s))`
   into Check 2's own invariance test directly (not as a separate,
   after-the-fact cross-reference against Check 3's failure list), so
   Check 2's reported number is the tight one on a single run. An earlier
   draft of this diagnosis measured Check 2 with K-relabeling alone,
   reporting 94.38%/281 failures — **that undercounted**; folding in the
   idempotence probe raises true non-invariance to 348/5,000, corrected
   throughout this document (data preserved in the script's
   `detected_via_idempotence` / `detected_via_random_relabeling` per-row tags
   so the two mechanisms below remain distinguishable in the corrected total).
3. **Idempotence** — `canonical(canonical(s)) == canonical(s)`. Reported
   separately per the task's requirement, and useful in its own right as
   the specific probe that catches Root cause 3 below — but not an
   independent finding from Check 2 (see the note above).
4. **Semantic structure parity** — for every Check-1 mismatch, reparse
   chematic's canonical output through RDKit and compare **actual structure**
   (formula, heavy-atom multiset, aromatic atom/bond counts, SSSR ring-size
   multiset, CIP stereocenter label multiset, bond E/Z multiset, isotope/
   charge/atom-map multiset) against RDKit's own canonicalization of the
   original input — never string comparison alone. Implemented as a
   pure-data decision tree (`MolFeatures` → `classify_real_diff`), unit-tested
   in isolation from RDKit object construction (`--self-test`) so the
   classification logic itself is verified, not just exercised.

Every Check-1/Check-4 residual is placed in exactly one of: `aromaticity_
kekulization`, `tetrahedral_parity`, `ez_direction`, `ring_closure_ordering`,
`bridged_fused_spiro`, `disconnected_fragment_ordering`, `isotope_charge_
atommap`, `symmetry_tie_breaking`, `writer_token_bug`, `unclassified`. No
silent drops — `check4_bucket_totals_check.accounted_for` asserts
`sum(bucket_counts) == check1_mismatch_count` (verified `true`).

## Results (5,000-mol corpus, K=8, seed=0)

| Check | Result | Note |
|---|---|---|
| 1. Exact string parity | **0.16%** (8/5,000) | Expected — different algorithms. Not a bug signal by itself; only an entry point into Check 4. |
| 2. Permutation invariance | **93.04%** (4,652/5,000); **348 real permutation-invariance bugs** | Includes the idempotence probe (see method note above) — **two distinct mechanisms**, see Root causes 1 and 3 below. All 348 are molecule-preserving (0 flip the encoded chemistry). |
| 3. Idempotence | **98.44%** (4,922/5,000); 78 failures | Matches `canonical_diff.py`'s fresh 98.44% exactly. A **subset** of Check 2's 348 (78 = the 67 idempotence-only failures + 11 also caught by random relabeling). All 78 independently confirmed (RDKit structural re-comparison) to be molecule-preserving (same structure before/after) — 0/78 show genuine structural drift. |
| 4. Semantic structure parity (applied to all 4,992 Check-1 mismatches) | **100% cosmetic — 0 real semantic bugs** | Verified directly: `sum(not r["semantically_identical"] for r in rows) == 0`, not inferred from bucket labels (a bug class can reach the "real diff" branch and still land in a bucket — e.g. `bridged_fused_spiro` — that also has a cosmetic-branch meaning; checked the boolean directly to avoid that trap). |

**On "cosmetic" vs "real bug" language**: every one of the 348 Check-2
failures is a **real permutation-invariance bug** by the task's own
definition, full stop — none is dismissed as cosmetic. What legitimately
varies is *where the eventual fix belongs* (a self-contained SMILES-writer
change for Root cause 1, vs. the shared aromaticity-perception engine for
Root cause 3) and whether the encoded *chemistry* is corrupted (it never is,
across all 348) — that second property is reported as "semantically
identical / molecule-preserving," not as "not a bug."

### Check 2 breakdown: two distinct mechanisms, by detection probe

| | Idempotence probe only | Random relabeling only (K=8) | Both probes | Total |
|---|---|---|---|---|
| Count | 67 | 270 | 11 | 348 |
| Has an E/Z marker in some output | 1/67 | 270/270 | 11/11 | 282/348 |

Root cause 1 (E/Z marker-selection instability, below) accounts for
270 + 11 = **281** failures (all E/Z-marker-bearing, 100%). Root cause 3
(aromaticity round-trip inconsistency, below) accounts for the 67
idempotence-only failures (66/67 have no E/Z marker at all — a genuinely
separate, non-stereo mechanism). The 11 "both" cases are E/Z-marker
molecules where the marker-selection instability itself also breaks the
idempotence fixed point (feeding chematic's own output back in re-triggers
a different marker choice) — a compound manifestation of Root cause 1, not
a third mechanism.

### Check 4 bucket breakdown (of the 4,992 Check-1 mismatches — all cosmetic)

| Bucket | Count | Share |
|---|---|---|
| `ring_closure_ordering` | 4,647 | 93.1% |
| `bridged_fused_spiro` | 276 | 5.5% |
| `symmetry_tie_breaking` | 37 | 0.7% |
| `unclassified` | 32 | 0.6% |
| `aromaticity_kekulization` | 0 | — |
| `tetrahedral_parity` | 0 | — |
| `ez_direction` | 0 | — |
| `disconnected_fragment_ordering` | 0 | — |
| `isotope_charge_atommap` | 0 | — |
| `writer_token_bug` | 0 | — |
| **Total** | **4,992** | **accounted_for: true** |

The 32 `unclassified` were inspected individually (not just counted): all 32
are acyclic, single-fragment, automorphism-untied molecules (e.g. `CNCC(O)CO`
→ `CNCC(CO)O`) where chematic and RDKit simply pick a different valid
start-atom/branch order — confirmed cosmetic, just not explained by any of
the three structural heuristics (`is_multi_fragment`, `has_symmetry_tie`,
`ring_relationship`). **Caveat on the cosmetic sub-buckets**: since virtually
every Check-1 mismatch is, at bottom, "two different valid canonicalization
algorithms disagreeing" (not a bug), the `ring_closure_ordering` /
`bridged_fused_spiro` / `symmetry_tie_breaking` split is a best-effort
**descriptive** heuristic over the molecule's topology (does it have a ring
/ a bridged-fused-spiro system / a tied canonical rank), not a proof that
the specific string divergence was *caused* by that trait — see the
docstring on `classify_cosmetic` in the script for the full caveat and the
priority order used (ring/fragment topology checked before the near-
ubiquitous symmetry-tie catch-all, so it doesn't swallow more specific
buckets).

## Root cause 1 (Check 2, dominant, 281/348 = 80.7% of the non-invariance found, 5.6% of the corpus, 54% of E/Z-bearing molecules): trisubstituted/tetrasubstituted double-bond marker-selection instability

**Mechanism, pinned to exact code.** `crates/chematic-smiles/src/canonical.rs`'s
`normalize_ez` (~L443–470) normalizes the *flip polarity* of an E/Z system's
directional bonds so the first-encountered bond in canonical write order is
always `Up` (`/`) — this part is correct and does not depend on input
spelling. But at a **trisubstituted or tetrasubstituted** stereo double-bond
carbon, SMILES only requires *one* of the ≥2 substituent bonds to carry an
explicit `/`/`\` marker (the other substituent's relative position is
implied). Which of the ≥2 candidate substituent bonds gets that mark is
decided by the **parser**, not re-derived by the canonical writer — the
writer's emission logic (`write_chain`/`dfs_mark`, ~L585–715) only checks
"does *this specific* bond already carry a directional `order`", inherited
straight from parse time. Two RDKit-valid respellings of the *identical*
molecule can validly choose to mark either substituent, so chematic parses
two different (but chemically-equivalent) `Molecule`s and its writer, having
no canonical rule for "which substituent should carry the mark," just
reproduces whichever choice the input made.

**Confirmed with direct evidence**, not inferred: `chematic_output_a` /
`_b` for `COc1ccc(/C=C(\C)C(=O)c2cc(OC)c(OC)c(OC)c2)cc1O` trace back to two
RDKit-generated variants that mark the acyl bond vs. the methyl bond,
respectively (see `validation/canonical_residual_fixtures.jsonl`,
`perm-inv-02-trisub-ez-chalcone`). RDKit's own `BondStereo` on both chematic
outputs is identical (`STEREOZ` for the imine fixture) — **the encoded
geometry never changes**, only the token choice.

**Fully quantified, 0 unclassified**: of the 281 corpus failures attributable
to this mechanism (270 caught by random relabeling alone + 11 also caught by
the idempotence probe — see the breakdown table above), **281/281 (100%)**
have at least one stereo double bond with ≥2 heavy-atom substituents on one
side; **0/281** show an actual E/Z value flip across the divergent outputs
(checked via RDKit `BondStereo` on each reparsed output, not string
comparison). This is a **writer canonicalization-stability bug** (violates
the canonicalizer's "one string per molecule" contract, and is squarely a
real permutation-invariance bug per the task's definition) — it does not,
however, corrupt the encoded chemistry.

**Independent of both parallel tracks.** This mechanism lives entirely in
the SMILES writer's own marker-selection logic; it does not touch CIP/
tetrahedral assignment (the 2D-stereo track's domain) or aromaticity
perception (the aromaticity track's domain) at all. See "What this diagnosis
does not do" below for why it is still not fixed in this round.

## Root cause 2 (Check 2, hand-constructed, 0/5,000 corpus occurrences but confirmed still live): bridged-bicyclic ring-closure ordering

> **CORRECTED, see "Update (C2)" at the top of this document.** Everything
> below this note is the ORIGINAL (incorrect) diagnosis, kept verbatim for
> the historical record. The two SMILES quoted immediately below are, per
> RDKit `MolToInchi`, two different molecules, not two spellings of one —
> "confirmed still live" here was never independently checked against a
> structural oracle. There is no bridged-bicyclic ring-closure-ordering bug.

Independently reconfirmed at commit `659baca`:
`chematic.from_smiles("C1CC2CCC1CC2").smiles` → `C12CCC(CC2)CC1`, while
`chematic.from_smiles("C1CCC2CC1CC2").smiles` → `C1C2CCC1CCC2` — two
spellings of bicyclo[2.2.2]octane (no stereo, no aromaticity, no E/Z at all)
give different chematic canonical output. This is **not** the same
mechanism as Root cause 1 (no directional bonds involved at all) — it is
purely about which ring atom the canonical DFS reaches first among tied
bridgehead candidates, already documented as a non-panicking known gap in
`crates/chematic-smiles/tests/canonical_robustness.rs`
(`bridged_bicyclic_canonical_gap_documentation`). Zero occurrences in the
5,000-mol ChEMBL corpus attributable to *this specific* mechanism (drug-like
molecules rarely contain plain undecorated aliphatic bridged bicyclics with
no E/Z and no aromatic ring at all): of the corpus's 232 bridged-ring
molecules, 13 fail the (corrected, unified) Check 2 — 10 own to Root cause 1
(E/Z marker present) and the remaining 3 own to Root cause 3 below (all 3
idempotence-only, i.e. invisible to K=8 random relabeling, on large
fused-aromatic macrolide/ellagitannin-like structures — the exact profile
Root cause 3 already predicts, confirmed individually, not assumed from the
aggregate count). **0 of the 232 own to Root cause 2 specifically** — this
mechanism is real (hand-confirmed above) but this particular corpus happens
not to sample it; it is not merely under-sampled evidence in favor of it,
the corpus's bridged-ring failures are fully and specifically explained by
the *other two* mechanisms. Pinned as fixture
`perm-inv-01-bridged-bicyclooctane`.

## Root cause 3 (Check 2's idempotence-only failures, 67/348 = 19.3% of the non-invariance found, 78 total idempotence failures once the 11 E/Z-overlap cases are included): aromaticity-perception round-trip inconsistency — already documented, reconfirmed, deferred

All 78 idempotence failures were checked (not assumed) to be
molecule-preserving: RDKit re-canonicalizes `canonical(s)` and
`canonical(canonical(s))` to the *identical* molecule in all 78/78 cases —
0 genuine structural drift. This matches `docs/rdkit_compat.md`'s
already-documented finding ("Canonical idempotency on large fused
polycyclics... Aromaticity-perception round-trip inconsistency... The
molecule is preserved (InChI invariant); only the representation differs").
This diagnosis does not re-derive that root cause (out of scope — it
belongs to the parallel aromaticity-diagnosis track) but does independently
reconfirm the failure rate, the "molecule-preserving, not structural drift"
property, and — new in this round — that it is specifically **the
idempotence probe, not K=8 random relabeling, that surfaces it**: 66/67 of
these failures have no E/Z marker at all, so they are not Root cause 1
recurring under a different name; they are a genuinely distinct mechanism
that only manifests when chematic's own canonical spelling is fed back
through itself (see the Check 2/3 relationship note in "Method" above for
why that is expected — chematic's own DFS traversal order is a needle in
the haystack for uniform random relabeling to hit by chance).

## What this diagnosis does not do

- **No production code is fixed or modified.** `crates/chematic-smiles/src/
  canonical.rs`'s `normalize_ez` and the bridged-ring DFS traversal are
  pinpointed above but not patched.
- **Per this diagnosis's mandate: production fixes to canonical SMILES wait
  until the parallel 2D-stereo implementation and aromaticity-diagnosis
  efforts elsewhere in this project land and stabilize**, since both are
  likely root causes for a chunk of the residual found here — concretely,
  Root cause 3 (the 78-case idempotency residual) **is** that chunk, sharing
  the exact aromaticity-perception mechanism the aromaticity track owns.
- **Clarifying which residuals that deferral actually covers**, so a future
  implementer does not assume everything here is blocked: Root cause 1 (the
  dominant 281-case Check-2 finding) and Root cause 2 (bridged-bicyclic
  ring-closure ordering) are **self-contained SMILES-writer bugs**,
  independent of both the 2D-stereo and aromaticity tracks — they touch
  `canonical.rs`'s own directional-bond marker selection and ring-closure
  DFS ordering, not CIP/tetrahedral assignment or Hückel perception. They
  are not fixed in this round either (diagnosis-only mandate for this
  branch), but their non-fix is not because they are blocked on the other
  two tracks — a future fix PR for them could proceed independently.
- **No re-litigation of the aromaticity round-trip mechanism itself** — that
  belongs to the parallel aromaticity-diagnosis effort; this doc only
  reconfirms the failure rate and its molecule-preserving (non-structural-
  drift) nature at the current commit.
- **No changes to any other agent's branch, PR, or files.**

## Reproduce

```bash
# Self-test (unit-level, no corpus needed — verifies the classification
# decision tree discriminates every bucket, plus one positive control: the
# known bicyclo[2.2.2]octane permutation-invariance bug must still reproduce)
.venv/bin/python scripts/canonical_residual_diagnosis.py --self-test

# Full corpus (reproducible: same --seed => byte-identical JSON output)
.venv/bin/python scripts/canonical_residual_diagnosis.py ~/Downloads/SMILES.csv --k 8 --seed 0
```

Writes `validation/results/canonical_residual_diagnosis.jsonl` (one row per
Check-1 mismatch, with bucket + detail) and
`validation/results/canonical_residual_diagnosis_summary.json` (headline
numbers, bucket examples, and the full permutation-invariance failure
sample). Frozen fixtures (hand-constructed + representative corpus-sampled,
each independently re-verified to reproduce): `validation/
canonical_residual_fixtures.jsonl`.
