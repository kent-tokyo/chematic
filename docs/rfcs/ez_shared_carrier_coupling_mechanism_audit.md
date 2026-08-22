# Wave 3 audit: general shared-carrier coupling-component residual (issue #149)

Status: **Diagnosis only. No production behavior change.** `Refs #149`.

## Scope

PR #351 (Wave 2D, merged 2026-08-20) implemented the ring-endocyclic-in-a-
ring-smaller-than-8-atoms exclusion in `compute_stereo_alkene_ends`
(`crates/chematic-smiles/src/canonical.rs`), closing all 18 pinned
`EZ_SHARED_CARRIER_FULLY_RESOLVED` fixtures. Its own commit message and
`docs/rfcs/ez_ring_constrained_residual_audit.md`'s verdict section are
explicit that this closes only **~10% (3 of 31) of the corpus's general
shared-carrier coupling-component population** — the other **~90% (28 of
31)** was reported as "a separate, still-unidentified mechanism," and issue
#149 was kept open on that basis.

This audit measures that remaining population directly, rather than
inheriting its size, shape, or residual status from the ~90% figure — which
was a **topological-presence count** (does this molecule have a coupling
component at all), not a **confirmed-permutation-invariance-failure count**
(does chematic's canonical output actually change across valid respellings
of it). No production code is touched or generalized by this audit. Does not
conflate this population with `ROADMAP.md` backlog item 6's separate
"346 abstained E/Z bonds" ledger (85 lost-in-canonicalization + 261
CarrierConflict, an unrelated, unchanged code path).

**Files created:**
`crates/chematic-smiles/examples/ez_shared_carrier_coupling_mechanism_audit.rs`
(new Rust example, public-API-only, three subcommands: `scan`/`axis1`/`axis2`),
`scripts/ez_shared_carrier_coupling_mechanism_diagnosis.py` (driver),
`validation/results/ez_shared_carrier_coupling_mechanism_audit.jsonl` +
`_summary.json` (this audit's raw output), this doc. **Files touched:** none
under `crates/*/src/**`.

## Why a new tool, not either existing example

`ez_shared_carrier_component_audit.rs` (Wave 2B) and
`ez_ring_constrained_residual_audit.rs` (Wave 2C) both predate PR #351 and
reimplement `stereo_alkene_end_nodes` **without** its ring-endocyclic gate.
Re-running the former directly against `scripts/descriptor_census_corpus.smi`
today still reports the pre-fix topology (31 components, confirmed by
direct re-run, unchanged). This audit's new example mirrors the **current**
`compute_stereo_alkene_ends` exactly (including PR #351's gate), so its
topology reflects `main` as it stands today.

## Method

1. **Provenance gate** — do the 18 pinned fixtures / the 2 never-corrupts
   SMILES appear verbatim in the committed 5,000-mol corpus
   (`scripts/descriptor_census_corpus.smi`)? Measured: **5/18** pinned
   fixtures appear in the corpus (the 5 hydrazone-imine ones, matching Wave
   2C's own finding); **0/2** never-corrupts SMILES appear. Fixture-derived
   and corpus-derived populations are never pooled below.
2. **Current-topology scan** (`scan` subcommand, ring-gate-aware) —
   measured **28 coupled (size ≥2) components** in the corpus today, **every
   one size exactly 2, every one shape "path," 0 cycles** — exactly matching
   PR #351's own "28 of 31" claim (31 pre-fix − 3 closed by the ring gate =
   28). This also gives a fresh, current re-confirmation of the branch-point
   proof (below).
3. **Axis 1 (RDKit relabeling)** — 16 seeded (`random.Random(seed)` for
   `seed in range(16)`, matching this project's own Check-2 convention),
   reproducible `Chem.RenumberAtoms` relabelings per coupled molecule (448
   total variants across the 28 components), fed through chematic's real
   `canonical_smiles()`. **0/28 molecules show more than one distinct
   canonical output across their 16 relabelings; 0/448 cross-spelling
   correspondence checks failed.**
4. **Axis 2 (single-end mark relocation)** — reimplements the private
   `alternate_ez_markings` test helper via public `Molecule::with_bond_order`.
   **Measured finding: this probe is structurally incapable of testing any
   genuinely coupled 2-node component.** Traced mechanistically on a real
   corpus example
   (`CCOC(=O)/C=C(C)/C(F)=C/C=C(C)/C=C/c1c(C)cc(C)c(Cl)c1C`, one of the 5
   acyclic components): the shared bond between the two ends is the one
   currently carrying the directional mark for *both*; relocating either
   end's mark strips that shared bond's mark, which the other end's own
   geometry reading depends on — `alternate_ez_markings`'s own
   geometry-preservation check (`geometry_fingerprint(&alt) != baseline_geo`)
   rejects every such move. Confirmed empirically too: **0 alternates
   generated for any of the 28 coupled components or the never-corrupts
   pair** (both candidates being `Aromatic`-stashed for 23/28 components is
   a second, independent reason axis 2 doesn't fire there — but the acyclic
   5/28, whose candidates are literal `Single`/`Up`, show the *same* zero
   result for the structural reason above, proving the limitation is
   general, not aromatic-representation-specific). Only the singleton
   negative-control fixture (`CC1=C2CC[C@H](/C=N/N=C(N)N)...`, a known-fine
   fixture sharing the ring-endocyclic *shape* but never a residual)
   produced any alternate at all (1), and it did not diverge — confirming
   the tool correctly finds nothing to flag on a known-good case, as
   expected, not that the tool is silently broken.
5. **RDKit stereogenicity oracle** (`Chem.FindPotentialStereo` +
   `Chem.AssignStereochemistry`, same index-correspondence-first discipline
   as Wave 2C's `crosscheck_row`) — **28/28 coupled components have BOTH
   ends independently, genuinely RDKit-`Specified`** (real stereocenters).
   None show the asymmetric "one fake ring-endocyclic + one genuine" shape
   PR #351's predicate already targets — as expected, since that shape is
   now fully excluded from `compute_stereo_alkene_ends`'s candidacy.
6. **Structural classification** (measured feature tuple, buckets named
   only once occupied — see Results).
7. **Calibration cross-check** on the permanent regression fixture
   `ez_carrier_shared_bond_between_two_stereo_systems_never_corrupts`
   (canonical.rs) — see "A stale doc-comment claim, found and disclosed"
   below.

## Results

### Structural classification

Two buckets, both sharing the same RDKit-stereogenicity signature:

| Bucket | Count | Ring? | Candidate-bond representation | Both ends RDKit-`Specified`? |
|---|---|---|---|---|
| Aromatic-ring-adjacent | 23/28 | has_ring | `Aromatic` (direction stashed via `Molecule::bond_direction`) | yes |
| Acyclic conjugated | 5/28 | no_ring | literal `Single`/`Up` | yes |

**One mechanism, not several.** Every one of the 28 components — in both
buckets — is the same underlying shape at the RDKit-stereogenicity level:
two independently, genuinely stereogenic double bonds sharing one physical
candidate bond. This is the identical shape already known to be handled
correctly by `resolve_component_jointly` for 5 of the 18 pinned fixtures
(the hydrazone-imine family, `hypothesis_holds = False` in Wave 2C's own
audit) and, per this audit's calibration check, for the never-corrupts
fixture too. The ring/aromatic-vs-acyclic split is a real, secondary
structural distinction (worth naming), but it does not correspond to a
different *failure* mechanism — no failure was found in either bucket.

### Divergence measurement: 0/28

Axis 1: 0/28 molecules diverge across 16 relabelings each (448 variants, 0
correspondence failures). Axis 2: structurally inapplicable to any of the
28 (see Method §4) — its 0-alternates result is not informative either way
for coupled pairs, only confirmatory for the singleton negative control.

### A stale doc-comment claim, found and disclosed

`ez_carrier_shared_bond_between_two_stereo_systems_never_corrupts`
(canonical.rs) pins two spellings of a real corpus-derived molecule (two
exocyclic imines on a four-membered ring sharing the ring-closure bond) with
a doc comment stating: *"This is a real corpus molecule that a fully
general carrier choice does NOT resolve to one canonical string (a known,
documented residual...)."* The test itself only asserts **no corruption**
(`ez_before == ez_after` for each spelling independently) — it does not
assert non-convergence.

This audit measured convergence directly: both pinned spellings, `OC(=O)
[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c1/c(c(c1O)O)=N/CCCCC` and `OC(=O)
[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c\1c(/c(c1O)O)=N/CCCCC`, produce
the **identical** canonical output
(`c3ncc(c(c3Cl)C(Nc1ccc(C[C@H](/N=c\2c(/c(c2O)O)=N/CCCCC)C(=O)O)cc1)=O)Cl`),
and 10 fresh RDKit relabelings of the first spelling all agree with it too
(0/10 divergent). **This claim appears stale** — most likely superseded by
PR #229's joint-component solver (`resolve_component_jointly`, designed
generally for exactly this "two independently stereogenic ends sharing a
bond" shape) landing after the test/comment were originally written,
without the comment being revisited. **Not fixed here**: even a
comment-only edit to `canonical.rs` is a change under `crates/*/src/**`,
outside this diagnosis-only PR's declared scope. Recommended as a tiny,
separate, low-risk follow-up (re-verify + update or retire the stale claim;
the test's own no-corruption assertion remains correct and unaffected
either way).

### Field-by-field findings (per the original request)

1. **Canonical rank per atom + relabeling consistency** — per-end ranks via
   `morgan_ranks`/`canonical_atom_order` (descriptive); the real invariance
   verdict is axis-1/axis-2 above (measured, not inferred from the
   topological ~90% figure).
2. **Candidate carrier bonds per end** — trivial, 2 per end by construction
   (`candidate_bonds` field in the `scan` JSONL).
3. **Coupling conflict graph, explicit per component** — emitted per
   `component` row (`members`, `edges: [{a,b,bond_idx}]`, `shape`).
4. **Marker constraint propagation order** — `resolve_component_jointly` is
   exhaustive `2^k` enumeration plus a rank-based minimal-deviation
   tie-break, **not** sequential constraint propagation, at any size — a
   structural fact about the algorithm, confirmed by direct source
   inspection, not a size-2 degeneracy. `propagation_order_applicable:
   false` on every emitted row. The one node order that *is* genuinely
   outcome-relevant is the tie-break comparison order (`tie_break_order`
   field — the rank-ascending sequence the real algorithm sorts by before
   enumerating).
5. **First branch point per component** — proven, not merely observed: every
   node has exactly 2 candidate substituents by the ambiguity precondition
   itself (`substituents(...).len() == 2`), so the coupling graph's maximum
   degree is 2 for any molecule — a branch point (degree ≥3) is
   **impossible**, not just absent from this corpus. Freshly re-confirmed
   empirically too: all 28 current components are size exactly 2, shape
   "path," 0 cycles.
6. **First RDKit/chematic divergence point** — three tiers:
   - *Tier 1 (well-defined)*: RDKit stereogenicity vs. chematic ring-gated
     candidacy for the same bond — measured above (28/28 agree: both
     genuinely stereogenic, correctly candidate on both sides).
   - *Tier 2 (well-defined output-level proxies)*: geometry-preservation
     (confirmed via the correspondence-trick reparse) and
     first-divergent-relabeling-index (N/A here — 0 divergence found).
   - *Tier 3, explicitly NEEDS-RESEARCH*: a true shared-intermediate-state
     comparison is not well-defined below tier 1 — chematic (rank-based
     joint enumeration) and RDKit (CIP-based canonical ranking) run
     unrelated algorithms with no shared intermediate representation, and
     RDKit's writer has no explicit per-substituent carrier-choice concept
     to compare against at all. Stated as a structural limitation, not
     forced into an answer.

## Verdict: NEEDS-RESEARCH, leaning GO (already-likely-resolved)

0/28 coupled components show canonical-output divergence under axis 1
(RDKit relabeling, K=16, 0 correspondence failures). Axis 2 cannot test
coupled pairs at all (a structural limitation, not a gap introduced by this
audit). The previously-cited never-corrupts calibration example — itself an
instance of this exact shape — also converges. Together this points toward
**"already resolved by the existing joint solver,"** contradicting PR #351's
own pessimistic "~90% unidentified mechanism" framing (a topological-
presence count, not a confirmed-failure count).

**This is not treated as proof.** K=16 relabeling is a sample, not
exhaustive (matching this project's own standing caveat on every prior
permutation-invariance check). Whether RDKit's own relabel-and-reserialize
process ever varies *which specific bond* carries the mark for a *shared*
bond specifically (as opposed to other, non-shared bonds it's also free to
remark) was not independently confirmed — this is the one open question
axis 1's "0 divergence" result cannot rule out.

**Recommended next step** (not a fix, and not done here): hand-construct 2–3
genuine alternate spellings that explicitly move the mark to the *other*
candidate bond by direct SMILES-text authorship (bypassing axis 2's
single-end-relocation limitation entirely, in the same style the
never-corrupts fixture's own `a`/`b` pair was originally hand-built) for a
sample of the 28 measured components, and check whether chematic still
converges. If it does, issue #149 could very plausibly be closed outright
(pending a decision on the stale never-corrupts doc comment); if it does
not, that would be the first genuine repro of the "other 90%" mechanism the
project has been looking for since PR #351.

## Reproduce

```bash
cargo run -p chematic-smiles --release --example ez_shared_carrier_coupling_mechanism_audit -- scan
cargo run -p chematic-smiles --release --example ez_shared_carrier_coupling_mechanism_audit -- scan scripts/descriptor_census_corpus.smi
.venv/bin/python3 scripts/ez_shared_carrier_coupling_mechanism_diagnosis.py --self-test
.venv/bin/python3 scripts/ez_shared_carrier_coupling_mechanism_diagnosis.py
```

Writes `validation/results/ez_shared_carrier_coupling_mechanism_audit.jsonl`
(per-component axis1/axis2/classification detail) and
`..._summary.json` (topology, axis summaries, structural buckets,
calibration check, verdict).
