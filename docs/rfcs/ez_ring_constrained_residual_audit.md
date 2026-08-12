# Wave 2C audit: ring-constrained E/Z residuals (issue #149)

Status: **Diagnosis only. No production behavior change.** `Refs #149`.

This is a follow-up to PR #229 (Wave 2B), which replaced the single-end
shared-bond abstain guard with a joint component solver
(`resolve_component_jointly` in `crates/chematic-smiles/src/canonical.rs`)
and split the 18 `EZ_SHARED_CANDIDATE_BOND_RESIDUALS` fixtures from
`docs/rfcs/canonical_smiles_residual_rfc.md` into 10
`EZ_SHARED_CARRIER_FULLY_RESOLVED` (fully permutation-invariant) and 8
`EZ_SHARED_CARRIER_RING_CONSTRAINED_RESIDUALS` (still not permutation-
invariant, but confirmed semantically safe — 0 corruption either way).

The doc comment directly above `EZ_SHARED_CARRIER_RING_CONSTRAINED_RESIDUALS`
proposes a hypothesis: every one of the 8 residuals' coupled components
includes an alkene end whose own C=X double bond is *endocyclic* in a 5- or
6-membered ring — real-world cis/trans fixed by the ring topology, not a
free stereochemical choice — which `compute_stereo_alkene_ends` has no gate
for. This audit tests that hypothesis empirically, classifies each of the 8
individually (never pooled), and measures the corpus-wide blast radius of
three candidate gating rules that could exclude such ends from candidacy.

**No production code was changed to produce this doc.** Everything below
comes from:

- `crates/chematic-smiles/examples/ez_ring_constrained_residual_audit.rs`
  (new, read-only Rust example — public `chematic_core`/`chematic_smiles`/
  `chematic_perception` API only)
- `scripts/ez_ring_constrained_residual_diagnosis.py` (new, drives the Rust
  example via `cargo run`, independently cross-checks every row against live
  RDKit 2026.03.3)
- `validation/results/ez_ring_constrained_residual_audit.jsonl` (1,387
  per-end rows, full `scripts/descriptor_census_corpus.smi` corpus)
- `validation/results/ez_ring_constrained_residual_audit_summary.json`
  (fixture classification + blast-radius table)

Reproduce:

```bash
.venv/bin/python3 scripts/ez_ring_constrained_residual_diagnosis.py --self-test
.venv/bin/python3 scripts/ez_ring_constrained_residual_diagnosis.py
```

## How `marker_placed` is measured without touching private API

`canonical.rs`'s `ez_marker` map is `pub(crate)`, unreachable from an
example. The Rust example instead composes the public
`chematic_smiles::canonical_atom_order` (a pure, relabeling-invariant
function of molecule structure — the same `winning_individualized_ranks`
`canonical_smiles` itself uses) computed on both the original molecule and a
reparse of its own canonical output, to recover atom correspondence across
the round trip with no molecule mutation. This is safe because
`initial_invariant` (canonical.rs line ~451) collapses `Single`/`Up`/`Down`/
`Dative` to the identical bond-order class for ranking purposes — rank
computation cannot see *which* candidate bond carries a mark, only that a
bond exists.

An earlier design tried tagging atoms via `atom_map` for correspondence and
was rejected: `canonical_partition.rs` folds `atom_map` into the initial
invariant, so tagging would have perturbed the very ranks/marker-choice this
audit is trying to observe.

The composition is verified per molecule, never assumed: every one of the
original molecule's bonds is replayed through the mapping and checked
against the reparsed molecule (same atom pair, same bond-order class).
**6 of 1,387 corpus rows (0.4%) fail this verification** — all 6 are
singleton (non-coupled) endocyclic ring-fusion alkene ends in a
gem-dimethyl-tetrahydroquinoline scaffold family, where a genuine
automorphism-orbit rank tie apparently picks a different winning branch
between the original molecule and the reparse of its own canonical output.
These rows have `marker_placed: null` and `correspondence_ok: false`, and
are excluded from every marker-placement-based count below (not guessed).

## RDKit oracle notes

RDKit 2026.03.3's `Chem.StereoSpecified` enum has only `Unspecified`,
`Specified`, `Unknown` — **no `NOT_POSSIBLE` value exists in this version**
(confirmed by direct enum dump). The operational equivalent, confirmed
mechanistically: `Chem.FindPotentialStereo(mol, cleanIt=False,
flagPossible=True)` simply omits a `Bond_Double` entry for a bond it
considers structurally incapable of real E/Z:

| SMILES | ring size | `FindPotentialStereo` entry |
|---|---|---|
| `CC1=C(C)CCCC1` (1,2-disub cyclohexene) | 6 | **absent** |
| `CC1=C(C)CCCCCC1` (1,2-disub cyclooctene) | 8 | present, `Unspecified` |

"Absent from the list" is the signal rules (b)/(c) below use. `Chem.
AssignStereochemistry(mol, cleanIt=True, force=True)` independently confirms
this: the excluded bonds come back `bond.GetStereo() == STEREONONE`, never
an assignable E/Z value.

RDKit atom-index correspondence (chematic parse index vs RDKit parse index,
same input string) was verified explicitly per row (element-symbol match on
both the end atom and its double-bond partner), never assumed:
**0 of 1,387 rows failed** this check.

## Classification of the 8 residual fixtures (individual, not pooled)

Every one of the 8 shows the SAME shape: a coupled pair where one end's own
double bond is endocyclic and RDKit independently confirms it is not a real
stereocenter (`STEREONONE`, absent from `FindPotentialStereo`), coupled via
the shared candidate bond to a genuinely, independently stereogenic partner
(real `STEREOE`/`STEREOZ`).

| # | SMILES (truncated) | culprit end (atom idx, element) | culprit ring size | RDKit verdict on culprit bond | coupled partner (atom idx) | RDKit verdict on partner bond |
|---|---|---|---|---|---|---|
| 1 | `CC1=C2CC[C@H](/C=N/N=C(N)N)...` | 1 (C), ring C1=C2 | 6 | not potential, `STEREONONE` | 16 (C), exocyclic C=N | potential, `STEREOE` |
| 2 | `CC1=C2CC[C@@H](/C=N/N=C(N)N)...` | 1 (C), ring C1=C2 | 6 | not potential, `STEREONONE` | 16 (C), exocyclic C=N | potential, `STEREOE` |
| 3 | `COC(=O)/C=C/[C@H]1CCC2=C(C)/C(=N/...)...` | 10 (C), ring C=C | 6 | not potential, `STEREONONE` | 12 (C), exocyclic C=N | potential, `STEREOE` |
| 4 | `CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccc(I)cc1` | 9 (C), ring C=N | 6 | not potential, `STEREONONE` | 11 (C), exocyclic C=C | potential, `STEREOZ` |
| 5 | `CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]` | 6 (C), ring C=N | 6 | not potential, `STEREONONE` | 5 (C), exocyclic C=C | potential, `STEREOZ` |
| 6 | `CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(C)c1` | 9 (C), ring C=N | 6 | not potential, `STEREONONE` | 11 (C), exocyclic C=C | potential, `STEREOZ` |
| 7 | `CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(OC)c1` | 9 (C), ring C=N | 6 | not potential, `STEREONONE` | 11 (C), exocyclic C=C | potential, `STEREOZ` |
| 8 | `CCO/C(O)=C(\C1=NCCN1)c1nnc(N)s1` | 6 (C), ring C=N (imidazoline) | **5** | not potential, `STEREONONE` | 5 (C), acyclic C=C | potential, `STEREOE` |

All 8/8 confirmed: `hypothesis_holds = True` for every one, individually
(never pooled — see `fixture_classification` in the summary JSON for the
full per-fixture detail this table summarizes). Fixture 8's culprit ring is
a **5-membered** imidazoline, and its coupled partner is fully acyclic (not
merely exocyclic-but-ring-adjacent like the others) — a structurally
distinct shape within the family, called out explicitly per this session's
rule against rounding differences away.

## The hypothesis is necessary but NOT sufficient — a real finding, not an artifact

The same "endocyclic-culprit + real-stereo-partner" shape, individually
RDKit-confirmed the same way as the 8 above, is **also present in 5 of the
10 `EZ_SHARED_CARRIER_FULLY_RESOLVED` fixtures**:

| SMILES (truncated) | culprit end | RDKit verdict | still fully resolved? |
|---|---|---|---|
| `CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(I)c1` | 9 (C), ring C=N, 6-ring | `STEREONONE`, not potential | yes |
| `CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1C(F)(F)F` | 6 (C), ring C=N, 6-ring | `STEREONONE`, not potential | yes |
| `CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccccc1OC` | 9 (C), ring C=N, 6-ring | `STEREONONE`, not potential | yes |
| `CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccc([N+](=O)[O-])cc1` | 6 (C), ring C=N, 6-ring | `STEREONONE`, not potential | yes |
| `CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1cccc(C(F)(F)F)c1` | 6 (C), ring C=N, 6-ring | `STEREONONE`, not potential | yes |

These 5 fixtures share the **exact same core scaffold** as 4 of the 8
residuals (`CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1...` / `CCO/C(O)=C1\C(C)=NC
(C)=C(C(=O)OC)C1...`), differing only in the pendant phenyl ring's
substituent identity and position (e.g. meta-iodophenyl → resolved vs.
para-iodophenyl → residual; para-nitrophenyl → resolved vs.
ortho-nitrophenyl → residual). No consistent ortho/meta/para symmetry
pattern was found across all pairs that explains the split (checked, not
assumed — the iodo pair and the nitro pair invert the same
ortho/meta/para axis relative to resolved/residual). **What
specifically differentiates these 5 from the 8 residuals was not
identified in this audit** — it is a genuinely separate, second mechanism
(likely an interaction with tie-breaking or automorphism elsewhere in the
molecule), out of scope for this audit's own hypothesis test, and named
here explicitly rather than rounded into a clean story.

The other 5 of the 10 `EZ_SHARED_CARRIER_FULLY_RESOLVED` fixtures (the
`.../N=c1\c(O)c(O)\c1=N/...` hydrazone-imine family) use a **completely
different coupling mechanism**: `hypothesis_holds = False` for all 5 — both
coupled ends are independently, genuinely stereogenic exocyclic imines
(both `STEREOZ`), with no endocyclic culprit at all. This is a different,
unrelated shared-carrier coupling shape that the current solver already
handles correctly, and is untouched by anything measured here.

**Practical consequence**: a gating rule that excludes endocyclic-small-ring
ends from `compute_stereo_alkene_ends` candidacy would change the current
marker-placement input for **13 of the 18 pinned fixtures** (5 resolved + 8
residual), not just the 8 residual ones. For the 5 already-resolved
fixtures, the excluded end currently carries a mark (`marker_placed: true`)
in every case — see the blast-radius table below,
`blast_radius_per_fixture_effect.a_bond_endocyclic_lt8.fully_resolved` in
the summary JSON. Whether those 5 would *remain* resolved after such a
change is not measured here (that requires implementing the rule, which is
out of scope) — but the exclusion is real, not merely theoretical.

## Full-corpus blast radius (5,000-molecule `scripts/descriptor_census_corpus.smi`, 1,387 stereo-alkene ends, 995 molecules, 62 coupled ends / 31 coupling components)

Each rule measured independently, never combined. "Excluded" = this rule
would remove the end from `compute_stereo_alkene_ends` candidacy.
"Excluded-and-marked" is an explicit **upper bound** on output change, not a
measurement of it — determining whether output would actually change
requires implementing the rule, which this diagnosis-only PR does not do.
The **confirmed** column counts only rows where chematic's own atom
correspondence check succeeded (`marker_placed` is `true`); the **incl.
unknown** column additionally credits the 6 rows (all endocyclic, all
excluded by every rule below) where correspondence failed and marker status
is genuinely unknown — reported as a range, not a single number, so an
unmeasured row is never silently folded into "not marked."

| Rule | Ends excluded | % of ends | Excluded & marked, confirmed (upper bound) | Excluded & marked, incl. unknown (true upper bound) | Excluded & coupled | Distinct molecules affected | Distinct molecules affected & marked |
|---|---|---|---|---|---|---|---|
| (a1) end atom in ring size < 7 | 901 | 65.0% | 62 | 68 | 52 | 660 | 62 |
| (a1) end atom in ring size < 8 | 959 | 69.1% | 62 | 68 | 52 | 709 | 62 |
| (a2) **double bond endocyclic** in ring size < 7 | 783 | 56.5% | 10 | 16 | 6 | 580 | 10 |
| (a2) **double bond endocyclic** in ring size < 8 | 837 | 60.3% | 10 | 16 | 6 | 625 | 10 |
| (b) RDKit: bond absent from `FindPotentialStereo` | 933 | 67.3% | 10 | 16 | 6 | 672 | 10 |
| (c) (b) restricted to endocyclic bonds | 837 | 60.3% | 10 | 16 | 6 | 625 | 10 |

All 6 marker-status-unknown rows are endocyclic 6-ring ends (the
gem-dimethyl-tetrahydroquinoline automorphism-tie family from the
correspondence-verification section above), so all 6 fall inside the
excluded set for every rule in this table — the "confirmed" and "incl.
unknown" columns for rules (a2)/(b)/(c) differ by exactly these 6 rows
(10 → 16), not a coincidence.

Key findings from this table, each verified directly (not inferred from
matching totals alone):

- **Rule (a1) ("atom merely sits in a small ring") is badly over-inclusive**:
  it excludes 901–959 ends (65–69% of the whole corpus!) because most ring
  carbons bearing an *exocyclic* double bond also get caught, even though
  their own double bond is genuinely free (exactly the atom-16 shape from
  fixture 1 — in a ring, but its own C=N points outward). This confirms the
  module doc comment's warning: ring *membership* of the atom is not the
  right predicate; ring *membership of the double bond itself* (endocyclic)
  is.
- **Rule (a2) ("double bond itself is endocyclic in a small ring") is the
  correctly-scoped structural predicate.** At threshold 8 it excludes 837
  ends (60.3%) vs. rule (a1)'s 959 ends (69.1%) at the same threshold — 122
  fewer ends, i.e. rule (a1)'s naive framing over-excludes by 122 ends
  (8.8 percentage points of the whole corpus) purely from the atom-vs-bond
  confusion.
- **Rule (a2) at threshold 8 and rule (c) (RDKit-confirmed, restricted to
  endocyclic bonds) agree on every single row of the corpus — 0 of 1,387
  mismatches**, checked row-by-row, not just at the aggregate-count level.
  This is strong, direct empirical support for "ring size < 8" as the
  concrete threshold (not merely "5- or 6-membered" as originally
  hypothesized in the doc comment — a 7-membered endocyclic ring is
  empirically indistinguishable from a 5-/6-membered one on this corpus,
  all excluded by RDKit the same way). **The corpus data has a real gap
  that this agreement does not close**, checked explicitly: endocyclic
  ring sizes observed in the corpus are `{4, 5, 6, 7, 16, 21, 24, 29}` —
  nothing between 8 and 15. The 897 endocyclic ends split as 837 in rings
  ≤7 (all RDKit-excluded) and 60 in rings ≥16 (all RDKit-*included*, 0
  exceptions), so the corpus independently confirms the *direction* on
  both sides of the threshold (small rings excluded, large macrocycles
  not) but does not, by itself, pin the exact cutover point anywhere
  between 8 and 15 — that specific boundary rests on the hand-constructed
  cyclooctene control above (ring size 8, RDKit includes it), not on
  corpus-observed data at size 8. "Ring size < 8" is the recommendation
  because it is the smallest threshold consistent with both the
  cyclooctene control and every corpus row; a value up to 15 would fit the
  corpus equally well and cannot be ruled out by this audit.
- **Rule (b) is genuinely broader than rule (c) — they do NOT coincide**:
  933 vs. 837, a difference of 96 ends (6.9% of the corpus). Spot-checked
  directly: every one of the 96 b-only exclusions is an *exocyclic* double
  bond RDKit excludes for an unrelated reason (identical substituents on
  both sides — e.g. guanidine `C(N)(N)=N-`, which has no real E/Z regardless
  of ring topology). Rule (c) correctly excludes these from its own count
  (they are not endocyclic), confirming (c) is a strictly narrower, more
  surgical measurement than (b), not a coincidental duplicate.
- **Only 6 of the corpus's 62 coupled ends (3 of 31 coupling components,
  verified by grouping on `(smiles, component_atom_idx_set)` rather than
  dividing the end-count in half)** are excluded under the recommended rule
  (a2/c at threshold 8). Unlike the pinned fixtures — where the coupled pair
  always has exactly ONE endocyclic end and one genuinely free end — these
  3 corpus components are a structurally different shape: fused,
  many-ring cage-like molecules (`COCCOCCOCCN1CC23C4=C5C6=C7...`, a
  crown-ether-substituted polycyclic aromatic hydrocarbon family) where
  **both** ends of the coupled pair are endocyclic in a 6-ring, contributing
  2 excluded ends per component. The remaining 28 components' coupling
  arises from some other mechanism entirely (plausibly, but not confirmed
  here, the "both ends genuinely stereogenic" shape seen in 5 of the pinned
  fixtures) — the ring-endocyclic mechanism explains all 8 pinned residual
  fixtures individually, but only a **minority (~10%)** of this particular
  corpus's general coupled-end population.

  The 18 pinned fixtures and the corpus's own coupling components are
  **not** fully disjoint sample sets, checked via an exact-string set
  intersection over all 18 (not a partial grep): **5 of the 18 fixture
  SMILES appear verbatim in the 5,000-molecule corpus** — specifically all
  5 of the "both-ends-genuinely-stereogenic" `EZ_SHARED_CARRIER_FULLY_
  RESOLVED` hydrazone-imine fixtures (`hypothesis_holds = False` in the
  classification above), each contributing 1 of the corpus's 31 coupling
  components. None of the 13 fixtures relevant to the endocyclic hypothesis
  (5 resolved-with-the-shape + all 8 residual) are present in the corpus —
  so the "minority (~10%)" finding is a genuine, independent measurement of
  the corpus's own population, not inflated or deflated by fixture overlap
  with the endocyclic mechanism specifically.

## Recommended production predicate (NOT implemented — specification only)

Add a guard to `compute_stereo_alkene_ends` (or equivalently to
`end_has_substituent`/the per-end loop) that excludes an end when:

```text
the atom's own C=X double bond is endocyclic — i.e. some SSSR ring (from
chematic_perception::find_sssr) contains BOTH atoms of the double bond
  AND
the smallest such ring has fewer than 8 atoms
```

Concretely: for a stereo-alkene end candidate `end` with double-bond partner
`partner`, compute `rings = find_sssr(mol).rings()`, and exclude `end` when
`rings.iter().any(|r| r.contains(&end) && r.contains(&partner) && r.len() <
8)`. This is exactly the `double_bond_endocyclic`/`double_bond_endocyclic_
ring_sizes` predicate the Rust example already computes for this audit,
generalized from ring size < 7 to < 8 per the row-level RDKit agreement
found above.

**Not "exclude all atoms in rings of size < 8"** (rule a1) — that is the
over-inclusive naive reading empirically ruled out above.

## Verdict: CONDITIONAL GO

The core hypothesis is confirmed, individually, for all 8 residual
fixtures, cross-checked against an independent RDKit oracle at the level of
the specific bond involved (not just aggregate counts). The recommended
predicate (ring size < 8, endocyclic) has row-level 100% agreement with
RDKit's own independent stereo-possibility judgment across the full
5,000-molecule corpus (0/1,387 mismatches) — unusually clean empirical
support, though (per the corpus-data-gap note above) the corpus itself only
pins the threshold's *direction* on both sides (≤7 excluded, ≥16 included),
not the exact cutover between 8 and 15; the specific value 8 additionally
relies on the hand-built cyclooctene control.

This is **not an unconditional GO** for three reasons, none of which
contradicts the hypothesis but all of which bound its scope honestly:

1. **The predicate is necessary, not sufficient.** 5 of the 10 already-fully-
   resolved fixtures share the identical endocyclic-culprit shape and would
   also have a currently-marked end excluded by this rule. A future
   implementation must re-verify those 5 stay permutation-invariant after
   the change (not merely re-verify the 8 residuals converge) — this audit
   does not measure that, by design (no production code was changed to
   produce it).
2. **The remaining 71% of corpus-general coupling components (28 of 31) are
   untouched by this predicate.** It closes the specific mechanism behind
   the 8 pinned fixtures, not the general "two stereo-alkenes share a
   carrier bond" coupling problem — issue #149 should stay open for the
   other mechanism(s) after this predicate ships, not just for these 8.
3. **A second, unidentified factor governs the 5-vs-8 split** among fixtures
   sharing the endocyclic shape. A correct implementation might still leave
   some of those 13 fixtures non-convergent for a different reason — this
   audit cannot rule that out without implementing and re-measuring.

**Recommendation**: proceed with a follow-up PR implementing the predicate
above, scoped explicitly to include a permutation-invariance re-check of
*all 18* pinned fixtures (not just the 8), not just the 10/8 split assumed
going in — and keep issue #149 open afterward pending the other ~90% of the
general corpus coupling population.
