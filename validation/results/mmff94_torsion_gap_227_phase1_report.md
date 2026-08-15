# Issue #227 Phase 1: MMFF94 torsion parameter gap — audit and resolution

Date: 2026-08-15. Corpus: 265 molecules
(`validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_{a,b}.json`).
Branch point: `c079926` (v0.16.0 release commit, `main`). RDKit oracle:
`rdkit==2026.03.4` (this session's actual installed version — see the
environment manifest below for why this does not confound the comparison).

Every number below names its producing tool and its denominator. `bonds`/
`angles`/`torsions`/`stretch_bend` "missing" figures are **type-only**
diagnostics from `mmff94_term_coverage_audit.rs` unless stated otherwise;
gate pass/fail figures name their exact producing policy call.

## 1 molecule excluded from every count

`force_field_unsupported_probe` (SMILES `[P](C)(C)(C)=C`,
`primary_category: "force_field_unsupported"` in its own manifest entry —
name and category confirm this is a deliberate fail-closed probe, not
investigated further) fails MMFF94 numeric typing (`atom 0 (Element(15))
was assigned MMFF94 numeric type 20 (CR4R), whose registry element is
Element(6)`) before torsion enumeration ever runs. Excluded from every count
in this report — an effective 264/265 typing-succeeded population — stated
here explicitly rather than silently narrowing every denominator below.

## Phase 1A — T0 audit (`mmff94_term_coverage_audit.rs`, pre-fix state)

| Metric | Value | Denominator |
|---|---|---|
| `torsions_missing` instances | 257 | across all torsion instances in 264/265 typing-succeeded molecules |
| `torsions_missing` molecules | 62 | / 264 typing-succeeded (Tier A=1, Tier B=61) |
| `present_at_different_classification = Some` | 254 | / 257 missing instances |
| `present_at_different_classification = None` (`table_gap`) | 3 | / 257 missing instances |
| `bonds_missing` (type-only) | 80 | across all bond instances / 264 |
| `angles_missing` (type-only) | 191 | across all angle instances / 264 |
| `stbn_type_only_missing` | 1865 | across all angle-triple instances / 264 (never gated by any `Mmff94BondAngleStrict` policy) |
| bond+angle gate-would-fail molecules | 14 | / 264 (this is `bond_angle_fully_covered() == false`, the type-only proxy the audit tool computes directly — NOT the same measurement as the `minimize_with_policy` gate rows below, which run the real production path) |

Classification of the 257 missing torsion instances against the directive's
required taxonomy:

- **Classification/routing-bug candidate** (`present_at_different_classification
  = Some`): 254/257. Root-caused below to a bond-order-source bug, not an
  equivalent-atom-type or wildcard-resolution gap in the lookup logic itself
  (chematic's existing exact+reverse+2-single-wildcard+double-wildcard+
  type0-generic chain was never the problem — it was never even being asked
  the right question).
- **Genuine table gap, resolvable via a specification-backed empirical
  rule's applicability condition**: 3/257 — but the applicable Halgren rule
  here is "omit the term entirely" (linear central atom), not a numeric
  formula. See `torsion_no_term_by_design` below.
- **Equivalent-atom-type-resolution gap**: 0 (investigated, falsified — see
  Phase 1B).
- **Wildcard-resolution gap**: 0 (the existing wildcard chain already
  reaches the correct row once fed the correct bond order).
- **Atom-typing-is-wrong**: 0 (types independently verified against the
  live oracle for all 257 instances — 100% match).
- **Input-chemistry-out-of-MMFF94-scope**: 0.
- **RDKit-itself-cannot-parameterize**: 0.
- **Cause-undetermined**: 0.
- **Duplicate/symmetric-torsion**: 2 of the 3 `table_gap` instances
  (`chembl_tier_b_0001`'s two nitrile-adjacent torsions, atoms
  `[19,20,21,22]` and `[23,20,21,22]`) share `j,k,l=(20,21,22)` and differ
  only in the terminal `i` atom (two symmetric aryl ring neighbors of the
  same ipso carbon) — genuinely distinct torsion instances, not a
  double-count, but structurally symmetric as flagged in the original brief.

## Phase 1B — oracle comparison

Environment manifest (see `validation/results/
pipeline_v2_vs_rdkit_environment_record.json` for the house-style full
manifest shape this reuses): `rdkit_version=2026.03.4` (differs from
`PROVENANCE.md`'s earlier 2026.03.3 pin; that prior file's own
`rdkit_version_isolation_control` entry already established this exact
version bump does not confound status/coverage transitions on this corpus),
`python_version=3.13.6`, `os_arch=aarch64-apple-darwin` (Apple M4, 10 cores),
`chematic_commit=c079926` base, corpus SHA unchanged since the 2026-08-06
baseline (`validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_{a,b}.json`,
byte-identical), no embedding/geometry involved (torsion parameter lookup is
purely topological — `GetMMFFTorsionParams`/`GetBondBetweenAtoms` need no
conformer), `mmff_variant=MMFF94` (not MMFF94s) on both sides.

For all 254 `present_at_different_classification = Some` instances (not a
sample — every one), and separately for all 3 `table_gap` instances:

1. **Two hypotheses tested and falsified** before the real cause was found
   (full detail: `scripts/mmff94_provenance/PROVENANCE.md`'s Torsion entry,
   tooling: `crates/chematic-3d/examples/mmff94_torsion_empirical_diagnostic_227.rs`):
   - A from-scratch Halgren empirical torsion rule (rule (b), "aromatic
     central bond" — structurally the case that matches all 254, since
     every one has `bond_order_jk == Aromatic`), derived from OpenBabel's
     public source comments citing Halgren Part V (the primary paper is
     paywalled). Predicts a uniform `(V1,V2,V3)=(0,6.0,0)` independent of
     the terminal atoms. **0/254 oracle matches.**
   - RDKit's real eqLevel terminal-atom substitution ladder (reusing
     Angle's own production `MMFF94_EQ_LEVEL` table). **0/254 additional
     table hits** beyond chematic's existing lookup chain.
2. **The real cause**: chematic's own, unmodified 926-row
   `MMFF94_TORSION_ENERGY` table already has a row matching the oracle's
   value exactly, at a DIFFERENT classification code than
   `torsion_type_for` computed — **254/254**, zero exceptions — because
   `torsion_type_for` was reading the wrong `BondOrder` for the central j-k
   bond (chematic's general perception: `Aromatic`; RDKit's real, Kekulized
   sanitizer output: `Single`/`Double` — oracle-confirmed
   `GetBondBetweenAtoms(j,k).GetIsAromatic() == False` for **254/254**).
   Where RDKit succeeds, this is squarely "chematic's classification input
   was wrong" — not an equivalent-type, wildcard, or empirical-rule
   insufficiency, and not RDKit itself falling back to anything (RDKit's own
   real table lookup succeeds directly at the correct code).
3. **The 3 `table_gap` instances**: oracle confirms `GetMMFFTorsionParams`
   returns `None` for all 3 too — genuinely unresolved on both sides, not a
   chematic-only gap. Root cause: each central atom (type 4 CSP or type 53
   `=N=`) has MMFF's `lin` flag; Halgren's real empirical-rule cascade omits
   the torsion term entirely for a linear central atom (rotating around a
   bond whose other end is a linear 180° center changes no real geometry).
   `torsion_no_term_by_design` implements this check; nothing here is left
   `unresolved` in the directive's sense — the mechanism is understood and
   oracle-confirmed, it just isn't a parameter to look up.

No case was left genuinely `unresolved` (the directive's stop-condition
category for "cannot determine why"). No empirical rule was guess-added to
cover an unclear case.

## Phase 1C — resolution order

**Implemented**: none of exact/equivalence/wildcard/empirical needed a new
tier. The existing exact→reverse→wildcard chain (`mmff94_torsion_energy`,
unmodified) already implements the directive's priority order correctly;
the bug was entirely upstream of it (wrong input, not wrong lookup logic).
Equivalence resolution (eqLevel ladder) was investigated (Phase 1B #2 above)
and found to contribute zero cases on this corpus, matching the same
"real mechanism, unexercised" verdict `PROVENANCE.md` already reached for
Bond/Angle's own eqLevel investigation before this PR — not built, to avoid
shipping dead code (ponytail: YAGNI). The empirical-rule tier was
investigated (Phase 1B #1) and its one falsifiable, central-bond-only
prediction failed against the live oracle on the population it was
supposed to cover — not implemented, per the explicit instruction not to
ship an unvalidated or already-falsified formula.

**What WAS implemented**: `assign_mmff94_numeric_types_with_view`
(`crates/chematic-ff/src/mmff94_numeric.rs`), fixing the bond-order input
`bond_type_for`/`angle_type_for`/`torsion_type_for`/`stretch_bend_type_for`
all consume, threaded through chematic-ff's 5 production energy/gradient
entry points and chematic-3d's `compute_mmff94_coverage`. Plus
`torsion_no_term_by_design`/`Mmff94Resolution::NoTermByDesign`, a
denominator correction (not a resolution) for the 3 linear-central-atom
cases.

## Phase 1D — diagnostic API

`Mmff94Resolution::NoTermByDesign` (`crates/chematic-ff/src/mmff94_energy/mod.rs`,
`pub`, matching the existing Bond/Angle Stage C `Mmff94Resolution` variants'
visibility). `Mmff94CoverageReport::torsions_no_term_by_design: usize`
(`crates/chematic-3d/src/minimize.rs`) is reachable from
`compute_mmff94_coverage`'s real production output and from
`mmff94_term_coverage_audit.rs`'s per-molecule/aggregate JSONL — a real
measurement, not an estimate. No `ParameterSource`/`TorsionResolutionTrace`
pair was added as a SEPARATE type: since no new resolution tier exists to
distinguish (the directive's own example enum's `EquivalentType`/`Wildcard`/
`EmpiricalRule` variants would all be permanently unreachable dead code on
this corpus), reusing the existing `Mmff94Resolution` enum with one new,
genuinely-reachable variant is the minimal faithful implementation.

## Phase 1E — tests

Added (see the PR's test-focused commit for the full list): both-sides/
neither-side/unknown-type-fails-closed unit tests for
`torsion_no_term_by_design`; the exact `chembl_tier_b_0001`/
`chembl_tier_b_0080` table_gap shapes verified to have no table row AND be
flagged `NoTermByDesign`; end-to-end resolution of caffeine's dione-ring
torsion through the real production call shape (`assign_mmff94_numeric_types_with_view`
→ `torsion_type_for` → `mmff94_torsion_energy`), matching the live oracle's
`GetMMFFTorsionParams` value exactly (not just the diagnostic tool);
reversal symmetry; a `Mmff94BondAngleStrict` +
`minimize_with_policy_gated(...,true,true)` integration test confirming
caffeine now passes the complete bonded-term gate directly, no UFF fallback;
`assign_mmff94_numeric_types_with_view` determinism (same input twice → same
types + same bond orders) and thin-wrapper equivalence with
`assign_mmff94_numeric_types`. **Atom-numbering permutation invariance**
(directive checklist item, reviewer follow-up): the determinism test above
only proves no hidden internal randomness on a FIXED atom order, not
renumbering-invariance — added
`caffeine_reperceived_bond_order_is_invariant_under_atom_renumbering` (32
`deterministic_permutation` relabelings, bond re-identified by unique
atom-type-pair content signature, not index) and
`benzene_reperceived_ring_bond_orders_are_invariant_under_atom_renumbering`
(the textbook genuine-Kekule-tie case) to check this directly. **Result:
32/32 renumberings identical on both — no renumbering-dependence bug
found**, a checked negative result, not assumed. Several items from the
directive's full checklist (equivalent-type resolution, wildcard
specificity, ambiguous wildcard, empirical-rule hit) are N/A given no new
resolution tier was implemented — already covered by this table's
pre-existing torsion tests where applicable (exact hit, reversed hit,
wildcard, ring torsion, permutation invariance of `torsion_type_for`
itself).

## Phase 1F — re-measurement

See `mmff94_torsion_gap_227_phase1_summary.json` (machine-readable,
same numbers) and `scripts/mmff94_provenance/PROVENANCE.md`'s Torsion entry
for the complete before/after with every number's producing tool named.
Headline: `torsions_missing` 257→0 instances (62→0 molecules);
`bonds_missing` 80→1; `angles_missing` 191→46 (side effects of the shared
root cause); `pipeline_v2_mmff94_strict`-equivalent bond+angle gate `Ok`
248→249/265 (+1); `complete_bonded_term_gate` `Ok` 187→249/265 (+62); zero
success→failure regressions on either gate, verified by a genuine
per-molecule-ID join (not count-level) against the pre-fix baseline.

## Files

- Machine-readable summary: `validation/results/mmff94_torsion_gap_227_phase1_summary.json`
- This report: `validation/results/mmff94_torsion_gap_227_phase1_report.md`
- Full T0 JSONL dump (pre-fix, all 265 molecules × all 5 term kinds),
  committed verbatim, matching the existing `_v0_13_0`-suffixed versioned-
  snapshot naming convention already in this directory:
  `validation/results/mmff94_coverage_227_term_audit_v0_16_0_prefix.jsonl`
  (+ its `_stderr.log`). Captured from `mmff94_term_coverage_audit.rs` at
  commit `c079926` (the v0.16.0 release commit, this PR's base) — reproduce
  independently via `git checkout c079926 -- crates/chematic-3d/examples/mmff94_term_coverage_audit.rs
  crates/chematic-ff && cargo run --release -p chematic-3d --example
  mmff94_term_coverage_audit` against a fresh checkout of that commit (not a
  stash trick on a working tree that has since changed).
- Investigation tooling (kept, documents the two falsified hypotheses):
  `crates/chematic-3d/examples/mmff94_torsion_empirical_diagnostic_227.rs`.
- Provenance: `scripts/mmff94_provenance/PROVENANCE.md` (Torsion entry).
