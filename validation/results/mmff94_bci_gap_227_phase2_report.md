# Issue #227 Phase 2: MMFF94 BCI partial-charge bug + 3-state `embed_pipeline_v2` quality re-measurement

Date: 2026-08-16. Corpus: 265 molecules
(`validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_{a,b}.json}`). Every
success-rate/RMSD/TFD/coverage/stereo/wall-time number below is measured
through the production `embed_pipeline_v2` entry point, arm
`chematic_pipeline_v2_mmff94_strict` (`pipeline_v2_mmff94_strict` in Phase
1's own naming) — **never** the coverage-gate-only tooling
(`mmff94_strict_gate_remeasure_227.rs`) Phase 1 used for its own headline
numbers. Where a number is Phase 1's own gate-tool number, it is named as
such explicitly.

## 0. Scope correction found during Step 1

The Phase 1 follow-up note ("`mmff94_charges_numeric` reads bond order from
the caller's original molecule, the same root-cause shape" as the
bond/angle/torsion/stretch-bend fix) was **half right**. Investigation
(`scripts/mmff94_provenance/PROVENANCE.md`'s Charges/BCI entry) found a
**compound** bug: a private, standalone `bond_type_for(order) -> u8` in
`mmff94_numeric.rs` encoded bond *multiplicity* (`Double->1, Triple->2,
Aromatic->4`), structurally unrelated to RDKit's real `getMMFFBondType`
(0 unless the bond is formally SINGLE *and* both atom types are
`sbmb`/`arom`-flagged) — a formula bug that fires even where the reperceived
and original bond order **agree** (e.g. every real C=O double bond), not
just where they differ. Fixed by deleting the wrong function and reusing
`crate::mmff94_minimizer::bond_type_for(ti, tj, order)`, the already-fixed,
oracle-validated function bond-stretch/angle/torsion/stretch-bend already
use, fed the reperceived (`assign_mmff94_numeric_types_with_view`) bond
order.

## 1. BCI oracle-check verdict: **FIXED** (real bug, not "already correct")

Source-level proof: RDKit's real `computeMMFFCharges`
(`AtomTyper.cpp:3071-3488`, pinned commit
`e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`) calls the identical
`getMMFFBondType(bond)` its own `getMMFFBondStretchParams` calls, on the
identical sanitized/Kekulized `mol`. There is no separate "charge bond type"
concept in RDKit's own algorithm.

Empirical falsification, full corpus (all 264 typing-succeeded molecules,
not a sample), live oracle `rdkit==2026.03.4`
(`scripts/mmff94_bci_charges_oracle_227.py`,
`crates/chematic-3d/examples/mmff94_bci_charges_dump_227.rs`):

| | Atoms mismatched | Molecules w/ any mismatch | mean \|Δ\| | p90 \|Δ\| | p99 \|Δ\| | max \|Δ\| |
|---|---|---|---|---|---|---|
| Before fix | 1,687 / 6,693 (25.2%) | 206 / 264 (78.0%) | 0.0189 e⁻ | 0.076 e⁻ | 0.239 e⁻ | 1.0 e⁻ |
| After fix | 67 / 6,693 (1.0%) | 11 / 264 (4.2%) | 0.00116 e⁻ | 0.0 e⁻ | 0.0144 e⁻ | 1.0 e⁻ (unrelated outlier) |

Per-atom join (before-state vs after-state, both against the same oracle
snapshot): **0 regressions** (0 atoms that matched the oracle exactly
before now mismatch), **1,620 improvements** (mismatched → exact match),
5,006 unchanged-match, 67 unchanged-mismatch (the residual, below).

**Residual, 67/6,693 atoms / 11/264 molecules, explicitly out of scope for
this fix**: the largest-magnitude residual molecules
(`chembl_tier_b_0080`/`_0159`/`_0161`) show chematic's and RDKit's MMFF atom
TYPES agreeing exactly at every mismatched atom, with a charge difference of
almost exactly 0.5 e⁻ — consistent with `mmff94_charges_numeric`'s
formal-charge/`fcadj` redistribution step for charge-separated species
(nitro/azide/charged-sulfoxide), NOT the bond-type BCI step this fix
addresses. Confirmed independent of this fix by construction (unmoved
before/after). Flagged as follow-up, not root-caused further here.

Full writeup, source citations, table-level sensitivity analysis:
`scripts/mmff94_provenance/PROVENANCE.md`'s Charges/BCI entry.

## 2. Three-state `pipeline_v2_mmff94_strict` measurement

Environment: `rustc 1.97.0`, `rdkit==2026.03.4` (matches the environment
record's pin — RDKit oracle rows (`pipeline_v2_vs_rdkit_rdkit_rows.jsonl`)
were **reused unchanged** across all 3 states, not regenerated, since RDKit
behavior depends only on the RDKit installation and the corpus, neither of
which changed between states). `aarch64-apple-darwin`, Apple M4, 10 cores.
Same tooling throughout: `crates/chematic-3d/examples/pipeline_v2_vs_rdkit_dump.rs`
(unmodified), run to completion in 3 separate, sequential (never
concurrent), machine-quiet git worktrees to avoid the CPU-contention/
timeout-boundary corruption failure mode `pipeline_v2_vs_rdkit_environment_record.json`
already documents for this exact tool.

- **State 1**: commit `c079926` (v0.16.0 release, pre-torsion-fix, pre-BCI-fix).
- **State 2**: commit `a2baac4` (current `main`, post-torsion-fix, pre-BCI-fix).
- **State 3**: this PR's branch tip (post-torsion-fix, post-BCI-fix).

| Metric | State 1 | State 2 | State 3 |
|---|---|---|---|
| `pipeline_v2_mmff94_strict` success | 240/265 (90.57%) | 241/265 (90.94%) | 241/265 (90.94%) |
| typed_failure | 25 | 24 | 24 |
| timeout | 0 | 0 | 0 |
| internal_error (crash) | 0 | 0 | 0 |
| RMSD mean / median / p75 / p90 / p95 (Å) | 1.698 / 1.531 / 2.536 / 3.339 / 3.909 | 1.685 / 1.531 / 2.494 / 3.339 / 3.909 | 1.685 / 1.504 / 2.493 / 3.456 / 3.876 |
| TFD mean / median / p75 / p90 / p95 | 0.2245 / 0.1873 / 0.3213 / 0.4694 / 0.5888 | 0.2233 / 0.1873 / 0.3213 / 0.4694 / 0.5888 | 0.2228 / 0.1797 / 0.3216 / 0.4779 / 0.5888 |
| Coverage @0.5 / 1.0 / 2.0 Å | 17.36% / 28.68% / 55.85% | 17.74% / 29.81% / 56.98% | 17.74% / 30.19% / 56.23% |
| Stereo declared / satisfied / violated | 146 / 83 / 63 | 146 / 83 / 63 | 146 / 82 / 64 |
| Stereo satisfaction rate | 56.85% | 56.85% | 56.16% |
| Wall time mean / p50 / p90 / p95 (ms) | 1923.9 / 1176 / 5125 / 6268 | 1843.4 / 1056 / 5084 / 6229 | 1827.9 / 1025 / 5051 / 6047 |

Best-of-1 vs best-of-10: **not measurable with existing tooling** —
`pipeline_v2_vs_rdkit_dump.rs` runs exactly one embedding-attempt sequence
per (molecule, arm); `MAX_ATTEMPTS=8` in `EmbedParameters` is
retry-on-failure, not a best-of-N conformer-selection loop. Building
best-of-N tooling is out of scope for this measurement PR per the
directive's own "use what's feasible" clause.

Wall time is flat-to-slightly-decreasing across all 3 states (no slowdown
from the BCI investigation or fix — the extra electrostatic-term work is
within noise).

Machine-readable per-state summaries:
`validation/results/mmff94_bci_gap_227_state{1,2,3}_report.json`.

## 3. Subset evaluations

### 3a. The 62 molecules Phase 1's torsion fix touched (State 1 → State 2)

Per-molecule join (`mmff94_bci_gap_227_transition_state1_to_state2_torsion62subset.json`):
0 regressions, 1 coverage improvement (caffeine, `typed_failure` →
`success`). RMSD 47/58 improved, 11/58 worsened; mean 1.156 → 1.115 Å,
median 1.210 → 1.128 Å; coverage@0.5/1.0/2.0 all increased (29.0%→30.6%,
38.7%→43.5%, 77.4%→82.3%). TFD roughly balanced by count (25 improved / 26
worsened) but flat-to-improved on aggregate (mean 0.130 → 0.127, median
unchanged at 0.114).

**Verdict: BOTH IMPROVED** — coverage AND geometry quality both moved in the
improving direction for this subset, not just coverage. This directly
answers Phase 2's central question: Phase 1's coverage-gate improvement
(87→17 failing on the gate tool, 62→0 `torsions_missing` molecules) is
accompanied by a real geometry-quality improvement on the production
pipeline for the same molecules, not offset by it.

### 3b. Molecules the BCI fix changed at least one charge on (State 2 → State 3, n=206)

Per-molecule join (`mmff94_bci_gap_227_transition_state2_to_state3_bciaffectedsubset.json`):
0 coverage regressions, 0 coverage improvements (expected — partial charges
do not gate `Mmff94BondAngleStrict`/`Mmff94WithUffFallback` coverage
eligibility, only the electrostatic energy term). RMSD 83/186 improved,
101/186 worsened, mean delta −0.0009 Å (essentially zero — a genuine
noise-level wash, not a directional shift; 12/186 molecules show a
meaningful >0.1 Å shift in either direction, consistent with a few cases
landing in a different local minimum, not a systematic quality change). TFD
92/186 improved, 84/186 worsened (also roughly balanced). One new stereo
violation: `chembl_tier_b_0082` (0 → 1 violation).

**Verdict: NO CHANGE** — coverage unaffected as structurally expected;
aggregate geometry quality is flat/noise-level in both directions, not a
quality-only improvement or a regression. The BCI fix is a real correctness
fix (Section 1) whose effect on final minimized geometry quality is small
and mixed-sign, exactly as expected for a fix to one term (electrostatics)
among several in the MMFF94 energy function — it changes which local
minimum the minimizer converges to for some molecules, not which basin is
reachable at all. The one new stereo violation is reported plainly, not
hidden: an isolated case, not part of a systemic pattern (0 regressions
among the other 205 molecules' coverage or stereo state).

## 4. Full-corpus transition table, State 1 → State 3 final

(`mmff94_bci_gap_227_transition_state1_to_state3_full.json`, genuine
per-molecule-ID join, matching `mmff94_torsion_gap_227_phase1_summary.json`'s
own `per_molecule_join_regressions`/`_improvements` methodology — not
aggregate-count arithmetic.)

| Transition | Count |
|---|---|
| success → success | 240 |
| typed_failure → success | 1 (caffeine) |
| typed_failure → typed_failure | 24 |
| success → typed_failure / timeout | **0** |
| timeout → anything | 0 (no timeouts in either state) |

- `per_molecule_join_regressions`: **0**
- `per_molecule_join_improvements`: **1** (caffeine)
- Pre-registered wall-clock-timeout-boundary jitter molecules
  (`chembl_tier_b_0166`/`_0114`/`_0117`, `atorvastatin_fragment`,
  `cholesterol` — see `pipeline_v2_vs_rdkit_3point_paired_diff_summary.json`'s
  own `known_jitter_molecules`/`byte_identical_verification`, which shows
  even a byte-identical re-run of the SAME commit flips these): **0 flips**
  this round — none of the numbers above are jitter-confounded.
- RMSD: 113/239 improved, 70/239 worsened; mean delta −0.0064 Å (net small
  improvement, driven by the torsion-fix subset above; median delta exactly
  0.0 since most of the corpus is untouched by either fix).
- TFD: 93/234 improved, 82/234 worsened.
- Stereo newly violated: 1 (`chembl_tier_b_0082`, same molecule as 3b —
  entered between State 2 and State 3, i.e. attributable to the BCI fix,
  not the torsion fix).

## 5. Stop-condition check

None triggered: RDKit-real-BCI-source was determined with high confidence
(source read + full-corpus oracle falsification, not a guess); success rate
did not worsen anywhere (240→241→241); no clear RMSD/TFD worsening
accompanies the success-rate improvement (Section 3a shows both improved
together); no molecule regressed from success to failure/timeout at any
transition (0 regressions, all 3 pairwise diffs); canonical SMILES output is
untouched by this PR (no canonicalization code was touched — confirmed by
`git diff --stat` scope, only `mmff94_numeric.rs`'s BCI function and new
measurement tooling changed); no unrelated refactor was needed; all 3
`embed_pipeline_v2` measurements completed to full 265-molecule corpus size
(no truncation).

## 6. Recommended next action toward "95 points"

The MMFF94 bond/angle/torsion/BCI classification-input bug class (issue
#227) is now closed for its 4 originally-identified sub-bugs (bond, angle,
torsion, BCI charges) — `pipeline_v2_mmff94_strict` coverage moved a total
of +1/265 across both Phase 1 and Phase 2 relative to the fresh State 1
baseline (240→241), smaller than the coverage-gate-tool's own +62/265
because `pipeline_v2_mmff94_strict` never gated torsion/stretch-bend
coverage to begin with (see Phase 1's CHANGELOG entry) — the real payoff of
both fixes is in geometry/charge *correctness* for molecules that already
passed, not in unlocking new passes. The next highest-leverage items,
un-touched by Phase 1 or 2:

1. **`DistanceGeometry::BoundsConstructionFailed`** (11/265, unchanged
   across all 3 states) is now `pipeline_v2_mmff94_strict`'s single largest
   failure-cause bucket (tied with `MissingParameters`) — an embedding-stage
   failure upstream of MMFF94 entirely, not touched by issue #227's scope.
   Root-causing this is likely the next-highest-value coverage lever.
2. The 67-atom / 11-molecule BCI residual (Section 1) — formal-charge/
   `fcadj` redistribution for charge-separated species — is small but real
   and independently root-causeable.
3. `mmff94_with_uff_fallback`'s and the `RepairAndVerify` arms' own 3-state
   numbers were not re-measured this round (Phase 2 scoped to
   `pipeline_v2_mmff94_strict` per the directive) — worth a follow-up pass
   if MMFF94-with-fallback coverage/quality is separately prioritized.
4. Stereo satisfaction rate (~56-57% across all 3 states, `Ignore` policy)
   remains RDKit's largest measured advantage on this corpus (RDKit: 100%)
   — `StereoPolicy::RepairAndVerify`'s own paired-arm numbers (not
   re-measured this round) are the existing lever for this gap.

## Files

- This report: `validation/results/mmff94_bci_gap_227_phase2_report.md`
- Machine-readable summary: `validation/results/mmff94_bci_gap_227_phase2_summary.json`
- Provenance / full BCI investigation writeup: `scripts/mmff94_provenance/PROVENANCE.md`
- Raw dumps: `mmff94_bci_gap_227_state{1,2,3}_*_chematic_rows.jsonl`,
  `_common_scored_paired.jsonl` (paired RMSD), `_tfd_paired.jsonl` (paired
  TFD), all in `validation/results/`.
- BCI charge investigation dumps: `mmff94_bci_charges_227_chematic_dump_{OLD_prefix,NEW_fixed}.jsonl`,
  `mmff94_bci_charges_227_rdkit_oracle.jsonl`,
  `mmff94_bci_charges_227_compare_{OLD,NEW}_summary.json`.
- Subset/transition JSONs: `mmff94_bci_gap_227_state{1,2}_torsion62_subset_report.json`,
  `mmff94_bci_gap_227_transition_state1_to_state2_torsion62subset.json`,
  `mmff94_bci_gap_227_transition_state2_to_state3_bciaffectedsubset.json`,
  `mmff94_bci_gap_227_transition_state1_to_state3_full.json`.
- Molecule ID lists: `mmff94_bci_gap_227_torsion62_subset_molecule_ids.json`
  (re-derived from the committed T0 JSONL dump
  `mmff94_coverage_227_term_audit_v0_16_0_prefix.jsonl`, matches Phase 1's
  reported 62 exactly), `mmff94_bci_gap_227_bci_affected_206_molecule_ids.json`.
- Tooling (new, all measurement-only, no production algorithm code):
  `crates/chematic-3d/examples/mmff94_bci_charges_dump_227.rs`,
  `scripts/mmff94_bci_charges_oracle_227.py`,
  `scripts/mmff94_bci_charges_compare_227.py`,
  `scripts/pipeline_v2_vs_rdkit_tfd_227.py`,
  `scripts/mmff94_bci_phase2_report_227.py`,
  `scripts/mmff94_bci_phase2_transition_227.py`; extension:
  `crates/chematic-3d/examples/pipeline_v2_vs_rdkit_common_scorer.rs` (new
  opt-in `--pair` mode, default behavior unchanged).
