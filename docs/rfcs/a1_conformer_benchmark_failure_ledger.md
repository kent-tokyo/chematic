# A1: Conformer Benchmark & Failure Ledger (diagnosis only — no implementation)

Status: draft, measurement round complete (see "Measured results," below;
the aggregation script itself could not complete — Finding 3). Scope:
`chematic-3d` (conformer generation) and
`chematic-ff` (force fields), per `docs/rfcs/openeye_materials_advantage_rfc.md`'s
Track A, item A1. **Diagnosis only — zero changes under `crates/*/src/**`.**
This is the RFC's own recommended next step after the OpenEye/materials
competitive audit (PR #369, merged).

## Goal

Re-audit the existing 265-molecule 3D corpus fresh against current `main`,
formalize success/RMSD/TFD/stereo-retention/runtime/failure-reason as
explicitly separate reported axes (not blended into one pass/fail number),
and re-verify — not assume — every prior finding this round's own research
cited, before any future A2 (Conformer Ensemble Core) design work builds on
top of them.

## Method

Reused the *existing* measurement infrastructure exactly as designed —
no new benchmark script was written:

- `cargo run --release -p chematic-3d --example pipeline_v2_vs_rdkit_dump`
  — chematic-side dump, 13 `pipeline_v2` arms + the legacy `etkdg` arm,
  against `validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_{a,b}.json`
  (65 + 200 = 265 molecules), fixed `EMBED_SEED = 20260801`,
  `MAX_ATTEMPTS = 8`.
- `scripts/pipeline_v2_vs_rdkit_oracle.py` — independent RDKit ETKDGv3
  oracle (4 arms incl. `etkdgv3_best_of_n`, N=10), reading the same
  manifests, never consuming chematic's own output.
- `cargo run --release -p chematic-3d --example
  pipeline_v2_vs_rdkit_common_scorer` — independent geometry+stereo
  scorer over both engines' already-saved coordinates; run once in its
  default flag-free form (refreshes the committed
  `..._common_scored_rows.jsonl`) and once with its existing, already
  opt-in `--pair` flag (uncommitted, ad hoc RMSD-vs-RDKit measurement —
  see "Measured results," below).
- `scripts/gen_pipeline_v2_vs_rdkit_report.py` — aggregation, intended to
  write `validation/results/pipeline_v2_vs_rdkit_aggregate.json` and
  update `docs/rfcs/pipeline_v2_vs_rdkit_etkdgv3_benchmark.md`. **Did not
  complete on this round's fresh data** — see Finding 3, below.

**Environment note**: the chematic-side dump run (`cargo run --release`)
took substantially longer than a quick sanity check — legitimately, not a
hang, since several arms retry up to 8 times with real MMFF94/UFF/DREIDING
minimization per attempt, and a subset of molecules sit near this
project's own documented `total_timeout_ms=20000` boundary
(`docs/rfcs/pipeline_v2_vs_rdkit_etkdgv3_benchmark.md`'s own prior notes
on `chembl_tier_b_0166` flip near that boundary). The run was interrupted
once by an out-of-process event unrelated to the benchmark itself and
restarted with `nohup`/`disown` for resilience.

## Findings confirmed so far (independent of the final aggregate numbers)

### 1. A previously-cited repro is stale — corrected, not just re-flagged

`docs/rfcs/openeye_materials_advantage_rfc.md` (PR #369, merged) cited
`docs/rfcs/etkdg_3d_gap_rfc.md`'s decane repro
(`worst_bond(m.conformer_ensemble(1, 0.0, 'dreiding', 0.0)[0])` ≈ 11.3 Å)
verbatim, without re-running it — exactly the kind of unverified citation
A1 exists to catch. Re-tested fresh (fresh `.venv` rebuild via
`maturin develop --release`, current `main`, commit `7a7dc814`):

- Decane, dreiding/mmff94 × noise 0.0/30.0, 20 attempts: **0/20 torn**
  (all landed ~1.5–1.6 Å, sound).
- Naphthalene (ring-containing, targets the ring-unaware-rotation
  mechanism), both force fields, noise 30.0, 10 attempts: **0/10 torn**.

Likely explanation: an unrelated ring-placement fix (PR #253, the
issue #277/#256 era) changed `dg::generate_coords`'s base geometry enough
that `etkdg.rs`'s under-iterated constraint-repair and ring-unaware
torsion rotation no longer trigger on *these specific examples* — an
incidental side effect, not a deliberate fix, and not established to
generalize to other molecules or ring shapes.

**A third, distinct mechanism is independently re-confirmed still live**:
MMFF94 silently contributes zero energy/gradient for atom-type pairs its
parameter tables don't cover, letting VdW repulsion push atoms apart
unopposed. On a halomethane stereocenter (`[C@H](F)(Cl)Br`) under
`conformer_ensemble(1, 0.0, 'mmff94', 0.0)`, central-carbon-to-halogen
distances land at 8.9–11.3 Å (real C–F/C–Cl/C–Br bonds: ~1.4–1.9 Å) —
verified twice independently (once via a dedicated bug-check subagent,
once via a direct per-atom distance check). Corrected in
`validation/openeye_advantage_fixtures.jsonl`'s `oe-01` row,
`docs/competitive/openeye_gap_matrix.md`'s Omega row, and
`openeye_advantage_holdout.jsonl`'s decane reference (commit `849ef60`).

**Takeaway for A2's design**: the public `conformer_ensemble()` API's
soundness defect is real and still live, but the demonstrating example
must be the halomethane/MMFF94 case, not decane — and the underlying
claim ("chematic has zero best-of-N arm, no seed control on this path")
is unaffected by this correction.

### 2. Six macrocyclic molecules are mis-categorized in the manifest, and the report generator can't see past that

Audited `validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_{a,b}.json`
directly (271 lines of Python, not assumed): of 265 molecules, **11 have
a ring ≥9 atoms** (`MACROCYCLE_MIN` in
`crates/chematic-3d/src/etkdg_knowledge/classify.rs`). Five are correctly
tagged `primary_category` as macrocycle-related (`cyclododecane`,
`crown_12_4`, `cyclooctadecane`, `cyclononane`, `macrolactam_12`) — but
**six ChEMBL molecules with 28–32-atom rings
(`chembl_tier_b_0009/0023/0028/0029/0030/0034`) are tagged
`primary_category: "drug_like"`**, not macrocycle. All six are the same
bis-pyridinium macrocycle family already tracked in `ROADMAP.md`'s
Backlog item 13 for an *unrelated* reason (issue #337's symmetrized-SSSR/
MMFF94-aromaticity residual) — a useful cross-reference, not a
coincidence worth re-investigating here.

Checked `scripts/gen_pipeline_v2_vs_rdkit_report.py` directly
(`grep -n "primary_category"`): the aggregation script groups results
**exclusively** by the manifest's own `primary_category` field, with no
independent ring-size-based re-classification. **This means every past
report this benchmark has ever produced has silently blended these six
molecules' macrocycle-class 3D-embedding performance into the general
"drug_like" average** — directly violating this RFC's own stated
acceptance criterion ("macrocycle accounted separately from normal
rings... never blended into the overall pass rate"). This is a real,
previously-undocumented methodology gap in an already-existing,
multiply-reused benchmark script — not a new defect in the underlying
`pipeline_v2` algorithm itself.

**Not fixed this round** (diagnosis only — re-categorizing manifest
entries or changing the report generator's grouping logic is itself a
methodology change that deserves its own small, explicit decision, not a
silent edit buried inside a "just re-running the numbers" round). Flagged
as Discovered Work for a tiny, clearly-scoped follow-up.

### 3. The aggregation script's own self-consistency check blocked report generation this round — and the trip is a third, independent symptom of Finding 2

`scripts/gen_pipeline_v2_vs_rdkit_report.py` has a hard `assert` (never
before tripped in this project's history, across the original Wave 1 run,
Priority 2, and Priority 2B) that fails the whole script if widening a
coverage gate turns a molecule's failure into a success with no
independently-verified explanation. On this round's fresh data it tripped:

```
AssertionError: chematic_pipeline_v2_mmff94_with_uff_fallback_complete_bonded_term_gated:
widening a coverage gate turned a failure into a success for
[{'name': 'chembl_tier_b_0030', 'earlier_status': 'timeout',
  'earlier_failure_cause': 'Timeout', 'earlier_elapsed_ms': 21946,
  'later_elapsed_ms': 18914,
  'why_unexplained': 'later-stage force_field_actual != UffOnly'}]
with no independently-verified timeout-rescue explanation
```

The check's recognized "legitimate explanation" shape is specifically a
UFF-fallback rescue (a stricter policy times out, a looser one falls back
to UFF and succeeds). That is not what happened here: on
`chembl_tier_b_0030`, the later (more-gated) arm succeeded via genuine
MMFF94 (`force_field_actual: Mmff94BondAngleStrict`, no fallback), not a
UFF rescue — an explanation shape the check does not know how to
recognize, so it correctly refuses to accept the flip silently.

Reading all 14 arms' rows for `chembl_tier_b_0030` directly: MMFF94
coverage is identical (0 missing bond/angle/torsion/oop terms) across
every gated variant, so the gates themselves are computational no-ops for
this molecule — yet its mmff94-family runtimes range from 6,460 ms (the
plain, least-gated `mmff94_strict` arm) up to 21,946 ms (the one arm that
timed out), with several gated/repair variants landing at 17.0–18.9 s in
between. Exactly one arm, one run, crossed the 20,000 ms
`total_timeout_ms` boundary. This spread across otherwise-equivalent
computations (same coverage, same converged bond length on every
successful arm: 1.5152690531578321 Å) is consistent with wall-clock
scheduling variance on a molecule that already sits close to the timeout
boundary across most of its mmff94-family arms — not a deterministic
function of which gate is enabled. This project has one prior, similarly
unexplained-but-disclosed boundary case (`chembl_tier_b_0166`, first seen
in Priority 2/2B), so a molecule in this corpus occasionally landing near
this exact wall-clock cutoff is not itself new; what is new is that this
round's flip happened to hit the specific comparison this assertion
checks.

**`chembl_tier_b_0030` is one of the six mis-tagged macrocycles from
Finding 2 above**, and it is measurably one of the slower-minimizing
molecules in the corpus across its mmff94-family arms. This is plausible
but not confirmed to be caused by its 28–32-atom ring specifically — this
round did not run a comparison against non-macrocyclic molecules of
similar heavy-atom count, so ring size is named as a correlation worth
noting, not asserted as the mechanism. Either way, this is a third,
independent way Finding 2's manifest mis-categorization has real
consequences: these six molecules are not just miscounted in the
"drug_like" bucket's average, they are disproportionately represented
among the corpus's timeout-boundary-adjacent cases, which a
category-blended report cannot surface.

**Per this round's own no-silent-scoring-changes discipline** (the same
reasoning as Finding 2's "not fixed this round"), the assertion was not
weakened, bypassed, or patched. Effect: `gen_pipeline_v2_vs_rdkit_report.py`
could not complete, so `docs/rfcs/pipeline_v2_vs_rdkit_etkdgv3_benchmark.md`
and `validation/results/pipeline_v2_vs_rdkit_aggregate.json` were **not**
regenerated this round — see "Files touched," below. The numbers in the
next section were instead computed directly from the two already-produced
JSONL files (`pipeline_v2_vs_rdkit_chematic_rows.jsonl`,
`..._rdkit_rows.jsonl`) plus one existing, unmodified Rust binary
(`pipeline_v2_vs_rdkit_common_scorer`, run both in its default flag-free
form and once with its existing opt-in `--pair` flag), via a throwaway,
uncommitted script — no scoring logic in any committed file was changed.

**Parked, not decided this round**: whether/how to extend
`_verify_timeout_rescue`'s recognized-explanation set to also cover
"succeeded within budget via the non-fallback path, on a molecule already
known to sit near the timeout boundary" is a real methodology decision
(it would need to define what independent evidence makes a flip
"explained" versus accepting these as permanently rare, disclosed,
report-blocking events). Not made this round; needs explicit
authorization before touching `scripts/gen_pipeline_v2_vs_rdkit_report.py`'s
scoring logic.

## Measured results (fresh 265-molecule re-run, 2026-08-23, commit range on this branch)

All numbers below are computed directly from the fresh
`pipeline_v2_vs_rdkit_chematic_rows.jsonl` / `..._rdkit_rows.jsonl` (both
265/265-molecule-complete, verified) plus one `--pair` invocation of the
existing `pipeline_v2_vs_rdkit_common_scorer` binary — not from
`gen_pipeline_v2_vs_rdkit_report.py`, which could not run (Finding 3).

- **`pipeline_v2_mmff94_strict`: 241/265 (90.9%) — CONFIRMED unchanged**
  from the previously-recorded v0.17.0 figure (commit `e2876bb`).
  v0.18.0's `assign_c_type` fix did not move this number on this corpus.
  Every one of the 241 successes is also `sound: true` (0 successes with
  unsound geometry on this arm).
- **Macrocycle subset reported separately, not blended** (the 11
  ring-≥9-atom molecules identified in Finding 2, all six mis-tagged ones
  included): **11/11 success, 11/11 sound** under `mmff94_strict` — this
  corpus's macrocycles do not show a lower `mmff94_strict` success rate
  than the rest of the corpus (230/254 = 90.6% for everything else). Note
  this cuts the other way from what Finding 2 might suggest at first
  glance: blending these six into "drug_like" very slightly *raises*
  `mmff94_strict`'s reported pass rate here (241/265 = 90.9% vs. 235/259 =
  90.7% with the six excluded), not lowers it — the real risk this
  mis-categorization hides is the timeout-boundary exposure in Finding 3,
  not a hidden pass-rate penalty, at least on this specific 265-molecule
  corpus.
- **Failure-cause classification, `mmff94_strict`**: all 24 typed
  failures resolve to exactly two `PipelineV2FailureCause` outer variants
  — 13 `ForceField(...)`, 11 `DistanceGeometry(...)` — summing to 24 with
  no unclassified residual, and 0 timeouts on this specific arm. The
  typed enum fully covers this arm's failure modes on this corpus.
  (Other, more heavily-gated arms do see timeouts — e.g. the 1 flip in
  Finding 3 — but a process-level wall-clock timeout is a separate
  escape valve, not itself a `PipelineV2FailureCause` variant.)
- **Stereo satisfaction, `mmff94_strict` successes**: 146 declared
  stereocenters across the 241 successful embeds, 82 satisfied (56.2%),
  64 violated (43.8%), 0 unevaluable.
- **Runtime, `mmff94_strict` successes**: median 1,301 ms, p90 5,570 ms,
  max 13,032 ms (n=241) — 0 molecules exceed 15,000 ms on this specific,
  least-gated arm (the near-timeout cases in Finding 3 only appear on
  more heavily-gated/repair-variant arms for the same molecule).
- **RMSD vs. RDKit, `chematic_pipeline_v2_mmff94_strict` vs.
  `rdkit_etkdgv3_mmff94`** (RDKit's own MMFF94 arm) — the **first time
  this specific 265-molecule benchmark has measured this number**;
  `docs/benchmark.md`/`docs/rdkit-migration.md`'s "no RMSD/TFD comparison
  exists" claim (corrected in PR #369) referred to a *different*,
  227-molecule BCI-residual corpus, not this one. Computed via the
  existing `pipeline_v2_vs_rdkit_common_scorer` binary's opt-in `--pair`
  flag (no committed-file or scoring-logic change): n=240 pairable
  molecules (1 of the 241 `mmff94_strict` successes has no RDKit-side
  MMFF94 success to pair against — not a join failure, RDKit itself
  didn't succeed on that molecule under its own MMFF94 arm), mean 1.678
  Å, median 1.503 Å, p90 3.456 Å, max 6.144 Å. Macrocycle subset (n=11):
  mean 1.610 Å, max 3.092 Å — comparable to, not worse than, the
  full-corpus figure. **TFD (torsion fingerprint deviation) is not
  computed anywhere in this project's codebase** (confirmed absent by
  direct grep) — it is not a gap introduced or left open by this round,
  it has simply never existed as a measured axis in chematic, on this or
  any other corpus.
- **Best-of-N parity, reconfirmed structurally**: chematic still has 0
  best-of-N arm across its 13 `pipeline_v2`/legacy arms; RDKit's
  `rdkit_etkdgv3_best_of_n` (N=10) remains the only multi-restart arm in
  this benchmark. Unchanged from the RFC's original claim — this round
  found no new information here, only confirmed the prior claim still
  holds structurally.

## Files touched by this round

- `docs/rfcs/a1_conformer_benchmark_failure_ledger.md` — this document.
- `validation/results/pipeline_v2_vs_rdkit_chematic_rows.jsonl` —
  regenerated fresh (`cargo run --release -p chematic-3d --example
  pipeline_v2_vs_rdkit_dump`, 265/265 molecules, 13 `pipeline_v2`/legacy
  arms + 1 diagnostic probe row, ~70 minutes; the first attempt was
  killed by something outside this session's control at 67% progress,
  restarted with `nohup`/`disown` for resilience and completed cleanly on
  the second attempt).
- `validation/results/pipeline_v2_vs_rdkit_rdkit_rows.jsonl` —
  regenerated fresh (`scripts/pipeline_v2_vs_rdkit_oracle.py`, 265/265
  molecules, 4 arms).
- `validation/results/pipeline_v2_vs_rdkit_common_scored_rows.jsonl` —
  regenerated fresh, default flag-free invocation (same schema as
  before, fresh data derived from the two files above).
- `validation/results/pipeline_v2_vs_rdkit_aggregate.json` and
  `docs/rfcs/pipeline_v2_vs_rdkit_etkdgv3_benchmark.md` — **NOT
  regenerated this round.** `scripts/gen_pipeline_v2_vs_rdkit_report.py`
  could not complete; see Finding 3. Both files remain at their prior
  (pre-A1) committed content in this PR — this ledger's "Measured
  results" section above is the source of truth for this round's fresh
  numbers, computed by a separate, uncommitted script, not by the blocked
  report generator.
- `validation/openeye_advantage_fixtures.jsonl`,
  `docs/competitive/openeye_gap_matrix.md`,
  `validation/openeye_advantage_holdout.jsonl` — corrected (commit
  `849ef60`, already pushed on this branch).

No changes under `crates/*/src/**` or `crates/*/examples/**` (both
`pipeline_v2_vs_rdkit_dump` and `pipeline_v2_vs_rdkit_common_scorer` were
run as-is, never edited).
