# A1: Conformer Benchmark & Failure Ledger (diagnosis only — no implementation)

Status: draft, in progress. Scope: `chematic-3d` (conformer generation) and
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
- `scripts/gen_pipeline_v2_vs_rdkit_report.py` — aggregation, writing
  `validation/results/pipeline_v2_vs_rdkit_aggregate.json` and updating
  `docs/rfcs/pipeline_v2_vs_rdkit_etkdgv3_benchmark.md`.

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

## Pending (this section will be filled in once the corpus re-run + oracle + aggregation complete)

- [ ] Confirm `pipeline_v2_mmff94_strict` still measures 241/265 (or
      report the real current number) after v0.18.0's `assign_c_type`
      fix, not previously re-validated against this exact corpus.
- [ ] RMSD/TFD/stereo-satisfaction/runtime reported as explicitly
      separate axes (already the aggregation script's own convention —
      confirm the fresh run reproduces this shape, not a regression to
      one blended number).
- [ ] Macrocycle-corpus subset (all 11, with the 6 mis-categorized ones
      manually re-grouped for this ledger's own reporting, independent of
      the aggregation script's current blind spot) reported as its own
      row.
- [ ] Best-of-N parity gap restated with the fresh numbers (chematic has
      zero best-of-N arm; RDKit's `etkdgv3_best_of_n` arm is the only
      multi-restart comparison point in this benchmark today).
- [ ] Full failure-reason classification: does `PipelineV2FailureCause`'s
      existing enum already cover 100% of observed failures on this
      fresh run, or does the ledger surface an unclassified residual.

## Files touched by this round

- `docs/rfcs/a1_conformer_benchmark_failure_ledger.md` — this document.
- `validation/results/pipeline_v2_vs_rdkit_chematic_rows.jsonl`,
  `..._rdkit_rows.jsonl`, `..._aggregate.json` — regenerated (existing
  files, same schema, fresh data).
- `docs/rfcs/pipeline_v2_vs_rdkit_etkdgv3_benchmark.md` — regenerated by
  `scripts/gen_pipeline_v2_vs_rdkit_report.py` (existing file, same
  schema, fresh data).
- `validation/openeye_advantage_fixtures.jsonl`,
  `docs/competitive/openeye_gap_matrix.md`,
  `validation/openeye_advantage_holdout.jsonl` — corrected (commit
  `849ef60`, already pushed on this branch).

No changes under `crates/*/src/**`.
