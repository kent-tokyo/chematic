# chematic Validation

Documented evidence that chematic's descriptors agree with industry-standard tools.

## Corpora

### Platinum coordination-chemistry compatibility benchmark (18-compound corpus)

Measurement-only survey of whether chematic can represent, parse,
round-trip, and canonicalize anticancer platinum(II)/platinum(IV)
coordination complexes (cisplatin, transplatin, carboplatin, oxaliplatin,
nedaplatin, lobaplatin, picoplatin, dicycloplatin, satraplatin, iproplatin,
tetraplatin, oxoplatin, plus charged/sulfur-donor/carbon-donor diversity
cases and 2 non-platinum generalization-gate cases) without silently
corrupting their coordination chemistry — **not** an anticancer-activity
project; no IC50/resistance/toxicity/pharmacokinetics prediction is
implemented or claimed anywhere in this benchmark. Not to be confused with
the unrelated third-party "Platinum Diverse Dataset" referenced in
`validation/manifests/dataset_provenance.json` (a drug-conformer dataset
that happens to share the word "platinum" in its name, nothing to do with
platinum coordination chemistry).

- **Files:** `validation/platinum/pt_corpus.jsonl` (corpus, with
  per-compound source/formula/charge/coordination-number/cis-trans
  provenance), `validation/results/platinum_baseline_chematic.jsonl`
  (unmodified-`main` baseline), `platinum_after_fix_chematic.jsonl`
  (post-fix), `platinum_attribution_{valence,mass}_only.jsonl` (per-fix
  isolated attribution runs), `platinum_rdkit_oracle.jsonl` (independent
  RDKit oracle, same corpus)
- **Reference tools:** RDKit 2026.03.3 (primary oracle), Open Babel 3.1.1
  (secondary, narrower comparison — see the full report for what was and
  wasn't checked)
- **How to regenerate:** `cargo run --release -p chematic-mol --example
  platinum_benchmark -- validation/results/platinum_after_fix_chematic.jsonl`
  + `python scripts/platinum_rdkit_oracle.py
  validation/results/platinum_rdkit_oracle.jsonl` + `python
  scripts/platinum_compare_chematic_rdkit.py`
- **Full report:** `validation/platinum/FEASIBILITY.md` — corpus
  provenance (including two independently-confirmed internally-
  inconsistent PubChem records), baseline, failure taxonomy, the
  cisplatin/transplatin killer benchmark, what was fixed vs. explicitly not
  fixed, before/after numbers with per-fix attribution, and the RDKit/Open
  Babel comparison.
- **Known findings from this benchmark, fixed:** dative-bond donor implicit
  hydrogen count (`chematic-core`), MOL V2000/V3000 bond-type-9 silent
  corruption (`chematic-mol`), periodic-table mass-data gap for every
  element outside a ~24-element hardcoded list (`chematic-chem`) — see
  `CHANGELOG.md`'s `[Unreleased]` section.
- **Known finding from this benchmark, measured but explicitly NOT fixed
  this round:** chematic has no square-planar stereo representation at all
  — cisplatin and transplatin currently canonicalize to the same identity.
  Reported as the top P0 gap for future work, not silently patched over.

### pipeline v2 vs RDKit ETKDGv3 independent 3D benchmark (Wave 1, 265-mol corpus)

Independent, measurement-only comparison of `chematic_3d::pipeline_v2::embed_pipeline_v2`
(9 force-field-policy/stereo-policy/gate-scope arms as of Priority 2's stretch-bend coverage
gate -- see `docs/rfcs/pipeline_v2_vs_rdkit_etkdgv3_benchmark.md` for the full current arm list)
+ the legacy `etkdg` entry point, against RDKit ETKDGv3
(raw/+UFF/+MMFF94/best-of-10), across a 65-molecule curated stress corpus (Tier A,
reusing the existing `pipeline_v2_integration_gate.rs` corpus) and a 200-molecule
deterministic ChEMBL subset (Tier B). Chematic side calls `embed_pipeline_v2` directly
from a Rust example (no Python binding overhead in the measurement); RDKit side is an
independent Python oracle reading the same corpus manifests -- neither feeds the other.
Atom mapping verified for 265/265 molecules (heavy-atom element-sequence match, not
assumed). No pipeline v2/force-field algorithm code was changed to produce these numbers;
historical numbers from earlier legacy-`etkdg`-only diagnoses are not reused as current.

Follow-up round: both engines' already-saved heavy-atom coordinates are additionally scored
by one common, independent scorer (`pipeline_v2_vs_rdkit_common_scorer.rs`) -- not chematic's
own internal `final_validation.sound`, which only ever existed for chematic's side. The same
scorer applies chematic's own `verify_stereo` judge to RDKit's geometry too, so stereo
preservation is compared with an identical judge rather than assumed from config flags.
Explicit per-arm denominators (`total`/`success`/`independently_sound`/`usable_coverage`)
replace the original "100% sound" framing, which implicitly meant "100% of successes," not
"100% of the corpus." Process-level (separate-process, 5-run) wall-clock timing supplements
the original in-process p50/p95/p99. A cyclopentane crash found in the original round was
scoped via ablation (`useSmallRingTorsions`/`enforceChirality`/5 seeds/3 stages) rather than
reported as a general RDKit failure.

- **Files:** `validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_{a,b}.json` (corpus
  manifests, with source/license/selection-rule/hash provenance), `validation/results/
  pipeline_v2_vs_rdkit_{chematic,rdkit}_rows.jsonl` (row counts grow with the arm matrix --
  see `pipeline_v2_vs_rdkit_environment_record.json` for the exact run this reflects;
  nothing silently dropped), `pipeline_v2_vs_rdkit_common_scored_rows.jsonl`
  (independent geometry + stereo judgment for both engines), `pipeline_v2_vs_rdkit_aggregate.json`,
  `pipeline_v2_vs_rdkit_process_level_perf.json`, `pipeline_v2_vs_rdkit_cyclopentane_ablation.jsonl`,
  `mmff94_coverage_227_term_audit{,_summary}.json(l)` (Priority 2 stretch-bend sub-classification),
  `*_config_snapshot.log`
- **Reference tool:** RDKit 2026.03.3 (`AllChem.ETKDGv3`)
- **How to regenerate:** `python scripts/gen_pipeline_v2_vs_rdkit_tier_a_manifest.py` +
  `python scripts/gen_pipeline_v2_vs_rdkit_tier_b_manifest.py` + `cargo run --release -p
  chematic-3d --example pipeline_v2_vs_rdkit_dump` + `python scripts/
  pipeline_v2_vs_rdkit_oracle.py` + `cargo run --release -p chematic-3d --example
  pipeline_v2_vs_rdkit_common_scorer` + `bash scripts/pipeline_v2_vs_rdkit_process_level_perf.sh`
  + `python scripts/pipeline_v2_vs_rdkit_cyclopentane_crash_ablation.py` + `python
  scripts/gen_pipeline_v2_vs_rdkit_report.py` (see each script's docstring for exact invocation)
- **Full report:** `docs/rfcs/pipeline_v2_vs_rdkit_etkdgv3_benchmark.md` (per-class/per-metric
  conclusions, auto-generated from the aggregate JSON -- no single overall win/loss claim)
- **Known issue filed from this benchmark:** MMFF94 parameter coverage gap
  ([#227](https://github.com/kent-tokyo/chematic/issues/227), separate from #185/#188)

#### v0.11.0 / pre-2B / v0.12.0 three-point paired comparison

Same harness (`pipeline_v2_vs_rdkit_dump.rs`), run unmodified at three historical
commits to see how chematic's own full-pipeline output changed release to release,
independent of RDKit: `ac52800` (v0.11.0 tag's library code, verified library-identical
to the `v0.11.0` tag via `git diff v0.11.0 ac52800 --stat -- crates/` -- only the
example script changed, adding RepairAndVerify arms), `c42627a` (PR #248+#249 merged,
immediately before Priority 2B's Dfsb production port), and `33eb2c3` (`v0.12.0` tag,
current release). Each commit is checked out wholesale (in an isolated `git worktree`
+ isolated `CARGO_TARGET_DIR` each, to avoid cross-contaminating build caches), so each
snapshot's arm matrix is whatever existed contemporaneously at that point -- newer arms
(the four stretch-bend/complete-bonded-term-gated arms, added in PR #249) are reported
only where they exist, not fabricated for earlier points. `scripts/pipeline_v2_vs_rdkit_oracle.py`
and both corpus manifests are confirmed byte-identical (via `git diff`) since the
v0.11.0 tag, so the RDKit oracle side was run once (RDKit 2026.03.3) rather than three
times -- this comparison is chematic-vs-itself across releases, not chematic-vs-RDKit
(see the main benchmark section above for that).

**v0.11.0 -> pre-2B**: 5/2121 common rows changed status, all `timeout<->success`/
`timeout<->typed_failure` flips on borderline-slow (~15-20s) molecules, consistent with
wall-clock timeout-boundary jitter already documented for this same molecule class in
Priority 2B's own measurement (`chembl_tier_b_0166` recurs here). No other movement --
matches expectation (PR #248 measurement-only, PR #249's gate defaults `false`
everywhere, `dg.rs` untouched in this range).

**pre-2B -> v0.12.0**: 168/3181 common rows changed status. The large majority are
`typed_failure -> success` on the four stretch-bend/complete-bonded-term-gated arms
(expected: Priority 2B's Dfsb production port resolved stretch-bend coverage from
2,107 missing instances to 0, see the Priority 2B changelog entry). A handful more
`timeout`-boundary flips on the same recurring borderline-slow molecule set as above
(`chembl_tier_b_0166`, `chembl_tier_b_0114`, `atorvastatin_fragment`).

**Finding, reported prominently rather than buried: `chematic_legacy_etkdg` regressed.**
This is the *only* arm in this benchmark that calls `dg::generate_coords` (confirmed:
`embed_pipeline_v2`/`distance_geometry_v2` never call it -- `grep generate_coords`
across `pipeline_v2.rs`/`distance_geometry_v2.rs` finds no hit -- so PR #253's
issue #185/#252 fix cannot and does not affect the production pipeline through this
benchmark). `sound` count: 265 -> 265 -> **248** across the three points -- 0 molecules
newly sound, 17 newly unsound, introduced entirely between pre-2B and v0.12.0 (i.e. by
PR #253). All 17 affected molecules already had high `bond_violation_rate_15pct`
(0.57-0.91) and nonzero `gross_clash_count` at *all three* snapshots -- these were
already marginal, messy ETKDG-refined geometries, not clean ones that broke outright.
PR #253 changed `dg::generate_coords`'s starting geometry (by design, to fix real
ring-placement bugs), which measurably shifted the downstream ETKDG local-refinement
outcome for these already-borderline molecules -- net negative on this specific metric
for this specific, non-production-path arm. Not investigated further here (would
require per-molecule geometry inspection, out of scope for a measurement-only PR) --
flagged as a candidate follow-up.

**MMFF term coverage** (missing-instance counts, summed across all mmff94-involving
arms present at each point): `stretch_bend_missing` 8,601 (pre-2B) -> **0** (v0.12.0),
matching the Priority 2B Dfsb port's own already-documented effect. `torsion_missing`
5,002 -> 5,371 (+7.4%) -- plausible but *not independently verified* explanation:
Dfsb resolving stretch-bend let more molecules proceed further into the pipeline (past
what used to be an early stretch-bend-driven failure) where their already-present
torsion gaps are now reached and counted for the first time, a counting-exposure
effect rather than a new torsion regression. Stated as a hypothesis, not a fact.

**Repair arms** (`mmff94_strict_repair`/`mmff94_with_uff_fallback_repair`):
`stereo_before_violations` constant at 64 across all three points as expected
(input-derived, policy-independent). `stereo_repaired_count` (52) and
`stereo_repair_failed_count` (12) are unchanged across all three points. The +/-1
variation in `final_stereo_violations` has not been causally attributed and is not
evidence of a measured repair-success-rate change.

- **Files:** `validation/results/pipeline_v2_vs_rdkit_v0_11_0_chematic_rows.jsonl`,
  `..._pre_2b_chematic_rows.jsonl`, `..._v0_12_0_chematic_rows.jsonl` (frozen per-point
  snapshots, distinct from the main benchmark's current-main-only rolling files above),
  `pipeline_v2_vs_rdkit_3point_paired_diff_summary.json` (full per-arm breakdown,
  status-change lists, the `chematic_legacy_etkdg` finding, MMFF coverage deltas)
- **Reproducibility:** v0.12.0 (candidate) side re-run twice. 3,173/3,181 rows
  (99.75%) byte-identical excluding `elapsed_ms`; the 8 differing rows all involve
  exactly the same 3 already-flagged recurring timeout-boundary molecules
  (`atorvastatin_fragment`, `chembl_tier_b_0114`, `chembl_tier_b_0166`) -- not claimed
  as fully deterministic, these 3 molecules' pass/fail is wall-clock-sensitive in this
  environment, consistent with every other timing-jitter observation in this report,
  not contradicting the byte-identical claim for the other 99.75%.
- **Known issue filed from this benchmark:** none yet -- the `chematic_legacy_etkdg`
  regression is reported here as measurement evidence; whether to file a tracking issue
  is a judgment call for the maintainer

#### v0.13.0 combined remeasurement (issue #227, PR #281 + PR #282 together)

Same harness (`pipeline_v2_vs_rdkit_dump.rs`), a fourth point added to the comparison
above. This is the **first time both MMFF94 production fixes are measured together**:
PR #281 (stretch-bend classification-key fix, `MMFF94_STBN` keyed by `stretch_bend_type`
not `angle_type`) and PR #282 (torsion classification fix, `torsion_type_for` keyed by
the j-k bond's real MMFF bond type, not atom-type membership). Both were individually
measured+merge-ready in their own PRs already; the repo owner's explicit condition for
proceeding to a v0.13.0 release is a clean, **one-time, non-iterated** combined
remeasurement confirming zero regressions. **Measured commit:
`2b608d3` (main's tip at RC freeze)** -- this also includes PR #268 (per-atom
stereocenter candidates) and PR #283 (atropisomer notation-invariance fix), neither of
which touches `chematic-ff`/MMFF94 code (`git diff c57cf58..2b608d3 --stat` -- 3 files,
all in `chematic-chem`/`chematic-perception`). Baseline: the same frozen
`pipeline_v2_vs_rdkit_v0_12_0_chematic_rows.jsonl` used by PR #281/#282 individually and
by the three-point comparison above.

**`mmff94_term_coverage_audit`, fresh on `2b608d3`:** StretchBend `routing_bug_candidate`
**180**/`table_gap` 1680, Torsion `routing_bug_candidate` **254**/`table_gap` 14 --
exact match to both PRs' individually-reported numbers. Zero interaction effect between
the two fixes on this corpus.

**Control arms (never call MMFF94 code -- `dreiding`/`no_ff`/`uff_only`/`legacy_etkdg`/
`ring_torsion_failclosed_probe`, 1,061 rows): byte-identical to the v0.12.0 baseline**
(excl. `elapsed_ms`), checked first as a harness/environment sanity gate before trusting
any MMFF94-arm number -- a moving control would have invalidated the whole comparison.

**Paired diff, all 3,181 common rows: zero `success -> failure`, zero `sound ->
unsound`, zero `non-timeout -> timeout`.** All 73 status changes are explained by two
already-documented mechanisms: **65** `typed_failure -> success` + the 2 corresponding
`timeout -> success` flips (chembl_tier_b_0166) in the `*_complete_bonded_term_gated`
arms (sound **85 -> 149** for `mmff94_strict_complete_bonded_term_gated`, **250 -> 252**
for the UFF-fallback variant -- a direct, deterministic consequence of the torsion
`routing_bug_candidate` drop, matching PR #282's own +64/+2 finding exactly), plus **6**
more `timeout -> success` and **2** `timeout -> typed_failure` flips, all on the 3
already-named recurring boundary-timeout molecules (`atorvastatin_fragment`,
`chembl_tier_b_0114`, `chembl_tier_b_0166`) behaving exactly as PR #281/#282 each
already documented for these same molecules -- zero new/unexplained molecules this run
(an earlier attempt at this same measurement, on a contended machine, *did* show several
unexplained novel timeout flips; isolating all six flagged molecules with an extended
120s `total_timeout_ms` found every one completing in 6.9-18.3s with `sound=true` and a
tiny, unchanged `distance_geometry_ms`, matching PR #281's own "CPU contention, not a
code effect" verdict -- that contended run was discarded, not reported as the measurement
of record, per the pre-committed run1/run2 policy below).

**RepairAndVerify** (`mmff94_strict_repair`/`mmff94_with_uff_fallback_repair`): repair
success rate, computed on the paired intersection of rows that reached a repair verdict
in both snapshots (not raw arm totals, which mix in bucket-membership effects from
molecules that reached repair in only one snapshot) -- **52 repaired / 12 failed / 0.8125
success rate, identical before and after, in both arms.** `final_stereo_violations` on
that same intersection: **10 -> 10** (strict), **12 -> 12** (UFF-fallback) -- zero
change. A naive arm-total reading shows +1 in each arm; that +1 is entirely
`chembl_tier_b_0166`'s already-covered `timeout -> typed_failure(FinalStereoViolation)`
bucket move (PR #282's own causally-verified finding), not a new violation on any row
comparable in both snapshots.

**E/Z + tetrahedral stereo, checked across all arms (paired on `final_stereo_violations`
non-null both sides), flagged prominently rather than buried: two molecules regress by
exactly 1 declared-stereocenter satisfaction each.** `chembl_tier_b_0126`
(`CC(=O)/C=C/CC1C(=O)N2[C@@H](C(=O)O)C(C)(C)S(=O)(=O)[C@@H]12`, violations 1 -> 2) and
`chembl_tier_b_0168` (the C12 epimer, violations 0 -> 1) -- both bicyclic beta-lactam
sulfones -- lose exactly one satisfied declared stereocenter, reproducibly, in all 4
un-gated/stretch-bend-gated MMFF94 arms (`mmff94_strict`,
`mmff94_strict_stretch_bend_gated`, `mmff94_with_uff_fallback`,
`mmff94_with_uff_fallback_stretch_bend_gated`); `sound` stays `True` in every case (not a
soundness regression). Reproduced identically across two independent commits measured in
this session (`c57cf58`, before PR #268/#283 landed, and `2b608d3`, after) -- since it
appears on both, it is caused by the stretch-bend+torsion classification fix itself
(present in both), not by the later, unrelated stereo/atropisomer commits. Arm-level
totals still net-improve (e.g. `mmff94_strict` corpus-wide `final_stereo_violations`:
49 -> 48) because more molecules gain a satisfied stereocenter than lose one -- but on a
strict per-molecule reading this does not cleanly clear a "zero regression" bar. Not
investigated further here (would need an RDKit-oracle re-check of which configuration is
actually correct for these two molecules, out of scope for this chematic-vs-itself
measurement) -- reported as a real, small (2/265 molecules), plausible consequence of the
corrected MMFF94 energy landscape converging to a different local minimum, not
hand-waved away.

**Coverage** (`stretch_bend_missing_count`/`torsion_missing_count`, paired
intersection): `stretch_bend_missing` stays **0** in every arm, both snapshots (Priority
2B already achieved full coverage). `torsion_missing`, `mmff94_strict`: **378 -> 0**;
`mmff94_with_uff_fallback`: **1,068 -> 248** -- consistent with the
`routing_bug_candidate` 1,107 -> 254 drop measured directly by the term-coverage audit
above.

**Reproducibility:** full corpus re-run twice on `2b608d3` (run1 is the measurement of
record for the paired diff above, decided before either run completed; run2 exists only
for this reproducibility check). **761/3,181 rows (23.9%) byte-identical in full
including `elapsed_ms`; 2,412/3,181 (75.8%) differ only in `elapsed_ms`-type fields
(coords/status/every other field bit-identical); 8/3,181 (0.25%) differ semantically.**
All 8 semantic differences involve exactly one molecule, `chembl_tier_b_0166`, across its
8 MMFF94-involving arms, flipping from run1's success/typed_failure (PR #282's own
isolation check put this molecule at ~19.1-19.3s, just under the 20s budget) to run2's
timeout -- zero other molecules differ. **8/3,181 = 0.25% unexplained-by-known-flake
rate matches the v0.12.0 three-point comparison's own 99.75% figure exactly.**

**Packaging / API smoke test:** `cargo build --release -p chematic-mcp -p chematic`
(the umbrella + MCP crates, re-exporting nearly everything, confirming the public API
surface compiles clean end to end after both fixes' breaking signature changes) --
PASS. `.venv/bin/maturin develop --release -m crates/chematic-py/Cargo.toml` +
`python3 -c "import chematic; chematic.from_smiles('c1ccccc1').mw"` -- PASS. WASM build
already covered by CI's own `test-wasm` job (not duplicated locally).

- **Files:** `validation/results/pipeline_v2_vs_rdkit_v0_13_0_chematic_rows.jsonl` (run1,
  measurement of record), `..._v0_13_0_run2_chematic_rows.jsonl` (run2, reproducibility
  only), `..._v0_13_0_paired_diff_summary.json` (full per-arm breakdown, all status
  changes, repair/stereo/coverage detail, reproducibility detail),
  `mmff94_coverage_227_term_audit_v0_13_0{,_summary}.json(l)` (fresh combined audit)
- **Not done here:** `docs/rfcs/pipeline_v2_vs_rdkit_etkdgv3_benchmark.md` was checked and
  deliberately NOT hand-edited -- it declares itself auto-generated from
  `pipeline_v2_vs_rdkit_aggregate.json`, which requires a fresh RDKit-oracle run +
  common-scorer run + process-level perf + cyclopentane ablation, none of which is in
  this measurement's scope (this is a chematic-vs-itself paired diff against a frozen
  baseline, not a chematic-vs-RDKit remeasurement); this section + `validation/results/
  *_v0_13_0_*` is this measurement's canonical home instead.

### MMFF94 strict-gate raw funnel remeasurement (issue #227, Priority 3 population)

A **low-level diagnostic harness**, not the production embedding entry point: calls
`minimize::minimize_with_policy(ForceFieldPolicy::Mmff94BondAngleStrict, ...)` directly
per molecule over the same 265-mol Tier A+B corpus, on top of `dg::generate_coords`
starting geometry. This is deliberately **not** `embed_pipeline_v2` (the production
entry point) and **not** issue #227's own posted reproduction path (which calls
`embed_pipeline_v2`, embedding via `distance_geometry_v2::embed_distance_geometry_v2_with_adjustments`
before minimizing -- a different, better starting geometry; see issue #252) -- and not
`mmff94_term_coverage_audit.rs`'s simplified bond/angle-only `Some`/`None` check either.
It exists to isolate the raw strict-minimization population directly on
`generate_coords` output, so that population's movement (e.g. across Priority 2B) can
be tracked in isolation from the full embedding pipeline. Fully deterministic
(`MinimizeConfig` has no RNG or wall-clock component, fixed `max_steps=200`) --
verified byte-identical across 2 back-to-back runs before being trusted for a
before/after diff.

Used to re-determine Priority 3's (Stage 1C, `MinimizationFailed` root-causing) target
population after Priority 2B's Dfsb stretch-bend production port (PR #250, merged as
`c92e075`). Finding: the pre/post-Priority-2B `MinimizationFailed` COUNT coincidentally
stayed at 28, but the molecule SET did not -- 14/28 (50%) churned (7 newly resolved to
`Ok`, 7 newly failing), confirmed deterministic via the byte-identical-rerun check on
both sides, not wall-clock jitter. `MissingParameters` (106) and `UnsupportedAtomType`
(1) sets are exactly unchanged, confirming Dfsb only perturbs energy/gradient among
molecules whose bond+angle parameters already fully resolve, and does not gate under
the legacy (non-gated) strict policy. Post-Priority-2B `MinimizationFailed` breaks down
19 `CatastrophicBondBlowup` / 9 `ExcessiveResidualForce` -- see summary file. Follow-up
diagnosis (issue #252) found this 28-molecule population is a `generate_coords`
starting-geometry artifact with no production impact -- a DIFFERENT, smaller population
than the 11-15 `typed_failure` molecules issue #227 itself flagged from the
`embed_pipeline_v2` funnel; the two are not the same measurement and are not known to
overlap.

- **Files:** `validation/results/mmff94_strict_gate_remeasure_227_rows.jsonl` (current
  main, one row per molecule, `MinimizationFailed` rows carry the full
  `MinimizationFailureDetail` -- `reason`/`converged`/`iterations`/
  `max_residual_force`/`worst_bond_length`/`distance_geometry_v2_retry_attempted`),
  `mmff94_strict_gate_remeasure_227_pre_priority2b_baseline.jsonl` (frozen snapshot at
  commit `c42627a`, NOT reproducible by running today's script against today's main --
  see its provenance note in the summary file), `mmff94_strict_gate_remeasure_227_summary.json`
  (churn analysis, set-identity checks, reason breakdown)
- **How to regenerate (current-main side only; the pre-Priority-2B baseline is frozen):**
  `cargo run --release -p chematic-3d --example mmff94_strict_gate_remeasure_227 >
  validation/results/mmff94_strict_gate_remeasure_227_rows.jsonl 2>
  validation/results/mmff94_strict_gate_remeasure_227_stderr.log`

#### Post-#253 re-measurement (issue #185/#252)

Same harness, re-run after PR #253 (merged as `c370cb3`) fixed three
`dg::generate_coords` atom-coincidence bugs (root/ring-vertex collision,
ring-fusion-order mismatch, new-island direct-bond anchor). Baseline is
`ada800d` (main's tip immediately before #253), frozen as
`mmff94_strict_gate_remeasure_227_pre_253_baseline.jsonl` (distinct file
from the pre-Priority-2B baseline above -- do not conflate the two).
`mmff94_strict_gate_remeasure_227_rows.jsonl` now reflects the post-#253
(candidate) side; `mmff94_strict_gate_remeasure_227_summary.json` was
regenerated to compare baseline vs. candidate instead of
pre/post-Priority-2B.

Finding: all 28 `MinimizationFailed` molecules (19 `CatastrophicBondBlowup`
+ 9 `ExcessiveResidualForce`) resolved to `Ok`, 0 remaining, 0 new
failures. `MissingParameters`/`UnsupportedAtomType` sets unchanged.
Candidate side re-run twice, byte-identical, before trusting this
comparison. None of the 28 resolved molecules contain a fused-polycyclic
aromatic ring system (checked by SMILES inspection) -- the anthracene-class
fused-ring-seam limitation found during PR #253's review (see issue
tracker) is not implicated in this population. Does not establish which
individual PR #253 sub-fix was causal per molecule (no ablation study
run) -- only that the tracked population is fully resolved with no
regression.

### Symmetric RMSD oracle check (Priority 4 groundwork, `rmsd_symmetric`)

`chematic_3d::conformer::rmsd_symmetric` (`crates/chematic-3d/src/conformer.rs`)
is a port of RDKit's `GetBestRMS` (`Code/GraphMol/MolAlign/AlignMolecules.cpp`):
symmetry-aware Kabsch RMSD via brute-force enumeration of the molecule's own
graph automorphisms (VF2 self-match, `uniquify: false`), not a fixed
atom-index correspondence. Needed before Priority 4's reference-conformer
benchmark, since a plain fixed-index RMSD is wrong on any molecule with
permutation-equivalent atoms (e.g. a terminal `-CF3`'s three fluorines).

- **Dump:** `cargo run --release -p chematic-3d --example
  rmsd_symmetric_oracle_dump > /tmp/rmsd_symmetric_oracle_dump.jsonl` -- 6
  cases (propane control; `-CF3`, neopentane, and benzene testing
  automorphism recovery on leaf atoms vs. ring atoms; acetate testing the
  known gap below; ibuprofen as a no-symmetry drug-like control), each a
  `dg::generate_coords` conformer paired with a fixed rigid
  rotation+translation of itself with two atoms' positions additionally
  swapped.
- **Oracle:** `.venv/bin/python scripts/rmsd_symmetric_oracle_check.py
  /tmp/rmsd_symmetric_oracle_dump.jsonl` -- independently recomputes each
  pair's RMSD via RDKit's `rdMolAlign.GetBestRMS` on the same coordinates.
- **Result:** 5/6 cases agree with RDKit to within 2e-6 Å (propane,
  `-CF3`-ethane, neopentane, benzene, ibuprofen). Benzene's nonzero residual
  (1.143095 Å, both engines) is correct, not noise: swapping two *para ring*
  atoms' positions is not a graph automorphism of a plain hexagon (unlike
  swapping two terminal leaf atoms on the same parent, as in `-CF3` or
  neopentane's methyls), so no relabelling can recover 0 -- confirmed by
  independent derivation, not just tool agreement.
- **Known, documented gap:** acetate disagrees (chematic 0.356 Å vs. RDKit
  0.0 Å). RDKit's `GetBestRMS` additionally runs
  `symmetrizeConjugatedTerminalGroups` before matching, treating a
  carboxylate's two oxygens as interchangeable despite their different
  formal bond orders; `rmsd_symmetric` does not port that preprocessing step
  (documented in its own doc comment). Expected, not a bug -- tracked as a
  follow-up, not blocking Priority 4's benchmark groundwork.

### 175-mol drug-like corpus

A curated set of 175 drug-like molecules covering common scaffolds (benzoic acid derivatives,
heterocycles, amino acids, steroids, macrolides). Used for per-descriptor regression testing.

- **File:** `scripts/rdkit_ref_properties.tsv` (175 rows)
- **Columns:** name, smiles, mw, logp, tpsa, hac, hbd, hba
- **Reference tool:** RDKit 2026.03.3
- **How to regenerate:** `python scripts/gen_rdkit_reference.py`

### 4,999-mol ChEMBL subset

A random sample from ChEMBL used for large-scale agreement testing on HBA, HBD, and aromatic ring count.

- **File:** external (not committed; requires download)
- **Reproduce:** `python scripts/bench5k.py ~/Downloads/SMILES.csv`

### Morgan/ECFP RDKit environment-parity diagnostic (5,000-mol corpus + 41 fixtures)

Locates, per molecule, the first stage at which chematic's Morgan/ECFP expansion diverges
from RDKit's (radius-0/1/2 invariants, redundant-environment suppression, sparse counts,
2048-bit folding, bitInfo). Diagnostic only -- production `ecfp4()`/`ecfp6()`/
`morgan_fp_counts()` are unchanged.

- **Files:** `ecfp_rdkit_environment_parity_manifest.json`, `_summary.json`,
  `_rows.jsonl` (41 edge-fixture molecules), `_first_divergence.tsv` (full 5,041-input run)
- **Reference tool:** RDKit 2026.03.3 (`rdFingerprintGenerator.GetMorganGenerator`)
- **How to regenerate:** `python scripts/gen_ecfp_rdkit_environment_oracle.py` +
  `cargo run -p chematic-fp --release --features diagnostics --example
  morgan_rdkit_environment_trace` + `python scripts/ecfp_rdkit_environment_parity.py`
  (see each script's docstring for exact invocation)

### Morgan/ECFP RDKit environment-suppression parity (Phase B, 5,041+4-input set)

Whether chematic's RDKit-equivalent redundant-environment suppression
(`crates/chematic-fp/src/morgan_environment.rs`, `SuppressRdkitRedundant` mode
-- additive/experimental; production `ecfp4()`/`ecfp6()`/`morgan_fp_counts()`
unchanged) emits the same set of `(atom_idx, radius)` environments, and the
same raw-identifier sparse-count *shape*, as RDKit's own default
(`includeRedundantEnvironments=False`) generator, on PR #120's original
5,041-input set plus 4 pinned representative-swap fixtures (5,045 total; see
`scripts/ecfp_rdkit_suppression_representative_swap_fixtures.csv`).

Implementation verified directly against RDKit's real C++ source: commit
[`0062b670640352ab63d6256be608615e87e1af53`](https://github.com/rdkit/rdkit/blob/0062b670640352ab63d6256be608615e87e1af53/Code/GraphMol/Fingerprints/MorganGenerator.cpp),
`MorganEnvGenerator<OutputType>::getEnvironments` -- a specific commit SHA,
not a mutable `master` reference.

**Results:**

| Metric | Result |
|---|---|
| Emitted `(atom_idx, radius)`-pair-set exact match | 5,032/5,045 (99.74%) |
| Raw-identifier sparse-count *shape* exact match (multiset of per-id emission counts) | 5,044/5,045 (99.98%) |
| `sparse_count_mismatch` fixtures (8, from the Phase A diagnostic) now shape-resolved | 8/8 |
| Tanimoto-vs-RDKit Pearson r, before (`ecfp4_rdkit_invariants`) → after (`ecfp4_rdkit_environment_experimental`) | 0.9479 → 0.9547 (Δ+0.0068, improved; n=300 sample, seed=42, 44,850 pairs -- non-gating reference) |
| Full-corpus (5,000 mol) wall time, baseline → suppression (median of 5 independent process runs) | 1.315s → 1.508s (1.146x) |
| Full-corpus peak RSS, baseline → suppression (median of 5 runs, `/usr/bin/time -l`) | 18.7 MiB → 19.3 MiB (1.030x) |

Pair-set mismatches: 13 of 5,045 (the same 9 residuals from the original
5,041-input run, plus their 4 pinned duplicates), all single-pair swaps at
the same radius -- **not** a claim that the swapped atoms are chemically
equivalent or near-equivalent, and **not** a claim that the two candidates
provably compute the identical cumulative bond environment (that would
require diagnosing raw bond-index-sets directly, which this validation
doesn't do). What's actually measured, precisely: two different atoms
produce the same *raw identifier*, and the selected representative differs
because RDKit and chematic currently order those candidates using different
hash values (FNV-1a vs RDKit's own hash never match by construction -- same
"not bit-compatible, partition/set-only" scope as every other RDKit-parity
mode in this crate). See the pinned fixtures for concrete cases: `CC(=O)NO`
(atoms 1 vs 3, not a symmetric pair), an isotope-labeled methyl pair, a
steroid-like fused-ring epoxide, and a large polycyclic aromatic -- each
verified to be *exactly* a 1-pair swap with total-emitted-count,
sparse-count shape, and unique-*raw-identifier*-count (deliberately not
called "unique bond-environment count" -- a raw identifier can in principle
be shared by two structurally different environments via hash collision, as
the pyridine case below demonstrates) all preserved; only which atom
represents one shared identifier differs. **These 4 fixtures, plus the 8
`sparse_count_mismatch` fixtures, plus every "both"-bucket mismatch anywhere
in the input, are hard GATES in `scripts/ecfp_rdkit_suppression_parity.py`
(nonzero exit on any regression) -- not just reported numbers.**

Sparse-count-shape mismatch: 1 of 5,045, a pair-set *exact match*
(`C1=CC=NC=C1`, Kekulé pyridine) whose count multiplicities still differ --
traced to accidental cross-radius hash collisions that differ between
FNV-1a and RDKit's hash for this molecule's structurally-symmetric ring
carbons, not a suppression-algorithm defect (the underlying emission
*decision*, i.e. which atoms survive at which radii, is provably identical
between the two implementations for this molecule).

- **File:** `ecfp_rdkit_suppression_parity_summary.json`,
  `ecfp_rdkit_suppression_tanimoto_summary.json`
- **Reference tool:** RDKit 2026.03.3 (same oracle rows as the Phase A
  diagnostic above -- `default` variant's `sparse_bit_info`/`sparse_counts`/`folded_on_bits`)
- **How to regenerate:** `cargo run -p chematic-fp --release --features
  diagnostics --example morgan_suppression_dump` +
  `cargo run -p chematic-fp --release --example morgan_suppression_tanimoto_dump`
  + `python scripts/gen_ecfp_rdkit_environment_oracle.py` +
  `python scripts/ecfp_rdkit_suppression_parity.py` +
  `python scripts/ecfp_rdkit_suppression_tanimoto.py` (see each script's
  docstring for exact invocation). Performance record (not a merge gate):
  `cargo run -p chematic-fp --release --example morgan_suppression_benchmark`.

### Morgan M4-A0: RDKit-exact raw-identifier hash port (diagnostic, 5,048-input set)

Diagnostic-only, source-verified port of RDKit's actual Morgan hash-combine
machinery (`crates/chematic-fp/src/rdkit_morgan_hash.rs`) -- unlike every
prior Morgan-parity mode in this crate, which only claims *partition*
agreement (which atoms are chemically equivalent) and *lifecycle* agreement
(who wins/dies under suppression), this compares actual 32-bit hash VALUES
against real RDKit, atom by atom and radius by radius. Not wired into any
production API at the time of this milestone (later promoted to production
by the "Phase B" section below); `ecfp4()`/`ecfp6()`/`ecfp4_rdkit_invariants()`/
`ecfp4_rdkit_environment_experimental()` are all unchanged (production
snapshot verified byte-identical, see Results).

Ported directly from RDKit's real C++ source, commit
[`8afba32ec539dcb2369bc84549d802aca3f7eb39`](https://github.com/rdkit/rdkit/blob/8afba32ec539dcb2369bc84549d802aca3f7eb39/Code/GraphMol/Fingerprints/MorganGenerator.cpp)
(the true resolution of tag `Release_2026_03_4`, independently verified via
the GitHub tags API this session -- see `THIRD_PARTY_NOTICES.md` for the
attribution and the note on two other, imprecise SHAs already in this
project's history under the same tag label, since fixed in place).

**Two aromaticity-preprocessing paths were compared, and their results are
NOT pooled into one number** -- an earlier draft of this section did pool
them (a Hueckel fallback silently substituted for 2 rows where the
RDKit-parity engine failed, then counted as an RDKit-parity "match"),
which is exactly the measurement accident
[`apply_aromaticity_rdkit_parity_experimental`]'s own `Result` contract
exists to prevent. Corrected: every row is tagged with an explicit
`aromaticity_status`, and the RDKit-parity exact-match rate is computed
ONLY over rows where that engine actually succeeded.

**Results (5,048-input set = 5,000-mol corpus + PR #120's 41 fixtures + PR
#123's 4 fixtures + `ecfp_rdkit_m4a0_hash_fixtures.csv`'s 3):**

**Update (2026-07-22, `fix/kekulize-charge-aware-k1` #141 + `fix/core-te-normal-valences`
#142, merge commit `97e1fd50a19703b0f9c8829df6bd219e85c8fd72`):** the table below is
regenerated against current `main`, not the original 2026-07-20 numbers. Two
chematic-core fixes changed `atom_must_be_matched`'s (`kekulization.rs`) charge-blind
lone-pair-donor rules and added a missing `Te` valence-table entry
(`element.rs::normal_valences`); one of the two RDKit-parity preprocessing failures below
-- pyridinium's `c1cc[nH+]cc1` -- now succeeds and is bit-exact, dropping the error count
from 2/5,048 to 1/5,048. Only the bridgehead-N purine-like ring remains. Verified via the
same scripts, same 5,048-input corpus, same options -- not re-derived from the old
numbers.

| Path | Metric | Result |
|---|---|---|
| Production Hueckel aromaticity | Radius-0 numeric exact match | 5,048/5,048 (100%) |
| Production Hueckel aromaticity | Full numeric exact match (radius 0-2, representative selection, sparse counts, folded bits, bitInfo) | 4,989/5,048 (98.83%) |
| RDKit-parity aromaticity (`apply_aromaticity_rdkit_parity_experimental`, no fallback) | Preprocessing succeeded | 5,047/5,048 |
| RDKit-parity aromaticity | Preprocessing failed (`KekulizationFailed`, pinned as a fixture -- see below) | 1/5,048 |
| RDKit-parity aromaticity | Full numeric exact match **among the 5,047 successful rows** | **5,047/5,047 (100%)** |
| RDKit-parity aromaticity | Non-exact among successful rows | 0 |
| Hueckel control on JUST the 1 error row (non-gating -- answers "does the OLD path agree with RDKit here", not "does RDKit-parity work here") | Exact match | 1/1 |
| PR #123's 9 unique representative-selection residuals | Resolve to `exact_match` under the RDKit-exact hash alone (no aromaticity-engine swap needed) | 9/9 |
| PR #123's Kekule-pyridine sparse-count-shape mismatch (`C1=CC=NC=C1`) | Resolves to `exact_match`; confirms the documented root cause (FNV-1a-specific hash collision, not a suppression defect) | resolved |
| Production API byte-identical (`ecfp_regression_snapshot`, before/after SHA-256) | confirmed | 0 change |
| Oracle regeneration determinism (`--verify-determinism`, full 5,048 input) | byte-identical across two runs | confirmed |
| Positive controls (radius-0/1 identifier, bond invariant, 32-bit-wrapping removal, representative swap, folded-bit, dropped row, duplicate row ID) | all correctly cause non-zero exit; reverted, never committed | 8/8 |

**Cross-referencing the 59 Hueckel-path residuals against the RDKit-parity
path, row by row (not just comparing aggregate counts):** all 59 had
RDKit-parity preprocessing succeed, and all 59 became `exact_match` under
it -- `resolved_by_rdkit_parity: 59, not_evaluable_due_to_aromaticity_error:
0, still_mismatching: 0`. The 1 remaining RDKit-parity error row does not
overlap with the 59 Hueckel residuals (it was already an exact match
under Hueckel).

The 59-row residual under production Hueckel aromaticity traces to ONE
mechanism: chematic's Hueckel-based aromaticity *perception* disagreeing
with RDKit's own aromaticity model on specific fused/macrocyclic ring
systems (e.g. `C[Si](C)(C)c1ccc(C2=Cc3ccccc3C3=NCCCN23)cc1`) -- not a hash
defect. `apply_aromaticity_rdkit_parity_experimental`
(`crates/chematic-perception/src/rdkit_parity.rs`, built for exactly this
kind of disagreement, in an earlier milestone) resolves it.

**The 1 remaining RDKit-parity preprocessing failure (the bridgehead-N
purine-like ring) is pinned as a permanent fixture**, not just recorded in
a JSON summary -- `scripts/ecfp_rdkit_m4a0_rdkit_parity_kekulization_gap_fixtures.csv`
(now 1 line; pyridinium's line was removed post-K1, not duplicated -- it
remains in the 5,048-input corpus itself, at its original position in
`ecfp_rdkit_edge_fixtures.csv`, now correctly classified as a success/
exact-match row rather than a gap fixture), plus
`chematic-perception::rdkit_parity::tests::production_api_reports_kekulize_failure_not_panic`
/ `production_api_does_not_mutate_input_on_failure` and
`chematic-fp::rdkit_morgan_ecfp4::tests::kekule_bridgehead_n_purine_reports_kekulization_failed_not_a_fallback_result`
/ `hueckel_fallback_would_be_detectable_if_silently_reintroduced` (both renamed
from their pre-K1, pyridinium-specific names -- same general contract, not a
pyridinium-specific fact) -- so a future kekulization fix is verified
against the actual engine returning `Ok`, not just a number changing. Not a
corpus/fixture duplicate (the purine SMILES still appears exactly once
across the whole 5,048-input set).

A real bug in this diagnostic's own trace logic was found and fixed during
this milestone (not a hash defect either): an early version shared one
`dead`-atom array between RDKit's two `includeRedundantEnvironments`
lifecycles, so an atom suppressed under the `default` lifecycle silently
stopped being *computed* under the `full` lifecycle too in later rounds --
caught by the full-corpus comparator (hand verification on a handful of
fixtures missed it, since it only checked value equality on entries present
on both sides, not entry *count*). Fixed by running two fully independent
passes and merging by `(atom, radius)` key; regression-pinned in
`rdkit_morgan_hash.rs`'s own test suite.

- **Files:** `ecfp_rdkit_raw_identifier_parity_summary.json` (production
  Hueckel aromaticity run), `ecfp_rdkit_raw_identifier_parity_aromaticity_variant_summary.json`
  (RDKit-parity aromaticity engine, full corpus, honest success/error
  denominators), `ecfp_rdkit_raw_identifier_parity_oracle_manifest.json`
- **Reference tool:** RDKit 2026.03.3 (`rdFingerprintGenerator.GetMorganGenerator`,
  `includeRedundantEnvironments` True and False variants; same pinned option
  set as every other Morgan-parity mode in this crate)
- **How to regenerate:** `python scripts/gen_ecfp_rdkit_environment_oracle.py`
  + `cargo run -p chematic-fp --release --features diagnostics --example
  rdkit_morgan_hash_dump` + `python scripts/ecfp_rdkit_raw_identifier_parity.py`
  (see each script's docstring for exact invocation). RDKit-parity-engine
  comparison (no Hueckel fallback): `cargo run -p chematic-fp --release
  --features diagnostics --example rdkit_morgan_hash_dump_aromaticity_variant`
  + `python scripts/ecfp_rdkit_raw_identifier_parity_aromaticity_variant.py`.

**Implemented as Phase B, same day (2026-07-20)** -- see the section
immediately below for the production API and its own full corpus results.

### Phase B: `rdkit_morgan_ecfp4_experimental` -- production, fallible, RDKit-bit-exact ECFP4

Promotes M4-A0's reference engine (`crates/chematic-fp/src/rdkit_morgan_hash.rs`)
to a real public API in `crates/chematic-fp/src/rdkit_morgan_ecfp4.rs`:
`pub fn rdkit_morgan_ecfp4_experimental(mol: &Molecule) -> Result<RdkitMorganEcfp4, RdkitMorganError>`.
Scope is intentionally narrow, matching exactly what M4-A0 verified numerically:

- **radius = 2 (ECFP4) only** -- not ECFP6/radius = 3; M4-A0 never compared radius 3
  against the oracle, so claiming bit-exactness there would be unverified.
- Uses `apply_aromaticity_rdkit_parity_experimental` internally as a fallible `Result`
  step -- **no Hueckel fallback anywhere in the public path.** No entry point accepts a
  pre-aromatized `Molecule` (would let a caller bypass the engine and silently lose the
  bit-exactness guarantee).
- `RdkitMorganError`: `Aromaticity(AromaticityError)` (kekulization/internal-invariant
  failure, wrapping `rdkit_parity.rs`'s own error), `UnsupportedBondOrder { bond_idx, order }`
  (a `BondOrder` with no real RDKit `Bond::BondType` counterpart -- only chematic's
  SMARTS-query-only variants, which cannot occur for SMILES-parsed input; confirmed via a
  programmatically-built `Molecule` test since it can't be reached via `parse()`),
  `InternalInvariantViolation { reason }`.
- One shared computation, not independent per-field loops: a single pass over RDKit's
  `includeRedundantEnvironments=false` ("default") lifecycle populates all four
  `RdkitMorganEcfp4` fields (`fingerprint`, `sparse_counts`, `raw_bit_info`,
  `folded_bit_info`) at once.

**Results (same 5,048-input M4-A0 corpus, fresh dump + comparison against the same RDKit
oracle rows, not re-derived from M4-A0's own numbers):**

**Update (2026-07-22, same `fix/kekulize-charge-aware-k1` #141 + `fix/core-te-normal-valences`
#142 as the M4-A0 update above):** regenerated against current `main`; pyridinium moves
from the error count into the success count, now bit-exact.

| Metric | Result |
|---|---|
| Preprocessing succeeded (`status: "success"`) | 5,047/5,048 |
| Preprocessing failed (`rdkit_parity_kekulization_failed`, the same pinned fixture as M4-A0) | 1/5,048 |
| Full exact match (default-lifecycle raw pairs, sparse counts, folded on-bits, folded bitInfo) among the 5,047 successful rows | **5,047/5,047 (100%)** |
| Hermetic equivalence to `rdkit_morgan_raw_trace`'s already-oracle-validated `raw_identifier_default` output, same already-aromatized molecule | confirmed (unit test, 4 representative fixtures) |
| Non-regression: `ecfp_regression_snapshot` (10 existing entry points: `ecfp4`, `ecfp6`, `ecfp` chiral, `ecfp_with_bitinfo`, `morgan_fp_counts`, `ecfp4_rdkit_invariants`, `ecfp6_rdkit_invariants`, `ecfp4_rdkit_environment_experimental`, `ecfp6_rdkit_environment_experimental`, `ecfp_with_bitinfo_rdkit_environment_experimental`), full 5,048-input corpus, SHA-256 before/after (git-worktree baseline at the pre-Phase-B commit) | byte-identical, 0 change |
| Unsupported-bond-order path (`BondOrder::QueryAny`, programmatically built -- cannot arise from `parse()`) | explicit `Err(UnsupportedBondOrder)`, confirmed by test |
| Positive control: a silently reintroduced Hueckel fallback would be numerically detectable (Hueckel perceives the bridgehead-N purine-like ring's real Hückel partition where the real path correctly errors -- fixture swapped from pyridinium's `c1cc[nH+]cc1` post-K1, since that molecule no longer fails) | confirmed by test |

**Performance vs. `ecfp4_rdkit_environment_experimental` baseline** (5 independent process
runs each, full 5,048-corpus, median wall time, `/usr/bin/time -l` for peak RSS -- not a
Criterion-registered benchmark, see `feedback_criterion_gate_pseudo_replication`):

| | Baseline | Candidate | Ratio |
|---|---|---|---|
| Median wall time (5 runs) | 4.862s | 9.734s | **2.00x** |
| Peak RSS | ~20.2 MB | ~20.4 MB | 1.01x |

The ~2x ratio is fully attributable, not an unexplained regression: the baseline reads
whatever aromatic flags are already on the input `Molecule` and never calls an aromaticity
engine, while the candidate performs its own kekulization + RDKit-parity aromaticity
perception on every call (a per-molecule breakdown shows preprocessing is 46-56% of the
candidate's time on aromatic-ring-heavy molecules like benzene/aspirin/a steroid-like
fused system, dropping to 19-20% on large acyclic alkanes with no rings to perceive).
Per the acceptance-gate policy (stop and explain, don't silently tune), this is reported
as measured rather than optimized against.

Any `BondOrder` this engine cannot map to a real RDKit `BondType` (verified against
`Bond.h` during M4-A0: SINGLE=1, DOUBLE=2, TRIPLE=3, QUADRUPLE=4, AROMATIC=12, DATIVE=17,
ZERO=21; only chematic's SMARTS-only `Query*` variants have no RDKit equivalent) is an
explicit `Err`, never an implicit/guessed mapping.

- **Files:** `crates/chematic-fp/src/rdkit_morgan_ecfp4.rs`,
  `crates/chematic-fp/examples/rdkit_morgan_ecfp4_dump.rs`,
  `crates/chematic-fp/examples/rdkit_morgan_ecfp4_benchmark.rs`,
  `scripts/ecfp_rdkit_morgan_ecfp4_parity.py`.
- **How to regenerate:** `cargo run -p chematic-fp --release --example rdkit_morgan_ecfp4_dump
  -- <SMILES.csv> <out.jsonl>` + `python scripts/ecfp_rdkit_morgan_ecfp4_parity.py --chematic
  <out.jsonl> --rdkit-oracle <gen_ecfp_rdkit_environment_oracle.py output>`. Self-test:
  `python scripts/ecfp_rdkit_morgan_ecfp4_parity.py --self-test`.

### Descriptor census — `crates/chematic-chem/src/descriptors.rs` (5,000-mol ChEMBL corpus)

Full census of all ~196 individually-named values across the 71 `pub fn` in `descriptors.rs`
specifically (not the whole `chematic-chem` crate -- see the RFC for the sibling-file scope
boundary). Diagnostic only, no production code changed. Root-causes 3 concrete defects
(`num_unspecified_stereocenters` massive over-counting, `num_hydrogens` double-counting on
bracket-H stereocenters, `molecular_weight` ignoring isotope labels) down to the exact line of
source with minimal reproducers, plus a real VF2/PAINS performance hang on symmetric
macrocycles (traced via `/usr/bin/sample`, attributed to `alerts.rs`/`drug_score.rs`, out of
this file's scope).

- **Files:** `scripts/descriptor_census.py`, `scripts/descriptor_census_corpus.smi` (5,000
  freshly-downloaded ChEMBL SMILES), `crates/chematic-chem/examples/descriptor_census_unbound.rs`
  (dumps the 5 functions with no Python/WASM/MCP binding, plus `bcut2d`/`carbon_types` which have
  no individual Python getter), `descriptor_census.json` / `descriptor_census_unbound.jsonl`
  (`validation/results/`).
- **Reference tool:** RDKit 2026.03.3.
- **Full writeup:** [`docs/rfcs/descriptor_census_rfc.md`](../docs/rfcs/descriptor_census_rfc.md).
- **How to regenerate:**
  ```bash
  cargo run -p chematic-chem --release --example descriptor_census_unbound \
      < scripts/descriptor_census_corpus.smi > validation/results/descriptor_census_unbound.jsonl
  .venv/bin/python scripts/descriptor_census.py \
      --corpus scripts/descriptor_census_corpus.smi \
      --unbound validation/results/descriptor_census_unbound.jsonl \
      --json validation/results/descriptor_census.json
  ```

### IO-1: SMILES table file I/O (`.smi`/`.smiles`/`.csv`/`.tsv`/`.txt`)

New streaming `SmilesRecordReader`/`SmilesRecordWriter` in
`crates/chematic-mol/src/smiles_table.rs`, built from a source-cited audit of
RDKit's `SmilesMolSupplier`/`SmilesWriter` (RDKit commit
`8afba32ec539dcb2369bc84549d802aca3f7eb39`, the true resolution of tag
`Release_2026_03_4`) rather than guessed behavior. `chematic.SmilesMolSupplier`/
`chematic.SmilesWriter` (Python, pyo3) match RDKit's constructor signatures;
`rdkit_compat.py`'s own wrapper classes of the same names were rewritten to
delegate to these (previously a separate, non-streaming, whole-file pure-Python
implementation that used no Rust parser at all).

**Oracle methodology (deliberately avoids the chematic/RDKit canonical-SMILES
divergence, which is a known, separately-tracked, unrelated issue):** never
compares chematic-canonical vs. RDKit-canonical SMILES strings directly.
Instead: (1) exact string equality of extracted `name`/property values
(pure tokenization, no chemistry); (2) each tool's *own* self-consistency
against the fixture's known ground-truth SMILES, canonicalized only within
that same tool — proves each tokenizer extracted the right substring
without ever comparing the two canonicalizers against each other.

**Results (235 rows, 8 scenarios covering space/tab/comma delimiters,
header/no-header, name/no-name column, extra properties, quoted CSV, blank
lines, comments, malformed SMILES, isotopes, charges, disconnected
fragments, stereochemistry):**

| Metric | Result |
|---|---|
| Status parity (success vs. unparseable, per row) | 235/235 (100%) |
| Known-malformed rows correctly rejected by both tools | 5/5 |
| Name/property exact-match (excluding 2 documented divergences below) | 100% |
| Chematic self-consistency vs. known ground truth | 230/230 (100%) |
| RDKit self-consistency vs. known ground truth (non-gating) | 230/230 (100%) |

**Two deliberate, documented divergences from RDKit found via this oracle
(not bugs — see `smiles_table.rs`'s module doc comment for the full
citation):**
1. `name_column=None` (RDKit's `nameColumn=-1`): RDKit falls back to the
   physical line number as `_Name`; chematic's `MoleculeRecord::name` is
   simply empty. Judged a low-value RDKit implementation detail not worth
   reproducing.
2. **RDKit's `SmilesMolSupplier` has no CSV-quote-awareness at all** for its
   comma-delimiter mode — confirmed via this oracle, not merely inferred: a
   quoted field like `"has, a comma"` is split into extra raw columns by
   RDKit's literal comma-splitting. Chematic implements a real RFC 4180
   *subset* (quoted fields, doubled-quote escaping, no multi-line quoted
   fields) — a genuine improvement, not a matched behavior.

CXSMILES is not recognized in the SMILES column (parsed via
`chematic_smiles::parse`, which has no CXSMILES support) — this matches
RDKit's *own* default for `SmilesMolSupplier`, which explicitly disables
CXSMILES for this entry point too.

**Performance** (10,000-record synthetic corpus, ~2% deliberately malformed rows to also
measure invalid-record recovery throughput; 5 independent process runs, `/usr/bin/time -l`):

| | chematic (`SmilesRecordReader`) | RDKit (`SmilesMolSupplier`, Python) |
|---|---|---|
| Records/sec (median of 5 runs) | ~137,000 | ~4,200 |
| Peak RSS | ~2.3 MB | ~45 MB |
| Success/error split | 9,800/200 | 9,800/200 (identical) |

The ~30x throughput difference is Python-interpreter-call overhead, not a controlled
same-language comparison — reported as reference/informational only, per the
"performance is never traded for correctness, and cross-language numbers aren't a gate"
policy. Both tools agree exactly on which 200/10,000 rows are malformed.

**Adversarial/fuzz-style coverage:** no `cargo-fuzz`/libfuzzer harness exists anywhere in
this workspace yet, and introducing that toolchain for one text-tokenizer module was judged
disproportionate — instead, 9 deterministic adversarial unit tests (empty input, truncated
input mid-record and mid-quote, a line exceeding `max_line_bytes`, a 500KB property value
within the limit, invalid-UTF-8 byte handling, 5,000-column rows, a 3,000-atom SMILES field,
and a 2,000-iteration seeded random-mutation corpus) assert only "no panic, no hang, no OOM" —
never a specific output — since malformed input must degrade to a clean `Err`, never worse.

- **Files:** `crates/chematic-mol/src/smiles_table.rs`, `crates/chematic-mol/examples/smiles_table_dump.rs`,
  `crates/chematic-mol/examples/smiles_table_benchmark.rs`,
  `scripts/gen_smiles_table_fixtures.py`, `scripts/gen_rdkit_smiles_table_oracle.py`,
  `scripts/smiles_table_io_parity.py`.
- **Reference tool:** RDKit 2026.03.3.
- **How to regenerate:** `python scripts/gen_smiles_table_fixtures.py --out-dir <dir> --corpus
  <SMILES.csv> --manifest-out <manifest.json>` + `cargo run -p chematic-mol --release --example
  smiles_table_dump -- <manifest.json> <dir> <out.jsonl>` + `python
  scripts/gen_rdkit_smiles_table_oracle.py --manifest <manifest.json> --fixtures-dir <dir> --out
  <oracle.jsonl>` + `python scripts/smiles_table_io_parity.py --chematic <out.jsonl> --rdkit-oracle
  <oracle.jsonl> --manifest <manifest.json>`. Self-test:
  `python scripts/smiles_table_io_parity.py --self-test`. The generated fixture files themselves
  are not committed (regenerable byte-for-byte from the corpus + script); only the generator
  scripts and this summary are.

### IO-2: Daylight TDT file I/O (`.tdt`)

New streaming `TdtRecordReader`/`TdtRecordWriter` in `crates/chematic-mol/src/tdt.rs`, reusing
`MoleculeRecord`'s (IO-1's shared record type) `coordinates_2d`/`coordinates_3d` fields (unused
by the SMILES-table format, exactly what TDT's `2D`/`3D` tags need). Built from a source-cited
audit of RDKit's `TDTMolSupplier`/`TDTWriter` (same pinned commit
`8afba32ec539dcb2369bc84549d802aca3f7eb39`).

**Four deliberate, documented divergences from RDKit — three fix real bugs the audit/oracle
found in RDKit itself, not guessed improvements:**
1. **RDKit's own coordinate-list parser drops the last atom's position** — its comma-tokenizer
   treats the token containing the trailing `;>` as "found the terminator" and never pushes that
   token's own numeric value. Confirmed against a live RDKit run (both via source tracing during
   the audit, and reproduced again by this PR's own oracle: the last atom comes back at
   `(0,0,0)` from real RDKit). Chematic parses the full list correctly.
2. **RDKit's `TDTWriter` hard-codes its name tag as `"NAME"` while `TDTMolSupplier`'s own
   `nameRecord` defaults to `""`** (no tag recognized) — a bare RDKit writer+reader round trip
   silently loses the molecule name by default, confirmed empirically. Chematic's reader/writer
   both default `name_tag` to `Some("NAME")`, so the round trip preserves the name out of the box.
3. **RDKit's malformed-tag recovery has an infinite-loop hazard** — the exception thrown for a
   missing `>` isn't caught inside `TDTMolSupplier::next()`'s own position-advance bookkeeping,
   so naively retrying re-throws on the same record forever (confirmed empirically; this PR's own
   RDKit oracle script had to switch to explicit index-based access, `sup[idx]`, to work around
   it and make progress). Chematic's reader always scans forward to the next record boundary
   internally before returning `Err`, so the next `Iterator::next()` call always makes progress.
4. **RDKit drops the final tag line when a file has no trailing newline — discovered via this
   PR's own oracle run, not predicted by the initial source audit.** `$SMI<CC>\nNAME<ethane>`
   with no trailing `\n` yields `_Name == ""` in real RDKit rather than `"ethane"`. Chematic's
   `BufRead::read_line`-based reader returns a final unterminated line correctly.

**Results (205 rows, 8 scenarios covering SMI/NAME/arbitrary-property records, `|` terminator,
empty property, repeated property (last-wins), unknown tags, malformed-tag recovery, EOF
mid-record, 2D/3D coordinate tags, isotopes, charges, disconnected fragments, stereochemistry):**

| Metric | Result |
|---|---|
| Status parity | 205/205 (100%) |
| Known-malformed row correctly rejected by both tools | 1/1 |
| Name/property/coordinate exact-match (excl. the 2 scenarios with documented divergences above) | 100% |
| Chematic self-consistency vs. known ground truth | 204/204 (100%) |
| RDKit self-consistency vs. known ground truth (non-gating) | 204/204 (100%) |

**Performance** (10,000-record synthetic corpus, ~2% deliberately malformed, 5 independent
process runs, `/usr/bin/time -l`):

| | chematic (`TdtRecordReader`) | RDKit (`TDTMolSupplier`, Python, index-based access) |
|---|---|---|
| Records/sec (median of 5 runs) | ~128,000 | ~5,100 |
| Peak RSS | ~2.4 MB | ~45 MB |
| Success/error split | 9,800/200 | 9,800/200 (identical) |

RDKit's own numbers use index-based access (`sup[idx]`), the workaround this PR's own oracle
script needed for RDKit's malformed-tag recovery hazard — an even less apples-to-apples
comparison than a plain iterator would be. Reported as informational cross-language reference
only, same policy as IO-1.

**Adversarial/fuzz-style coverage:** 11 deterministic adversarial unit tests (empty/truncated
input, an oversized line, a 500KB property value, invalid-UTF-8 bytes, 5,000-property records, a
3,000-atom SMILES field, a malformed/never-terminated coordinate list, a 2,000-iteration seeded
random-mutation corpus) — same no-cargo-fuzz-harness-exists rationale as IO-1. **One of these
tests caught a real bug before this PR shipped**: an oversized-line or otherwise-mid-record read
error left the reader's position at a leftover fragment (e.g. the record's own `|` terminator,
already consumed into a since-rejected oversized buffer) that got misinterpreted as the start of
a phantom next record. Fixed by centralizing recovery-to-next-record-boundary at every error exit
from record-body parsing, not just the malformed-tag case that was tested first.

- **Files:** `crates/chematic-mol/src/tdt.rs`, `crates/chematic-mol/examples/tdt_dump.rs`,
  `crates/chematic-mol/examples/tdt_benchmark.rs`,
  `scripts/gen_tdt_fixtures.py`, `scripts/gen_rdkit_tdt_oracle.py`, `scripts/tdt_io_parity.py`.
- **Reference tool:** RDKit 2026.03.3.
- **How to regenerate:** `python scripts/gen_tdt_fixtures.py --out-dir <dir> --corpus <SMILES.csv>
  --manifest-out <manifest.json>` + `cargo run -p chematic-mol --release --example tdt_dump --
  <manifest.json> <dir> <out.jsonl>` + `python scripts/gen_rdkit_tdt_oracle.py --manifest
  <manifest.json> --fixtures-dir <dir> --out <oracle.jsonl>` + `python scripts/tdt_io_parity.py
  --chematic <out.jsonl> --rdkit-oracle <oracle.jsonl> --manifest <manifest.json>`. Self-test:
  `python scripts/tdt_io_parity.py --self-test`.

### IO-3: ChemAxon Marvin file I/O (`.mrv`)

New `parse_mrv`/`write_mrv` in `crates/chematic-mol/src/mrv.rs`, reusing `MoleculeRecord` (IO-1's
shared record type). A purpose-built, dependency-free, nesting-aware XML tokenizer is used
instead of `cml.rs`'s line-by-line flat scanner (too weak for MRV's nested elements and child
text content, e.g. `<bondStereo>W</bondStereo>`) and instead of adding an external XML crate
dependency (whose default DTD/entity posture would need independent verification against this
module's own mandated security limits anyway). Built from a source-cited audit of RDKit's
`MolFromMrvBlock`/`MolToMrvBlock` (same pinned commit `8afba32ec539dcb2369bc84549d802aca3f7eb39`).

**Scope, deliberately bounded, not an incomplete draft:** `molecule`/`atomArray`/`bondArray`, atom
IDs, element/charge/isotope/atom-map, bond order (single/double/triple/aromatic), 2D/3D
coordinates, wedge/dash stereo. S-groups, polymers, reactions, multicenter bonds, query
atoms/bonds, R-groups, enhanced stereo groups, and embedded/compressed data are explicitly out of
scope and return `MrvError::UnsupportedFeature` rather than a silent partial parse — RDKit itself
*does* support several of these; chematic's port does not, by explicit scope decision.

**Security (chematic-only — RDKit's own MRV parser passes no `xml_parser_flags` to Boost's
`read_xml` at all, confirmed via this session's source audit, so "RDKit does it" cannot justify
skipping any of this):** `<!DOCTYPE`/`<!ENTITY` rejected outright (verified against a
billion-laughs-style entity chain and an XXE local-file-reference payload, both rejected before
any chemistry is interpreted), input byte-size limit, element-nesting-depth limit,
attribute-value-length limit, duplicate atom-ID detection (a real robustness improvement over
RDKit, which has none), missing bond-endpoint detection.

**Oracle methodology — corrected mid-session after an initial, flawed design.** The first attempt
compared chematic's own canonical SMILES against chematic's own canonical SMILES of the fixture's
known ground-truth SMILES (both canonicalized by chematic itself) and showed 149/206 "mismatches."
Investigating one case showed the two strings differed only in which ring-closure digit was
assigned to which ring — a canonicalizer-enumeration artifact, not a structural bug — but
building the *actual* fix (never compare two different canonicalizers' output strings directly;
re-canonicalize both sides through the *same* tool, RDKit) revealed that most of the original
149 were **not** ring-closure artifacts at all, but a real, separate, non-MRV-specific finding
(below). The corrected, final methodology, per fixture:

- **A** = RDKit reads the original `.mrv` fixture directly (`Chem.MolFromMrvBlock`).
- chematic reads the same fixture and emits isomeric SMILES; RDKit re-parses it → **B**.
- chematic's own MRV write of the parsed record is read by RDKit → **B′** (chematic-write →
  RDKit-read leg).
- `Chem.MolToSmiles(A)` vs `Chem.MolToSmiles(B)` / `Chem.MolToSmiles(B′)` — both sides
  canonicalized by RDKit — is the gate. A 0-atom RDKit "success" is never counted as a match
  (RDKit's own MRV parser can return an empty, error-free `RWMol` for some malformed input).
  InChIKey is logged as an auxiliary cross-check only, never the deciding signal.
- On mismatch: atom/bond count, element/charge/isotope multisets, bond-order histogram, aromatic
  atom count, fragment count, and — critically — the actual assigned stereocenter/E-Z **labels**
  (not just a raw stereocenter *count*, which can't distinguish "assigned R" from "unassigned")
  are compared independently to classify the cause.

**Results (206 RDKit-generated fixtures — 160 drawn from the M4-A0 corpus + hand-picked coverage
for acyclic/aromatic/fused-ring/isotope/charge/radical/tetrahedral-stereo/E-Z-stereo/atom-map/2D-
and-3D-coordinates/disconnected-fragments):**

| Metric | Result |
|---|---|
| Parse success | 206/206 (100%) |
| chematic read→write→read round trip (chematic-only, no RDKit) | 206/206 (100%) |
| Phase 1 (RDKit-read vs. RDKit-reread-of-chematic's-SMILES) exact match | 125/206 |
| Phase 1 unexplained (real structural diff) | **0** |
| Phase 2 (chematic-write → RDKit-read) exact match | 203/206 |
| Phase 2 unexplained (real structural diff) | **0** |

Every mismatch is accounted for — zero unexplained residuals. All 81 phase-1 and 3 phase-2
mismatches fall into one of three documented, non-gating buckets:

1. **69 cases (phase 1 only): tetrahedral/E-Z stereo is lost when converting a wedge-bond-encoded
   MRV molecule to SMILES.** chematic has a reader for wedge/dash bond direction + 2D coordinates
   → CIP R/S label (`chematic_perception::stereo2d::apply_stereo_from_2d`), but no converter from
   that representation into `Atom.chirality` (the field `chematic_smiles::write` actually reads).
   **MRV's own read→write→read round trip is fully correct and unaffected** (phase 2 matches for
   all of these) — only conversion to a *different* format (SMILES) loses the assignment. MOL
   V2000 has the identical limitation (same wedge-bond representation, same missing converter).
   Pre-existing, cross-cutting, **not specific to MRV** — flagged for a dedicated follow-up.
2. **9 cases: a pre-existing `chematic-smiles` writer bug, discovered incidentally via this
   oracle work.** `chematic_smiles::write()` omits the implicit-H count for a bracket atom forced
   by isotope/charge/atom-map when `Atom.hydrogen_count` is `None` (implicit/inferred — how every
   non-SMILES-parser format reader in the workspace builds atoms) rather than `Some(n)` (explicit,
   as the SMILES parser itself always sets). Confirmed via a minimal repro with **no MRV
   involvement at all**: `Atom::new(N)` + `charge=1` writes `[N+]` instead of `[NH4+]`. Affects
   every non-SMILES format reader (MOL/SDF/CML/CDXML/MOL2/PDBQT/etc.), not just MRV. **MRV's own
   write path is unaffected** (confirmed: 0 of these 9 recur in the phase-2, MRV-native write
   comparison — `write_mrv` writes the isotope/charge/atom-map attributes directly, never routing
   through the buggy SMILES-writer bracket logic). Flagged for a dedicated follow-up fix in
   `chematic-smiles`, out of scope for this PR.
3. **3 cases (both phases): radical-electron information loss on write.** `chematic_core::Atom`
   has no radical-electron slot — a radical atom (e.g. methyl radical, dioxygen biradical) gets an
   extra implicit H when RDKit re-reads chematic's write output. Documented, deliberate — the same
   convention as `mol2000.rs`'s pre-existing "doublet radical — treated as neutral."

**Dedicated kekulize/stereo option checks** (independent of the 206-fixture pool, since the main
pool always uses the default `kekulize=false, include_stereo=true`): `kekulize=True` writes
alternating single/double bonds (no `order="A"` token) and RDKit reads it back to the same
canonical structure as `kekulize=False`'s `order="A"` form (2/2 pass). `include_stereo=False`
correctly drops a tetrahedral wedge/dash assignment on RDKit re-read while preserving connectivity
(2/2 pass); `include_stereo` has **no effect on E/Z double-bond stereo** by design — E/Z is
geometry-derived from the 2D coordinates chematic always writes, so RDKit perceives it regardless
of the option (2/2 pass, correctly not expecting a drop).

**Adversarial/fuzz-style coverage:** 32 unit/adversarial tests total, including: empty input,
deeply-nested elements (depth limit), oversized attribute, oversized input, a billion-laughs-style
entity chain, an XXE local-file-reference payload, truncated input, non-ASCII multi-byte content
adjacent to byte-index slicing, excessive atom count (20,000 atoms), NaN/Infinity coordinate
values, and a 500-iteration seeded random-mutation corpus — none panic, hang, or OOM. Duplicate
atom IDs, unknown atom references, and unknown element symbols are covered by dedicated
non-adversarial unit tests. "Excessive property count" has no distinct MRV analogue (no arbitrary
key/value property list the way SDF/TDT have); it folds into the existing attribute-count/nesting
limits. "Huge line" has no distinct MRV analogue either (the parser is a whole-string
recursive-descent tokenizer, not line-based); the analogous protection is the overall
input-size/attribute-length limits, already covered.

**Performance** (206-document fixture pool, informational only): ~18,000 documents/sec. MRV is
one document per molecule (unlike SMILES-table/TDT's one-file-many-records shape), so this
reports parse throughput over a directory of individual `.mrv` files rather than a single
10k-record file.

- **Files:** `crates/chematic-mol/src/mrv.rs`, `crates/chematic-mol/examples/mrv_dump.rs`,
  `crates/chematic-mol/examples/mrv_kekulize_stereo_dump.rs`,
  `crates/chematic-mol/examples/mrv_benchmark.rs`, `scripts/gen_rdkit_mrv_oracle.py`,
  `scripts/mrv_io_parity.py`, `scripts/mrv_kekulize_stereo_check.py`.
- **Reference tool:** RDKit 2026.03.3, pinned commit `8afba32ec539dcb2369bc84549d802aca3f7eb39`.
- **How to regenerate:** `python scripts/gen_rdkit_mrv_oracle.py --corpus <SMILES.csv> --out-dir
  <fixtures_dir> --manifest-out <manifest.json>` + `cargo run -p chematic-mol --release --example
  mrv_dump -- <manifest.json> <fixtures_dir> <out.jsonl> <written_dir>` + `python
  scripts/mrv_io_parity.py --chematic <out.jsonl> --manifest <manifest.json> --fixtures-dir
  <fixtures_dir> --written-dir <written_dir> --summary-out <out.json>`. Self-tests:
  `python scripts/mrv_io_parity.py --self-test` and `python
  scripts/mrv_kekulize_stereo_check.py --self-test`.

### Wave 2C: ring-constrained E/Z residuals (issue #149, 5,000-mol corpus + 18 pinned fixtures)

Audit-only follow-up to PR #229 (Wave 2B), which split the 18 pinned
`EZ_SHARED_CANDIDATE_BOND_RESIDUALS` fixtures (`crates/chematic-smiles/src/
canonical.rs`) into 10 fully-resolved and 8 still-residual. Tests the
hypothesis that every residual's coupled component includes an alkene end
whose own double bond is endocyclic in a small ring -- confirmed
individually for all 8 (never pooled) via a live RDKit oracle, but found
necessary-not-sufficient (the same shape occurs in 5 of the 10 resolved
fixtures too). Measures the corpus-wide blast radius of 3 candidate
production gating rules, each independently, never combined. No production
code changed.

- **Files:** `validation/results/ez_ring_constrained_residual_audit.jsonl`
  (1,387 per-end rows, full corpus), `ez_ring_constrained_residual_audit_
  summary.json` (fixture classification + blast-radius table)
- **Reference tool:** RDKit 2026.03.3 (`Chem.FindPotentialStereo`,
  `Chem.AssignStereochemistry`)
- **How to regenerate:** `cargo run -p chematic-smiles --release --example
  ez_ring_constrained_residual_audit -- scripts/descriptor_census_corpus.smi`
  + `python scripts/ez_ring_constrained_residual_diagnosis.py`. Self-test:
  `python scripts/ez_ring_constrained_residual_diagnosis.py --self-test`.
- **Full report:** `docs/rfcs/ez_ring_constrained_residual_audit.md` (per-fixture
  classification, blast-radius table, recommended predicate,
  CONDITIONAL GO verdict)

### Wave 3: general shared-carrier coupling-component residual (issue #149, 5,000-mol corpus + 18 pinned fixtures + 2 never-corrupts fixtures)

Follow-up to PR #351 (Wave 2D), which closed only ~10% (3 of 31) of the
corpus's general shared-carrier coupling-component population per its own
commit message -- the other ~90% (28 of 31) was reported as "a separate,
still-unidentified mechanism." This audit measures that remaining
population directly instead of assuming its residual status from that
topological-presence figure. Result: 28 coupled components confirmed today
(all size-2 paths, 0 cycles), 0/28 confirmed permutation-invariance
failures under 16 RDKit relabelings each; a second probe (mark relocation)
found structurally incapable of testing coupled pairs at all. All 28 share
one RDKit-confirmed shape (both ends genuinely stereogenic). A calibration
check found the permanent never-corrupts regression fixture's own "does not
converge" doc claim likely stale. Verdict: NEEDS-RESEARCH, leaning GO
(already-likely-resolved) -- not proof. No production code changed.

- **Files:** `validation/results/ez_shared_carrier_coupling_mechanism_audit.jsonl`
  (per-component axis1/axis2/classification detail),
  `ez_shared_carrier_coupling_mechanism_audit_summary.json` (topology, axis
  summaries, structural buckets, calibration check, verdict)
- **Reference tool:** RDKit 2026.03.4 (`Chem.RenumberAtoms`,
  `Chem.FindPotentialStereo`, `Chem.AssignStereochemistry`)
- **How to regenerate:** `.venv/bin/python3 scripts/ez_shared_carrier_
  coupling_mechanism_diagnosis.py` (drives `cargo run -p chematic-smiles
  --release --example ez_shared_carrier_coupling_mechanism_audit`
  internally). Self-test: `.venv/bin/python3 scripts/ez_shared_carrier_
  coupling_mechanism_diagnosis.py --self-test`.
- **Full report:** `docs/rfcs/ez_shared_carrier_coupling_mechanism_audit.md`
  (method, per-bucket structural classification, field-by-field findings,
  NEEDS-RESEARCH/leaning-GO verdict)

### Explainable Standardization Phase 1 -- acceptance fixtures + holdout (implemented)

`crates/chematic-chem/src/standardize.rs` had no formal design/RFC prior to
this pass. Before drafting the RFC, three concrete defects in the
already-shipped code were confirmed empirically (a throwaway example run
against `main` at commit `1ac442a`, deleted after use, not committed): (a)
tied-fragment-size selection was spelling-order-dependent, not
atom-order-invariant (`largest_fragment` kept a different fragment for
`"CCC.CCN"` vs `"CCN.CCC"`); (b) fragment "size" counted raw atom count
including explicit hydrogens, not heavy atoms; (c) the named `SaltCatalog`'s
"ammonium" SMARTS entry false-positived on real organic cations (e.g.
choline), previously masked rather than fixed by an unrelated size
comparison. All three are now fixed and pinned as named regression tests.
Implementing the fixtures as tests also caught two more issues (see the
RFC's section 8): fragment extraction was silently corrupting
stereocenters (a `stereo_neighbor_order` remap bug shared with the legacy
`remove_salts_with_catalog` path, fixed for both at once), and one holdout
fixture's own note had an arithmetic error (phosphoric acid is 5 heavy
atoms, not 4).

- **Files:** `validation/standardization_phase1_fixtures.jsonl` (34 rows),
  `validation/standardization_phase1_holdout.jsonl` (10 rows, held out from
  rule design). Categories: simple salts, zwitterions, hydrates,
  organometallics, multi-organic fragments (incl. a genuine 2-API cocrystal
  case), isotope-containing compounds, equal-size ties (incl. 3 named
  regression fixtures), charged fragments, and edge cases (empty molecule,
  all-inorganic input, duplicate fragments, a carbon-containing solvate the
  default policy deliberately does not strip).
- **Implementation:** `crates/chematic-chem/src/standardize.rs` --
  `FragmentPolicy`, `select_fragment`, `TransformationRecord`/
  `FragmentRecord`/`FragmentSnapshot`/`FragmentDecision`, re-exported from
  `chematic-chem`'s crate root. `largest_fragment`/`remove_salts` are now
  thin wrappers around `select_fragment` with the default policy --
  `SaltCatalog`/`remove_salts_with_catalog` remain unchanged, opt-in-only.
- **How to regenerate:** all 44 fixtures are hand-transcribed as `#[test]`
  functions (`phase1_*`) in `standardize.rs`'s test module, not a
  JSONL-driven runtime harness -- `cargo test -p chematic-chem --lib
  standardize::tests::phase1`.
- **Full report:** `docs/rfcs/explainable_standardization_phase1_rfc.md`
  (status-quo defect audit, fragment-policy design, audit-log data model,
  section 8: implementation deviations, bugs found, and disclosed gaps)

## Summary results

See [rdkit/README.md](rdkit/README.md) for per-descriptor breakdowns.

| Metric | Corpus | Agreement |
|--------|--------|-----------|
| HBA / HBD / ARC | 4,999 mol | **100%** |
| MW, HAC | 175 mol | **100%** |
| TPSA | 175 mol (drug-like) | **100%** (±0.1 Å²) |
| TPSA | 4,999 mol | 93.3% (±0.1 Å²) |
| LogP | 175 mol | ~99% (±0.3) |

## Methodology

- Reference values are generated by RDKit Python API (rdkit-sys ≥ 2024.x)
- chematic values are computed via `chematic.from_smiles(smi).descriptors()`
- Agreement = fraction of molecules within the stated tolerance
- TPSA uses Ertl (2000) SMARTS-based approach in both tools
- Scripts are deterministic and pinned to RDKit 2026.03.3

## How to run

```bash
# Fast regression on 175-mol in-repo corpus (no download required)
pip install chematic rdkit pandas
python scripts/rdkit_benchmark.py

# Large-scale agreement on 5k ChEMBL subset
python scripts/bench5k.py ~/Downloads/SMILES.csv
python scripts/bench5k.py ~/Downloads/SMILES.csv --detail   # show mismatches

# Per-molecule LogP (Crippen) mismatch breakdown, beyond bench5k's aggregate number
python scripts/analyze_logp_mismatches.py ~/Downloads/SMILES.csv --tolerance 0.01
```
