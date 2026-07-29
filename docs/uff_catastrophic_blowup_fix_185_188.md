# UFF catastrophic-blowup fix: issues #185 and #188

Companion to `docs/uff_robustness_diagnosis_185_188.md` (diagnostic pass,
no fix). This document covers the actual fix: what was tried, what was
measured, what shipped, and what didn't.

## Recap of the diagnosis

`dg::generate_coords`'s starting geometries have enormous, unrelieved van
der Waals clash energy for essentially every molecule (not just the ones
that end up unsound). `chematic_ff::uff::minimize_uff`'s plain fixed-step
steepest descent eventually relieves this for every molecule tested, but
how many iterations a specific molecule's clash-relief trajectory needs is
not predicted by any single starting-geometry metric (worst bond length,
starting energy) — and can be 1-3+ orders of magnitude larger than the
default 200-iteration budget. `embed_distance_geometry_v2`'s geometry does
not have this problem: every molecule tested reaches a sound geometry within
budget from that start.

## Ablation: six candidates measured independently

Before writing any fix, six candidate approaches were implemented and
measured against the same 11 blow-up molecules (both #185's and #188's
named cases) plus 4 negative-control molecules, all from identical starting
geometries:

| Candidate | Description | Result |
|---|---|---|
| A | Raise `max_iter` to 20,000 (algorithm unchanged) | 7/11 sound; hexane/penicillin_core/testosterone/cholesterol still blow up; up to 58s wall-clock on the worst case |
| B | Fix the step-growth condition (`energy - new_energy < prev_energy * 1e-7`, where `prev_energy` starts at `f64::MAX`, making growth almost unconditional early on) | Same 7/11 as A, but ~2x fewer iterations where it succeeds — a real, low-risk bug fix, but doesn't rescue anything A doesn't already rescue |
| C | Backtracking line search (rejected steps don't consume the outer iteration budget) | 9/11 sound, cheaper than A/B where it works — but testosterone/cholesterol still blow up even at 20,000 evaluations |
| D | Adaptive stall-detection (extend budget unless worst-bond progress stalls) | **Unsafe as implemented**: the specific stall heuristic (check every 500 iterations, require >1% improvement) triggers false early termination — naphthalene needs 4,539 iterations but the check gives up at 500, well before the real, slow-crawl-then-relieve trajectory turns around |
| E | Cheap vdW-only pre-relaxation before the standard minimizer | **Dangerous**: with no bond/angle term to hold chain atoms together, this explodes flexible molecules — cholesterol's worst bond reached 3,054,035 Å |
| F | Detect a catastrophic first-attempt failure and retry from `embed_distance_geometry_v2`'s geometry instead | **11/11 sound**, at the unchanged 200-iteration budget, cheaply — the only candidate that fully closes the gap |

Full outcome matrix and raw numbers: see the PR description / `scratchpad`
ablation logs (not committed — this doc is the permanent record of the
methodology and result).

## What shipped

**F only.** B, C, D, and E are not included in this fix:
- B is a real, valid bug, but doesn't fix anything F doesn't already fix on
  its own, and touching `minimize_uff`'s core algorithm (used by every other
  caller, and pinned by several exact-iteration-count regression tests)
  is a larger, separate change with its own risk. Left for a future,
  independent, narrowly-scoped perf PR — not mixed into this one.
- C, D, E are all rejected outright by the ablation (D is actively unsafe;
  E is actively dangerous for flexible chains; C alone doesn't close the gap
  on the two hardest cases).

### Design

`run_uff_bridge` (`crates/chematic-3d/src/minimize.rs`) tries the
caller-provided coordinates first, exactly as before — every existing
caller whose first attempt already succeeds sees zero behavior change. Only
when that first attempt fails soundness with `CatastrophicBondBlowup` or
`NonFiniteCoordinates` specifically does it retry once from
`embed_distance_geometry_v2`'s geometry. This is a **post-hoc retry on the
actual outcome**, not a pre-minimization heuristic: raw starting energy
cannot predict which molecules need the rescue (anthracene's raw
`generate_coords` energy is ~5 orders of magnitude worse than naphthalene's,
yet anthracene never needs it).

Never silent: `PolicyMinimizeResult::starting_geometry` (`Option<UffStartingGeometry>`,
`Some` whenever UFF actually ran) discloses `AsProvided` vs.
`ReplacedWithDistanceGeometryV2` on every success. A rescue attempt that
doesn't help is disclosed too, via
`MinimizationFailureDetail::distance_geometry_v2_retry_attempted`, and the
ORIGINAL failure (not a new, possibly-confusing one from the retry) is what
gets returned.

### A second real bug found during implementation: the rescue can silently invert declared stereochemistry

`embed_distance_geometry_v2`'s default parameters do not enforce declared
chirality (`EmbedParameters::enforce_chirality` defaults to `false`), and
setting it `true` does not *fix* chirality — it makes embedding itself
refuse any molecule with declared stereo outright (a "fail closed" design,
per that module's own doc, not a correction mechanism). Measured directly
via `verify_stereo` (`chematic_3d::stereo_constraints`) on the corpus's
chiral blow-up molecules: penicillin_core's rescue happens to preserve all
3 declared stereocenters, but testosterone (2/6 violated) and cholesterol
(3/8 violated) do not.

Shipping a bond-length-sound geometry that silently inverts declared
stereochemistry would be a **worse** outcome than the honest failure it
replaced. `rescue_with_distance_geometry_v2` therefore requires
`verify_stereo(...).is_fully_satisfied()` in addition to
`check_minimization_soundness` before accepting a rescue — if the rescued
geometry violates a declared stereocenter or E/Z bond, the rescue is refused
and the original failure returned. This check is scoped to the rescue path
only: the pre-existing, unrelated fact that this bridge's first-attempt path
has never verified stereo either is a known, disclosed, out-of-scope gap
(`examples/cf_integration_smoke_test.rs`'s own closing note), not one this
fix introduces or is scoped to close.

## Measured result (58-molecule corpus, `examples/cf_integration_smoke_test.rs`)

`UffOnly` end-to-end success, legacy `generate_coords` source:

| | before this fix | after this fix |
|---|---|---|
| Sound | 41/58 | **53/58** |
| Blown up | 17/58 | **5/58** |

The 5 remaining failures — **ibuprofen_S, naproxen_S, testosterone,
cholesterol, atorvastatin_fragment** — are exactly the declared-stereo
molecules whose rescue geometry doesn't preserve stereochemistry. All 5
already succeed when fed `embed_distance_geometry_v2`'s geometry directly
(the `embedder-fed` column) — this is specifically a legacy-start-plus-
chirality gap, tracked as issue #210, not fixed here.

No regression is possible in this count: `rescue_with_distance_geometry_v2`
is only ever invoked from the branch where the first attempt already failed
with `CatastrophicBondBlowup`/`NonFiniteCoordinates` — the first-attempt
code path itself is completely unchanged, so any molecule that succeeded
before still succeeds identically (same code, same inputs), and any molecule
reported as a failure after this fix was, by construction, already a failure
before it.

`Mmff94BondAngleStrict` (never touches UFF) is unaffected: 32/58 before and
after, confirmed by the same corpus run. `Mmff94WithUffFallback` moves from
41/58 to 55/58 (same shared mechanism, since its fallback path is exactly
`run_uff_bridge`).

## Performance (measured, not estimated)

Per-molecule wall-clock, `UffOnly` + `generate_coords`, 33-molecule subset of
the corpus, split by which path each molecule took:

| Path | n | p50 | p95 | p99/max |
|---|---|---|---|---|
| As-provided (no rescue triggered) | 19 | 16.68ms | 81.63ms | 109.09ms |
| Rescued (successful) | 12 | 169.60ms | 276.56ms | 326.09ms |
| Rescue attempted, still failed (stereo-refused) | 2 | 1182.07ms | 1182.07ms | 1182.07ms |

Molecules whose first attempt already succeeds (the overwhelming majority)
pay **zero** extra cost — the rescue path is never entered. Only the
previously-failing molecules pay a one-time extra cost (one failed attempt
+ one embed + one retry minimization), which is the correct tradeoff: they
were hard failures before, at any cost.

## Test plan

- `crates/chematic-ff/src/uff.rs`: `uff_energy_breakdown_total_matches_uff_total_energy`,
  `traced_minimizer_matches_untraced_minimizer_byte_for_byte` — pin the new
  `uff_energy_breakdown`/`minimize_uff_with_trace`/`uff_worst_bond_length`
  diagnostic instrumentation (added during this fix's own ablation work, to
  measure per-term energy and per-iteration trajectories across the six
  candidates) as byte-identical/numerically-agreeing with the existing,
  unmodified `uff_total_energy`/`minimize_uff` — additive only, zero risk to
  either published function.
- `crates/chematic-3d/src/minimize.rs`:
  - `mmff94_with_uff_fallback_reports_typed_failure_when_fallback_itself_is_unsound`
    and the renamed `uff_only_succeeds_for_hexane_from_generate_coords_via_distance_geometry_v2_rescue`
    — rewritten from their pre-fix "pins the bug" form to "pins the fix,"
    since that's exactly what changed.
  - `uff_only_reports_as_provided_when_no_rescue_was_needed` — negative
    control: a molecule whose first attempt already succeeds must report
    `AsProvided`, never a spurious rescue.
  - `uff_only_refuses_a_rescue_that_would_violate_declared_stereochemistry` —
    cholesterol's rescue is measured to violate stereo; pins that this
    surfaces as the original typed failure, `distance_geometry_v2_retry_attempted: true`,
    never a silent `Ok`.
  - `chematic_ff_own_uff_minimizer_blows_up_naphthalene_independent_of_this_bridge`
    — unchanged; this fix lives entirely in `run_uff_bridge`, not in
    chematic-ff, so chematic-ff's own minimizer still reproduces the
    original blow-up in isolation, as it should.
- `examples/cf_integration_smoke_test.rs` — full 58-molecule corpus gate,
  re-run to produce the before/after table above.

## Not fixed here (scope)

- `minimize_uff`'s step-growth bug (candidate B) — real, but independent of
  this fix's actual gap-closing mechanism; left for a future, narrowly-scoped
  perf PR.
- The stereo-blind rescue gap (5/58 remaining) — issue #210.
- `pipeline_v2`'s own embedding path already uses `embed_distance_geometry_v2`
  directly and was never affected by this bug in the first place; this fix
  only changes `run_uff_bridge`'s behavior for callers that hand it a
  `generate_coords`-sourced (or otherwise severely clashed) starting
  geometry.
