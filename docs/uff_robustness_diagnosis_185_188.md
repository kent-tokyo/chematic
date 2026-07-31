# UFF robustness diagnosis: issues #185 and #188

Diagnostic pass only — no fix in this document or its companion PR. Regression
fixtures pinning the measurements below live in
`crates/chematic-3d/src/minimize.rs`'s `mod tests` (search for "Issues #185 /
#188 diagnosis").

## Question

Both issues report `ForceFieldPolicy::UffOnly` blowing the worst bond length
past `MAX_SANE_BOND_LENGTH` (3.0 Å) for some molecules and not others:
- **#185**: naphthalene blows up; anthracene, a larger fused-ring aromatic,
  does not.
- **#188**: 8 non-aromatic molecules (hexane, decane, triethylene glycol,
  hexanediol, hexadecane, penicillin core, testosterone, cholesterol) blow up
  under the same policy.

The task instruction governing this investigation explicitly required *not*
assuming #185 and #188 share a root cause. They were investigated
independently below; the shared mechanism stated in the conclusion is a
measured result, not a starting assumption.

## Method

Three starting-geometry sources exist in `chematic-3d`:
- `dg::generate_coords` — legacy, deterministic, rule-based ring/chain
  template placer.
- `distance_geometry_v2::embed_distance_geometry_v2` — stochastic distance
  geometry + classical MDS + bounds refinement. No force-field minimization.
- `pipeline_v2::embed_pipeline_v2` — a 12-stage pipeline whose embedding stage
  calls `embed_distance_geometry_v2_with_adjustments` (confirmed by reading
  `pipeline_v2.rs:637`), i.e. the same embedder as above with optional
  pair-bound adjustments layered in. **Not independently re-measured
  numerically in this pass** — traced to the same underlying embedder rather
  than measured as a third distinct starting point; stated here as an
  explicit gap, not silently assumed equivalent.

For each molecule, both sources' raw (pre-minimization) worst bond length were
measured, then `chematic_ff::uff::{assign_uff_types, minimize_uff}` was run
directly (bypassing this crate's policy bridge, matching the existing
`chematic_ff_own_uff_minimizer_blows_up_naphthalene_independent_of_this_bridge`
test's isolation) at varying iteration budgets, recording energy and worst
bond length at each cutoff.

`MinimizeConfig::max_steps` defaults to **200** (`minimize.rs:72`), threaded
into `run_uff_bridge` at both of its call sites (`minimize.rs:1998,2045`) —
this is the actual production iteration budget `ForceFieldPolicy::UffOnly`
and `Mmff94WithUffFallback` use, and the budget both issues' reports
implicitly measured against.

## Finding 1: raw starting geometry is not the discriminator

| molecule | `generate_coords` worst bond (Å) | `embed_distance_geometry_v2` worst bond (Å) |
|---|---|---|
| naphthalene | 2.2644 | 1.7503 |
| anthracene | **3.7346** | 1.7600 |
| quinoline | 2.2644 | 1.7143 |
| hexane | 1.5400 | 1.5900 |
| penicillin_core | 2.6791 | 1.9416 |

Naphthalene's raw geometry (2.26 Å) is *better* than anthracene's (3.73 Å) —
already above the 3.0 Å "sane" threshold before any minimization runs. This
rules out a dg.rs ring-placement defect specific to naphthalene's fused-ring
shape: the molecule that ends up sound (anthracene) starts from the worse raw
geometry of the two.

## Finding 2: raw UFF energy is not the discriminator either

| molecule | UFF energy from `generate_coords` (kcal/mol) | UFF energy from `embed_distance_geometry_v2` (kcal/mol) |
|---|---|---|
| naphthalene | 125,557.70 | 109.06 |
| anthracene | 27,452,810,971.38 | 174.99 |
| quinoline | 108,433.33 | 106.84 |
| hexane | 15,488,929.48 | 3.78 |
| penicillin_core | 158,309,885,781.48 | 205.21 |

`generate_coords` produces enormous starting energies (dominated by
unrelieved van der Waals overlap, not bond stretch) for *every* molecule
tested — not only the ones that end up blowing up. Anthracene's starting
energy (2.7×10¹⁰) is five orders of magnitude worse than naphthalene's
(1.3×10⁵), yet anthracene is the one that converges to a sound geometry at
the default 200-step budget. Starting energy alone does not predict outcome.

## Finding 3: full trajectories — under-budgeted, not divergent

Energy and worst bond length at increasing `minimize_uff` iteration cutoffs,
same `generate_coords`-seeded start throughout (deterministic — no
randomness in either `generate_coords` or the minimizer):

**naphthalene** (default budget: 200 steps)

| max_iter | energy | worst bond (Å) | converged |
|---|---|---|---|
| 200 (default) | 4,727.84 | 4.7358 | false |
| 400 | 1,803.72 | 3.4355 | false |
| 800 | 775.51 | 2.5053 | false |
| 2,000 | 662.62 | 2.2114 | false |
| 5,000 | 660.58 | 2.1876 | **true (at iter 4,539)** |

**anthracene** (converges comfortably inside the default budget)

| max_iter | energy | worst bond (Å) | converged |
|---|---|---|---|
| 8 | 6,685,535,898 | 4,300.82 | false |
| 20 | 19,456,533 | 169.28 | false |
| 100 | 7,970.28 | 4.21 | false |
| 200 (default) | 24.74 | 1.4651 | false |
| 400+ | 24.60 | 1.4688 | **true (at iter 262)** |

**hexane** (#188 class; needs far more steps than naphthalene)

| max_iter | energy | worst bond (Å) |
|---|---|---|
| 200 (default) | 5,159,575.89 | 165.79 |
| 5,000 | 710,035.72 | 45.72 |
| 20,000 | 144,108.51 | 23.80 |
| 100,000 | 244.66 | **2.04 (sound)** |

Anthracene's own trajectory transiently blows one bond out to 4,300 Å at
iteration 8 before fully recovering by iteration 262 — the same qualitative
"clash-relief overshoot then recover" pattern as naphthalene and hexane, just
compressed into far fewer steps. All three molecules follow the same
mechanism; they differ only in how many steps their specific clash-relief
trajectory needs. Given enough iterations, every molecule checked in this
pass (naphthalene, anthracene, hexane) reaches a sound, low-energy geometry —
none diverges permanently. This is measured directly (`cargo test -p
chematic-3d --lib -- --ignored uff_direct`, both tests pass) and pinned as a
regression fixture, not left as an inference.

## Conclusion

**Shared mechanism across #185 and #188**: `generate_coords` reliably
produces starting geometries with severe, unrelieved van der Waals clashes
(reflected in enormous starting UFF energies) for essentially any molecule.
`minimize_uff`'s plain fixed-step steepest descent (no conjugate-gradient or
quasi-Newton acceleration) does relieve this and reach a sound geometry given
enough iterations — but how many iterations a given molecule's specific
clash-relief trajectory needs is not predicted by any single starting metric
measured here (worst bond, total energy), and can be 1–3+ orders of magnitude
larger than the 200-step default budget `ForceFieldPolicy::UffOnly` and
`Mmff94WithUffFallback` use. Molecules whose trajectory happens to settle
within 200 steps (anthracene) report success; molecules that need more
(naphthalene: ~4,500; hexane: >20,000, sound by 100,000) are correctly
rejected as unsound by `check_minimization_soundness`'s bond-length gate —
which is doing its job — but the underlying cause is an iteration budget
mismatch, not a broken algorithm or a `dg.rs`-specific defect.

This is a measured conclusion, not an assumption made before investigating:
Findings 1–2 above independently ruled out a naphthalene/anthracene
structural difference and a starting-energy-magnitude explanation before
Finding 3 established the shared mechanism.

**A candidate mechanism for *why* the required iteration count varies so
widely** (documented as a hypothesis for the fix PR to investigate, not
proven or fixed here): `minimize_uff`'s step-size growth rule
(`crates/chematic-ff/src/uff.rs`, inside `minimize_uff`) grows the step by
1.2× whenever `energy - new_energy < prev_energy * 1e-7`. Early in a
trajectory `prev_energy` can itself be enormous (10⁵–10¹⁰ range here), making
that relative threshold enormous too — so the step can grow on almost any
early iteration where the *absolute* energy drop is nowhere near converged,
not only near an actual minimum. This is a plausible source of the erratic,
molecule-specific step trajectories observed above, but was not isolated or
proven as the mechanism in this pass — flagging it for the fix PR rather than
patching it here (no step-clamping/symptom-masking, per task instruction).

## Not fixed here (scope)

- No change to `MinimizeConfig::max_steps`'s default, `minimize_uff`'s
  algorithm, or any other existing default/behavior.
- No new `chematic-ff` public API (an energy-breakdown or termination-reason
  type would help future diagnosis, e.g. distinguishing "gradient converged"
  from "step size collapsed to the 1e-8 floor" from "budget exhausted" — none
  of `UffMinimizeResult`'s three possible return paths are currently
  distinguishable from outside `minimize_uff` except by re-deriving them from
  `iterations`/`converged`). Deferred to the fix PR, since it would be a
  public-API change on a published crate and this program's own instructions
  require explicit authorization for that, not implicit inclusion in a
  diagnostic PR.
- Per-term (bond/angle/vdW) UFF energy breakdown, line-search step size at a
  given iteration, and non-finite-onset point were not obtained — none are
  observable through `chematic-ff`'s current public API, and reproducing them
  by re-deriving the internal formulas externally would duplicate
  (and risk drifting from) `uff_total_energy`'s real implementation. Same
  deferral as above.
- `embed_pipeline_v2`'s pre-force-field coordinates were traced to the same
  underlying embedder by reading source, not independently re-measured
  numerically (see Method).

## Two items spun out, not part of #185/#188

- **`chematic_3d::generate_and_minimize_uff` (`lib.rs:162`, doc comment
  "Generate 3D coordinates and minimize using UFF force field") does not run
  chematic-ff's UFF.** It calls `generate_coords` then this crate's own
  `minimize::minimize_uff` (`minimize.rs:88` — itself honestly documented as
  "identical to calling `minimize(mol, coords)`"), which delegates to the
  generic `minimize()` DREIDING-default hand-rolled harmonic force field,
  never `chematic_ff::uff::minimize_uff`. The inner function's own doc
  comment discloses this; the outer, public `generate_and_minimize_uff`'s
  doc comment does not, and its name claims a force field it never runs.
  Filed as its own issue (not fixed here — a public function rename/behavior
  change needs its own authorization).
  - **Resolved in issue #204's fix PR.** `#[deprecated]`, with an honest doc
    comment. The fix investigation also refined the framing above: the
    "DREIDING-default" phrasing here undersells it — `minimize_with_config`'s
    dispatch only special-cases `ForceField::MMFF94`, so `ForceField::UFF`
    and `ForceField::DREIDING` are indistinguishable on that path and both
    fall through to the same generic, untyped, element-pair harmonic engine.
    That engine is a **third** implementation, distinct from both real UFF
    and this crate's own typed DREIDING engine (`minimize_dreiding`, which
    `generate_and_minimize_dreiding` actually calls) — confirmed by a test
    that initially asserted DREIDING-equivalence and failed. See
    `CHANGELOG.md`'s `[Unreleased]` entry for the corrected description.
- **Issue #185's reported "~481.27 Å" worst-bond number could not be
  reproduced.** Naphthalene's worst bond length from `generate_coords` at
  the default 200-step `UffOnly` budget, measured independently twice — once
  in this diagnostic pass (pinned as a regression test above) and once
  already in the previously-merged
  `mmff94_with_uff_fallback_reports_typed_failure_when_fallback_itself_is_unsound`
  test's own docstring ("1.43 Å -> 4.74 Å") — is **4.7358 Å**, not ~481 Å,
  and the full 200-step trajectory (Finding 3) never passes anywhere near
  481 Å at any point. 481.27 is plausibly a `max_residual_force` value
  (kcal/mol/Å) rather than a bond length in Å: values in that range are
  exactly what this bridge's own code documents for *other* blown-up
  molecules' residual force (`MAX_SANE_RESIDUAL_FORCE`'s doc comment cites
  quinoline at 337.99). This is a plausible reconciliation supported by two
  independent measurements agreeing with each other and disagreeing with the
  issue, not a confirmed trace of the original "one rerun" that produced
  481.27 — stated with that caveat so a future reader does not chase a 481 Å
  bond that neither measurement here reproduces.
