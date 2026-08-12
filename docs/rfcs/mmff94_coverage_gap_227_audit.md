# MMFF94 strict coverage gap — issue #227 Phase 1A audit

Audit-only. No production algorithm code changed. Baseline: `main` @
`382a5f33791646b758e1f9d4fd4859b91211d9de` (v0.10.0), reconfirmed current
at audit start — no rebase needed. Branch: `audit/mmff94-coverage-227`.

## TL;DR

`ForceFieldPolicy::Mmff94BondAngleStrict` fails 221/265 (83.4%) of the Wave 1
corpus today (fresh remeasurement — the historical issue #227 number, 216,
was measured on an earlier `main` and is not reused). The root cause is
singular and well-characterized: **`assign_mmff94_numeric_types`'s aromatic
atom typing (`aromatic_c_type` / `aromatic_n_type` in
`crates/chematic-ff/src/mmff94_numeric.rs`) assigns numeric MMFF94 type IDs
that do not match the real Halgren/RDKit-verified numbering**, while
chematic's own parameter tables (`MMFF94_ANGLE_ENERGY` /
`MMFF94_BOND_ENERGY` / `MMFF94_TORSION_ENERGY` / `MMFF94_STBN`) are
correctly populated under the *real* numbering. This single root cause
explains **98.5%** of all missing terms across every energy-term kind
(bond/angle/torsion/stretch-bend/out-of-plane) by volume.

That is a real, singular root cause — but it does **not** clear the bar for
a same-session Phase 1B production fix. A concrete negative-result
experiment (below) shows a partial, context-blind correction resolves only
7.7% of molecule-level gate failures and causes a real regression. The fix
needs comprehensive, ring-context-aware re-derivation of C+N+O+S aromatic
typing together, verified against RDKit **energy parity** (not just
coverage `Some`/`None`), because some currently-"successful" molecules are
silently computing energies from semantically wrong-but-present parameter
rows (furan, demonstrated below). That is out of scope for this session per
the program's own stop condition ("multiple independent root causes /
correctness verification not yet available → audit-only Ready PR, stop").

## Reproduction

```bash
# Faithful remeasurement (exact production entry point, same methodology as issue #227)
cargo run --release -p chematic-3d --example mmff94_strict_gate_remeasure_227

# Term-level per-molecule JSONL dump (bond/angle/torsion/oop/stretch-bend/vdW/charge)
cargo run --release -p chematic-3d --example mmff94_term_coverage_audit \
  > validation/results/mmff94_coverage_227_term_audit.jsonl \
  2> validation/results/mmff94_coverage_227_stderr.log

# chematic's own numeric MMFF94 types, per atom, whole corpus
cargo run --release -p chematic-3d --example mmff94_numeric_type_dump \
  > validation/results/mmff94_chematic_numeric_types.jsonl

# RDKit oracle: real MMFF94 atom types + energy/FF construction, whole corpus
.venv/bin/python scripts/mmff94_rdkit_type_oracle.py \
  > validation/results/mmff94_rdkit_type_oracle.jsonl

# Negative-result probe: does a naive constant-swap fix work? (it doesn't — see below)
cargo run --release -p chematic-3d --example mmff94_typing_fix_simulation
```

Environment: rustc 1.97.0 (2026-07-07), cargo 1.97.0, RDKit 2026.03.3.
Corpus: same Wave 1 265-molecule corpus as PR #226/issue #227 — Tier A
sha256 `6a478ea0f5d4ef067a4d1739e77a7209e8f76ecaa837e3f487c723dd6f465d6b`,
Tier B sha256 `b3cde3fedcc68391ba3d3cbae228acd3057cadea6e6fc17499592b23bdc7550a`
— **not** a new corpus.

## Fresh baseline (today's main, not the historical 216)

| Metric | Count |
|---|---|
| Total corpus | 265 |
| `Mmff94BondAngleStrict` success (Ok) | 44 |
| `Err(MissingParameters)` ("unsupported") | 221 |
| `Err(UnsupportedAtomType)` | 0 |
| `Err(MinimizationFailed)` ("typed_failure") | 0 |

The historical issue text (38 success / 216 unsupported / 11 typed_failure)
was measured on an earlier main; today's main has moved (PR #226, #228,
#229, #230, #231, #234 etc. merged since). Both measurements agree on the
*shape* of the problem (>80% unsupported); the exact counts differ because
the program's own rule is "regenerate fresh, never reuse historical
numbers."

## Term-level missing-parameter counts (independent bond+angle recheck)

| Term kind | Total evaluated | Missing | Gates `Mmff94BondAngleStrict`? |
|---|---|---|---|
| Bond | — | 808 | Yes |
| Angle | — | 5,782 | Yes |
| Stretch-Bend | — | 6,900 | **No — never checked at all** |
| Torsion | — | 7,322 | No (measured/reported, not gated) |
| Out-of-plane | — | 67 | No (measured/reported, not gated) |
| van der Waals | — | 0 | No |
| Charge (whole-molecule) | — | 0 | No |

221/265 molecules fail because of Bond and/or Angle misses (the only two
term kinds `Mmff94BondAngleStrict`'s gate actually checks — confirmed by
reading `run_mmff94_bridge`/`compute_mmff94_coverage` in
`crates/chematic-3d/src/minimize.rs`). vdW and charge coverage is 100%
clean for this corpus — not a contributor.

## Root cause: atom typing wrong, not a missing-parameter data gap

### The benzene contradiction, resolved

Issue #227 flagged a real contradiction: `mmff94_angle_params(C_Aromatic,
C_Aromatic, C_Aromatic)` passes its own unit test, yet real benzene fails
the pipeline's coverage check for the same physical angle. Traced across
all 7 requested stages:

| Stage | What happens | Type/key involved |
|---|---|---|
| 1. Parser | `c1ccccc1` parses, all 6 atoms get `atom.aromatic = true` | — |
| 2. Aromaticity perception | Already aromatic from the SMILES lowercase atoms; no Kekulé pass needed | — |
| 3. Atom typing | `assign_mmff94_numeric_types` → `aromatic_c_type` sees a 6-ring, returns **63** | numeric type 63 |
| 4. Coverage report (`compute_mmff94_coverage`) | For each ring angle: `angle_type_for(...)` computes **2** (both flanking bonds classify `bond_type=1`, `bt_sum=2`, no 3/4-ring); looks up `mmff94_angle_energy(2, 63, 63, 63)` → `None`; falls back to `mmff94_angle_energy(0, 63, 63, 63)` (type≠0 fallback) → **also `None`** → reported missing | key `(2, 63, 63, 63)`, fallback key `(0, 63, 63, 63)` |
| 5. FF construction (`run_mmff94_bridge`) | `coverage.has_gate_failure()` is true → refuses with `Err(MissingParameters(...))` before any minimizer is built | — |
| 6. Energy evaluation | Never reached (gate refused first) | — |
| 7. Gradient evaluation | Never reached | — |

**Two structurally separate MMFF94 angle implementations exist in this
codebase** (same pattern CLAUDE.md already documents for ECFP4):

- `mmff94_angle_params(MMFF94Type, MMFF94Type, MMFF94Type)` —
  enum-keyed, a small hand-written `match` in `chematic-ff/src/mmff94_params.rs`
  that happens to have a hardcoded `(C_Aromatic, C_Aromatic, C_Aromatic)`
  arm. It has a real caller (`crates/chematic-3d/src/minimize.rs:821`,
  `angle_energy_mmff94` → `minimize_mmff94_with_config`), a legacy/simpler
  gradient-descent path — **not** used by `pipeline_v2`'s
  `Mmff94BondAngleStrict`.
- `mmff94_angle_energy(u8, u8, u8, u8)` — numeric-type-keyed, backed by the
  2,245-row `MMFF94_ANGLE_ENERGY` table in
  `chematic-ff/src/mmff94_energy/angle.rs`. This is what
  `compute_mmff94_coverage` and the real MMFF94 minimizer
  (`chematic_ff::mmff94_minimizer`) actually call.

The passing unit test (`test_angle_params_aromatic`) exercises the first
function. Real benzene, through `pipeline_v2`, goes through the second. The
test was never covering the failing path — it is not "wrong," it simply
tests different, unrelated code.

### Why the angle lookup misses: wrong numeric type, not missing data

Checked directly: does `MMFF94_ANGLE_ENERGY` have a row for `(*, 63, 63,
63)` at *any* angle type 0–8? **No** — confirmed by exhaustive scan. Table
sortedness/dedup also confirmed clean (`sort -c`, no duplicate keys) — this
is not a binary-search-on-unsorted-data bug either.

But cross-checked against an RDKit oracle (`AllChem.MMFFGetMoleculeProperties`,
RDKit 2026.03.3 — a real, independent, standards-compliant MMFF94
implementation) on the same 265-molecule corpus:

| Molecule | Atom | chematic numeric type | RDKit numeric type |
|---|---|---|---|
| benzene | ring C (×6) | 63 | **37** |
| pyridine | ring C | 63 | 37 |
| pyridine | ring N | (chematic's own scheme differs too) | 38 |
| furan | alpha C (adjacent to O) | 37 | **63** |
| furan | beta C | 38 | **64** |
| furan | O | 6 (plain ether O) | **59** |

`chematic`'s own table *does* have a correctly-populated row at the real
numbering: `(0, 37, 37, 37, ka=0.669, theta0=119.977°)` — chemically exactly
benzene's ring angle (θ≈120°). **The parameter table is right. The atom
typer assigns the wrong ID.** `aromatic_c_type` in `mmff94_numeric.rs`
(lines 813–837) literally returns `63` for 6-ring aromatic carbon and
`37`/`38` for 5-ring alpha/beta — the reverse of the real Halgren/RDKit
scheme.

### Scope: this explains the overwhelming majority of missing terms

Cross-referencing every missing term (bond/angle/torsion/stretch-bend/oop,
20,879 rows total across the corpus) against the RDKit oracle's per-atom
ground truth, classified exclusively:

| Root-cause bucket | Count | % |
|---|---:|---:|
| **`atom_typing_wrong`** (≥1 involved atom's chematic numeric type ≠ RDKit's) | 20,559 | 98.5% |
| `other_classified` — torsion-type classification bug (correctly-typed atoms, row exists at a different `torsion_type`) | 119 | 0.6% |
| `other_classified` — stretch-bend never gated at all (see below) | 116 | 0.6% |
| `other_classified` — stretch-bend never gated + also a routing candidate | 65 | 0.3% |
| `parameter_table_genuinely_absent` (correctly-typed atoms, no row at any classification code) | 20 | 0.1% |
| **unclassified** | **0** | **0%** |

98.5% one root cause; the residual 1.5% is real but small and belongs to
distinct, separately-scoped follow-ups (below) — not folded into the
primary root cause.

## Why this does NOT clear the same-session Phase 1B bar

The program's Phase 1A→1B gate requires, among other things, that the fix
be safe, existing-parameter-only, and have a *measured* blast radius. All
three were tested directly, not assumed:

**Experiment**: `mmff94_typing_fix_simulation.rs` (diagnostic-only, marked
known-wrong in its own header — not a proposed fix) applies the simplest
possible correction — a context-blind 3-constant swap, carbon only (`63↔37`,
`38→64`), no ring-size/element-together awareness — and re-runs the
bond+angle gate.

| | Count |
|---|---|
| Molecules failing before | 221 |
| Molecules failing after | 205 |
| **Flipped fail → pass** | **17 (7.7% of failures)** |
| **Regressions (pass → fail)** | **1 — `furan`** |

7.7% is nowhere near the ≥80% same-session bar, and the `furan` regression
is a real, demonstrated violation of the merge gate's "previously
successful molecules: 0 regressions" requirement — from a *partial*
correction, which is the most anyone could safely attempt without a full
C+N+O+S re-derivation.

### The `furan` regression is itself the headline finding, not a footnote

`furan` passes `Mmff94BondAngleStrict` **today**, with chematic's current
(wrong) types `[C:38, C:38, C:37, O:6, C:37]` (RDKit's real types:
`[64, 64, 63, 59, 63]`). Directly querying the real energy functions with
furan's actual (wrong) types:

```
bond 0-1 (C-C, both type 38) -> Some(BondEnergyParams { kb: 5.002, r0: 1.246 })
```

Type 38 in the real Halgren numbering is pyridine-type nitrogen, not
carbon — `r0 = 1.246 Å` is not a plausible furan aromatic C–C bond length
(real furan C–C ≈ 1.36–1.44 Å; this looks like an N-involved table row).
**`furan` "succeeds" today by silently landing on a semantically wrong but
numerically present parameter row.** This means the true scope of issue
#227 is larger than "216 (or 221) unsupported": an unknown subset of the 44
currently-"successful" molecules may also be computing physically wrong
MMFF94 energies, never caught because `compute_mmff94_coverage` only checks
`Some`/`None`, never whether the returned row is the *semantically correct*
one for that atom's real chemistry. A correctness gate (RDKit energy
parity), not just a coverage-count gate, is required for any future fix PR
to be trustworthy.

## Residual buckets — explicitly out of scope for this audit

Named individually per the program's "no unclassified" rule, not folded
into the primary root cause:

- **Stretch-bend is never gated at all.** `compute_mmff94_coverage` (the
  function backing `Mmff94BondAngleStrict`'s refusal check) enumerates
  bond/angle/torsion/oop only — it has no stretch-bend loop.
  `mmff94_energy_breakdown`'s `stretch_bend_energy` silently contributes
  `0.0` for any missing term (`if let Some(...) = mmff94_stbn(...) { ... }`,
  no `else` branch) — never a typed failure, never counted in the 216/221.
  6,900 missing stretch-bend terms exist in this corpus today, invisible to
  every current metric.
- **Torsion classification candidate bug** (119/7,322 = 1.6% of torsion
  misses, on *correctly-typed* atoms, where a row exists at a different
  `torsion_type` than `torsion_type_for` computes). Real, but: (a) small
  relative to the corpus, (b) torsion coverage doesn't gate
  `Mmff94BondAngleStrict` at all, so fixing it would not move the 221
  number. Independent follow-up candidate, not evidence for the primary
  root cause.
- **Bond lookup has zero fallback mechanism** — unlike angle (type→0
  fallback), torsion (7-tier fallback chain), oop (6 orderings + 4
  wildcard tiers), and stretch-bend (type→0 + generic wildcard),
  `mmff94_bond_energy` is a single direct `binary_search` with no fallback
  at all. An inconsistency across term kinds worth noting for any future
  bond-table work, not itself the dominant driver (bonds are only 808/20,879
  = 3.9% of missing-term volume).
- **Aromatic oxygen/sulfur typing** (furan O: chematic assigns plain-ether
  type 6, RDKit's real type is 59) and the messier nitrogen mismatches
  (chematic type 67 maps to *six different* RDKit types — 38/10/9/58/40/54
  — depending on context; chematic type 8 maps to RDKit type 40 in 113/193
  cases) are real but are genuine **classification-logic** gaps (chematic's
  nitrogen/oxygen typer conflates distinct chemical environments the real
  spec differentiates), not simple constant relabeling like the carbon
  case. These need their own careful, separately-reviewed re-derivation.
- **vdW and charge coverage**: 0 missing across the whole corpus — clean,
  no action needed.

## Recommendation

**Do not attempt the fix in this session.** File a separate, dedicated
follow-up issue/PR scoped as:

1. Re-derive the full aromatic C+N+O+S numeric MMFF94 type assignment
   against a citable primary source — RDKit's own
   `Code/ForceField/MMFF/Params.cpp` is the natural choice (chematic already
   cites this exact file for the PBCI table in `mmff94_numeric.rs:38`,
   establishing precedent for this provenance pattern) — not a guess, not
   back-solved from RDKit's *output* on this corpus.
2. Fix carbon, nitrogen, oxygen, and sulfur aromatic typing **together in
   one PR** — the `furan` regression proves a carbon-only partial fix
   creates new inconsistent states.
3. Gate on RDKit **energy parity** for the currently-passing molecules, not
   just coverage `Some`/`None` — the `furan` finding shows coverage success
   today does not imply energetic correctness.
4. Expect the "44 successful" baseline to legitimately shift (up or down)
   once typing is corrected — restate the merge gate in terms of
   correctness, not raw pass-count, before that PR is judged.
5. Separately evaluate whether to add a stretch-bend coverage check to
   `compute_mmff94_coverage` (a distinct, smaller decision — currently by
   design out of `Mmff94BondAngleStrict`'s stated bond+angle-only scope,
   but the *silent* zero-contribution behavior for uncovered terms is worth
   an explicit decision either way).

## Artifacts

- `crates/chematic-3d/examples/mmff94_term_coverage_audit.rs` — term-level
  JSONL dump tool
- `crates/chematic-3d/examples/mmff94_numeric_type_dump.rs` — chematic's own
  per-atom numeric types, whole corpus
- `crates/chematic-3d/examples/mmff94_strict_gate_remeasure_227.rs` —
  faithful production-API remeasurement
- `crates/chematic-3d/examples/mmff94_typing_fix_simulation.rs` —
  diagnostic-only negative-result probe (explicitly marked not-a-fix)
- `scripts/mmff94_rdkit_type_oracle.py` — RDKit oracle
- `validation/results/mmff94_coverage_227_term_audit.jsonl` — 20,879
  missing-term rows + 265 molecule-summary rows
- `validation/results/mmff94_chematic_numeric_types.jsonl` — chematic's
  numeric types per atom, whole corpus
- `validation/results/mmff94_rdkit_type_oracle.jsonl` — RDKit oracle output,
  whole corpus
- `validation/results/mmff94_coverage_227_root_cause_classification.json` —
  exclusive classification counts
- `validation/results/mmff94_coverage_227_aggregate.json` — consolidated
  summary (this document's numbers, machine-readable)
