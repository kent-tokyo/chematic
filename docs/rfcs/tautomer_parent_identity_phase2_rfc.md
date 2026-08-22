# RFC: Tautomer & Parent Identity (ROADMAP.md Phase 2, v0.20.0)

Status: draft (RFC + acceptance fixtures only; no production code changes this round)

## 1. Why now — audit of existing code

`standardize.rs` and `tautomer.rs` were never audited end-to-end before this
RFC, the same situation Phase 1 found in the fragment-selection code. Both
modules already implement most of the primitives Phase 2's spec names as "new"
scope. The real work is auditing what exists, not building from scratch. Five
concrete findings, all empirically verified against `main` (`99401a6`):

### 1.1 Confirmed defect: canonical tautomer is not invariant across aromatic
lactam/lactim pairs

`canonical_tautomer` is supposed to return the same output regardless of which
tautomer was fed in. Tested pairs of the *same* real molecule (formula-checked
via RDKit 2026.03.4):

| pair | input A → output | input B → output | self-consistent? |
|---|---|---|---|
| 2-pyridone | `O=c1cccc[nH]1` → `[nH]1ccccc1=O` | `Oc1ccccn1` → `c1nc(O)ccc1` | **no** |
| 4-pyridone | `O=c1cc[nH]cc1` → `[nH]1ccc(=O)cc1` | `Oc1ccncc1` → `c1c(O)ccnc1` | **no** |
| cytosine amino/imino | `Nc1cc[nH]c(=O)n1` → unchanged | `Nc1ccnc(O)n1` → unchanged | **no** |
| uracil | `O=c1cc[nH]c(=O)[nH]1` → unchanged | `Oc1ccnc(O)n1` → unchanged | **no** |
| guanine-like purine | `Nc1nc2[nH]cnc2c(=O)[nH]1` → unchanged | `Nc1nc2[nH]cnc2c(O)n1` → unchanged | **no** |

By contrast, every non-aromatic-exocyclic case tested is already
self-consistent: keto-enol (acetylacetone), amide-iminol (acetamide),
guanidine, nitroso-oxime, and ring-*internal* NH shifts (imidazole, pyrazole,
tetrazole) all converge to one output regardless of input form.

**Root cause, confirmed by inspecting bond orders after parse:** ring bonds in
an aromatic system are stored as `BondOrder::Aromatic`, never Kekulized. Every
`TautomerRule`'s `donor_bridge_order`/`bridge_acceptor_order` is matched via
`BondOrderMatch::{Single,Double}`, and `BondOrderMatch::Double::matches` only
accepts `BondOrder::Double` — never `Aromatic`. So no 1,3- or 1,5-shift rule
can ever fire across an aromatic ring bond. The separate mechanism that *does*
handle aromatic systems, `find_direct_aromatic_matches` /
`transfer_hydrogen_aromatic`, only moves the mobile H between ring atoms of
one aromatic system (e.g. pyrazole N1↔N3) — it never considers an exocyclic,
singly-bonded substituent (the lactam/lactim `=O`/`-OH`) as a place the mobile
H can come from or go to. This is a real, well-scoped gap, not a fixture
artifact: it silently affects the canonical identity of pyridone-, cytosine-,
uracil-, and purine-class heterocycles — precisely the "heterocycle" and
"drug-like" corpora ROADMAP.md names for this phase, and the textbook test
case for RDKit's own `TautomerEnumerator`.

This is squarely Phase 2's headline deliverable, not a tangential nice-to-have
— the design in §4.4 below targets it directly, but implementation is
deferred to the next round (see §7).

### 1.2 Confirmed non-bug: no repeat of Phase 1's stereocenter-corruption class

Phase 1's biggest finding was fragment extraction silently dropping
`stereo_neighbor_order`/`bond_directions`/`stereo_groups` by rebuilding via a
fresh `MoleculeBuilder` without copying those side tables. `tautomer.rs`'s
core transforms (`transfer_hydrogen`, `transfer_hydrogen_aromatic`) also
rebuild via `MoleculeBuilder::new()`, so the same class of bug was plausible
here. It already isn't: both functions explicitly call
`copy_stereo_groups_from`, `copy_stereo_from`, and `copy_bond_directions_from`
after rebuilding, and this is pinned by existing regression tests
(`transfer_hydrogen_aromatic_preserves_remote_stereo_neighbor_order`,
`transfer_hydrogen_preserves_stereo_groups_and_bond_directions`, and others in
the same file). Verified directly: atom order and count are 1:1 preserved by
both functions (only `hydrogen_count` and specific bond orders change), which
is exactly the precondition those `copy_*_from` calls need to be correct. No
fix needed here — recorded so this round doesn't manufacture a finding where
none exists.

### 1.3 Confirmed gap: budget exhaustion is silent

`canonical_tautomer_with_config` runs at most `config.max_iter` outer
iterations (default 16) and simply breaks when the budget runs out or no rule
fires — there is no signal to the caller that the result may be
input-order-dependent because the budget was hit rather than because the
molecule converged. `enumerate_tautomers_with_config` behaves the same way
for `max_tautomers` (default 32). This is not hypothetical: the file already
carries a permanent, `#[ignore]`d regression test,
`test_max_iter_default_diverges_on_many_independent_sites`, which constructs a
25-independent-site "comb" molecule and proves the default `max_iter=16`
produces different canonical output depending on arm-build order — with a
companion test, `test_max_iter_1000_resolves_the_divergence`, proving this is
pure budget exhaustion (raising the limit fixes it), not a deeper algorithm
bug. No real-molecule instance of this has been found, but the mechanism is
real and the API gives a caller no way to detect it happened. This is exactly
the gap `StandardizationLimits` + a result-state enum (§4.1–4.2) exists to
close.

### 1.4 Confirmed gap: 5 `StandardizationStep` variants are unwired, but their
bodies mostly already exist

`StandardizationStep` has 10 variants; `StandardizationPipeline::run()` only
invokes 5 (`NeutralizeCharges`, `LargestFragment`, `ZwitterionNormalization`,
`RemoveExplicitHydrogens`, `CanonicalTautomer`). The other 5 —
`NormalizeGroups`, `FragmentParent`, `ChargeParent`, `IsotopeParent`,
`StereoParent` — exist only as enum cases plus one naming unit test
(`parent_variant_step_names_distinct`, which asserts `.as_str()` strings and
nothing else). Checking what stands behind each name:

| variant | wired in `run()`? | underlying function | function status |
|---|---|---|---|
| `NormalizeGroups` | no | `normalize_groups` (nitro/azide/sulfoxide) | exists, tested, just not called from `run()` |
| `FragmentParent` | no (but `LargestFragment` covers the same job under a different name) | `select_fragment` / `largest_fragment` (Phase 1) | exists, tested |
| `ChargeParent` | no (but `NeutralizeCharges` covers the same job) | `neutralize_charges` | exists, tested |
| `IsotopeParent` | no | `remove_isotopes` | exists, tested, unused outside its own tests |
| `StereoParent` | no | `remove_stereo` | exists, tested, unused outside its own tests |

Mechanically re-verified this round: `neutralize_charges`, `remove_isotopes`,
`remove_stereo`, and `largest_fragment` all still produce correct output on
representative inputs (ammonium acetate → acetic acid + ammonia; a
deuterated/¹³C-labeled molecule → unlabeled; alanine → destereoed; an HCl
amine salt → the free amine). So the gap isn't "these transforms don't exist,"
it's "there is no `Parent`-shaped concept" — RDKit's `Parent` functions are a
distinct idea from a mutating pipeline stage: an idempotent, order-independent
reduction along *one axis* of variability, meant to be used as a grouping/dedup
key, computable independently of (and without mutating) any other axis. §4.3
proposes exposing that concept directly on top of the existing primitives,
rather than wiring these 5 names into the mutating `run()` pipeline as more
stages (which would just make `run()` non-idempotent along more axes at once,
not give callers an addressable per-axis identity).

### 1.5 Confirmed gap: no tautomer rule/scoring extensibility

ROADMAP.md's implementation list includes "custom normalization rule / custom
tautomer scoring." Today `TautomerRule` (the 15-entry rule struct) and
`tautomer_score` (the O-H>N-H>S-H + aromaticity scoring function used to pick
among direct-aromatic-shift candidates) are both private, and
`TautomerConfig::enabled_rules` can only toggle *which* of the 15 built-in
rules run, by index — it cannot add a new rule or override the scoring
function. This is genuinely unimplemented scope, not merely unwired.

## 2. Goals / non-goals

**Goals this phase:**
- Make `canonical_tautomer` invariant across input tautomer spelling for the
  aromatic lactam/lactim class (§1.1), by generalizing the existing aromatic
  H-shift mechanism rather than adding a parallel one.
- Give callers a way to know when a result was budget-limited
  (`StandardizationLimits` + a result-state enum), instead of a silently
  possibly-order-dependent answer.
- Expose `FragmentParent`/`ChargeParent`/`IsotopeParent`/`StereoParent`/
  `TautomerParent`/`SuperParent` as an explicit, addressable "Parent identity"
  concept built on the existing primitives, each returning a lightweight audit
  record (mirroring Phase 1's `TransformationRecord`).
- A minimal custom-scoring hook for tautomer selection (exact shape decided in
  the implementation round — see §4.5).
- Rust/Python/WASM parity is explicitly **deferred** until the Rust core here
  is stable, mirroring Phase 1's own sequencing (ROADMAP.md's "直近の開始順").

**Non-goals this phase (explicitly out of scope):**
- 3D conformer generation, MMFF94, symmetrized SSSR, canonical-SMILES residual
  work (issue #149) — untouched, unrelated axes.
- A general SMARTS-based or fully user-authorable tautomer rule *language* —
  §4.5 proposes the smallest hook that satisfies the roadmap line item, not a
  plugin system.
- Rewriting `max_iter`'s pathological-comb divergence into a *fix* — the goal
  is making exhaustion **visible**, not eliminating the fundamental
  budget/completeness tradeoff (raising `max_iter` arbitrarily high is not
  free; see §4.1).
- Any change to `StandardizationReport`/`StandardizationStepReport`/
  `MoleculeSnapshot`/`PipelineStatus` (Phase 1's public types) — Parent
  identity and the new result-state enum are additive, orthogonal types, same
  resolution Phase 1's §4.5 landed on.

## 3. Confirmed non-goals validated by audit

Sections 1.2 confirms the stereo-preservation machinery in `tautomer.rs`
needs no repair. No further audit-driven scope changes were found beyond what
§1 already lists.

## 4. Design

### 4.1 `StandardizationLimits`

```rust
pub struct StandardizationLimits {
    pub max_restarts: usize,
    pub max_transforms: usize,
    pub max_tautomers: usize,
    pub timeout_ms: Option<u64>,
}
```

`max_transforms` and `max_tautomers` are direct generalizations of
`TautomerConfig::max_iter`/`max_tautomers` (§ note: the existing
`TautomerConfig` fields are the literal precedent — this struct doesn't
reinvent them, it lifts them to the standardization-pipeline level so
`FragmentParent`/`ChargeParent`/etc. share one limits type instead of each
growing its own). `max_restarts` is new: it bounds how many times the overall
Parent/canonicalization computation may restart from scratch (relevant once a
custom scorer or rule set is pluggable and could in principle loop against
itself — see §4.5). All three are **deterministic** budgets: same input, same
limits ⇒ same result, always — this is what makes `MaxTransformsReached`/
`MaxTautomersReached` states reproducible.

`timeout_ms` is a different kind of bound and is called out explicitly rather
than silently accepted at face value: a wall-clock timeout depends on
machine speed and load, so a result state of `TimedOut` is **not**
deterministic across environments, and cannot honestly carry the "same
canonical tautomer regardless of input" guarantee ROADMAP.md's own exit
criterion demands. Recommendation: keep `timeout_ms` as an optional escape
hatch (`None` by default, deterministic budgets alone govern normal
operation), and document `TimedOut` results as explicitly outside the
determinism guarantee — a safety valve for pathological inputs in a
service context, not a knob normal callers should reach for. This is a design
decision made now, not left open, per Phase 1's own "decide, don't defer"
precedent (its §8.5).

### 4.2 Result-state enum

```rust
pub enum ParentComputationStatus {
    Completed,
    MaxTransformsReached,
    MaxTautomersReached,
    TimedOut,
    Canceled,
    Abstained(String),
    InvalidInput(String),
}
```

This is intentionally a **new, orthogonal type**, not an extension of the
existing `PipelineStatus` (`Unchanged`/`Modified`/`CompletedWithWarnings`).
`PipelineStatus` answers "did the molecule change"; `ParentComputationStatus`
answers "did the computation reach a definite answer, and if not, why not" —
different questions, and conflating them would mean every existing
`StandardizationReport` consumer suddenly has to handle
`MaxTransformsReached` even though today's pipeline can't produce it. This
mirrors how Phase 1 resolved its own public-API-compatibility question (its
§8.1): add a new type, don't extend a shipped one. `Abstained`/`InvalidInput`
carry a reason string, matching `TransformationRecord::abstained`'s existing
`Option<String>` shape from Phase 1 rather than inventing a second
error-string convention.

### 4.3 Parent-generation functions

A **Parent** is an idempotent, deterministic reduction of one axis of
molecular variability, meant to be used as a grouping/dedup key — computing
it never depends on, or mutates, any other axis. This is the concept RDKit's
`MolStandardize` names `FragmentParent`/`ChargeParent`/`IsotopeParent`/
`StereoParent`/`TautomerParent`/`SuperParent`. Proposed public functions,
each thin wrapper over an existing, already-verified primitive (§1.4's table)
plus a lightweight audit record in the same shape as Phase 1's
`TransformationRecord`:

```rust
pub fn fragment_parent(mol: &Molecule) -> (Molecule, TransformationRecord); // = select_fragment
pub fn charge_parent(mol: &Molecule) -> (Molecule, TransformationRecord);   // = neutralize_charges + audit
pub fn isotope_parent(mol: &Molecule) -> (Molecule, TransformationRecord);  // = remove_isotopes + audit
pub fn stereo_parent(mol: &Molecule) -> (Molecule, TransformationRecord);   // = remove_stereo + audit
pub fn tautomer_parent(mol: &Molecule, limits: &StandardizationLimits)
    -> (Molecule, ParentComputationStatus, TautomerAuditRecord);
pub fn super_parent(mol: &Molecule, limits: &StandardizationLimits)
    -> (Molecule, ParentComputationStatus, Vec<TransformationRecord>); // composes all of the above, in a fixed order
```

`fragment_parent`/`charge_parent`/`isotope_parent`/`stereo_parent` return
`TransformationRecord` unchanged from Phase 1 (no new fields needed — these
are simple mechanical transforms, not multi-candidate selections).
`tautomer_parent` needs its own record shape, since "which tautomer won and
why" has no Phase 1 analog:

```rust
pub struct TautomerAuditRecord {
    pub selected: MoleculeSnapshot,
    pub candidate_count: usize,
    pub score_breakdown: Vec<(String, i32)>, // human-readable rule/criterion -> contribution
    pub applied_transforms: Vec<String>,      // rule names applied, in order
    pub lost_stereo: bool,
    pub affected_isotopes: bool,
}
```

`super_parent` runs the other five in a fixed, documented order (fragment →
charge → isotope → stereo → tautomer, matching `StandardizationPipeline::run`'s
existing ordering rationale of "neutralize before fragment-select" applied
consistently) and returns every intermediate `TransformationRecord`, so a
caller can see exactly which axis contributed which change — never a single
blended "it changed" boolean, following Phase 1's own audit-log discipline.

`StandardizationStep::{NormalizeGroups, FragmentParent, ChargeParent,
IsotopeParent, StereoParent}` stay as enum cases (their `.as_str()` names are
already public API, changing them would be a breaking rename) but are
reinterpreted as **labels a caller can report which Parent function ran**,
not as new stages silently added to `StandardizationPipeline::run()`. Wiring
them into the mutating pipeline was considered and rejected: `run()` already
has a specific, deliberate stage order for producing one standardized
molecule; a `Parent` is a different molecule per axis, not one more mutation
in that same sequence.

### 4.4 Fixing the aromatic lactam/lactim gap (design only, not implemented
this round)

Generalize the existing `find_direct_aromatic_matches` /
`transfer_hydrogen_aromatic` mechanism (which already treats "the mobile H
can be at any of several ring positions in this aromatic system") to also
admit one more kind of position: an exocyclic atom that is singly bonded to a
ring atom and currently holds the mobile H (or, symmetrically, a ring atom
that could receive it), toggling that one exocyclic bond's order between
`Single` and `Double` as the H arrives or leaves — without touching any
in-ring bond order, since those are correctly left as `Aromatic` throughout
(confirmed: converting `Oc1ccccn1` to `O=c1cccc[nH]1` changes exactly one
bond, the exocyclic C–O bond, from `Single` to `Double`; every ring bond,
including the one adjacent to the nitrogen that gains the H, stays
`Aromatic`). This is a natural generalization of the existing
ring-internal-shift code path, not a new algorithm — but it is real,
non-trivial work (identifying the exocyclic-substituent candidates correctly,
without misfiring on substituents that aren't part of a lactam/lactim system,
e.g. a plain phenol), so it's scoped as the implementation round's primary
deliverable rather than done inline in this RFC.

### 4.5 Custom rule/scoring hook

Recommendation: make `TautomerRule` and `tautomer_score` public
(`#[non_exhaustive]` on `TautomerRule` to keep future field additions
non-breaking), and add two new `TautomerConfig` fields:

```rust
pub extra_rules: Vec<TautomerRule>,
pub scorer: Option<fn(&Molecule) -> i32>,
```

This is the smallest hook that satisfies "custom normalization rule / custom
tautomer scoring" — it lets a caller add rules or override scoring without a
trait, registry, or plugin system. Deferred to the implementation round for
exact field types (a `fn` pointer vs. a boxed closure trade-off depends on
whether WASM/Python bindings need to pass a scorer across the FFI boundary,
which isn't yet decided — see §6).

## 5. Acceptance fixtures

Two files, mirroring Phase 1's split: a main design-driving set and a
held-out generalization set never used to shape the rules above.

- `validation/tautomer_parent_identity_phase2_fixtures.jsonl`
- `validation/tautomer_parent_identity_phase2_holdout.jsonl`

Two row shapes, since "tautomer identity" and "Parent generation" are
verified differently — a tautomer fixture asserts **self-consistency**
(`canonical_tautomer` of every listed variant must produce one identical
canonical SMILES; which specific tautomer wins is chematic's own choice,
informed by its own scoring, and is *not* required to match RDKit's
preferred form — this matches ROADMAP.md's own definition, "consistent
representation independent of input tautomer... not thermodynamic
stability"). A parent fixture asserts a **fixed expected output**, since
`charge_parent`/`isotope_parent`/`stereo_parent`/`fragment_parent` are
simple, unambiguous mechanical transforms with exactly one correct answer.

```json
{"id": "tp2-05-2-pyridone", "kind": "tautomer_self_consistency", "category": "heterocycle_lactam_lactim", "variants": ["O=c1cccc[nH]1", "Oc1ccccn1"], "currently_passing": false, "note": "confirmed failing on main (99401a6): outputs differ ([nH]1ccccc1=O vs c1nc(O)ccc1). Root cause and fix design in RFC section 1.1/4.4."}
{"id": "tp2-14-charge-parent-ammonium-acetate", "kind": "parent_generation", "parent": "charge_parent", "input": "CC(=O)[O-].[NH4+]", "expected_output": "CC(O)=O.N", "note": "mechanical, already correct via neutralize_charges (RFC section 1.4)."}
```

Categories covered in the main set (24 rows): keto-enol, amide-iminol,
guanidine/oxime (non-aromatic, already-passing controls), ring-internal NH
shift (imidazole/pyrazole/tetrazole, already-passing controls), aromatic
lactam/lactim (2-pyridone/4-pyridone/cytosine/uracil/guanine-like, the
confirmed-failing class), a zwitterion-interaction check (does the full
`StandardizationPipeline` — which runs `ZwitterionNormalization` before
`CanonicalTautomer` — converge to one output regardless of whether the input
was drawn as a zwitterion or a neutral form), a metal-adjacent check (does
`disconnect_metals` + `LargestFragment` running before tautomer
canonicalization correctly isolate the organic ligand's tautomer identity
regardless of how the original complex was drawn), and all five Parent
functions (`charge_parent`/`isotope_parent`/`stereo_parent`/
`fragment_parent`/`super_parent`).

The holdout set (5 rows) generalizes each category without ever having its
expected answer used to shape a rule: a sixth lactam/lactim pair
(hypoxanthine) beyond the five studied above, a fused-ring generalization of
the ring-internal-shift control (benzimidazole), a longer-chain zwitterion
generalization (GABA), a combined negative control for all four mechanical
Parent functions on one already-reduced input (toluene), and a companion
citation of the existing `max_iter=1000` regression test proving the
limit-exhaustion divergence is pure budget exhaustion, not a masked second
bug.

## 6. Open questions

1. **`extra_rules`/`scorer` FFI shape** — a `fn` pointer is trivial in Rust
   but doesn't cross the Python/WASM boundary; a boxed `dyn Fn` closure does
   but adds an indirection cost on a hot path. Since bindings are explicitly
   deferred (§2), this is left for the round that actually adds them, not
   decided speculatively now.
2. **Does the aromatic-lactam/lactim fix ever need to *dearomatize* a ring
   bond**, or is toggling only the exocyclic bond always sufficient? All five
   confirmed cases only required the exocyclic bond to change — but this
   needs re-checking against a larger corpus during implementation, not
   assumed to generalize from 5 examples.
3. **Should `max_restarts` default to 1** (no restarts, matching today's
   effective behavior) or something larger? No molecule found this round
   requires more than one pass; default conservatively to 1 unless the
   implementation round finds a real case that needs more.

## 7. What ships this round

RFC (this document) + the two fixture files (§5) + `ROADMAP.md`/
`validation/README.md` pointers. **No changes under `crates/*/src/**` this
round** — matching Phase 1's own RFC-round discipline (`standardize.rs`,
`tautomer.rs` are read-only inputs to this document, not touched). Opened as
a draft PR; implementation (§4.1–4.5, plus the §4.4 fix) proceeds only after
separate explicit authorization, the same rhythm Phase 1 followed
(RFC/fixtures round → "次roundへ進んで" → Rust core round).
