# RFC: Tautomer & Parent Identity (ROADMAP.md Phase 2, v0.20.0)

Status: draft, revision round 2A (RFC + acceptance fixtures only; no
production code changes this round)

> **Revision (2026-08-22, round 2A):** the original draft defined
> `charge_parent` as a bare wrapper over `neutralize_charges` (leaving
> multiple fragments in the output), left `super_parent`'s stage order
> self-contradictory, left `StandardizationLimits`'s counters undefined,
> returned untyped tuples/`String` reasons from the Parent API, and shipped
> one broken fixture (`tp2-04`, an accidental same-structure "pair"). All
> five are fixed below, and a sixth finding — a second, mechanistically
> distinct tautomer non-convergence defect (nitroso/oxime, §1.6) — surfaced
> while fixing that fixture. §4 and §5 are substantially rewritten; §1's
> audit findings are additive (nothing in the original audit was wrong, one
> more finding is added). See the end of §7 for the round split this
> revision also formalizes (2A/2B/2C/2D).

## 1. Why now — audit of existing code

`standardize.rs` and `tautomer.rs` were never audited end-to-end before this
RFC, the same situation Phase 1 found in the fragment-selection code. Both
modules already implement most of the primitives Phase 2's spec names as "new"
scope. The real work is auditing what exists, not building from scratch. Six
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
| hypoxanthine (holdout) | `O=c1[nH]cnc2[nH]cnc12` → unchanged | `Oc1ncnc2[nH]cnc12` → unchanged | **no** |

By contrast, every non-aromatic-exocyclic case tested is already
self-consistent: keto-enol (acetylacetone), amide-iminol (acetamide),
guanidine, and ring-*internal* NH shifts (imidazole, pyrazole, tetrazole,
benzimidazole) all converge to one output regardless of input form. (Nitroso/
oxime is *not* in this passing list — see the new §1.6 below.)

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
uracil-, purine-, and hypoxanthine-class heterocycles — precisely the
"heterocycle" and "drug-like" corpora ROADMAP.md names for this phase, and the
textbook test case for RDKit's own `TautomerEnumerator`.

This is squarely Phase 2's headline deliverable, not a tangential nice-to-have
— the design in §4.4 below targets it directly, but implementation is
deferred to round 2C (see §7).

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
the gap `TautomerLimits` + a result-state enum (§4.1–4.2) exists to close.

**Counting semantics, precisely defined (re-reading the existing loop code
directly, not assumed):**
- `max_iter`/`max_transforms` counts **outer loop iterations in which a
  transform was actually applied** — `canonical_tautomer_with_config`'s loop
  body tries each forward-preferred rule's *first* match, and counts one
  iteration only when `changed` becomes `true` (i.e., an accepted, not
  previously-seen fingerprint). A rule that matches but produces an
  already-seen fingerprint does not consume budget; it simply falls through
  to the next rule in the same iteration.
- `max_tautomers` counts **distinct enumerated structures added to the result
  set** — `enumerate_tautomers_with_config`'s `result.len() < config.max_tautomers`
  check gates both the rule-based (1,3-/1,5-shift) and the direct-aromatic-shift
  branches from the *same* counter (confirmed by reading both branches' guard
  conditions in the same `while` loop) — they share one budget, not two
  separate ones. Fingerprint-duplicate candidates are rejected before
  incrementing `result.len()`, so duplicates never consume budget either.
- 1,3-shift and 1,5-shift rules are **not weighted differently** — each
  successful application counts as exactly one transform, regardless of
  `path_len`.

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
| `ChargeParent` | no | *none yet* — see §4.3, this is **not** simply `neutralize_charges` | `neutralize_charges` exists but is not, by itself, `ChargeParent` (revised understanding, see below) |
| `IsotopeParent` | no | `remove_isotopes` | exists, tested, unused outside its own tests |
| `StereoParent` | no | `remove_stereo` | exists, tested, unused outside its own tests |

Mechanically re-verified this round: `neutralize_charges`, `remove_isotopes`,
`remove_stereo`, and `largest_fragment` all still produce correct output on
representative inputs (ammonium acetate → acetic acid + ammonia [as
`neutralize_charges` alone — **not** the `ChargeParent`/`charge_parent`
result, see §4.3's revision]; a deuterated/¹³C-labeled molecule → unlabeled;
alanine → destereoed; an HCl amine salt → the free amine). So the gap isn't
"these transforms don't exist," it's "there is no `Parent`-shaped concept" —
RDKit's `Parent` functions are a distinct idea from a mutating pipeline stage:
an idempotent, order-independent reduction to *one representative structure*,
meant to be used as a grouping/dedup key. §4.3 proposes exposing that concept
directly on top of the existing primitives — with `ChargeParent` corrected
this round to mean "select the parent fragment, then neutralize it," not "
neutralize every fragment and keep them all" (§4.3's revision explains why the
original draft's `charge_parent = neutralize_charges` was wrong).

### 1.5 Confirmed gap: no tautomer rule/scoring extensibility

ROADMAP.md's implementation list includes "custom normalization rule / custom
tautomer scoring." Today `TautomerRule` (the 15-entry rule struct) and
`tautomer_score` (the O-H>N-H>S-H + aromaticity scoring function used to pick
among direct-aromatic-shift candidates) are both private, and
`TautomerConfig::enabled_rules` can only toggle *which* of the 15 built-in
rules run, by index — it cannot add a new rule or override the scoring
function. This is genuinely unimplemented scope, not merely unwired.

### 1.6 Confirmed defect (distinct mechanism): nitroso/oxime interconversion
has no forward-applied rule

Discovered while correcting a broken fixture (the original `tp2-04` paired
`CC=NO` with `CC(=NO)` — literally the same structure with an added
redundant parenthesis, not two different tautomers; see §5's correction).
The real nitroso/oxime pair — `CCN=O` (nitrosoethane) and `CC=NO` (acetaldehyde
oxime), RDKit-formula-verified as the same molecule, C2H5NO — is **also not
self-consistent**: `canonical_tautomer("CCN=O")` stays `C(C)N=O` and
`canonical_tautomer("CC=NO")` stays `ON=CC`. Root cause, confirmed by reading
the rule table directly, is different from §1.1's aromaticity issue: the one
rule shaped to handle this, `"1,3-C-to-O-any-bridge"` (`donor_elem: 6`,
`bridge_elem: None`, `acceptor_elem: 8`, forward pattern "C-H + X=O → C=X +
O-H"), is marked `prefer_forward: false`, and `canonical_tautomer`'s iteration
loop only ever applies rules where `prefer_forward == true`
(`active_rules(config).into_iter().filter(|r| r.prefer_forward)`). So this
rule never fires from either direction in the main loop. This is very likely
a deliberate trade-off rather than an oversight: because `bridge_elem: None`
matches *any* bridge atom, applying this rule forward would also fire on a
plain ketone's own alpha C-H (e.g. acetone, where bridge=C is "any"),
converting it into its enol form — the opposite of the intended "prefer keto"
behavior that the dedicated `keto-enol` rule already enforces correctly.
Disabling this generic rule's forward direction avoids that regression, at
the cost of leaving less-common heteroatom-bridge cases like nitroso/oxime
(bridge=N) unconverted. Recorded as a second, distinct confirmed defect —
not folded into §1.1's fix, since the mechanism (rule generality/specificity
trade-off, not aromatic bond-order matching) is unrelated. Scoping its fix is
left as an explicit open question (§6) rather than assumed to be solved by
§4.4's aromatic-shift work.

## 2. Goals / non-goals

**Goals this phase:**
- Make `canonical_tautomer` invariant across input tautomer spelling for the
  aromatic lactam/lactim class (§1.1), by generalizing the existing aromatic
  H-shift mechanism rather than adding a parallel one.
- Give callers a way to know when a result was budget-limited
  (`TautomerLimits` + a result-state enum), instead of a silently
  possibly-order-dependent answer.
- Expose `FragmentParent`/`ChargeParent`/`IsotopeParent`/`StereoParent`/
  `TautomerParent`/`SuperParent` as an explicit, addressable "Parent identity"
  concept built on the existing primitives, each returning a typed audit
  record (§4.3).
- A declarative, cross-language custom-scoring hook for tautomer selection
  (§4.5).
- Rust/Python/WASM parity is explicitly **deferred** until the Rust core here
  is stable, mirroring Phase 1's own sequencing (ROADMAP.md's "直近の開始順").

**Non-goals this phase (explicitly out of scope):**
- 3D conformer generation, MMFF94, symmetrized SSSR, canonical-SMILES residual
  work (issue #149) — untouched, unrelated axes.
- A general SMARTS-based or fully user-authorable tautomer rule *language* —
  §4.5 proposes a declarative scoring config, not a plugin system.
- Rewriting `max_iter`'s pathological-comb divergence into a *fix* — the goal
  is making exhaustion **visible**, not eliminating the fundamental
  budget/completeness tradeoff (raising `max_iter` arbitrarily high is not
  free; see §4.1).
- Fixing §1.6's nitroso/oxime gap this phase — flagged, not committed to a
  fix design yet (§6).
- Any change to `StandardizationReport`/`StandardizationStepReport`/
  `MoleculeSnapshot`/`PipelineStatus` (Phase 1's public types) — Parent
  identity and the new result-state enum are additive, orthogonal types, same
  resolution Phase 1's §4.5 landed on.
- A real cancellation mechanism (token/callback) — since none is designed
  this phase, `Canceled` is **removed** from the result-state enum (§4.2)
  rather than kept as a state nothing can ever produce.
- `max_restarts` — removed from this phase's limits type entirely (§4.1); it
  belongs, if anywhere, to a future Normalization-phase limits type governing
  `StandardizationPipeline`-level restart behavior, which is not designed
  here.

## 3. Confirmed non-goals validated by audit

§1.2 confirms the stereo-preservation machinery in `tautomer.rs` needs no
repair.

## 4. Design

### 4.1 `TautomerLimits`

```rust
#[non_exhaustive]
pub struct TautomerLimits {
    pub max_transforms: usize,
    pub max_tautomers: usize,
    pub timeout_ms: Option<u64>,
}
```

Renamed from the original draft's `StandardizationLimits` and **scoped to
tautomer computation only** — `max_restarts` is removed (see §2's non-goals;
it doesn't have a defined meaning at this level, and inventing one just to
fill the field would be exactly the kind of un-evidenced field this revision
round exists to catch). `max_transforms`/`max_tautomers` are direct
generalizations of `TautomerConfig::max_iter`/`max_tautomers` — the existing
fields are the literal precedent, counted exactly as §1.3 defines. Both are
**deterministic** budgets: same input, same limits ⇒ same result, always —
this is what makes `MaxTransformsReached`/`MaxTautomersReached` states
reproducible.

`timeout_ms` is a different kind of bound and is called out explicitly rather
than silently accepted at face value: a wall-clock timeout depends on
machine speed and load, so a result state of `TimedOut` is **not**
deterministic across environments, and cannot honestly carry the "same
canonical tautomer regardless of input" guarantee ROADMAP.md's own exit
criterion demands. Recommendation: keep `timeout_ms` as an optional escape
hatch (`None` by default, deterministic budgets alone govern normal
operation), and document `TimedOut` results as explicitly outside the
determinism guarantee.

### 4.2 Result-state enum

```rust
#[non_exhaustive]
pub enum ParentComputationStatus {
    Completed,
    MaxTransformsReached,
    MaxTautomersReached,
    TimedOut,
    Abstained(AbstainReason),
    InvalidInput(InvalidInputReason),
}

#[non_exhaustive]
pub enum AbstainReason {
    NoConfidentOrganicParent,
    AmbiguousFragmentSelection,
}

#[non_exhaustive]
pub enum InvalidInputReason {
    EmptyMolecule,
    UnparsableStructure,
}
```

`Canceled` is **removed** from the original draft (§2): no cancellation
mechanism (token or callback) is designed this phase, and a state nothing can
ever produce is worse than no state at all — it would look like a supported
feature to any caller who pattern-matches on it. If a real cancellation
mechanism is added in a later round, `Canceled` returns then, alongside the
actual token/callback API it depends on.

`Abstained`/`InvalidInput` now carry a **typed reason enum**, not a `String`
— `#[non_exhaustive]` so new reasons can be added without a breaking change,
matching Phase 1's own `abstained: Option<String>` precedent in spirit (a
reason must always be inspectable) but fixing the original draft's mistake of
using a free-text `String` for something Python/WASM callers need to branch
on programmatically.

This is intentionally a **new, orthogonal type**, not an extension of the
existing `PipelineStatus` (`Unchanged`/`Modified`/`CompletedWithWarnings`).
`PipelineStatus` answers "did the molecule change"; `ParentComputationStatus`
answers "did the computation reach a definite answer, and if not, why not" —
different questions, and conflating them would mean every existing
`StandardizationReport` consumer suddenly has to handle
`MaxTransformsReached` even though today's pipeline can't produce it.

### 4.3 Parent-generation functions

A **Parent** is an idempotent, deterministic reduction of one axis of
molecular variability to **one representative structure**, meant to be used
as a grouping/dedup key — computing it never depends on, or mutates, any
other axis. This is the concept RDKit's `MolStandardize` names
`FragmentParent`/`ChargeParent`/`IsotopeParent`/`StereoParent`/
`TautomerParent`/`SuperParent`.

**Revision: `charge_parent` is not `neutralize_charges`.** The original draft
defined `charge_parent(mol) = neutralize_charges(mol)`, which — on an input
like `CC(=O)[O-].[NH4+]` — neutralizes *both* fragments and returns
`CC(O)=O.N`, a 2-fragment molecule. That directly contradicts this section's
own "one representative structure" definition: a Parent must select a single
answer, not leave the input's fragment-selection ambiguity unresolved.
`neutralize_charges` stays exactly as it is — a low-level, all-fragments
mechanical transform, useful in its own right (e.g. inside
`StandardizationPipeline::run`'s existing `NeutralizeCharges` stage, which
*should* neutralize every fragment before fragment selection runs). But
`charge_parent` the *Parent function* is defined as **select the fragment
parent first, then neutralize that one fragment**:

```rust
pub fn fragment_parent(mol: &Molecule) -> (Molecule, TransformationRecord); // = select_fragment, unchanged from Phase 1
pub fn charge_parent(mol: &Molecule) -> (Molecule, TransformationRecord) {
    let (parent, mut record) = fragment_parent(mol);
    let neutralized = neutralize_charges(&parent);
    // record's `after` snapshot and any additional charge-neutralization
    // warnings are updated to reflect `neutralized`, not `parent`.
    ...
}
pub fn isotope_parent(mol: &Molecule) -> (Molecule, TransformationRecord);  // = remove_isotopes + audit
pub fn stereo_parent(mol: &Molecule) -> (Molecule, TransformationRecord);   // = remove_stereo + audit
pub fn tautomer_parent(mol: &Molecule, limits: &TautomerLimits) -> ParentResult;
pub fn super_parent(mol: &Molecule, limits: &TautomerLimits) -> ParentResult;
```

Verified this round: `select_fragment` on `CC(=O)[O-].[NH4+]` already keeps
the acetate fragment (4 heavy atoms, has carbon) over the ammonium fragment
(1 heavy atom, no carbon) under Phase 1's existing ranking — so
`charge_parent("CC(=O)[O-].[NH4+]")` = `fragment_parent` (→ `CC(=O)[O-]`)
then `neutralize_charges` (→ `CC(O)=O`), a **single-fragment** result. §5's
`tp2-17` fixture is corrected to this expected output.

**Typed return value, not a bare tuple.** The original draft returned
`(Molecule, ParentComputationStatus, TautomerAuditRecord)` for
`tautomer_parent`/`super_parent`. Tuples of 3+ heterogeneous fields invite
positional mix-ups and can't grow a field without breaking every call site.
Revised:

```rust
#[non_exhaustive]
pub struct ParentResult {
    pub molecule: Molecule,
    pub status: ParentComputationStatus,
    pub audit: ParentAudit,
}

#[non_exhaustive]
pub enum ParentAudit {
    /// fragment_parent / charge_parent / isotope_parent / stereo_parent
    Transformation(TransformationRecord),
    /// tautomer_parent
    Tautomer(TautomerAuditRecord),
    /// super_parent: one entry per stage, in the fixed order below
    Composed(Vec<ParentAudit>),
}
```

`fragment_parent`/`charge_parent`/`isotope_parent`/`stereo_parent` keep
returning `(Molecule, TransformationRecord)` directly (Phase 1's shape,
unchanged — these are simple mechanical transforms, no result-state question
applies to them: they cannot time out, run out of a transform budget, or
abstain). Only `tautomer_parent` and `super_parent`, which *do* have a
result-state question, return `ParentResult`.

**Typed audit record**, replacing the original draft's `(String, i32)`/
`String` fields:

```rust
#[non_exhaustive]
pub struct TautomerAuditRecord {
    pub selected: MoleculeSnapshot,
    pub candidate_count: usize,
    pub score_breakdown: Vec<ScoreContribution>,
    pub applied_transforms: Vec<AppliedTransform>,
    /// Atoms whose stereo descriptor was removed or invalidated by the
    /// winning transform sequence (empty = none lost).
    pub lost_stereo: Vec<AtomIdx>,
    /// Atoms whose isotope label was affected (moved, or its host atom's
    /// bonding changed) by the winning transform sequence.
    pub affected_isotopes: Vec<AtomIdx>,
}

#[non_exhaustive]
pub struct ScoreContribution {
    pub term: TautomerScoreTerm,
    pub value: i32,
}

#[non_exhaustive]
pub enum TautomerScoreTerm {
    AromaticRing,
    HeteroatomHydrogen { element: u8 }, // atomic number: O=8, N=7, S=16 today
}

#[non_exhaustive]
pub struct AppliedTransform {
    pub rule_id: TautomerRuleId,
    pub affected_atoms: Vec<AtomIdx>,
    pub affected_bonds: Vec<BondIdx>,
}

#[non_exhaustive]
pub enum TautomerRuleId {
    KetoEnol,
    AmideIminol,
    IminolAmide,
    ImineEnamine,
    // ... one variant per existing named TautomerRule, plus:
    AromaticExocyclicShift, // the new §4.4 mechanism
}
```

`lost_stereo`/`affected_isotopes` are now atom-level lists, not booleans —
"did anything get lost" is a strictly weaker answer than "which atoms," and
the strictly-more-informative shape costs nothing extra to compute (the
information already exists at the point each transform is applied).

**`super_parent`'s stage order, fixed and pinned.** The original draft stated
the order as "fragment → charge → isotope → stereo → tautomer" but then
justified it by citing `StandardizationPipeline::run`'s "neutralize before
fragment-select" comment — the reverse order, and a self-contradiction. The
order is now fixed to exactly:

```
fragment_parent → charge_parent → isotope_parent → stereo_parent → tautomer_parent
```

with no appeal to the mutating pipeline's own (different, and differently
motivated) stage order — `super_parent` selects the representative fragment
*first* precisely because every subsequent Parent step should operate on one
fragment, not on a not-yet-resolved multi-fragment molecule. `ParentAudit::Composed`
carries one entry per stage in this exact order, so a caller can inspect
every intermediate result, not just the final molecule. §5's `tp2-23` fixture
is corrected to pin all five intermediate snapshots (input, after each of
the five stages), not just the final output, so a future reordering would be
caught by a fixture diff instead of silently changing behavior.

`StandardizationStep::{NormalizeGroups, FragmentParent, ChargeParent,
IsotopeParent, StereoParent}` stay as enum cases (their `.as_str()` names are
already public API) but are reinterpreted as **labels a caller can report
which Parent function ran**, not as new stages silently added to
`StandardizationPipeline::run()`. Wiring them into the mutating pipeline was
considered and rejected: `run()` already has a specific, deliberate stage
order for producing one standardized molecule; a `Parent` is a different
molecule per axis, not one more mutation in that same sequence.

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
without misfiring on substituents that aren't part of a lactam/lactim system
— §5's new negative controls, phenol/anisole/aniline/pyridine N-oxide, exist
precisely to catch over-eager matching here), so it's scoped as round 2C
rather than done inline in this RFC.

### 4.5 Custom rule/scoring hook

**Revision:** the original draft proposed `scorer: Option<fn(&Molecule) ->
i32>` on `TautomerConfig`. A bare `fn` pointer works in Rust but cannot cross
the Python/WASM boundary, and making `TautomerRule` merely
`#[non_exhaustive]`-public still wouldn't let an external caller construct
one via struct-literal syntax (a `#[non_exhaustive]` struct requires a
constructor function, which wasn't proposed). Revised recommendation: ship a
**declarative, serializable** configuration this phase, deferring arbitrary
callback logic:

```rust
#[non_exhaustive]
pub struct TautomerScoringConfig {
    pub aromatic_ring_weight: i32,
    pub carbonyl_weight: i32,
    pub hetero_hydrogen_weights: Vec<(u8, i32)>, // (atomic number, weight)
    pub substructure_terms: Vec<ScoreTerm>,      // exact ScoreTerm shape: implementation round
}
```

This is representable identically across Rust/Python/WASM (plain data, no
function pointers), which an arbitrary-callback API cannot be without a
separate per-language binding layer. An arbitrary-Rust-closure escape hatch
can be added later, if ever needed, as a clearly-named Rust-only function
(e.g. `canonical_tautomer_with_rust_scorer`) that does not appear in the
Python/WASM-facing surface — not folded into the cross-language
`TautomerConfig`/`TautomerScoringConfig` types.

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
preferred form). A parent fixture asserts a **fixed expected output** (or,
for `super_parent`, a fixed sequence of intermediate outputs, per §4.3's
revision).

```json
{"id": "tp2-05-2-pyridone", "kind": "tautomer_self_consistency", "category": "heterocycle_lactam_lactim", "variants": ["O=c1cccc[nH]1", "Oc1ccccn1"], "currently_passing": false, "note": "confirmed failing on main (99401a6). Root cause and fix design in RFC section 1.1/4.4."}
{"id": "tp2-17-charge-parent-ammonium-acetate", "kind": "parent_generation", "parent": "charge_parent", "input": "CC(=O)[O-].[NH4+]", "expected_output": "CC(O)=O", "note": "charge_parent = fragment_parent (keeps acetate, 4 heavy atoms, over ammonium, 1) then neutralize_charges -- single-fragment result, per section 4.3's revision. NOT the same as the low-level neutralize_charges(mol), which would keep both fragments."}
```

**Corrections applied this round:**
- `tp2-04-nitroso-oxime`'s variants were `["CC=NO", "CC(=NO)"]` — the same
  structure written twice with a redundant parenthesis, not a real tautomer
  pair. Corrected to `["CCN=O", "CC=NO"]` (nitrosoethane / acetaldehyde
  oxime, RDKit-formula-verified as C2H5NO both). This pair is **not**
  currently self-consistent (§1.6) — `currently_passing` corrected from an
  implicit "true" to explicit `false`, with the distinct root cause noted.
- `tp2-17`'s expected output corrected from the two-fragment
  `"CC(O)=O.N"` to the single-fragment `"CC(O)=O"` (§4.3's revision).
- `tp2-23-super-parent-composed` now pins all five intermediate snapshots
  (input → after fragment_parent → after charge_parent → after
  isotope_parent → after stereo_parent → after tautomer_parent), not just
  the final result, and the input molecule was re-verified end to end for
  this exact stage order.
- The `metal_adjacent` category is renamed `disconnected_metal_ion_interaction`
  throughout both files and this RFC — the two fixtures in this category
  (`tp2-15`, `tp2-16`) co-mingle an organic ligand with a bare metal ion via
  `.`-disconnected SMILES, they do **not** contain an explicit coordinate/
  dative bond, so the original name overstated what's being tested.

**New negative/metamorphic fixtures added this round** (guarding against
round 2C's aromatic-shift generalization over-firing, and checking two
general properties any tautomer-canonicalization function should have):
- `tp2-25-phenol-no-dearomatization` — phenol (`Oc1ccccc1`) must stay
  aromatic phenol, never get pushed toward a cyclohexadienone form.
- `tp2-26-anisole-no-mobile-h` — anisole (`COc1ccccc1`) has no mobile H at
  all anywhere in the molecule; must be a strict no-op.
- `tp2-27-aniline-no-exocyclic-imine` — aniline (`Nc1ccccc1`) must stay
  aromatic aniline, never get pushed toward an exocyclic-imine
  (cyclohexadiene-imine) form.
- `tp2-28-pyridine-n-oxide-not-lactam-lactim` — pyridine N-oxide
  (`[O-][n+]1ccccc1`) has an exocyclic O on an aromatic ring but is not a
  lactam/lactim system (the O carries the formal negative charge, not a
  mobile H); must not be reinterpreted as one.
- `tp2-29-idempotence` — `canonical_tautomer(canonical_tautomer(x))` must
  equal `canonical_tautomer(x)` for every already-passing fixture; spot
  checked on acetylacetone (confirmed holding).
- `tp2-30-uracil-atom-order-reorder` — the *same* uracil tautomer (not a
  different one), written via two different ring-traversal orders
  (`O=c1cc[nH]c(=O)[nH]1` vs. `O=c1[nH]c(=O)[nH]cc1`), must already produce
  one canonical output today — confirmed holding (`same=true`) even before
  round 2C's fix, isolating that §1.1's defect is specifically about
  *which tautomer* was input, not general atom-order sensitivity.

All four negative controls (phenol/anisole/aniline/pyridine N-oxide) and the
idempotence/reorder checks are currently passing on `main` — they exist to
stay passing through round 2C's implementation, not because they're failing
today.

Categories now covered in the main set (30 rows after the corrections and
additions above): non-aromatic tautomer controls, ring-internal NH-shift
controls, aromatic lactam/lactim (the confirmed-failing class),
`disconnected_metal_ion_interaction`, zwitterion-interaction, all Parent
functions, and the new negative/metamorphic controls.

The holdout set (5 rows, unchanged from the original draft) generalizes each
category without ever having its expected answer used to shape a rule:
hypoxanthine (a sixth lactam/lactim case), benzimidazole (a fused-ring
ring-internal-shift generalization), a longer-chain zwitterion (GABA), a
combined negative control for all four mechanical Parent functions on one
already-reduced input (toluene), and a companion citation of the existing
`max_iter=1000` regression test.

## 6. Open questions

1. **`TautomerScoringConfig`'s exact `substructure_terms`/`ScoreTerm` shape**
   — sketched but not fully specified (§4.5); decided in the round that
   implements it (2D), not speculatively now.
2. **Does the aromatic-lactam/lactim fix ever need to *dearomatize* a ring
   bond**, or is toggling only the exocyclic bond always sufficient? All six
   confirmed cases (§1.1) only required the exocyclic bond to change — but
   this needs re-checking against a larger corpus during round 2C, not
   assumed to generalize from 6 examples.
3. **Is §1.6 (nitroso/oxime) one mechanism with §1.1, or genuinely separate?**
   Current evidence says separate (different root cause: rule-generality
   trade-off vs. aromatic bond-order matching) — but no fix design is
   committed to yet for §1.6; a future round should decide whether fixing it
   needs a narrower bridge-elem-restricted variant of the existing
   `"1,3-C-to-O-any-bridge"`/`"1,3-C-to-N-any-bridge"` rules (e.g., allow
   forward application only when `bridge_elem` is restricted to N, not
   "any"), which would not risk the acetone-enolization regression that
   originally justified disabling the general rule.
4. **`max_transforms`/`max_tautomers` defaults** — kept at 16/32 (unchanged
   from `TautomerConfig`'s existing defaults) unless round 2B's
   implementation finds a reason to change them; no such reason found this
   round.

## 7. What ships this round

RFC (this document, including this revision) + the two fixture files (§5) +
`ROADMAP.md`/`validation/README.md` pointers. **No changes under
`crates/*/src/**` this round** — matching Phase 1's own RFC-round discipline.
PR #362 stays in **draft**; no merge, version bump, or release this round.

Per this revision, implementation is split into four rounds, each awaiting
separate explicit authorization (mirroring Phase 1's RFC → "次roundへ進んで"
→ Rust core rhythm):

- **Round 2A (this RFC revision):** design fixes only, no production code.
- **Round 2B:** `ParentResult`/`ParentComputationStatus`/`TautomerLimits`,
  `fragment_parent`/`charge_parent`/`isotope_parent`/`stereo_parent`/
  `super_parent`, limit-exhaustion signaling. **No** aromatic lactam/lactim
  fix in this round.
- **Round 2C:** the §4.4 aromatic lactam/lactim fix, validated against the
  negative controls added in §5 (phenol/anisole/aniline/pyridine N-oxide,
  idempotence, atom-order reorder) plus all confirmed-failing fixtures
  turning green.
- **Round 2D:** `TautomerScoringConfig` and the rule-customization surface.

Round order is fixed as 2B → 2C → 2D (Parent API and budget visibility
before the harder aromatic-shift algorithm work, scoring customization
last, since it depends on the audit-record shape 2B establishes).
