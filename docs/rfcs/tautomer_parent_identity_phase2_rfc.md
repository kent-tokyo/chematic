# RFC: Tautomer & Parent Identity (ROADMAP.md Phase 2, v0.20.0)

Status: draft, round 2C-1 (mechanism fixation only; no production code
changes yet this round)

> **Revision (2026-08-22, round 2C-1):** before touching
> `crates/*/src/**`, fixed the exact structural condition under which the
> §1.1 aromatic lactam/lactim shift applies, and re-verified the negative
> controls against two cases the original list didn't cover. New: §1.7 (a
> second, distinct, out-of-scope non-convergence found by this round's own
> rigor: cytosine's ring-N-H position is itself non-unique, unrelated to the
> lactam/lactim shift); §4.4a (the per-molecule mechanism table, the
> validity condition for the exocyclic donor/acceptor, and why the fix must
> be a directional step rather than feed the existing score-ranked pool —
> the latter would select the wrong tautomer, measured, not assumed). §5
> gains 5 new fixtures. See §4.4a for detail.
>
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

### 1.7 Confirmed defect (distinct mechanism, out of scope for round 2C):
cytosine's ring-N-H position is itself non-unique

Found while fixing round 2C-1's mechanism table (below), not by design —
the same discipline that found §1.6. Two spellings of cytosine's *keto*
tautomer, differing only in which of the two ring nitrogens flanking the
carbonyl carbon carries the mobile H (both real, chemically distinct: this
is the literature's own N1-H/N3-H cytosine ambiguity, not an artifact of
SMILES-writing), do not converge on `main`:

```
canonical_tautomer("Nc1cc[nH]c(=O)n1")  -> c1[nH]c(nc(N)c1)=O    (H on N4)
canonical_tautomer("Nc1ccnc(=O)[nH]1")  -> c1(=O)nccc([nH]1)N    (H on N7)
```

Both ring nitrogens (`AtomIdx(4)` and `AtomIdx(7)` in the first spelling)
are directly ring-bonded to the carbonyl carbon (`AtomIdx(5)`); the ring is
not otherwise symmetric (only one of the two flanks the exocyclic-amino
carbon too), so these are genuinely different tautomers, not a relabeling of
one. Root cause: `find_direct_aromatic_matches`/`enumerate_direct_aromatic_forms`
(the pre-existing ring-internal-only mechanism that already handles
imidazole/pyrazole/tetrazole/benzimidazole correctly) does not generate this
particular hop for cytosine's ring at all — from either starting spelling,
so neither converges toward the other. This is a **pre-existing gap in the
ring-internal mechanism itself**, unrelated to §1.1's exocyclic donor/acceptor
gap: round 2C's fix only moves H between a ring nitrogen and an *exocyclic*
oxygen, never between two ring nitrogens, so it neither causes nor repairs
this. Left unfixed and out of scope for round 2C, matching §1.6's precedent
— broadening this round to also fix ring-internal N-position selection would
mix two independent mechanisms and muddy which change fixed what. Round
2C-1's own acceptance fixtures for cytosine (§5) hold the ring-N-H position
fixed (always `AtomIdx(4)`'s position, matching the design-driving SMILES
already in `tp2-05`'s sibling rows) specifically to avoid conflating this
gap with the one round 2C fixes; genuine atom-permutation-invariance checks
for cytosine use a re-traversal that preserves which physical nitrogen holds
the H, not a re-traversal that happens to swap it.

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

### 4.4a Round 2C-1: mechanism fixation (structural conditions + negative
controls, fixed before any implementation)

**Per-molecule mechanism table.** For each design-driving molecule, `bridge`
is the ring atom carrying the exocyclic substituent; `donor`/`acceptor` are
named for the keto→enol direction (donor loses the H, acceptor gains it);
`ring path (bridge→acceptor)` is the shortest walk along that one SSSR ring's
own bonds. Atom indices are from parsing the keto SMILES exactly as written
(order-dependent labels, not a structural property — the mechanism itself
must not depend on them; see permutation-invariance acceptance criteria
below).

| molecule | keto SMILES | bridge | exocyclic acceptor (enol side) | ring-N donor (keto side) | ring path, bridge→donor | ring bonds touched |
|---|---|---|---|---|---|---|
| 2-pyridone | `O=c1cccc[nH]1` | C(1) | O(0) | N(6) | 1 (ortho) | none — all stay `Aromatic` |
| 4-pyridone | `O=c1cc[nH]cc1` | C(1) | O(0) | N(4) | 3 (para) | none — all stay `Aromatic` |
| cytosine | `Nc1cc[nH]c(=O)n1` | C(5) | O(6) | N(4) | 1 (ortho) | none — all stay `Aromatic`; exocyclic amino (N0–C1) is a spectator, never touched (§1.7) |
| uracil | `O=c1cc[nH]c(=O)[nH]1` | C(1) **and** C(5) | O(0) **and** O(6) | N(7) **and** N(4) | 1 (ortho), both sites | none — both sites shift together; see order-invariance note below |
| guanine-class purine | `Nc1nc2[nH]cnc2c(=O)[nH]1` | C(8) | O(9) | N(10) | 1 (ortho) | none — all stay `Aromatic`; the 5-ring (imidazole-type) and its own N–H are untouched |

The one structural fact every row shares, confirmed by direct execution, not
assumed: only the exocyclic bridge–acceptor bond changes order
(`Double`↔`Single`); **every ring-internal bond, including the one directly
joining bridge and donor, stays `Aromatic` in both forms.** The H moves, the
ring's own bond-order labels do not.

**Validity condition (not a heuristic — a structural requirement).** The
donor must be an aromatic ring nitrogen reachable from `bridge` *along that
same SSSR ring* at an **odd** bond distance (1 or, for a 6-ring, 3). This is
not a preference tie-break; it is the condition for a real alternating
single/double (Kekulé) path to exist between bridge and donor at all — at
even distance (meta, distance 2 in a 6-ring) no such path exists, so there
is no neutral lactam tautomer to draw, full stop. Confirmed by direct
execution on the discriminating case this analysis specifically requires:
**3-hydroxypyridine** (`Oc1cccnc1`, donor at ring distance **2**, meta) —
`canonical_tautomer` is a no-op today (`c1ncccc1O`) and must stay one; the
fix's own scope predicate excludes it structurally, not by a score that
happens to lose.

**Why an O-only acceptor scope, not O-or-N.** Two negative-control molecules
were checked specifically to probe whether an N-type exocyclic acceptor
(amino/imino) should also be in scope: **4-aminopyridine** (`Nc1ccncc1`,
para, odd distance — same shape as 4-pyridone) and **2-aminopyridine**
(`Nc1ccccn1`, ortho, odd distance — same shape as 2-pyridone). Both are real
amino/imino tautomer pairs by the same distance-parity condition above, and
both are currently no-ops (`c1(ccncc1)N`, `c1nc(N)ccc1`). They are
deliberately **not** brought into scope this round: every confirmed-broken
molecule in §1.1 is an O-type (carbonyl/hydroxyl) shift — cytosine's and
guanine's own exocyclic amino groups are spectators in their broken pairs,
never the reacting atom (table above) — so there is no evidenced defect to
fix on the N-acceptor side, and amino/imino tautomer preference is not
simply "the same rule with N instead of O" (aminopyridines favor the amino,
single-bond form; lactams favor the carbonyl, double-bond form — the
opposite bonded-H pattern), which would need its own separate mechanism
audit, not an extension of this one. Scoping the acceptor element to O keeps
these two molecules correctly excluded by construction (no O to react with),
not by a scoring judgment call.

**Why this must be a directional step, not fed into the existing
score-ranked pool.** `enumerate_direct_aromatic_forms`'s candidates are
ranked by `tautomer_score` (`score_breakdown`, descending) with canonical
SMILES as tiebreak — correct for the ring-internal-only shifts it already
handles (imidazole/pyrazole/tetrazole/benzimidazole), where the competing
forms tie on heteroatom-H weight. It is measurably wrong for this
mechanism: computed directly from `score_breakdown`'s existing weights
(O-H=100, N-H=50, aromatic-ring bonus=1000, both forms fully aromatic under
chematic's model) for 2-pyridone, the **enol** side scores **1100** (one
O-H + aromatic bonus) against the **keto** side's **1050** (one N-H +
aromatic bonus) — sorted descending, the existing pool would select the
chemically minor lactim form, silently inverting the fix. This is why §4.4's
new mechanism must be applied as a **directional step** (added to
`canonical_tautomer_with_config`'s `prefer_forward` loop, the same
architecture the other 42 rules already use to hard-code amide-over-imidic-
acid preference, converging to a fixed point *before* the score-ranked pool
ever runs) rather than as a new candidate generator feeding
`enumerate_direct_aromatic_forms`. It must be added to **both**
`canonical_tautomer_with_config`'s loop and `tautomer_parent`'s equivalent
loop, or the two functions diverge on the same input; `tautomer_parent`
must count each application toward `transforms_applied` so
`MaxTransformsReached` does not under-count.

**Order-invariance requirement, made explicit for the multi-site case.**
Uracil has two independent qualifying sites; guanine's/cytosine's bridge
carbon is flanked by two ring nitrogens (only one of which is the correct
donor per the table above). The implementation must enumerate *all*
qualifying (bridge, acceptor) pairs each pass and select deterministically
(canonical-SMILES minimum among the resulting candidates, mirroring the
existing tiebreak's own reasoning at line ~1126) — never `.first()`/`[0]`
over an unordered neighbor list, which would make the result track input
atom-index order rather than structure (the same order-dependence class
this project has already audited and fixed elsewhere).

**Negative controls, fixed set (all confirmed no-ops on `main` today, must
stay no-ops):**

| control | SMILES | why excluded |
|---|---|---|
| phenol | `Oc1ccccc1` | no ring heteroatom acceptor at all |
| anisole | `COc1ccccc1` | exocyclic O has no H to donate |
| aniline | `Nc1ccccc1` | no ring heteroatom acceptor (all-carbon ring) |
| pyridine N-oxide | `[O-][n+]1ccccc1` | neither side has an H to move (O is anionic, ring N has none) |
| simple amide (acetamide) | `CC(N)=O` | not aromatic — mechanism is ring-gated |
| 3-hydroxypyridine | `Oc1cccnc1` | donor at even (meta) ring distance — no valid Kekulé path (see above) |
| 4-aminopyridine | `Nc1ccncc1` | acceptor is N, not O — out of scope this round (see above) |
| 2-aminopyridine | `Nc1ccccn1` | acceptor is N, not O — out of scope this round (see above) |
| ring-internal NH shift (imidazole, pyrazole, tetrazole, benzimidazole) | e.g. `c1cc[nH]n1` | already-passing existing mechanism; must be untouched, not merely unaffected |
| isotope-bearing (new, §5) | `[18O]=c1cccc[nH]1` / `[18OH]c1ccccn1` | not a control against misfire — a positive case pinning that the isotope label must follow the O atom through the shift, not be dropped or misplaced |
| remote stereocenter (new, §5) | `O=c1c([C@@H](F)Cl)ccc[nH]1` / `Oc1c([C@@H](F)Cl)cccn1` | not a control against misfire — a positive case pinning that a stereocenter uninvolved in the shift must survive unchanged |

**Acceptance criteria this fixes (2C-2/2C-3), restated precisely:**
- All 5 design pairs above converge to one canonical form each.
- The 8 no-op negative controls stay no-ops (byte-identical canonical
  SMILES before/after the round 2C-2 change).
- The isotope and remote-stereocenter cases converge with the label/center
  provably preserved (atom-level check — isotope value and `Chirality` on
  the *specific* untouched atom, not just "canonical SMILES looks plausible").
- Atom-permutation invariance and idempotence hold for every converging
  case, checked with a re-traversal that preserves which physical ring
  atoms are structurally donor/acceptor (per §1.7's caveat for cytosine —
  never a re-traversal that would additionally exercise the separate,
  out-of-scope ring-internal N-position gap).
- Hypoxanthine (holdout) is checked only in round 2C-3, after the above is
  fixed and frozen — never used to adjust the table or predicate above.

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

**New fixtures added in round 2C-1** (§4.4a's mechanism fixation):
- `tp2-31-3-hydroxypyridine-meta-not-lactam` — 3-hydroxypyridine
  (`Oc1cccnc1`), donor at even (meta) ring distance from any exocyclic
  acceptor site; must stay a no-op — the discriminating check that the
  fix's scope predicate is a real distance condition, not a coincidence of
  the 5 design molecules all being ortho/para.
- `tp2-32-4-aminopyridine-not-in-scope` / `tp2-33-2-aminopyridine-not-in-scope`
  — both real amino/imino tautomer pairs at valid (odd) ring distance, but
  with an N acceptor, not O; must stay no-ops this round (§4.4a's "why
  O-only" argument) — currently passing, must keep passing.
- `tp2-34-2-pyridone-isotope-preserved` — `["[18O]=c1cccc[nH]1",
  "[18OH]c1ccccn1"]`, `currently_passing: false`; round 2C-2 must converge
  this pair **and** the `18O` isotope label must remain on the same
  physical oxygen atom in the result, not merely produce a plausible-looking
  canonical string.
- `tp2-35-2-pyridone-remote-stereocenter-preserved` —
  `["O=c1c([C@@H](F)Cl)ccc[nH]1", "Oc1c([C@@H](F)Cl)cccn1"]`,
  `currently_passing: false`; the stereocenter is on a ring substituent
  uninvolved in the shift and must survive unchanged.

Categories now covered in the main set (35 rows after round 2A's and round
2C-1's corrections and additions): non-aromatic tautomer controls,
ring-internal NH-shift controls, aromatic lactam/lactim (the
confirmed-failing class), `disconnected_metal_ion_interaction`,
zwitterion-interaction, all Parent
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
   bond**, or is toggling only the exocyclic bond always sufficient? Round
   2C-1 (§4.4a) answers this structurally for the odd-ring-distance case: no
   ring bond changes order, by construction — the H moves, the ring's own
   `Aromatic` labels don't. Still open: whether corpus molecules outside the
   5 design pairs + hypoxanthine hit some other exocyclic-shift shape this
   table doesn't cover; deferred to round 2C-3's holdout/audit pass, not
   assumed to generalize from 6 examples alone.
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
- **Round 2C:** the §4.4/§4.4a aromatic lactam/lactim fix, split into three
  sub-steps per the same audit-before-code discipline:
  - **2C-1 (this update):** mechanism fixation only — §4.4a's per-molecule
    table, validity condition, negative-control set (including two new
    discriminating checks, 3-hydroxypyridine and the two aminopyridines,
    and the isotope/stereocenter positive cases), and the finding that the
    fix must be a directional step, not fed into the existing score-ranked
    pool. §1.7 (cytosine's separate, out-of-scope ring-N-H ambiguity) also
    surfaced this round. No changes under `crates/*/src/**` yet.
  - **2C-2:** the minimal directional-step implementation against 2C-1's
    fixed conditions.
  - **2C-3:** hypoxanthine holdout validation (checked only now, never used
    to shape 2C-1/2C-2) plus the fused-purine, multi-carbonyl, and
    candidate-order-invariance checks §4.4a's acceptance criteria list.
  - Round 2C as a whole ends as a **draft PR** — not merged, not marked
    ready — per explicit instruction; nitroso/oxime (§1.6) and
    `TautomerScoringConfig` (§4.5) are not in scope for any of 2C-1/2C-2/2C-3.
- **Round 2D:** `TautomerScoringConfig` and the rule-customization surface.

Round order is fixed as 2B → 2C → 2D (Parent API and budget visibility
before the harder aromatic-shift algorithm work, scoring customization
last, since it depends on the audit-record shape 2B establishes).
