# RFC: Explainable Molecule Standardization — Phase 1 (Fragment Policy + Audit Log)

**Status**: Draft — RFC + acceptance fixtures only, no production code this round.
**Refs**: `ROADMAP.md`'s "chematicを100点へ到達させるROADMAP", Phase 1 (v0.19.0).
**Scope of this PR**: this document + `validation/standardization_phase1_fixtures.jsonl`
+ `validation/standardization_phase1_holdout.jsonl`. Zero changes under `crates/*/src/**`.
Production implementation of what's proposed here is explicitly deferred to a later,
separately-authorized round.

## 1. Why this, why now

The roadmap's stated Phase 1 goal is largest-fragment selection, salt/solvent removal,
and a structured transformation audit log, with three hard requirements: every
transformation must be deterministic and atom-order-independent; the audit trail must
return *why* a fragment was kept/removed, not just before/after; and the fragment
classification must NOT start from a large empirical/historical salt-name list — it
must start from a small, principled structural policy plus acceptance fixtures.

This is **not a greenfield design**. `crates/chematic-chem/src/standardize.rs` already
exists (1,872 lines), has no RFC, is not mentioned in `ROADMAP.md` prior to this
document, and evolved through a series of ad hoc "BUG #N" fix commits rather than a
top-down design. It already implements most of the requested surface — largest-fragment
selection, salt removal, a `StandardizationPipeline`/`StandardizationReport` audit-log
skeleton, zwitterion handling, charge neutralization, isotope/stereo stripping, group
normalization, uncharging/reionization, metal disconnection. Two design choices already
in that file are precisely what the roadmap says the new work must avoid, and one
correctness gap blocks the "atom-order-independent" requirement outright. All three were
verified empirically against `main` (commit `1ac442a`) before writing this RFC, not
assumed from reading the source alone.

### 1.1 Confirmed defects in the current implementation

**(a) Tied-fragment-size tie-break is spelling-dependent, not atom-order-invariant.**
`connected_components` (standardize.rs:128) sorts components by `Reverse(len())`. Rust's
sort is stable, so components of equal length keep their BFS-discovery order, which is
their atom-index order, which is input-spelling order.

```
largest_fragment(parse("CCC.CCN")) -> keeps CCC (propane)
largest_fragment(parse("CCN.CCC")) -> keeps CCN (ethylamine)
```

Same two fragments, same molecule, different spelling, different kept fragment. This is
exactly the failure mode this project spent issue #149 (four waves) establishing:
"deterministic for one string" is not the same guarantee as "invariant under
relabeling." Any fragment-selection policy that can decide two equal-size candidates
either needs an *intrinsic* tie key computed from the fragment itself (canonical SMILES,
canonical atom-rank sum, or similar) or must explicitly abstain and say so — never rely
on discovery order.

**(b) Fragment "size" counts every atom in the graph, including explicit hydrogens, not
heavy atoms.** `remove_salts_with_catalog` (standardize.rs:228) ranks fragments by
`component.len()` — raw atom count.

```
largest_fragment(parse("CCC.[H]C([H])([H])[H]"))
  -> keeps the 5-atom explicit-H methane over the 3-heavy-atom propane
```

A caller who's just parsed a file with explicit hydrogens on one fragment and not the
other gets a different, chemically wrong "largest fragment" purely because of
Hydrogen-notation style, not molecular size. `RemoveExplicitHydrogens` already exists as
a *later* pipeline stage (after `LargestFragment`), so this isn't hypothetical input —
it is production order.

**(c) The named `SaltCatalog`'s "ammonium" entry false-positives on real organic
cations, not just counterions.** The catalog (standardize.rs:95) matches
`[#7+;H0,H1,H2,H3]` — essentially any protonated or quaternary nitrogen — against a
whole fragment. Choline (`C[N+](C)(C)CCO`), a legitimate organic cation and common
pharmaceutical counterion-*partner*, matches this pattern:

```
SaltCatalog::default().is_salt(parse("C[N+](C)(C)CCO"))  -> true
```

In `C[N+](C)(C)CCO.[Cl-]` this makes *both* fragments classify as salt, so
`largest_non_salt` stays `None` and the function silently falls back to
`components[0]` — which happens to still be choline today only because choline (7
heavy atoms) dominates chloride (1 atom) in the `Reverse(len)` fragment ordering
upstream. The catalog's classification is wrong; the correct-looking output is masked by
an unrelated size comparison, not produced by correct reasoning. A different molecule
pairing a quaternary-ammonium organic cation with a same-size or larger counterion would
misfire visibly.

This is the general shape of the problem with a named-list catalog: each entry is a
plausible-sounding SMARTS for *some* real salt, but SMARTS-matches-anywhere-in-fragment
has no way to distinguish "this fragment IS a small inorganic/simple-organic counterion"
from "this fragment happens to CONTAIN a substructure common in counterions." A 25-entry
list catches named salts the list's author thought of and produces false positives on
anything else sharing a matched substructure. Growing the list (the natural next step
under the current design) makes both problems worse, not better — which is why the
roadmap's instruction is to not extend it, and this RFC does not extend it either.

**(d) The current audit-log data model cannot answer "which fragment, and why."**
`MoleculeSnapshot` (standardize.rs) carries only `atoms: usize`, `bonds: usize`, and a
hash, per pipeline stage. `StandardizationStepReport` wraps that with `enabled`/
`changed`/`before`/`after`/warnings. There is no way, from a `StandardizationReport`
today, to answer: which specific fragment was removed, what was its formula/SMILES, what
rule fired, why the kept fragment was preferred over it, or (when nothing clearly wins)
that the pipeline abstained and why. This is the literal gap the roadmap's "explainable"
requirement targets.

**(e) `disconnect_metals` runs unconditionally in `StandardizationPipeline::run()` before
any tracked step**, so `report.input.hash` can already differ from `steps[0].before.hash`
with nothing in the report explaining why. Metal disconnection is explicitly out of scope
for this round (per the roadmap's exclusion list); this RFC does not fix it, but the new
audit-log model must not silently inherit this hole — see §4.4.

### 1.2 What Phase 1 actually is, given this

Not "build largest-fragment selection, salt removal, and an audit log from scratch."
It is: **replace the classification mechanism** (named-list matching → structural
policy), **fix the two invariance bugs** in fragment ranking (spelling-dependent
tie-break, atom-count-vs-heavy-atom-count), and **redesign the audit-log data model** to
carry fragment-level rationale. `SaltCatalog` itself is not deleted by this RFC — see
§4.2 for why it's demoted to an explicit, opt-in, disclosed-as-non-default layer rather
than removed outright.

## 2. Goals / Non-goals

**Goals** (this RFC + a later, separately-authorized implementation round):
- A general, structural (non-named-list) fragment classification policy for
  largest-fragment selection and salt/solvent removal.
- Deterministic, atom-order-independent behavior for every decision, including ties.
- A per-fragment audit-log data model: which fragments were present, which were kept,
  which were removed and why, before/after structure, warnings, and explicit abstain
  reasons — not just an aggregate atom/bond-count/hash per stage.
- Small (~30-50 item) acceptance fixtures exercising the policy, plus a smaller holdout
  set not used to design the rules (used only to sanity-check generalization).

**Non-goals for this round** (explicitly out of scope, per the roadmap and prior
instruction):
- Python/WASM bindings for any of this.
- Tautomer canonicalization, uncharger, metal disconnection (all pre-exist in
  `standardize.rs` today, untouched by this RFC beyond the disclosure in §1.1(e)).
- Any production Rust implementation — this round stops at RFC + fixtures, draft PR.
- Full 5,000-molecule corpus remeasurement, 3D work, MMFF94 work, symmetrized-SSSR work.
- Deleting `SaltCatalog`/`remove_salts_with_catalog` outright (back-compat question,
  see §4.2 and §6).

## 3. Design: Fragment Policy

### 3.1 Ranking key (fixes §1.1(a) and §1.1(b))

Replace `component.len()` (raw atom count, spelling-order tie-break) with an explicit,
two-part sort key computed per fragment:

```rust
pub struct FragmentPolicy {
    /// Rank fragments by heavy-atom count, not total atom count (fixes 1.1(b)).
    pub count_heavy_atoms_only: bool,       // default: true
    /// Prefer a fragment containing carbon when heavy-atom counts tie.
    pub prefer_organic: bool,               // default: true
    /// Never strip isotopically-labeled fragments as salt/solvent, even if small.
    pub preserve_isotopes: bool,            // default: true
    /// If no fragment is confidently classifiable, keep the counterion rather
    /// than guess (abstain the CLASSIFICATION, not the whole pipeline step).
    pub preserve_counterion_if_required: bool, // default: false
}
```

Sort key, applied in order, each an intrinsic property of the fragment itself (never
input position):

1. `heavy_atom_count` (descending) if `count_heavy_atoms_only`, else `atom_count`.
2. `has_carbon` (descending, i.e. `true` before `false`) if `prefer_organic`.
3. **Tie-break**: the fragment's own canonical SMILES (lexicographic ascending), computed
   via the existing `chematic_smiles` canonicalizer. This is intrinsic to the fragment's
   structure, not to where it appeared in the input string, closing §1.1(a) — verified
   design requirement: `rank(fragment_from("CCC.CCN"))` and `rank(fragment_from(
   "CCN.CCC"))` must select the same molecule regardless of which literal component list
   order they came from.

If, after all three keys, two fragments are still exactly tied (i.e. they are the same
molecule appearing twice, or two genuinely indistinguishable structures), the policy
records an explicit `Abstained { reason: GenuineTie }` decision rather than picking
one arbitrarily — see §4.3.

### 3.2 Salt/solvent classification (replaces catalog-first, fixes §1.1(c))

Default classification is **purely structural**, evaluated per fragment, with no named
substance list:

- **Always-strip set** (intentionally tiny, chemically fundamental, justified
  individually, not "a list of known salts"): water (`O` with 0 or 2 implicit/explicit
  H, no other heavy atom), and monatomic ions of Na⁺/K⁺/Li⁺/Cl⁻/Br⁻/I⁻/F⁻ — single-atom
  fragments with a nonzero formal charge and a fixed, well-known element identity. This
  set is small enough to enumerate and defend one at a time, unlike a 25-entry organic
  salt-name catalog.
- **Structural salt/solvent heuristic** (generalizes the existing `is_salt_fragment`,
  keeping its shape but making the boundary explicit and testable): a fragment with no
  carbon AND at most N heavy atoms (N is a named constant, not a magic number scattered
  inline) is classified `Salt`; everything else is classified `Kept` by default.
- **Never** classify a fragment as salt purely because it contains a matched
  substructure (e.g. "has a charged nitrogen somewhere") — that is exactly the
  mechanism behind §1.1(c)'s false positive. Any SMARTS-pattern-based catalog is
  reclassified as an **opt-in, explicitly-named, non-default layer** (§4.2), never part
  of the default decision path.

This intentionally draws a narrower "safe to always strip" boundary than
`SaltCatalog`'s 25 entries. Named organic salts/solvates (citrate, mesylate, tosylate,
DMSO, ethanol-of-crystallation, etc.) that don't fit the always-strip set or the
no-carbon/≤N-heavy-atom heuristic are **kept by default and flagged**, not silently
guessed at — consistent with the "don't start with a large empirical list" instruction.
Expanding coverage later is a matter of widening the *structural* heuristic (e.g. a
principled small-polyol/small-carboxylic-acid rule with its own fixtures) or explicitly
opting into the legacy catalog, never appending more named patterns to a matched-anywhere
list.

### 3.3 Determinism / atom-order-invariance requirement

Every decision in §3.1 and §3.2 must be a pure function of the *fragment's own*
structure (its multiset of atoms/bonds/charges plus, where used, its canonical form) —
never of its position, atom-index range, or discovery order in the parent molecule. The
acceptance fixtures in §5 include explicit alternate-spelling pairs for every
tie-relevant case specifically to test this, not just to test the "obvious" spelling.

## 4. Design: Audit Log

### 4.1 Data model

```rust
pub struct FragmentSnapshot {
    pub atom_count: usize,
    pub heavy_atom_count: usize,
    pub formula: String,            // e.g. "C5H14NO+"
    pub canonical_smiles: String,
}

pub enum FragmentDecision {
    Kept { rank_key: String },                     // human-readable rank explanation
    Removed { rule_id: String, reason: String },   // e.g. "always_strip_monatomic_ion"
    Abstained { reason: String },                  // e.g. "genuine_tie", "no_confident_classification"
}

pub struct FragmentRecord {
    pub snapshot: FragmentSnapshot,
    pub decision: FragmentDecision,
}

pub struct TransformationRecord {
    pub rule_id: String,
    pub rule_version: u32,
    pub fragments: Vec<FragmentRecord>,   // every input fragment, not just the removed ones
    pub before: MoleculeSnapshot,
    pub after: MoleculeSnapshot,
    pub warnings: Vec<StandardizationWarning>,
}
```

This directly answers the roadmap's six required fields per transformation: input
fragments (`fragments[].snapshot`), adoption/exclusion rationale
(`fragments[].decision`), before/after structure (`before`/`after`, already present
today), warnings (already present today), and abstain reasons
(`FragmentDecision::Abstained`).

### 4.2 `SaltCatalog`'s disposition

Not deleted. Demoted: kept as an explicitly-named, opt-in *supplementary* layer
(`FragmentPolicy` gains no field enabling it by default; a caller who wants the legacy
25-entry behavior asks for it by name, e.g. `remove_salts_with_catalog(&mol, &catalog)`
stays exactly as it is today, unchanged, for callers who explicitly reach for it). The
new default path (`standardize`/`largest_fragment` per §3) does not consult it. This
avoids a silent behavior change for any existing direct caller of
`remove_salts_with_catalog` while making the *default* pipeline path structural, matching
the instruction to not build the new work on top of the named list.

### 4.3 Determinism of the audit log itself

`FragmentRecord`s within a `TransformationRecord` are sorted by the same canonical-form
tie-break key as §3.1, not by discovery order — so the *audit log's own field order* is
also atom-order-invariant, not just the kept/removed decision.

### 4.4 Existing gaps this data model does not silently inherit

Per §1.1(e), `disconnect_metals` currently runs before any tracked step. This RFC does
not fix metal disconnection (out of scope), but requires that whatever implementation
follows either (a) wrap it in its own `TransformationRecord` (even a no-op one, `rule_id:
"metal_disconnect_v1"`, empty `fragments` when nothing was disconnected) so
`report.input.hash` always equals some step's `before.hash`, or (b) explicitly document
in the report that pre-report mutation occurred and why. Silently-unexplained hash drift
is exactly the "not explainable" failure mode Phase 1 exists to close — it must not
reappear in the new model for a stage this round doesn't touch.

### 4.5 Public API compatibility

`MoleculeSnapshot`, `StandardizationStepReport`, `StandardizationReport`, and
`StandardizationWarning` are already `pub`, re-exported from `chematic_chem`'s crate
root, and (per a repo-wide grep performed for this RFC) have **zero internal consumers**
outside `chematic-chem` itself — `chematic-py` and `chematic-wasm` only call
`standardize()`/`StandardizationPipeline::new(StandardizeOptions{..})` for the
transformation itself, never inspect the report's fields. Internal blast radius of
changing these types is therefore zero. External blast radius is not zero:
`chematic-chem` is published independently on crates.io, so any outside caller
constructing these structs via struct literal or exhaustively destructuring them would
see a breaking change from added fields. Three options, **not decided by this RFC** —
left for explicit approval before implementation:

1. Add `transformations: Vec<TransformationRecord>` to `StandardizationReport` as a new
   field, mark both it and `StandardizationStepReport`/`MoleculeSnapshot`
   `#[non_exhaustive]` going forward. Smallest diff; still a semver-minor-at-best,
   semver-major-if-strict change (adding a field to a public, non-`#[non_exhaustive]`
   struct is technically breaking).
2. Introduce a parallel, new, fully-`#[non_exhaustive]` report type (e.g.
   `ExplainableStandardizationReport`) alongside the existing one, leaving
   `StandardizationReport` frozen. No breaking change at all; two report shapes to
   maintain going forward.
3. Treat this as the deliberate start of a `chematic-chem` v0.19 semver-minor (or
   -major, if the project's versioning policy calls added-but-non-breaking-in-practice
   fields major) bump, documented in `CHANGELOG.md` under "Breaking changes."

Recommendation for the next round: option 1 with `#[non_exhaustive]` — smallest real
diff, zero known internal breakage, and consistent with this project's past practice of
purely-additive minor releases (see `ROADMAP.md`'s v0.15.0-v0.18.0 "Shipped" sections,
all described as "purely additive... no breaking changes"). Final choice is a product
decision for the implementation round, not this RFC.

## 5. Acceptance fixtures

`validation/standardization_phase1_fixtures.jsonl` (design corpus, rules may be tuned
against these) and `validation/standardization_phase1_holdout.jsonl` (held out, used only
to sanity-check that the shipped policy generalizes — never used to pick thresholds).
Categories, matching the roadmap's list: simple salts, zwitterions, hydrates,
organometallics, multi-organic fragments, isotope-containing compounds, equal-size tie
cases (each with an alternate-spelling pair to test invariance), charged fragments, plus
the three confirmed-defect cases from §1.1 as named regression fixtures.

Row shape (JSONL, one object per line, following `validation/canonical_residual_
fixtures.jsonl`'s established convention of a loose, id/category/input/note shape rather
than a rigid schema):

```json
{"id": "std-p1-01-simple-inorganic-salt", "category": "simple_salt",
 "input": "CC(=O)O.[Na+]", "expected_kept_fragment": "CC(=O)O",
 "expected_removed": ["[Na+]"], "rationale": "monatomic_always_strip_ion",
 "note": "sodium acetate; Na+ is a monatomic always-strip ion, acetic acid is kept as the organic parent"}
{"id": "std-p1-14-equal-size-tie-alt-spelling", "category": "equal_size_tie",
 "input": "CCC.CCN", "alt_spelling": "CCN.CCC",
 "expected_behavior": "same kept fragment regardless of spelling (intrinsic tie-break, canonical-SMILES-ordered)",
 "note": "regression fixture for confirmed defect 1.1(a): pre-RFC largest_fragment() gives a DIFFERENT answer for these two spellings on main today"}
{"id": "std-p1-15-explicit-h-vs-heavy-atom", "category": "equal_size_tie",
 "input": "CCC.[H]C([H])([H])[H]",
 "expected_kept_fragment": "CCC",
 "note": "regression fixture for confirmed defect 1.1(b): pre-RFC largest_fragment() keeps the 5-atom explicit-H methane over 3-heavy-atom propane on main today"}
{"id": "std-p1-16-choline-chloride-catalog-false-positive", "category": "charged_fragment",
 "input": "C[N+](C)(C)CCO.[Cl-]",
 "expected_kept_fragment": "C[N+](C)(C)CCO",
 "expected_removed": ["[Cl-]"],
 "rationale": "chloride is a monatomic always-strip ion; choline is organic (has carbon) and heavy-atom-dominant, kept regardless of any nitrogen-charge pattern",
 "note": "regression fixture for confirmed defect 1.1(c): the legacy SaltCatalog's ammonium SMARTS also flags choline itself as salt; this fixture asserts the NEW structural policy keeps choline for the right reason (size+organic), not by masked catalog-fallback luck"}
```

Full fixture set is written directly into the JSONL files alongside this RFC (see
`validation/standardization_phase1_fixtures.jsonl` / `..._holdout.jsonl`), not
duplicated here.

## 6. Open questions for the implementation round

- §4.5's public-API compatibility strategy (recommendation given, not decided).
- Exact value of "N" in the no-carbon/≤N-heavy-atom structural salt heuristic (currently
  4 in the existing `is_salt_fragment`; the RFC does not re-derive this from first
  principles, it is inherited pending fixture-driven validation in the next round).
- Whether `preserve_counterion_if_required` (§3.1) is needed for v0.19.0 or can wait —
  no fixture in this round currently exercises it; flagged as a plausible near-term
  follow-up, not required to close Phase 1's stated three features.

## 7. What ships this round vs. later

This round: this RFC, the two fixture files, and this document's disclosure of the three
confirmed defects in the existing implementation. Nothing under `crates/*/src/**`
changes. The next round (separately authorized) implements §3/§4 in
`chematic-chem`, resolves §4.5's open compatibility question, and adds the Rust-side
tests that make the acceptance fixtures executable — Python/WASM bindings follow only
after the Rust core is stable, per the roadmap's explicit ordering.
