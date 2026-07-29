# ADR: automorphism-orbit pruning for VF2 (issue #139)

Status: proposed, not yet implemented. Written before any pruning code, per
this program's instruction to compare approaches before coding.

## Context

`chematic-smarts::match_vf2`'s backtracking search (`match_recursive`,
`crates/chematic-smarts/src/match_vf2.rs:333`) tries every target atom as a
candidate for the next unmapped query atom:

```rust
for t in 0..ctx.mol.atom_count() {
    ...
    if !eval_atom_query(&query.atoms[q_next].query, t_idx, ctx) { continue; }
    if !bonds_compatible(q_next, t_idx, mapping, query, ctx) { continue; }
    mapping.insert(q_next, t_idx);
    match_recursive(query, ctx, mapping, results, max);
    mapping.remove(&q_next);
}
```

When several target atoms are in the same local automorphism orbit of the
target molecule (e.g. a tert-butyl group's three methyl carbons, or a
gem-dimethyl pair), every one of them independently satisfies the same atom
and bond constraints and leads to an isomorphic remainder of the search —
but this loop tries all of them, each triggering a full recursive
sub-search. Issue #139's concrete case (`tert_butyl_B(1)` on a
di-tert-butylphenol scaffold) resolves to `BudgetExhausted` at the shipped
1,000,000-visit budget purely from this redundant re-exploration, folding to
a false-positive PAINS flag.

PR #193 (`crates/chematic-smiles/src/canonical_automorphism.rs`,
`canonical_partition.rs`, `canonical_search.rs`) solved an analogous but not
identical problem for SMILES canonicalization: pruning redundant
individualization branches when a Morgan-rank-tied cell is a genuine
automorphism orbit of the *molecule against itself*. It is used here as
design inspiration only (the exact-refinement + backtracking bijection +
verified-mapping structure), per this task's explicit instruction not to add
a `chematic-smarts` → `chematic-smiles` dependency.

**Why VF2's problem is not the same shape as canonicalization's**, and why
literally reusing PR #193's code is not a drop-in fit even where a
dependency edge already exists (`chematic-smarts` already depends on
`chematic-smiles` per this crate's layer table):

- Canonicalization searches for one graph's automorphism group against
  *itself*, to prune individualization order when generating a canonical
  string. VF2 searches for embeddings of a *query* graph into a *target*
  graph; the relevant symmetry to exploit is the target's own local
  automorphism group, used to prune *candidate choices* during backtracking
  — a different role in a different search.
- PR #193's `VertexColor` is deliberately keyed to `CanonicalWriter`-visible
  attributes: explicit-bracket-H state (`hydrogen_count.is_some()`), atom-map
  (verbatim, since canonical SMILES preserves it), and a stereo-pinning rule
  tuned to avoid corrupting encoded chirality in a *written string*. None of
  this vocabulary maps cleanly onto VF2's actual need — a SMARTS query atom
  can be a wildcard (`*`), a recursive `$(...)` pattern, or a primitive
  combination with no notion of "bracket-H state" or "atom map" at all.
  Reusing `VertexColor` as-is would require either loosening it (risking the
  recently-hardened, three-round-independently-reviewed canonicalization
  engine for an unrelated consumer) or wrapping it with an adapter mapping
  SMARTS semantics onto SMILES-writer semantics — itself a source of subtle
  mismatches, not a simplification.

## Decision drivers

- Issue #139's acceptance criteria are narrow: fix the tert-butyl false
  positive at the shipped budget, no slowdown on the normal corpus, preserve
  the existing `MatchOutcome` API and fold policy exactly. Not "build a
  general graph-automorphism library."
- No new/deepened dependency between `chematic-smarts` and
  `chematic-smiles`'s internals (explicit task constraint).
- Do not destabilize `chematic-smiles`'s canonicalization engine, which
  underwent three independent review rounds (correctness, false-prune,
  performance) as part of PR #193 — any change to its internals for an
  unrelated consumer's benefit re-opens that surface.
- Layer table (`CLAUDE.md`): `chematic-core` has zero dependencies;
  `chematic-smiles` → `core`; `chematic-smarts` → `core, perception, smiles`.
  A module placed in `chematic-core` is reachable by both `chematic-smiles`
  and `chematic-smarts` today with no new edge at all.

## Options considered

### A — Shared: reuse chematic-smiles's existing automorphism engine

Expose (`pub`, not `pub(crate)`) `canonical_automorphism.rs`'s bijection
search and call it from `chematic-smarts` across the dependency edge that
already exists in `Cargo.toml`.

Rejected. The edge existing in `Cargo.toml` does not make this the same as
"no new dependency" in the sense this task's constraint means: it would
create new *internal-API* reliance on a module whose coloring scheme is
purpose-built for SMILES-writer output, not SMARTS query matching (see
Context above). It also risks the canonicalization engine's stability for a
consumer it was never designed for, and either direction of adaptation
(loosen `VertexColor`, or wrap it) trades a small amount of new code for a
correctness risk in code that's already been hardened through real review
rounds.

### B — Narrow reimplementation: VF2-specific pruning inside chematic-smarts

Implement orbit-based candidate pruning directly in `match_vf2.rs` (or a
small private sibling module within `chematic-smarts`), tailored to VF2's
actual data model:

- Compute the target molecule's exact automorphism group once per
  `find_matches`/`find_matches_with_config` call, using a small backtracking
  bijection search over a VF2-specific coloring: element, isotope, formal
  charge, aromatic flag, and `BondOrder` on incident bonds — the attributes
  VF2's own atom/bond compatibility checks already use, nothing borrowed
  from `chematic-smiles`'s writer-output vocabulary.
- **The group licensing a skip at recursion depth *k* is the pointwise
  stabilizer of the atoms already fixed by the current partial mapping, not
  the global automorphism group computed once at the top.** This is not a
  simplification to defer — using the global group at every depth is
  outright wrong: in di-tert-butylphenol, the two tert-butyl groups are
  related by a global automorphism, but once the partial mapping has fixed
  an atom inside one of them, that automorphism no longer fixes the mapped
  atom and is no longer valid for pruning at that point in the search. The
  stabilizer must shrink (be recomputed or incrementally maintained) as the
  mapping grows, or the pruning silently discards real embeddings — a
  correctness regression, not a speedup, and one a "does the tert-butyl
  fixture reach `NotFound`" test alone would not catch.
- At each `for t in 0..ctx.mol.atom_count()` step in `match_recursive`, skip
  a candidate `t` only if it is in the same orbit, *under the current
  depth's stabilizer*, as a `t'` already tried and exhausted for the same
  `q_next` at this recursion level.
- Orbit computation is not free and must not be paid unconditionally: gate
  it on the target actually having a nontrivial stabilizer at that depth (or
  compute lazily on the first tie), so the 480-PAINS-pattern-against-a-
  zero-symmetry-molecule case pays nothing extra. The issue's own acceptance
  criterion ("no slowdown on the normal corpus") depends on this.
- No new module boundary reusable by other crates; all new code stays
  inside `chematic-smarts`, no new dependency anywhere.

Smallest footprint, most directly targeted at the reported failure mode,
lowest review surface, zero risk to `chematic-smiles`. Cost: if a second
consumer for exact graph automorphism appears later (e.g. an MCS
improvement), this code isn't reusable as-is and would need extracting then
— judged acceptable since no second consumer exists today (YAGNI).

### C — New graph module: generic, reusable, decoupled from chematic-smiles

Build a standalone exact-colored-graph-automorphism/orbit module from
scratch (not reusing chematic-smiles's code, just its structural ideas),
placed either inside `chematic-smarts` as a new file or in `chematic-core`
(reachable by both `smiles` and `smarts` with no new edge, since both
already depend on it) so a future `chematic-smiles` refactor could adopt it
without a new dependency — not required or attempted now.

Produces a genuinely reusable primitive, but at meaningfully higher upfront
cost: its own exact-refinement algorithm, its own backtracking search, its
own from-scratch correctness test suite (small-graph exhaustive + fuzz,
mirroring PR #193's verification depth) — all to serve a single consumer
today. Premature generalization for a problem with no second consumer yet.

## Decision

**Option B — narrow reimplementation**, scoped inside `chematic-smarts`.
Rationale: issue #139's acceptance criteria are narrow and specific; Option
A is structurally disqualified by both the explicit task constraint and its
own adaptation risk; Option C's generality has no second consumer to justify
its cost today. If a future issue needs the same primitive elsewhere,
extracting Option B's code into a shared module at that point (effectively
migrating toward C) is a small, well-motivated follow-up — not blocked by
anything this decision does now.

## Consequences

- New code lives entirely in `chematic-smarts` (likely a new
  `crates/chematic-smarts/src/automorphism.rs`, called from
  `match_vf2.rs`). No `Cargo.toml` changes.
- `MatchOutcome`'s three-way API and conservative fold policy are untouched
  — pruning must only change which redundant branches are skipped, never
  which leaf mappings are found. **This is the property to prove, not
  assert**, and it only holds if a skip is licensed by the pointwise
  stabilizer of the currently-mapped atoms at that exact recursion depth,
  not the target's global automorphism group (see Option B above) — a wrong
  implementation here silently drops real embeddings, which no
  budget/timing test would catch.
- Verification plan for the implementation PR, in order of how likely each
  is to catch the stabilizer-vs-global-group mistake:
  1. **All-embeddings equality, not boolean agreement**: for every
     molecule/pattern pair in the existing PAINS/Brenk corpus,
     `find_matches`'s full *result set* (not just `MatchOutcome` or a match
     count) must be identical before and after pruning. A `NotFound`-only
     check on the `tert_butyl_B(1)` fixture would pass even if pruning
     silently dropped unrelated embeddings elsewhere — this is the check
     that specifically distinguishes "pruned correctly" from "found fewer
     matches."
  2. A small-graph exhaustive suite (mirroring PR #193's `n<=5`/`n<=8`
     approach, adapted to two-graph query-into-target matching rather than
     single-graph self-automorphism): every embedding found by the pruned
     search equals the brute-force set for hand-constructed symmetric
     targets at multiple recursion depths, specifically including cases
     where a candidate becomes distinguishable only after part of the
     mapping is already fixed (the di-tert-butylphenol shape itself is one
     instance of this general case).
  3. The `tert_butyl_B(1)` / di-tert-butylphenol fixture resolves to
     `NotFound` at the shipped budget (not a raised budget) — necessary but
     not sufficient on its own, per point 1.
  4. A benchmark confirming no slowdown on non-symmetric patterns/molecules
     (orbit computation must be gated/lazy, per Option B's note on cost, so
     the common zero-symmetry case pays nothing extra).
