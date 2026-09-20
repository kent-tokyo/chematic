# Canonical E/Z encoding plan

## Scope

This note defines the remaining implementation boundary for canonical SMILES
Issue #503 and its coupled E/Z predecessor, #149.  It does not claim that
canonical identity is complete.  Until the acceptance conditions below pass,
`canonical_smiles_stable_key()` must continue to refuse the affected aromatic
direction-stash structures.

## Established evidence

The source-provenance-pinned 1,024-relabeling audit has four divergent coupled
components out of 28, with no cross-correspondence failure.  Three current
aromatic-stash families remain intentionally fail-closed.  Their input
spellings and geometry-preservation tests live in
`crates/chematic-smiles/src/canonical.rs`.

Two bounded writer searches are now rejected:

1. PR #597 exhausts every candidate-bond DFS-priority assignment for the
   coupled component, reparses every candidate, and retains only candidates
   that preserve the input E/Z geometry.  The equivalent spellings share no
   output.
2. PR #599 additionally exhausts every candidate ring-marker-side assignment
   together with every such DFS priority.  Their geometry-preserving output
   sets still have an empty intersection.

These results rule out selecting a raw atom index, a raw bond index, a local
DFS preference, or one ring-marker-side rule as the canonical tie-breaker.
They do not prove that SMILES cannot encode the molecules; they bound only the
current writer's local degrees of freedom.

## Current implementation stage

`CanonicalWriter` now captures rank-fixed E/Z geometry facts before marker
carrier resolution.  When a complete two-ended fact exists, the existing
component solver consumes that extracted relation rather than re-reading a raw
marker from the candidate it is considering.  The spelling-invariance gate
proves this extraction agrees across the three current residual families.

This is only the first stage.  It does not yet construct output slots, solve a
whole component, or make the residual stable key admissible.  Incomplete or
unspecified E/Z inputs retain the previous read path so this refactor does not
silently alter their output contract.

A test-only canonical-DFS inventory now also shows that the three residual
families expose the same eligible tree and ring-token output slots across
their equivalent spellings.  This localizes the remaining work to assigning
the extracted geometry relations to those slots as one component plan; it is
not a justification to pick a slot by parse order.

## Required architecture

The implementation must separate chemical geometry from its SMILES spelling.

1. **Extract geometry facts.** For each specified stereogenic double bond,
   create a representation-independent fact before canonical traversal.  It
   identifies rank-stable reference substituents at both ends and their
   same-side/opposite-side relation.  Aromatic direction stashes and literal
   directional bonds are input encodings of that fact, not canonical output
   locations.
2. **Build the canonical skeleton.** Determine the canonical atom order, DFS
   tree, branches, ring edges, and ring-digit occurrences without consuming
   directional markers.  Each eligible directional token position is an
   explicit output slot with its traversal orientation.
3. **Solve one component plan.** Choose output slots and `/`/`\` polarity for
   the full coupled component simultaneously.  Constraints must preserve every
   extracted geometry fact, respect SMILES ring-token syntax, and prevent a
   shared slot from receiving contradictory assignments.  The solution is
   selected only from canonical-rank and output-slot keys; it must never use
   parse-time atom or bond indices.
4. **Validate before adoption.** Reparse every candidate plan and compare the
   extracted geometry facts, not merely text stability.  Choose the
   lexicographically minimal valid canonical serialization only after this
   semantic check.  An empty, over-budget, or genuinely tied solution set is a
   typed/observable fail-closed outcome for stable-key callers.

The plan may use a bounded component solver, but its cap must be measured and
documented.  A cap must reject the whole component; it must never apply a
partial marker plan.

## Non-goals

- Do not use a CIP label as a substitute for the lower-level E/Z geometry
  fact.  CIP priority remains a separate correctness and adjudication scope.
- Do not change the public writer merely to make a diagnostic pass.  The
  stable-key fail-closed boundary is safer than choosing an unproven output.
- Do not expand 3D embedding breadth before this canonical/stereo boundary is
  closed.

## Acceptance criteria

The production implementation closes this item only when all of the following
are demonstrated:

- all four audited components converge across atom-order and equivalent SMILES
  spelling permutations;
- canonicalize → reparse preserves every extracted E/Z fact and existing
  tetrahedral stereo facts;
- canonicalization is idempotent and the stable-key API accepts only converged
  inputs;
- existing shared-carrier, automorphism, aromatic-stash, and ring-closure
  regressions pass; and
- Rust, Python, Node, and WASM binding contracts expose the same result and
  refusal boundary.
