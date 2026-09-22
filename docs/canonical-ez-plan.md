# Canonical E/Z encoding plan

## Scope

This note records the bounded canonical-SMILES implementation for Issue #503
and its coupled E/Z predecessor, #149. It does not claim that canonical
identity is complete. The semantic complete-slot planner resolves the four
formerly divergent corpus components and admits stable keys for the proven
aromatic direction-stash path. Other multi-E/Z structures remain fail-closed
until separately proven.

## Established evidence

The historical source-provenance-pinned 1,024-relabeling audit had four
divergent coupled components out of 28, with no cross-correspondence failure.
After the complete-slot planner was adopted, a clean-main rerun at commit
`2fe0cdd7` passes all 28 components across 1,024 seeded relabelings each, with
zero divergent outputs and zero cross-correspondence failures. Evidence:

- `validation/results/ez_shared_carrier_coupling_mechanism_audit_1024_2026-09-22.jsonl`
- `validation/results/ez_shared_carrier_coupling_mechanism_audit_summary_1024_2026-09-22.json`

The former 4/28 baseline and rejected 7/28 close-side experiment remain
preserved as historical evidence. The current gate does not reuse those
results as proof of the adopted implementation.

Earlier bounded searches established why a broader representation was needed:

1. PR #597 exhausted candidate-bond DFS priorities, reparsed every candidate,
   and found no common geometry-preserving output for the historical residuals.
2. PR #599 added every candidate ring-marker side and still found no common
   output.
3. The complete-slot probe fixed the canonical DFS skeleton and included every
   raw directional carrier as plain, `/`, or `\`, plus both legal occurrences
   of every ring token. This produced a common semantic output.
4. Re-ranking every geometry-preserving candidate retained that common output.

These results rule out raw atom indices, raw bond indices, local DFS preference,
or one ring-marker-side rule as canonical tie-breakers. The adopted selector
uses the lexicographic minimum of semantically valid serializations.

## Current implementation

`CanonicalWriter` captures rank-fixed E/Z geometry facts before marker-carrier
resolution. When a complete two-ended fact exists, the component solver uses
that relation rather than re-reading a raw marker from each candidate.

The production planner constructs the canonical output-slot universe, searches
at most eight slots and 65,536 candidates, reparses every candidate, checks a
non-recursive Morgan-rank E/Z signature, and selects the lexicographic minimum.
It is restricted to aromatic direction-stash input. Over-budget, empty, or
semantically invalid candidate sets retain the previous fail-closed path.

The permanent regression suite includes all four components that diverged in
the historical 1,024-relabeling audit. Existing shared-carrier, automorphism,
ring-closure, E/Z semantic-reparse, and tetrahedral-preservation tests remain
mandatory. The bounded corpus result closes Issue #149's measured residual; it
does not remove the conservative refusal boundary for unmeasured coupled shapes.

## Required architecture

The implementation keeps chemical geometry separate from its SMILES spelling:

1. **Extract geometry facts.** Identify rank-stable reference substituents and
   their same-side/opposite-side relation before canonical traversal.
2. **Build the canonical skeleton.** Determine atom order, DFS tree, branches,
   ring edges, and explicit output slots without consuming input markers.
3. **Solve one component plan.** Search coupled substituent candidates and raw
   directional carriers together, never using parse-time indices as a winner.
4. **Validate before adoption.** Reparse each candidate and compare geometry
   facts before selecting the lexicographic minimum.

The component cap rejects the whole component. It must never apply a partial
marker plan.

## Non-goals

- Do not use a CIP label as a substitute for the lower-level E/Z geometry fact.
- Do not expand the proven aromatic-stash scope merely to make a diagnostic pass.
- Do not expand 3D feature breadth before the remaining stereo/identity exits.

## Acceptance criteria

The bounded Issue #149 item is complete when:

- all 28 measured coupled components, including the four historical residuals,
  converge across 1,024 atom-order relabelings;
- the four historical residuals also retain explicit equivalent-spelling tests;
- canonicalize → reparse preserves E/Z and tetrahedral stereo facts;
- canonicalization is idempotent and stable keys admit only proven paths;
- shared-carrier, automorphism, aromatic-stash, and ring-closure regressions pass;
  and
- binding contracts preserve the same result or refusal boundary.

Broader A2 work, including phosphorus CIP adjudication and independent gold
evaluation, remains separate from Issue #149.
