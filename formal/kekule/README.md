# Kekulé assignment: Lean 4 specification and differential test (pilot)

Status: **experimental pilot**. Not part of the workspace, not a required CI
check, not shipped in any package. It exists to measure whether a Lean 4
specification pays for itself on one chematic kernel before doing more.

## What is here

| Path | What |
| --- | --- |
| `lean/Kekule/Spec.lean` | The specification `IsKekule E S D` (double bonds are aromatic bonds, no atom gets two, every must-match atom gets one) and an executable checker, with `checkKekule_iff`: the checker decides exactly the specification. |
| `lean/Kekule/Search.lean` | A deliberately simple backtracking reference `search`, with `search_sound`, `search_complete` and `search_none_iff`: it returns `none` **iff** no Kekulé assignment exists. |
| `lean/Main.lean` | `kekule-check`: reads cases on stdin and judges each Rust result. |
| `harness/` | Standalone Rust binary (own `[workspace]`) that runs `kekulize` / `kekulize_with` and exports inputs and results. |
| `run.sh` | Builds both and runs random graphs plus the repository's SMILES corpora. |

All theorems are proved with no `sorry`; `#print axioms` reports only Lean's
standard `propext`, `Classical.choice`, `Quot.sound`. No Mathlib dependency.

## What is and is not covered

The model covers the **graph part** of `chematic-core/src/kekulization.rs`:
given the aromatic bonds and the must-match set, is the Rust answer a valid
Kekulé assignment, and does Rust return `Err` only when none exists? This
checks the augmenting-path passes, the two-order retry and the Edmonds blossom
fallback together.

It does **not** model the chemistry: which atoms must be matched
(`atom_must_be_matched`) and the ring-fusion / bridgehead N/P exclusion
heuristics are RDKit-parity rules with no independent mathematical
specification. The harness exports the must-set Rust computed and a superset of
the atoms those heuristics may drop; results that use them are reported as
`pass_relaxed` (no full assignment existed) or `heuristic_deviation` (a full
assignment existed but Rust used the relaxed one). Neither fails the run.

Verdicts that fail the run: `INVALID_OUTPUT` (Rust returned something that
breaks the specification), `COMPLETENESS_BUG` (Rust returned `Err` although
Lean found an assignment), `BAD_INPUT`.

## Running

```sh
formal/kekule/run.sh                 # needs lake (elan) and cargo
KEKULE_RANDOM_CASES=50000 formal/kekule/run.sh   # quicker
```

## Pilot results (2026-11-28, Lean 4.29.1)

* 2.2 M random graphs (G(n,p) and fused-ring systems with chords, arbitrary
  must-sets) and about 33 k molecules from 7 repository corpora plus the
  `fuzz/corpus/smiles` seeds: **no `INVALID_OUTPUT` and no
  `COMPLETENESS_BUG`**. Kekulé-form inputs were aromatized first, so the
  aromatize -> kekulize round trip never failed either.
* `heuristic_deviation` appears for one molecule (in two corpora),
  `O=C1N2N(C(=O)C3=C1C=CC=C3)C(=O)C4=CC=CC=C4C2=O`: after aromatization a
  full assignment would need an N=N double bond, and the ring-fusion N
  heuristic keeps the input's N-N single bond instead. Chemically that is the
  better answer; the spec flags it only because the must-set says otherwise.
* Mutation check: disabling the blossom fallback is caught immediately (17
  completeness bugs per 200 k random cases). Two other mutants (dropping the
  reverse-order retry, dropping one blossom path-marking call) were not
  detected, because the blossom pass masks them; they are equivalent at the
  output level on these inputs.
* Side finding: the SMILES parser accepts a ring closure from an atom to
  itself (`C11`, `c11`, `n11n11`) and builds a self-loop bond; `kekulize`'s
  N/P exclusion then counts the loop twice. Found via a
  `heuristic_deviation` on a fuzz seed.
* Upgrading Lean 4.25.0 -> 4.29.1 needed no proof changes; only two
  deprecated `String.trim` calls in `Main.lean`.
