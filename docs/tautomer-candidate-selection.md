# Checked selection from supplied tautomer candidates

RDKit [#9601](https://github.com/rdkit/rdkit/issues/9601) tracks wrapper
contracts, temporary-object lifetimes and an empty-input crash in
`TautomerEnumerator.PickCanonical`. Chematic already provided enumeration and
scoring, but had no public supplied-candidate selection API.

## Python

```python
import chematic

candidates = (chematic.from_smiles(s) for s in ["C=C(C)O", "CC(C)=O"])
selected = chematic.rdkit_pick_tautomer(candidates, max_candidates=1000)
assert selected.rdkit_smiles == "CC(C)=O"
```

The function returns a new `Mol` with the selected source's graph, atom order,
labels, stereo metadata, title and typed properties. Properties are independent
of the source. The underlying immutable graph can share storage safely.
Each yielded `Mol` is owned before requesting the next item, so temporary
objects from generators and lazy sequences cannot become dangling pointers.

Selection uses the existing RDKit 2026.03.1 port's default tautomer score, then
lexicographically smallest RDKit-port canonical SMILES. Exact ties keep the
first source occurrence. No candidates are generated, and no changes or stereo
cleanup are applied to the winner. A singleton enol therefore stays that enol.

An empty candidate collection raises `ValueError`; a zero-atom `Mol` is still
a valid candidate. Unsupported chemistry and oversized candidate graphs also
raise `ValueError` with the zero-based candidate index. Non-Mol elements raise
`TypeError`, and exceptions from the source iterator propagate. There is no
partial winner on any error.

`max_candidates` defaults to 1000 and must be in 1..4096. Duplicate candidates
count toward the limit. At most one additional item is requested to detect
excess input. Each candidate is limited to 1024 atoms and 4096 bonds. Python
integer values outside the signed 64-bit argument range raise `OverflowError`
at conversion. These are input bounds, not a wall-clock timeout guarantee for
individual graph operations.

## Rust

```rust
use chematic_smiles::{parse, rdkit_pick_tautomer};
let enol = parse("C=C(C)O").unwrap();
let keto = parse("CC(C)=O").unwrap();
let pick = rdkit_pick_tautomer([&enol, &keto], 1000).unwrap();
assert_eq!(pick.index, 1);
assert_eq!(pick.score, 5);
assert_eq!(pick.smiles, "CC(C)=O");
```

The native function accepts an iterator of molecule references and returns
`RdkitTautomerPick { index, score, smiles }`. `RdkitTautomerPickError` distinguishes
empty input, invalid limits, count overflow, graph-size overflow and indexed
RDKit-port errors. The caller retains ownership of the original molecules.

## Scope and validation

[Custom integer score callbacks](scored-tautomer-selection.md) are now
available through the separate `rdkit_pick_tautomer_scored` Python function.

This selects from supplied molecules; it does not verify that they belong to
one tautomer family or find the most stable tautomer under specific physical
conditions. This default-score native API has no custom scoring callback or direct
RDKit object support.
The existing `canonical_tautomer`, `rdkit_canonical_tautomer` and enumeration
APIs retain their behavior. WASM and the RDKit compatibility shim do not expose
this new API.

Newer upstream implementations can prefer the number of specified stereo
centers before lexical ties. This API deliberately uses the existing port's
lexical policy and does not claim equivalence for mixed specified/unspecified
stereo ties, or reproduce upstream's final `AssignStereochemistry` step.

Reproduction and independent reference checks are in
`validation/rdkit_issues/tautomer_pick/probe.py`; local results are recorded in
`validation/results/tautomer-pick-2026-10-17/`. This is source validation,
not a package release or an upstream RDKit modification.

Local checks passed five new Rust tests, 34 new plus nine existing Python
tests, and 1,050 matching supplied-candidate selections against RDKit
2025.09.3 and 2026.09.1. The separately recorded mixed-stereo tie boundary is
not included in that parity count. Both upstream versions crashed on empty
candidate lists in isolated subprocesses; chematic raises `ValueError`.
