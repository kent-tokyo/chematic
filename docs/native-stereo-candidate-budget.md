# Native carbon-stereo candidate budget

This local source fix uses the aminoadamantane and inositol diagnostics in
[RDKit #9701](https://github.com/rdkit/rdkit/issues/9701) to improve native
stereo enumeration. The chematic baseline already canonicalized the two
reported inositol spellings to the same stable string and enumerated nine
assignments. It returned an empty list for aminoadamantane because CH2 atoms
incorrectly consumed its six-candidate budget. This fix removes that failure.

## Use

```python
import chematic

mol = chematic.from_smiles("NC12CC3CC(CC(C3)C1)C2")
forms = mol.enumerate_stereoisomers()
assert len(forms) == 4
assert len({form.smiles for form in forms}) == 4

inositol = chematic.from_smiles("OC1C(O)C(O)C(O)C(O)C1O")
assert len(inositol.enumerate_stereoisomers()) == 9
```

Rust continues to export
`chematic_chem::enumerate_stereoisomers(&Molecule) -> Vec<Molecule>`.
Python's existing method keeps its signature. Both retain the existing
six-candidate/64-combination boundary and return an empty collection above it.
This limit counts eligible assignment sites, not the `num_stereocenters`
descriptor: some graph-symmetric sites become distinguishable only after
neighbouring stereo is assigned.

## Candidate selection and preserved data

Candidates are unassigned, nonaromatic, non-wildcard carbon atoms with four
single-bond/implicit-H ligand slots. Single-bond directional carriers are
accepted; multiple, dative, zero-order and query bonds do not qualify.

- Two or more non-node hydrogens rule out a site immediately.
- Repeated identical terminal substituents also rule it out. A terminal
  signature compares element, isotope, charge, aromaticity and resolved
  hydrogen count. Atom maps and nonchemical tags cannot make identical
  terminal ligands stereogenic.
- An implicit ordinary H is compared with explicit terminal ordinary H.
  Deuterium stays distinct, so H/D-dependent tetrahedral assignments survive.
- Graph-symmetric ring paths are retained. Excluding all graph-tied ligands
  before assigning the surrounding centers would discard pseudoasymmetric
  possibilities in the diagnostic rings.

This is conservative local screening, not a complete potential-stereocenter
classifier. Equivalent larger substituent subgraphs may still consume the
budget. Canonical SMILES remains the deduplication key, including the input's
atom-map labels; the candidate change does not make maps chemically stereogenic.

Enumeration now clones the source graph and edits only the selected centers'
relative parity, hydrogen declaration and recorded neighbour frame. That
preserves indexed atoms/bonds, atom maps/isotopes, tags, R labels, existing
stereo frames, enhanced atom/bond groups, directional stashes and ring-closure
metadata. New assignments invalidate stored CIP annotations, including at
fixed centers where Rule 5 priorities may change; callers can recalculate
those labels. Pure copies keep their annotations. Derived caches start empty.
`Clockwise` and `CounterClockwise` describe parity relative to the recorded
frame, not absolute R/S labels.

Python also copies the molecule name and typed property store independently
to every returned form. Editing a form's properties leaves the source and
other forms intact. The no-candidate path uses the same metadata-preserving
copy behavior.

The native method enumerates carbon assignments only. It does not enumerate
unassigned E/Z, heteroatom or coordination stereo, and it does not test 3D
embeddability. The RDKit-compatible enumeration APIs remain separate. In
particular, four aminoadamantane assignments here match a `tryEmbedding=False`
reference; this is not a claim that four independently isolable cage
stereoisomers exist. No general canonicalization or RDKit-wide fix is claimed.

## Pinned diagnostic evidence

| Diagnostic | Native before | Native after | RDKit legacy | RDKit modern |
| --- | ---: | ---: | ---: | ---: |
| Aminoadamantane | 0 | 4 | 4 | 5 |
| Inositol | 9 | 9 | 9 | 13 |

RDKit 2026.09.1 ran with `onlyUnassigned=True`, `unique=True`,
`maxIsomers=1024`, `tryEmbedding=False`. The primary reference uses legacy
stereo perception; the modern counts reproduce the issue's duplicate-count
symptom on this host. Every native output is independently read by RDKit
under legacy perception, and its InChI is compared with the complete reference
set. This verifies membership, not just the total count.

The independent recipe covers 12 structures, including explicit-H
aminoadamantane, isotope-dependent centers, identical terminal ligands, a long
methylene-rich chain, a meso control and a specified center next to an
unassigned one. Identity, reverse, one cyclic shift and two seeded random
atom permutations give 60 source numberings and 160 compared assignments.
All InChI sets match, all native SMILES round trips are fixed points, and
indexed graph sizes/bonds remain unchanged.

The source hash of the original enumerator, module/source hashes, fixture
provenance, license, minimization, exact oracle settings, expected InChI sets
and normalized case hash are recorded in
`validation/results/rdkit-stereo-candidates-2026-10-15.json`.
The test recipe uses minimal abstract molecular graphs and original project
code under `MIT OR Apache-2.0`; no third-party corpus, source or attachment
was copied. The temporary runtime was recreated during validation, with the
same RDKit version pinned before final results were generated.

```sh
cargo test -p chematic-chem -p chematic-smiles --lib --test stereo_candidate_budget --offline
python -m pytest crates/chematic-py/tests/test_stereo_candidate_budget.py crates/chematic-py/tests/test_mol.py crates/chematic-py/tests/test_rdkit_stereoisomers.py
python scripts/validation/validate_stereo_candidate_budget.py
```

Rust passed 1,347 tests (948 chemical-library, 386 SMILES-library, 13 new
contracts), with two pre-existing ignored tests. Python passed 80 (17 new,
63 existing). Clippy with warnings denied, formatting, static bindings,
product naming and strict MkDocs passed. The refactor-budget checker passed
after removing the resolved size exception for the shortened enumerator;
function allowances decreased from 88 to 87. CI registers the integration
test explicitly.

The implementation is in an isolated source copy retaining preceding local
candidates. The portable patch contains only this fix relative to that
snapshot. It is uncommitted and unpublished; these results concern local source.
