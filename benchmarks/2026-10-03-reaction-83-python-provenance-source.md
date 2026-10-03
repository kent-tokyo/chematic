# Reaction 83: checked Python atom-origin and template-map source gate

This is an **unpublished source-candidate** result on branch
`codex/reaction-83-binding-provenance` from main `8bd18bad`. An editable
development-profile Python extension was compared with pinned RDKit 2026.03.6
on the existing 57+26 fixtures. The fixture SHA-256 values are
`65b446de1ce66b7d9c355cf8b6bbbfe61d670c955e02931e8e8b8e564c558d3b`
and `4e4ce3ebcea4bcea5398f109554074d97e9f252293fa906493ae3a61d40725a1`.

| Classification | Rows |
| --- | ---: |
| Graph, atom origin and template map all agree | 76 |
| Typed unsupported: reactant tetrahedral semantics | 3 |
| Diagnosed product-valence refusal | 1 |
| Invalid in both engines | 3 |

The source gate is `scripts/reaction_python_checked_provenance_gate.py`.
It reads `run_smirks_checked(..., rdkit_compat=True)`, reorders each product's
metadata using `Mol.smiles_with_atom_order()`, then compares graph, origin and
template-map labels separately with RDKit. It fails on any change to the above
outcome counts. The local JSON report SHA-256 was
`036c85cdc02ea778a16f24aa2178f665070933a141315834ef21562f1838d7d6`.
This hash describes the local development extension, **not** a release wheel.

The WASM checked JSON now returns the same two metadata arrays, aligned to
canonical-SMILES parse order. A native Rust test checks that reordering, and
the Node 83-row adapter test checks metadata shape, source-index bounds and
map-number bounds. The Node test is not itself a full RDKit provenance oracle.
CI now runs the Python oracle gate on rebuilt Linux and macOS release-profile
source wheels. The gate records the wheel SHA-256 and requires the installed
native extension bytes to match the wheel contents; those CI results remain
pending until the updated branch is run there.

The three typed unsupported rows and the valence refusal are not counted as
exact parity. Raw reaction embedding counts need not match. This corpus does
not establish broad SMIRKS compatibility, yield or selectivity accuracy, and
none of the published v1.0.31 artifacts have been remeasured here.
