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

The WASM checked JSON returns the same two metadata arrays, aligned to
canonical-SMILES parse order. A native Rust test checks that reordering; the
Node 83-row adapter test checks metadata shape and index bounds. It is not a
full RDKit provenance oracle.

PR [#741](https://github.com/kent-tokyo/chematic/pull/741) at `d3b27f1c`
passed [CI run 37111467272](https://github.com/kent-tokyo/chematic/actions/runs/37111467272).
Linux Python job `111170094572` and macOS arm64 job `111170094540` built
release-profile **source wheels** and installed them before the oracle gate.
Both wheel reports classify the rows as 76/3/1/3 above; their normalized row
arrays have the same SHA-256,
`f2c7b9715add53d75a8b700a7abb18e10317a34cdabbefe323bdbb3bbe8fbc59`.
The gate checked that the loaded native extension matched the wheel entry byte
for byte. Linux wheel SHA-256:
`9a3a9835b8ae90bd10b3c04c679670c5da3bce8bbec35cf8b0d349b8ae6a98ad`;
macOS wheel SHA-256:
`138be16072f5b450b0acfb161f96a1c844fad07e7df1320daf45c0cf2bb4e584`.
WASM job `111170094590` passed a release-mode Node build and the 83-row Node
test. These CI wheels were rebuilt from the PR; they are not registry artifacts.

The three typed unsupported rows and the valence refusal are not counted as
exact parity. Raw reaction embedding counts need not match. This corpus does
not establish broad SMIRKS compatibility, yield or selectivity accuracy, and
none of the published v1.0.31 artifacts have been remeasured here.
