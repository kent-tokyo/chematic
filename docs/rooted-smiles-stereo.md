# Stereo-checked rooted SMILES

RDKit [#9646](https://github.com/rdkit/rdkit/issues/9646) reports a configuration
flip when a constructed implicit-H bridgehead is written first. The reduced
case with tags only on the two bridgeheads reproduces in the pinned RDKit
2026.09.1 wheel: original mapped atoms 3/5 have CIP labels S/R; rooting at
atom index 2 gives R/R after parsing, and rooting at index 4 gives S/S.
Converting the source centers to explicit bracket-H counts changes this
upstream behavior, so that is a separate control rather than a substitute
for the implicit-H case.

## APIs

```python
import chematic

mol = chematic.from_smiles("C1=C[C@H]2CC[C@@H]1O2")
text, order = mol.smiles_rooted_checked(2)
# order[parsed_index] is the corresponding original graph atom index.
assert order[0] == 2
```

Rust: `chematic_smiles::write_rooted_checked(&mol, root)` returns
`Result<(String, Vec<AtomIdx>), RootedSmilesError>`. The new error type is
non-exhaustive. Python maps refusal to `ValueError` (negative integer
extraction raises `OverflowError`). WASM exposes
`mol.smiles_rooted_checked_json(root)`, returning JSON with `smiles` and
`atom_order` keys; non-integral, negative and non-finite roots reject.

Output is non-canonical and begins at the requested root. Remaining
disconnected components follow. The correspondence is the actual DFS output
order, not a canonical rank or the intermediate atom permutation. Source
atoms, tags, stereo frames and cached properties are not mutated.

## Verified stereo frames

The checked path requires every tetrahedral tag to have a recorded four-slot
reference frame: either four graph neighbors with no implicit H, or three
neighbors and exactly one inferred/bracket H represented by
`STEREO_H_SENTINEL`. Duplicate, foreign or missing ligand slots reject.
Constructed graphs must declare the frame to which their chirality tag
refers; adjacency order is never silently treated as that declaration.
The original implicit-H count can remain inferred.

The implementation remaps the recorded metadata with the existing checked
atom permutation, then reuses the native writer's ligand-parity remapping.
It maps the writer's output atom indices back to the original graph indices.
It does not introduce a special bridgehead parity rule or promote H counts
to make the reproduction disappear.

## Scope

- At most 256 graph atoms, checked before recursive traversal or permutation
  allocation. The bound also limits recursion in the reused native writer.
- Ordinary tetrahedral configurations, isotopes, charges, atom maps,
  disconnected components and the existing native E/Z/dative spelling.
- Axial/non-tetrahedral/lone-pair stereo, enhanced stereo groups and zero/query
  bond orders reject. Arbitrary source and drawing metadata are not serialized.
- No stereocenter discovery, redundant-stereo cleanup, canonicalization or
  geometric-embedding guarantee. Some formal bridgehead assignments in the
  test matrix need not have realizable 3D conformations.

This is an opt-in native API. Existing `write`, canonical SMILES and the
separately versioned RDKit-compatibility writer retain their behavior.
Matching a reference configuration does not imply byte-identical RDKit text
or support for every stereochemical class.

## Reproduction and evidence

The original constructed fixture is
`crates/chematic-smiles/tests/support/rooted_bridgehead.rs`. The native suite
combines two H modes, four formal bridgehead-tag assignments, sixteen edge
insertion orders and seven roots: 896 constructed round trips. It verifies
the four-ligand frames against the original atom correspondence, plus
independent invalid-input and scope controls.

```sh
cargo test -p chematic-smiles --test rooted_stereo
cargo run -p chematic-smiles --example rooted_stereo_probe > /tmp/rooted-cases.json
# Provide a freshly compiled chematic extension and RDKit 2026.09.1 on PYTHONPATH.
python3 scripts/validation/validate_rooted_stereo.py --cases /tmp/rooted-cases.json
```

The oracle runner constructs the reference graph independently, assigns
indexed CIP labels with RDKit's accurate `rdCIPLabeler`, and compares the
native and Python outputs. Optional `--wasm-module` compares a real Node/WASM
build with Python exactly. The compressed packet retains every reference
and candidate output, input parameter tuple, correspondence and CIP result;
the summary records source/runtime hashes and upstream mismatch counts.
These are source-candidate checks, not published-package results.

Local results on base `a12702eb` plus this change (2026-10-16):

| Check | Result |
| --- | --- |
| Rust contracts | 5 tests passed, including all 896 constructed round trips |
| Python contracts | 15 tests passed |
| Native/Python indexed CIP vs constructed reference | 896/896 each matched |
| Real Node/WASM vs Python text and correspondence | 896/896 exactly matched |
| RDKit 2026.09.1 implicit-H rooted control | 320/448 matched, 128/448 inverted |
| RDKit 2026.09.1 explicit-H rooted control | 448/448 matched |

Warnings-denied Clippy (three libraries, new test and probe), rustfmt, refactor
budget, binding inventory, shared manifest, stub syntax and strict MkDocs
passed. The 128 upstream failures belong to this fixed constructed matrix;
they are not a failure rate for general molecules or all RDKit versions.
