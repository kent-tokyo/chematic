# Proof-backed removal of one redundant tetrahedral declaration

RDKit [#9656](https://github.com/rdkit/rdkit/issues/9656) reports a bridged
molecule whose central tetrahedral declaration is globally redundant:
`C=C1C[C@H]2CC(=C)C[C@@H](C1)[C@H]2C`.
This opt-in API removes one selected declaration only after finding an exact
atom mapping that proves its inversion describes the same configured graph.
Parsing and writing defaults are unchanged.

```python
import chematic

source = chematic.from_smiles("C=C1C[C@H]2CC(=C)C[C@@H](C1)[C@H]2C")
report = source.remove_redundant_tetrahedral_stereo(10)
assert report["removed"]
cleaned = report["molecule"]
assert report["witness"][10] == 10
```

For this input, atom 10 is removed and atoms 3 and 8 retain their declarations.
Selecting either bridgehead individually retains it. Atom indices refer to
the supplied molecule, before any canonical serialization.

## Proof and output contract

The operation clones the source, inverts only the selected CW/CCW declaration,
and searches for a complete graph isomorphism from the original to that copy.
The mapping must fix the selected atom and preserve:

- element, isotope, charge, resolved hydrogen count, aromaticity and degree;
- atom-map numbers and nonchemical atom tags;
- adjacency and chemical bond orders;
- every other specified tetrahedral configuration through its complete
  four-ligand frame and permutation parity;
- every specified E/Z fact, including reference-substituent changes and
  reversed double-bond endpoints.

A successful witness permutes the selected atom's ligand frame oddly.
Canonical-SMILES equality and CIP labels are never used as proof.
The output clears only that atom's chirality tag and stored ligand frame,
and clears all stored atom CIP annotations for recomputation. All other
indexed atoms, bonds, tags and configurations are preserved. The output has
cold derived caches; source metadata and caches are unchanged. Python typed
molecule properties are copied independently.

If no witness exists, the API returns an unchanged independent copy with
`removed=False` and `witness=None`. This is conservative retention, not a
classification of the atom as stereogenic. For example, distinct atom-map
numbers can block a mapping. Explicit and implicit graph hydrogens occupy
different proof slots, so `[H][C@H](F)Cl` is deliberately retained.

Errors, including search-budget exhaustion, return no partial cleanup.
This API does not determine 3D feasibility, provide general stereo cleaning,
or certify auxiliary CIP labels. Removing a redundant declaration can affect
auxiliary labels, as discussed in #9656; recompute any required annotations.

## API surfaces

Python `Mol.remove_redundant_tetrahedral_stereo(atom_idx, *, max_visits=1000000)`
returns:

| Field | Meaning |
| --- | --- |
| `molecule` | Independent output molecule |
| `removed` | True only after a complete proof |
| `atom_index` | Selected source atom index |
| `search_visits` | Candidate feasibility probes consumed |
| `witness` | Original-to-inverted atom permutation, or `None` |

Rust exposes `chematic_chem::remove_redundant_tetrahedral_stereo(&molecule,
AtomIdx(index), max_visits) -> Result<TetrahedralStereoCleanup,
TetrahedralCleanupError>`. The result's target field is `atom: AtomIdx` and
its witness is `Option<Vec<AtomIdx>>`.

The WASM API requires an explicit budget:

```javascript
const source = wasm.parse_smiles("F[C@H](F)Cl");
const result = wasm.remove_redundant_tetrahedral_stereo(source, 1, 1000000);
const report = JSON.parse(result.report_json());
const cleaned = result.molecule;
console.assert(report.removed && report.witness[1] === 1);
cleaned.free();
result.free();
source.free();
```

WASM getters expose `molecule`, `removed`, `atom_index` and `search_visits`;
`report_json()` includes the witness. Each `molecule` access returns an
independent handle. Python and JavaScript reject boolean, fractional, text,
negative and out-of-32-bit-range indices or budgets instead of coercing them.

## Bounds and unsupported inputs

Inputs are limited to **64 atoms and 128 bonds**. The accepted visit budget is
**1 through 10,000,000**, with a Python default of 1,000,000. Exhaustion is an
error, not a negative equivalence decision.

The proof profile rejects enhanced stereo groups anywhere in the input,
wildcard/R groups, allene and other non-tetrahedral assignments, incomplete
or nonunique ligand frames, and unspecified locally tetrahedral sites. It
also rejects query/dative/zero-order bonds, loops, parallel edges and
conflicting, stashed or unaccounted literal directions. These checks apply
globally so another configuration cannot silently be ignored.

## Validation

The 18 fixtures include the issue's published molecular examples and original
small controls for identical/distinct ligands, isotope/map/charge/bond-order
distinctions, chiral branches, aromatic branches and E/Z branches. Fixture
provenance and expected decisions are stored in
`validation/rdkit_issues/tetrahedral-cleanup/cases.json`; no upstream attachment
or coordinates were copied.

The dated packet in `validation/results/tetrahedral-cleanup-2026-10-16/`
records source/runtime hashes and compressed comparisons:

- **144 native decisions**, covering both mirrored inputs and four atom orders;
- **144 Python decisions** and **36 actual Node/WASM decisions** agreeing with
  the corresponding native results;
- **64 witnesses** independently checked for bijection, fixed target, exact
  atom/bond attributes and full-frame parity at every specified center;
- RDKit **2026.09.1** independently checking **144 atom/bond graphs** and
  supplying **24 geometric facts** checked against the witness mappings;
- zero mismatches.

The reference verifies graphs and geometric facts; RDKit's cleanup output is
not the decision oracle. Separate tests passed 9 new Rust contracts and 16
existing shared-proof contracts, 60 new Python tests and 35 existing
shared-proof tests, and the Node/WASM fixture and input-boundary suite.
These are bounded local validation results, not performance or release claims.

Regenerate the proof packet with built native, Python and WASM modules:

```sh
cargo build -p chematic-chem --example tetrahedral_cleanup_dump
python scripts/validation/validate_tetrahedral_cleanup.py \
  --native-dump target/debug/examples/tetrahedral_cleanup_dump \
  --wasm-module /path/to/node-package/chematic_wasm.js
```

The Python environment must contain the matching candidate extension and
RDKit 2026.09.1. Node is required for the WASM comparisons.
