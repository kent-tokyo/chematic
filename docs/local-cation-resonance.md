# Local cation resonance

The unreleased source API enumerates local cation resonance forms for
`A-C=B+ ↔ A+=C-B`, with non-aromatic N/O/S termini and neutral carbon.
It implements the oxygen and sulfur examples from
[RDKit #9458](https://github.com/rdkit/rdkit/issues/9458), independently of
RDKit's implementation.

```python
import chematic

mol = chematic.from_smiles("CSC(C)=[N+](C)C")
report = mol.local_cation_resonance(max_forms=32, max_transforms=1024)
print(report["status"])             # complete_within_scope
print(sorted(m.smiles for m in report["forms"]))
```

The result contains two forms: one with the positive charge on iminium N,
the other with the positive charge on sulfonium S. Starting from either
form reaches the same pair. The oxygen analog behaves the same way.

Rust callers use:

```rust
use chematic_chem::{enumerate_local_cation_resonance, LocalResonanceLimits};
use chematic_smiles::parse;

let mol = parse("CSC(C)=[N+](C)C").unwrap();
let report = enumerate_local_cation_resonance(&mol, &LocalResonanceLimits::default()).unwrap();
assert_eq!(report.forms.len(), 2);
```

## Scope and preservation

The neutral terminal has total bond-order valence 3 for N or 2 for O/S;
the positive terminal has one more bond-order unit. The shift swaps one
single and one double bond and transfers +1 between the termini. All
hydrogens stay on their original atoms, including explicit hydrogen nodes
and isotope-labeled hydrogens. Atom/bond indices, connectivity, element
identities, isotope labels, atom maps, and total formal charge are retained.

The search keeps the starting form first, then visits states breadth first.
Deduplication uses charge and bond-order assignments on the **indexed**
graph. Symmetry-equivalent placements on different atoms remain separate;
guanidinium therefore has three indexed forms even though their canonical
SMILES can be identical. The engine does not need canonicalization to
terminate or deduplicate its states.

This API does not enumerate aromatic/Kekulé resonance, anions, charge-separated
forms, or tautomers, and does not rank forms by energy. Unsupported valence
patterns have no eligible shifts. Completion describes this local rule's
reachable states, not general resonance coverage or the validity of unrelated
atoms in the input.

Stereo-bearing and query inputs are rejected. Stereo is not removed
implicitly; callers who want a stereo-free projection must request that
projection separately before enumeration.

## Limits and status

| Limit | Default | Accepted range |
| --- | ---: | ---: |
| `max_forms` (includes input) | 32 | 1–1024 |
| `max_transforms` (includes transitions to previously seen states) | 1024 | 1–100000 |
| `max_atoms` | 256 | 1–1024 |

Python reports `scope = "local_nos_cation_shift_v1"`, a `forms` list of
`Mol` objects, and `transforms_attempted`. Status is one of:

- `complete_within_scope`: the local-rule search exhausted its reachable states.
- `form_limit_reached`: another distinct indexed form exceeded `max_forms`.
- `transform_limit_reached`: another transition exceeded `max_transforms`.

A search that finishes exactly at the form cap reports completion. Truncated
searches retain valid partial results but do not promise independence from
the starting form or traversal order. Invalid limits and unsupported inputs
return typed Rust errors; the Python binding maps these errors to `ValueError`.
Python argument conversion can also reject unrepresentable integer inputs.

## Source validation

The issue fixture is
`validation/rdkit_issues/resonance/local_cation_pairs.smi`.
On the locally installed RDKit **2026.09.1**, the issue's flags produced two
oxygen forms from either starting form, but only the original sulfur form.
This establishes the specific upstream behavior, not an RDKit-wide comparison.

The new Rust suite exercises both pairs, independent sites, symmetric indexed
placements, isotope/H/map conservation, atom reordering, exact and truncated
budgets, invalid config, stereo/query rejection, and negative controls.
The Python suite exercises the public report, error and budget behavior.

Focused validation passed nine Rust tests and 25 Python tests (13 for the new
API plus the existing RDKit-tautomer and stereoisomer suites). Clippy with
warnings denied passed for `chematic-chem` and `chematic-py`; changed Rust
files passed formatting checks.
The complete default-feature `chematic-chem` suite subsequently passed
1,135 tests with zero failures and five pre-existing ignored tests.

An independent RDKit sanitization/canonicalization cross-check validated all
eight returned forms across the four issue inputs, including expected pairs,
formula and total charge. The machine-readable result, RDKit version and
local extension hash are in
[`local-cation-resonance-2026-10-14.json`](../validation/results/local-cation-resonance-2026-10-14.json).
This used a local debug Python extension without native InChI support.

```sh
CARGO_INCREMENTAL=0 cargo test -p chematic-chem --test local_cation_resonance --offline
python -m pytest crates/chematic-py/tests/test_local_cation_resonance.py -q
```

These are local source APIs and validation; no package release is included.
