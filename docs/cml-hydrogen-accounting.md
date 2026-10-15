# CML hydrogen accounting

The v1.1.4 reader/writer preserves explicit and non-node
hydrogens without double counting. The change addresses the failure class in
[Open Babel #2655](https://github.com/openbabel/openbabel/issues/2655).

```python
import chematic

mol = chematic.from_cml(cml_text)
print(mol.formula)  # CH2 for a carbon bonded to two explicit H atoms
output = mol.to_cml()
assert chematic.from_cml_strict(output).formula == mol.formula
```

## Reader and writer contract

CML `hydrogenCount` describes all hydrogens attached to an atom, including
explicit H atom neighbors. This meaning is documented in Table 4 of
[Murray-Rust and Rzepa's CML paper](https://www.ch.ic.ac.uk/rzepa/ijc/search_channels/cml.pdf).
The rule for omitted counts follows the
[Open Babel CML format documentation](https://openbabel.org/docs/FileFormats/Chemical_Markup_Language.html).

- If `hydrogenCount` is present, subtract the explicit H neighbors and store
  the remainder as the core's non-node H count.
- If the attribute is absent and explicit H neighbors exist, do not infer
  any additional H on that atom.
- If both the attribute and explicit H neighbors are absent, retain the
  existing valence-based inference.
- H isotopes are explicit H neighbors too; their atom/isotope labels remain
  in the graph.
- Every output atom receives the total `hydrogenCount`, including zero.
  This preserves partial explicit-H graphs such as methane represented by
  one C node, one H node, and three non-node hydrogens.

For example, a carbon with two explicit H neighbors is CH2 when no attribute
is supplied, CH3 with `hydrogenCount="3"`, and CH4 with `hydrogenCount="4"`.
The attribute does not add three or four more hydrogens to those explicit nodes.

Invalid integers, negative counts, totals below the explicit-neighbor count,
and remainders above the core's u8 capacity return
`CmlError::MalformedXml` with an `invalid CML hydrogenCount` message. Python
exposes these as `ValueError`.
The total itself can exceed 255 when explicit nodes account for the excess.
This is a representation check, not a new valence sanitizer.

The same hydrogen rules apply to `parse_cml` and `parse_cml_strict`, and thus
to Python `from_cml` and `from_cml_strict`. Structural strictness remains a
separate choice. This change does not add spin-multiplicity, stereo, or 3D
coordinate interpretation.

## Compatibility

Older chematic CML output wrote only the core's non-node H count when that
field was set. In files also containing explicit H nodes, this was not the
CML total-H meaning. Such existing files may need regeneration from the
original molecule or correction of their total counts. The reader now
interprets the standard total, rather than guessing which producer created
an ambiguous value. A total smaller than the explicit-neighbor count is
rejected rather than repaired silently.

## Reproduction and validation

The minimized graph is retained as
`validation/openbabel_issues/cml/explicit-ch2.cml`.
The original issue's public attachment was also fetched and reproduced on
the locally installed **Open Babel 3.2.1**. Both that attachment and the
minimized graph were read as methane (`C`, CH4) by Open Babel.

Before the fix, seven of eight new Rust tests failed: explicit CH2,
total-count reconciliation, isotope handling, malformed/contradictory counts,
large representable totals, atom reordering, and writer count semantics.
The existing implicit-only inference control passed before and after.

Validation passed:

- 988 tests across the complete default-feature `chematic-mol` suite,
  including eight new hydrogen regressions; no failures or ignored tests.
- 46 selected Python tests: 29 new hydrogen checks, four existing CML contract
  checks, and 13 local-cation-resonance regressions.
- Clippy with warnings denied for all `chematic-mol` targets, and formatting
  checks for the changed Rust files.

A four-case live cross-read checked absent, 2, 3, and 4 total-H attributes.
chematic produced the expected CH2/CH2/CH3/CH4 compositions, independently
checked by RDKit **2026.09.1**. Open Babel also read every chematic-generated
CML output with the expected formula: explicitly writing the correct total
allows its reader to retain CH2. Version, source-lane and extension/fixture
hashes are recorded in
[`openbabel-cml-hydrogen-2026-10-14.json`](https://github.com/kent-tokyo/chematic/blob/main/validation/results/openbabel-cml-hydrogen-2026-10-14.json).

```sh
CARGO_INCREMENTAL=0 cargo test -p chematic-mol --offline --quiet
CARGO_INCREMENTAL=0 cargo clippy -p chematic-mol --all-targets --offline --no-deps -- -D warnings
python -m pytest crates/chematic-py/tests/test_cml_hydrogen_contract.py -q
```

Python validation used a local debug extension without native InChI support.
These are source checks; no package release or upstream Open Babel change
is included.
