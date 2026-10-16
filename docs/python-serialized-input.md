# Python serialized text: str and UTF-8 bytes

RDKit [#9601](https://github.com/rdkit/rdkit/issues/9601) distinguishes serialized
input from filesystem paths and identifiers. chematic implements the
`str | bytes` boundary for these thirteen Python inputs:

| Function | Existing keyword |
| --- | --- |
| `from_smiles`, `is_valid_smiles` | `smiles` |
| `from_cxsmiles` | `s` |
| `is_valid_smarts` | `smarts` |
| `from_mol_block`, `from_mol_block_with_coords`, `from_mol_block_with_diagnostics` | `block` |
| `from_mol_v3000`, `from_mol_v3000_with_coords`, `from_mol_v3000_with_diagnostics` | `block` |
| `from_mol2` | `mol2_str` |
| `from_cml`, `from_cml_strict` | `cml_str` |

```python
import chematic

mol = chematic.from_smiles(b"C[C@H](F)Cl")
assert mol.smiles == chematic.from_smiles("C[C@H](F)Cl").smiles
assert chematic.is_valid_smarts(b"[#6]-[#8]")
mol, metadata = chematic.from_cxsmiles("CO |$炭素;酸素$|".encode("utf-8"))
assert metadata["atom_labels"] == ["炭素", "酸素"]
```

Only the Python argument extraction changes. Both types feed the same native
parser and retain its return shape, indexed graph, metadata, coordinates,
diagnostics, parse limits and chemical policies. Existing argument names stay
valid. Immutable bytes and strings are borrowed for the native call rather
than converted through a new Rust `String`. This is not a promise of an
allocation-free parser or a performance improvement.

## Errors and scope

- Invalid UTF-8 bytes raise `UnicodeDecodeError` before parsing, including for
  the boolean validity functions. The error retains the original bytes and
  the failing byte range. There is no replacement decoding or encoding guess.
- Valid UTF-8 containing malformed chemistry follows the original parser's
  behavior: existing typed `ChematicInputError`, `ValueError`, diagnostics or
  `False`, as appropriate. NUL bytes are passed to the parser and never used
  as a C-string terminator.
- Other objects raise `TypeError` at the call boundary. Mutable `bytearray`,
  `memoryview`, `PathLike`, numbers and arbitrary conversion objects are not
  coerced. `__str__`, `__bytes__` and `__fspath__` are not invoked.
- Native `str`/`bytes` subclasses are accepted by reading their immutable
  storage. A Python string containing a lone surrogate retains the original
  string-extraction error behavior.

Identifiers such as property keys, file path APIs, batch inputs and unlisted
format readers retain their existing contracts. Binary CDX, arbitrary binary
formats and encoding detection are outside this text-only boundary. No Rust
chemistry API, WASM API, global configuration or parsing default is changed.
This implements one subset of #9601, not its complete wrapper checklist.

## Validation

`validation/rdkit_issues/python-serialized-input/cases.json` contains original
small SMILES/SMARTS/CXSMILES controls and baseline-generated CO format blocks.
The fixtures include isotope/tetrahedral stereo, UTF-8 labels and coordinate/
diagnostic return forms. No upstream source or attachments were copied.

`crates/chematic-py/tests/test_serialized_input.py` verifies positional and
existing-keyword equality, exact decoder ranges, unsupported types, conversion
hook avoidance, subclasses, malformed-text error preservation and stub types.
The type stub advertises `str | bytes` for exactly the selected thirteen inputs.

Local evidence under `validation/results/serialized-input-2026-10-16/` records
baseline/candidate extension hashes and the native string/bytes comparisons.
The baseline is a retained local candidate, not a published-package claim.
Hosted CI is separate from these local checks.

On the local Python 3.13 native candidate, 448 tests passed (263 new boundary
tests and 185 existing parser/error tests). The isolated comparison confirmed
all thirteen baseline byte rejections, preserved string outputs, matching
candidate byte outputs, 39 malformed-text outcomes and thirteen lone-surrogate
errors. These checks cover the named fixtures, not arbitrary chemistry.

To regenerate the comparison with two already-built packages:

```sh
python validation/rdkit_issues/python-serialized-input/compare.py \
  --baseline /path/to/baseline-package --candidate /path/to/candidate-package \
  --output /tmp/serialized-input-comparison.json
```

```sh
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=/path/to/candidate-package \
  python -m pytest crates/chematic-py/tests/test_serialized_input.py \
    crates/chematic-py/tests/test_input_errors.py -q
```
