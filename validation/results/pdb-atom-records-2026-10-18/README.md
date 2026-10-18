# Typed PDB atom attributes — 2026-10-18

Related open issue: https://github.com/openbabel/openbabel/issues/2708

This implements the PDB attribute-access slice, not a full MDAnalysis converter.
The existing native PDB/protein data already retain some attributes; the general
Python atom-selection view lacks occupancy/B factor. The new immutable typed
read view exposes all validated source ATOM/HETATM rows with those attributes,
source indices/lines and complete original text. Unknown metadata is retained;
no rewritten selection text is exported.

- Base: `f99c7b8ab33098a2e1fd9fd568aee370bf79ab38`.
- `pytest.log`: 242 tests pass, including 54 new cases and all 720 permutations
  of six atoms. Existing selection and compatibility regressions also pass.
- `oracle.py`, `oracle-summary.json`, `oracle.log`: independent Gemmi 0.7.3 PDB
  parser checks 96 cases / 2592 atoms. Coordinates agree exactly. Max occupancy
  difference 2.8611e-8, max B-factor difference 3.6622e-6, consistent with Gemmi's
  float32 fields versus source-derived Python float64. Both readers retain
  modeled chain/residue/atom IDs, altloc, insertion codes and ATOM/HETATM kinds.
- `attributes.pdb`: original two-atom fixed-column control. Significant trailing
  spaces are retained, with scoped .gitattributes whitespace rules.
- `openbabel.stdout`, `openbabel.stderr`: local Open Babel 3.2.1 rewrite keeps
  occupancies 0.5/0.25 but drops altloc A/B and emits B factors 0/0 rather than
  23.45/45.67. This is one observed read/write path, not an all-version/API claim.
- `typing-valid.log`, `typing-invalid.log`: strict helper/public consumer pass;
  negative consumer exits 1 with four expected argument errors. Legacy imported
  modules/stubs are silent; package-wide strict typing is not claimed.
- `docs-and-typing.log`, `mkdocs.log`: Python 3.9 source grammar, executed guide,
  strict build, rendered guide and three public API anchors pass.
- `source-gates.log`: source stubs, manifest, binding inventory, unsafe surface,
  refactor budget and naming gates pass. No wheel was provided.
- `native-source-comparison.log`: strict PDB/selection source and selection
  binding paths match the verified old extension's source revision. No fresh
  current-main native build or whole-extension equivalence claim.
- `provenance.json`: source, runtime fixture, reference and evidence hashes.

Bounds: 8 MiB ASCII input, 50000 atom records, 250000 physical lines and 1024
bytes per native line. Callers may only lower these limits. Validation uses the
existing strict single-model coordinate-selection profile with an empty internal
selection. Native ANISOU/TER/CONECT checks remain; hybrid-36 and multi-model
sources refuse. All original records stay available in source_text.

Blank/absent occupancy and B-factor fields return None, while explicit zero
returns 0.0. Truncated nonblank, malformed and nonfinite numeric fields refuse.
Raw element/charge slices are not chemical assignments. No altloc selection,
bond/conformer inference, physical occupancy assessment, graph editing,
MDAnalysis/trajectory ownership integration, real-structure corpus, performance
or published-wheel validation is claimed.

## Reproduce

Add a copy of this revision's Python package with a compatible native extension
to PYTHONPATH. For the oracle also add an independent Gemmi 0.7.3 installation.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.13 -m pytest \
  crates/chematic-py/tests/test_pdb_atom_records.py \
  crates/chematic-py/tests/test_pdb_selection.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 python3.13 validation/results/pdb-atom-records-2026-10-18/oracle.py
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-pdb-records-site
```

The existing selection tests use `validation/openbabel_issues/pdb-selection/`.
Oracle requires `/opt/homebrew/bin/obabel`, the local path on this host. Tests
require neither Gemmi nor Open Babel. Strict mypy uses silent imported modules
and stubs, strict=True, ignore_missing_imports=True, and MYPYPATH pointing to
the candidate package. Check the helper with typing_valid.py, then check the
negative consumer separately and expect four errors.
