# Source-preserving PDB coordinate merging — 2026-10-18

Related open issue: https://github.com/openbabel/openbabel/issues/2181

The native strict coordinate selector already preserves source atom identity;
no source-merging API exists. This adds Python assembly with explicit chain
namespaces, consecutive atom/TER serials, associated ANISOU/CONECT updates and
source-boundary TER. Source strings, atom mappings and excluded entry metadata
remain in immutable sidecars. Native selection and molecular writers are unchanged.

- Base: `f99c7b8ab33098a2e1fd9fd568aee370bf79ab38`.
- `pytest.log`: 231 tests pass, 43 new plus 26 existing selection and 162
  compatibility cases. Includes actual 32-source/50000-atom/250000-record caps,
  serial 99999 success/100000 refusal, lazy exceptions, mapping collisions,
  metadata consent, source identity, reference multiplicity and immutability.
- `oracle.py`, `oracle-summary.json`, `oracle.log`: independent Gemmi 0.7.3 reads
  96 assemblies / 382 sources / 1338 atoms / 382 ANISOU records / 1146 CONECT
  target references. Identity, coordinates, numeric fields, tensors and explicit
  references agree after declared serial/chain mappings; coordinate error is zero.
  Raw ATOM/HETATM/ANISOU fields also agree outside serial/chain columns.
- `protein.pdb`, `ligand.pdb`, `merged.pdb`: original three-atom fragments,
  not a complete protein or physical model. No upstream attachment/source code
  is copied. Fixed-column trailing blanks are retained via scoped .gitattributes.
- `openbabel.pdb`, `openbabel.stderr`, `openbabel-summary.json`: Open Babel 3.2.1
  `-j` joins the same control, changing two ATOM kinds, all three residue names
  and chains, with zero TER. chematic emits two source-boundary TER records and
  retains record kinds/names. This is one CLI observation; the original Windows
  GUI and linked structure are not tested, and no all-version/API claim is made.
- `typing-valid.log`, `typing-invalid.log`, `mypy.ini`: strict helper/public
  consumer pass; four expected negative errors. Legacy imports/stubs are silent;
  package-wide strict typing is not claimed.
- `docs-and-typing.log`, `mkdocs.log`: Python 3.9 grammar, executed guide, strict
  docs build and four rendered API anchors pass.
- `source-gates.log`: source stubs, manifest, binding inventory, unsafe surface,
  refactor budget and canonical names pass. No wheel supplied.
- `native-source-comparison.log`: native selection module/reader/writer, strict
  3D PDB source and Python binding match the verified fixture revision. No fresh
  current-main native build or whole-extension equivalence claim.
- `provenance.json`: source/runtime/reference/evidence hashes and issue snapshot.

Only implicit single-model sources are accepted, with disjoint chain namespaces
after explicit mapping. Every source must contain atoms. Unsupported whole-entry
metadata refuses unless explicit exclusion consent is supplied; supported record
errors, hybrid-36 and model wrappers always refuse. Numeric fields, atom/element
names and ATOM/HETATM classification are never inferred or regenerated.

Hard caps: 32 sources, 8 MiB/source, 16 MiB aggregate input/output, 50000 atoms,
250000 physical input/output records, 1024 bytes/line and five-column serials.
Missing final source TER is inserted, including HET-only sources. Existing TER
positions are retained; missing internal chain boundaries are not reconstructed.
All CONECT rows move after coordinates, retaining target order/multiplicity.
Output normalizes LF, blank lines and CONECT spacing; complete source stays
available. No file writing, coordinate alignment, new cross-source bonds,
hydrogen/element/residue perception, geometry quality, crystallographic assembly,
full PDB deposition, performance or published-wheel validation is claimed.

## Reproduce

Put this revision's Python package with a compatible native extension on
PYTHONPATH, plus independent Gemmi 0.7.3 for the oracle. Tests do not need Gemmi
or Open Babel; selection tests need validation/openbabel_issues/pdb-selection/.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.13 -m pytest \
  crates/chematic-py/tests/test_pdb_coordinate_merge.py \
  crates/chematic-py/tests/test_pdb_selection.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 python3.13 validation/results/pdb-coordinate-merge-2026-10-18/oracle.py
python3.13 -m mypy --config-file validation/results/pdb-coordinate-merge-2026-10-18/mypy.ini \
  crates/chematic-py/python/chematic/_pdb_merge.py \
  validation/results/pdb-coordinate-merge-2026-10-18/typing_valid.py
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-pdb-merge-site
```

The oracle expects /opt/homebrew/bin/obabel, this host's Open Babel executable.
For typing set MYPYPATH to the candidate package parent. Check the negative
consumer separately, expecting exit 1 with four errors.
