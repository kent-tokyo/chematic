# XBS scene writer local evidence — 2026-10-17

Baseline `35b932f9` has no public `write_xbs`, `XbsSpecies` or `XbsBondRule`.
Exact Python files were verified against Git and imported with the reused
extension. `provenance.json` records hashes, commands and upstream metadata.

- `pytest.log`: 209 passed, zero failures/skips (47 new / 162 existing).
- `crossread.json`: hofmann 0.21.1 parsed 96 synthetic scenes, 5,344 atom
  rows, 288 styles and 96 distance rules, plus a separate 10,000-atom scene.
  Coordinates matched binary64 bytes, including negative zero/subnormals.
  No rendering, bond detection or original XBS GUI was executed.
- The Open Babel 3.2.1 diagnostic retains original input/output and confirms
  `bs` is a different format and no `xbs` identifier is listed. This is not
  claimed as a defect in its `bs` implementation.
- `typing.json`: mypy 1.8.0 strict helper passes, valid public consumer
  passes, four invalid calls yield four expected errors. Consumers use
  `--follow-imports=silent --ignore-missing-imports`; package-wide strict
  cleanliness is not implied. The isolated helper copy is unchanged.
- `source-gates.json`: refactor budget, Python source stubs, native binding
  inventory, cross-binding manifest and unsafe gates pass. Native function
  inventory remains 213 registered / 317 WASM exports.
- `mkdocs.log`: strict documentation build passes (existing unlisted pages
  produce informational messages). `product-name.log`: canonical spelling.

The implementation is standard-library Python. The native extension was
reused solely for import/existing tests. The oracle was installed separately
with `--no-deps`; it is not a runtime dependency or vendored implementation.
The generated fixtures are original project data under MIT OR Apache-2.0.
No fresh native build, published wheel/release, completed hosted CI or speed
measurement is represented by these records.
