# CIF-to-XSF geometry evidence — 2026-10-17

Baseline `35b932f9` has no `write_xsf_checked` / `PeriodicStructure.to_xsf`.
Python files were verified against Git; runtime probes use the reused native
artifact. Full hashes, source URLs and commands are in `provenance.json`.

- `pytest.log`: 198 passed, zero failures/skips: 17 new / 19 crystal / 162
  compatibility tests. Covers geometry/no wrapping, source preservation,
  metadata consent, occupancy/symmetry rejection, overflow and the site cap.
- `crossread.json`: ASE 3.26.0 read 60 cases / 1,360 atom rows with cell and
  Cartesian arrays matching native getter bytes. Additional CIF-pipeline
  and 118-element checks passed. Open Babel 3.2.1 read 48 cases; XYZ output
  error at most 4.984001341412636e-6 A, within 5.1e-6 rounding tolerance.
  Open Babel's coordinate-only XYZ output is not a cell-parity check.
- Open Babel 3.2.1 lists `xsf` as read-only; CIF-to-XSF attempts fail with
  `Invalid output format`. Original probe input is checked in.
- `typing.json`, `strict-mypy.ini` and logs: mypy 1.8.0 targeted strict
  helper passes, valid consumer passes, four invalid cases yield four errors.
  A direct strict run found 81 pre-existing generic-annotation diagnostics
  in the legacy root stub. The targeted check silences imported-stub
  diagnostics (`follow_imports_for_stubs=True`) and does not establish
  package-wide strict cleanliness. Existing native symmetry properties and
  the from_cif expansion flag now have accurate public stub signatures.
- `source-gates.json`: refactor/source-stub/native-inventory/cross-binding/
  unsafe checks pass. No new native function (213 Python / 317 WASM remain).
- `mkdocs.log`: strict build passed; existing unlisted pages informational.
  `product-name.log`: canonical spelling passed.

The helper calls existing native crystal getters, including NumPy arrays;
it is not independent of the native extension. The extension is reused,
not freshly built. Source/native runtime crystal behavior is unchanged.
XSF exports approved geometry only; no label/metadata round trip, GUI run,
occupancy collapse, symmetry search, conventional/primitive reduction,
publication or performance result is represented by these records.
The XCrySDen source host timed out; its distribution documentation mirror
was fetched and hashed. Implementation and synthetic input are original
project work under MIT OR Apache-2.0.
