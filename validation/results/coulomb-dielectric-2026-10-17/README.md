# Explicit-charge CoulombDielectric follow-up — 2026-10-17

Extends PR #995 from `5473b7de`; baseline main `35b932f9`.
The pre-follow-up Python package was extracted from Git and imported with
its reused native extension: `Coulomb` exists, `CoulombDielectric` does not.
Source, reference and extension hashes are in `provenance.json`.

- `pytest.log`: 284 passed, zero failures/skips: 59 dielectric, 63 vacuum,
  162 existing compatibility tests. Analytical tests cover all depth-table
  transitions, threshold-excluded/zero-charge neighbors, exact 4 A boundaries,
  equal permittivities, signs, source copying, repeated/concurrent evaluation
  and the 1,024-site cap checked before indexing or quadratic work.
- `rdkit-comparison.json`: RDKit 2026.09.1, 3,168 constructors, 154,080
  ordinary evaluations and 72 exact/adjacent neighbor-boundary checks pass
  2e-12 relative / 2e-9 kJ/mol absolute tolerance. Maximum absolute error
  9.313225746154785e-10 kJ/mol. Separate 85,536 energy-threshold diagnostics
  have zero differences above tolerance. This is synthetic-profile evidence,
  not a universal boundary-parity claim. Earlier vacuum evidence retains
  its 108 diagnostic differences and is pinned to the previous commit.
- `typing.json`: mypy 1.8.0, isolated two-helper package strict check passes,
  valid public consumer passes, four invalid calls yield four expected errors.
  Consumers use `--follow-imports=silent --ignore-missing-imports`; legacy
  package-wide strict checking is not claimed.
- `source-gates.json` and logs: refactor, source stub, native inventory,
  cross-binding manifest and unsafe gates pass. Inventory remains 213 Python
  functions / 317 WASM exports; this feature is Python-only.
- `mkdocs.log`: strict build passes (existing unlisted pages are informational).
  `product-name.log`: canonical product spelling passes.

The model uses an explicit empirical neighbor-depth table and pinned legacy
RDKit prefactor. It is not a general dielectric solver or new calibration.
Both helper modules use the standard library. The old native extension is
reused for package import/old tests only, not freshly built. No release,
performance result or completed hosted CI is represented by this record.
The reproduction command and selected upstream scope are in provenance;
`docs/coulomb-dielectric-field.md` states the full checked profile.
