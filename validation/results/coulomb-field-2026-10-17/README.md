# Typed explicit-charge Coulomb field — local evidence

Date: 2026-10-17. Baseline: `35b932f9` (full SHA in `provenance.json`).
The exact baseline Python files were checked against Git and imported with
the reused native extension: `hasattr(rdkit_compat, "rdMIF")` was false.

- `pytest.log`: 225 passed, zero failures/skips (63 new, 162 existing).
- `rdkit-comparison.json`: 864 seeded constructors and 120,960 ordinary-point
  scalar comparisons passed against RDKit 2026.09.1. Absolute error at most
  4.656612873077393e-10 kJ/mol, tolerance 2e-12 relative / 2e-9 absolute.
- Separate rounded-boundary diagnostics: 88,128 evaluations, 108 differences
  above tolerance, each within eight ULP of a site boundary. These are
  **not passing parity cases**. The first comparison initially failed on one
  such boundary; the final harness separates these diagnostics and still
  fails for ordinary-point discrepancies or differences away from boundaries.
- `typing.json` and logs: mypy 1.8.0, strict helper clean, public valid
  consumer clean, four invalid consumer calls give four expected errors.
  Consumers use `--follow-imports=silent --ignore-missing-imports`; this
  does not establish strict cleanliness of the legacy package. The helper
  was copied unchanged to a standalone file to isolate its strict check.
- `source-gates.json` and logs: refactor budget, Python source stub gate,
  native binding inventory, cross-binding manifest and unsafe boundary pass.
- `mkdocs.log`: strict documentation build passed; unlisted existing pages
  are reported as informational messages.
- `product-name.log`: tracked public text uses canonical product spelling.

`provenance.json` records upstream metadata, source URL/hash, candidate-file
hashes, commands and the extension hash/build commit. The extension is reused
for package import and old tests, not freshly built. The new helper itself
uses only standard-library Python and has no native calls. Existing native
binding inventory is unchanged. No wheel, release, new performance result or
completed hosted CI is represented by this local record.

Reproduction requires a built chematic package carrying the candidate Python
files plus RDKit 2026.09.1 and NumPy on PYTHONPATH. Run
`validation/validate_coulomb_field.py --output /tmp/coulomb-comparison.json`.
The synthetic case generator is seeded with 9601 and records a canonical
case-parameter SHA-256. API contracts and deliberate refusals are in
`docs/coulomb-interaction-field.md`.
