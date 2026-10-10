# Numerical and format boundary coverage

This follows #810 without changing coverage scope or exclusions. The processed Linux Codecov report determines the percentage.

- ETKDG angle interval and distance restraints: compare known energies to analytic values, analytic gradients to central differences away from interval boundaries, zero gradients at the boundaries, rigid-motion invariance, and guarded coincident-coordinate behavior. UFF inversion contributions use planar and nonplanar geometries for the same energy/gradient check. An empty field is already minimized and leaves coordinates unchanged. Collinear inversion planes are outside the finite-difference contract because their plane normal is undefined.
- Both deterministic-math adapters: test f32/f64 trigonometric identities, inverse functions and atan2 quadrants, exp/log inversion, powers, negative cube roots, overflow/underflow-resistant hypotenuses, IEEE domain errors, and signed zero.
- InChI layer updates: sparse, permuted external atom numbering updates the intended atom's charge/isotope and preserves unmodified attributes, element identities, hydrogen counts, and graph connectivity. These helper contracts do not claim that the pure-Rust parser supports every standard InChI layer encoding.
- RXN document export: collect every unsupported field (multiple steps, document/step provenance, conditions, coefficients), preserve the original document on failure, and round-trip atom-map identities once the unrepresentable fields are removed.
- Depiction geometry: verify grid membership/insertion order, closest-pair tie-breaking, nonbonded stacking counts, reflection involution and rigid-distance preservation, and intersection symmetry with touching/collinear segments excluded.
- Aromaticity traces: charged carbon and heterocycle examples report explicit electron contributions and reasons, including the strict model's exclusion of phosphorus. The trace still reports every ring atom when one is ineligible.
- Diagnostics: stereo ligand errors, overvalenced atom errors, extension lookup/failure propagation, square-planar MOL export, and invalid ensemble RMSD configuration retain the relevant identity/value in their messages.
- MOL2: expand each of the 58 existing Tripos records to all eight combinations of sanitization, hydrogen removal, and Corina cleanup. Compare acceptance, canonical SMILES where sanitized, atom/bond counts, atomic numbers, charges, and finite coordinate counts against exactly RDKit 2026.03.1 (464 cases, 460 accepted by RDKit). The JSON is compact to keep the extra reference data small; input blocks are stored only once per record.

Targeted contract tests and workspace all-target Clippy pass. The existing full coverage command passed 5,320 tests with nine ignored, including native InChI and the loopback HTTP transport tests. After the geometry/trace/diagnostic additions, all seven changed libraries were remeasured with `--no-clean`: 2,218 passed, five ignored. The final report uses the original exclusions and covers 340 files. Its local result is 95.12% (152,045/159,850 lines), not a processed Codecov claim.

Regenerate the MOL2 reference with `python scripts/generate_mol2_boundary_fixtures.py` under RDKit 2026.03.1; the second regeneration matched byte for byte (SHA256 `6d7d9f2d5e894d726e9f8b0b68b05e429b006aeeadeb85a0ed3e019cd84f2248`).

Coverage measures execution. These checks establish the named numerical and boundary contracts, not complete chemistry correctness or general RDKit equivalence.
