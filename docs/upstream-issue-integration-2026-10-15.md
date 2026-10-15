# Upstream issue candidate integration — 2026-10-15

This change integrates the previously isolated RDKit, Open Babel and COSMolKit
source candidates into current main, based on `15150957`. It retains the
subsequent v1.1.4 public-error compatibility fixes, coverage workflow, pinned
ring-version contract and rustls update.

The integrated APIs cover checked atom permutations and enhanced stereo groups,
query-H folding and correlated alternatives; bounded graph paths, atom-radius
environments, Layered fingerprints and native distance matrices; corrected
fingerprint target-density folding; indexed SDF selections; mmCIF/PDB identity,
PQR/PDBQT preservation, GRO, TURBOMOLE, RDFile, GZMAT, EXYZ, FHI-aims, MOPAC,
DMol3 and the strict molecular ChemDoodle JSON profile.

All production feature sources match the validated final candidate. Three
configuration/document files merge current-main changes: CHANGELOG, the CI
integration matrix and narrow whitespace attributes for archived patches.
Archived patch bytes and historical oracle records remain unchanged.

The dated per-feature reports describe their original isolated source lanes.
Their isolation, commit and runtime fields are historical evidence, not claims
about the current integration or a published package. The new integration
manifest is `validation/results/upstream-issue-integration-2026-10-15.json`.
It records source hashes and current-main preservation separately.

Local formatting, binding inventory, refactor-budget, release-documentation,
strict MkDocs and staged-diff checks pass. A fresh integrated Python extension
passes 1032 tests across 27 added suites and two existing adapter/preservation
suites. All 27 added Rust integration targets also pass (303 tests). The first full-workspace local test
attempt was interrupted during an existing long CIP oracle and is not counted
as a complete pass. GitHub CI must pass the protected main checks, all added
integration targets and binding/coverage jobs before merging. No package
publication or upstream RDKit/Open Babel change is included.

A Windows source-wheel run exposed platform-default decoding in new stub-reader
tests. Explicit UTF-8 reads correct that test portability issue; production
feature code is unchanged. The related local suites pass 354 tests.
