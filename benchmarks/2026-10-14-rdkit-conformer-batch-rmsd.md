# RDKit 2026.09.1 batch conformer RMSD follow-up

Date: 2026-10-14

## Scope

RDKit 2026.09.1 added `GetAllConformerBestRMSToRef`, which returns the
best-fit RMSD for every probe/reference conformer pair. chematic already had
the single-pair RDKit-compatible alignment engine but no batch API.

The source candidate adds:

- Rust `rdkit_all_conformer_best_rms` and
  `rdkit_all_conformer_best_rms_to_ref`;
- Python `Mol.rdkit_all_conformer_best_rms` and
  `Mol.rdkit_all_conformer_best_rms_to_reference`;
- matching WASM JSON methods.

Both paths enumerate atom mappings once and reuse them. The same-molecule
result order is `(1,0), (2,0), (2,1), ...`; the cross-molecule order is
reference-major, matching RDKit's 2026.09.1 implementation.

## Provenance

- base: `4be080328d666af9288fb3f5bfc256751efac9ed` (`origin/main`)
- implementation: `1691ee5d1814a939c3c97ab2d1572df553bd3d77`
- official comparator source: RDKit tag `Release_2026_09_1`, PR #9558
- local executable comparator: RDKit Python `2025.09.3`
- Python: CPython 3.13.6, macOS arm64
- candidate wheel SHA-256:
  `2e8b24b028feb0cfeda0a776bd001d5527ea6a090acb242d7848ddc0a526e758`

The local RDKit version predates `GetAllConformerBestRMSToRef`. Its existing
`GetAllConformerBestRMS` was used directly for the same-molecule lane. The
cross-molecule lane was checked against the Cartesian sequence of RDKit
`GetBestRMS` calls in the result order defined by the 2026.09.1 source.

## Results

The three same-molecule pair values are bit-identical to RDKit. The four
cross-molecule values differ by at most `3.3306690738754696e-16`. The focused
Python regression suite has 3 passing tests; the complete `chematic-smiles`
and `chematic-wasm` suites have 380 and 391 passing tests respectively, with
one unrelated ignored test in `chematic-smiles`. Clippy passes for the Rust,
Python and WASM binding crates.

## Boundary

This record establishes API ordering, atom-map reuse, explicit error handling,
and focused numeric agreement. It is not a broad conformer-quality result, a
published-package result, or a claim that all RDKit 2026.09.1 features are
implemented. ETKDGv4, internal-coordinate embedding, synthon-space shape
search, and nanobind are separate projects with substantially different scope.
