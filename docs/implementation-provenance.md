# Implementation provenance and FTO boundary

This document records source and implementation boundaries. It is not legal
advice and does not constitute a patent non-infringement or freedom-to-operate
opinion.

## Independent implementations

- SMILES, MOL/SDF, XYZ, CDXML, RXN, crystal, and binding contracts are
  implemented from the applicable public format behavior and chematic's own
  typed models.
- MMFF94 uses the published Halgren model and cited parameter tables. The
  implementation is independent; the presence of published equations or
  parameters is not treated as patent clearance.
- ETKDG-style distance geometry and torsion knowledge are implemented as
  independent, bounded heuristics. chematic does not ship RDKit ML weights,
  RDKit binaries, or copied source for the default 3D path.
- 3D descriptors identify their literature or standard source in module
  documentation where applicable. No paper text, figures, datasets, or
  trained artifacts are distributed as implementation material.
- Core descriptor provenance and tolerance semantics are checked in at
  `validation/cross_binding_contract.json`; the same fixtures are consumed by
  Rust, Python, and Node/WASM contract tests. Standard format support is
  documented in `docs/format-capabilities.md` and does not imply byte-level
  or proprietary compatibility.
- Core ECFP4 and MACCS fingerprints use the same manifest for packed-byte
  shape, bit order, configuration, sparse/count semantics, and implementation
  provenance. This freezes the binding boundary without claiming RDKit value
  parity; held-out parity and explanations are separate contracts.

## Explicit exclusions

Spectrophores and the Issue #464 geometry-aware spectral fingerprint are not
shipped because their FTO status has not been independently cleared. A
statement that a patent is expired must not be added without a jurisdiction,
patent number, maintenance/expiry evidence, and review record.

## Release rule

Before adding a research-derived algorithm or data to a public release, record
the exact source, license, implementation boundary, patent/FTO review status,
and whether the feature is standard, optional, or excluded. Unreviewed
research features remain private or are rejected from the release surface.
