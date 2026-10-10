# Ranking and diagnostic coverage contracts

This follow-up extends the tests merged in #809 (processed PR coverage 94.49%). The coverage workflow and measured file scope are unchanged. The next processed Linux report, not local line-index estimates, determines the percentage; this change alone does not claim 96%.

- Raw canonical rankings: 256 RDKit 2026.03.1 cases cover all six raw double-bond stereo states, alternative controlling substituents, partial fragments, chirality, isotope/map labels, and custom atom symbols. Every active rank is compared exactly. The generator documents a pinned Python wrapper argument-forwarding quirk: `includeChiralPresence=True` selects the C++ `includeAtomMaps=true` contract that the Rust raw fragment API implements. Unselected fragment atoms are masked to match Python's output convention.
- Directional canonical SMILES: the fixture grows from 444 to 646 rows with closed conjugated 12- and 16-atom macrocycles and deterministic E/Z patterns. Serialized SMILES are reparsed before reference canonicalization, so the expected result describes the serialized input, including any stereo RDKit discarded while writing it.
- Public APIs: test bond tokens/valence contracts, coordination mirror involutions, fingerprint alias configuration and vector cosine identities, precomputed ADMET scores, fallible graph descriptor resource errors, analytic carbon-distance graphs, and empty/partial-coordinate depiction with atom metadata and export dimensions.
- Diagnostics: test nucleic-acid error codes/paths and atomic failed edits, CIP comparison traces and edge parentage, and force-field setup errors and symmetry-invariant missing-parameter citations. Raw InChI cleanup tests preserve atom/connectivity identity and compare selected unusual N/S/Cl valence states to explicit product graphs; these are contract examples, not a general chemical correctness claim.
- Legacy BCI charges: test charge conservation and bond-storage-orientation invariance on representative functional groups with and without explicit H. These approximate legacy models are separate from the RDKit MMFF numeric model.
- CML: round-trip every supported bond type and explicit isotope/charge/H annotations; check malformed XML and non-finite coordinate rejection.
- PubChem: extract the existing request into a private lookup closure, keeping the timeout and bounded read unchanged. Deterministic tests check UTF-8 URL encoding, success payloads, malformed/missing response fields, request failure propagation, and argument validation before lookup. Tests make no external requests.
- The ring-order introsort is checked against Rust's integer sorter for ascending, descending, equal-key, sawtooth, organ-pipe, and deterministic pseudo-random input, including the heap fallback.

Regenerate both external fixtures with exactly RDKit 2026.03.1:

```sh
python scripts/generate_ranking_boundary_fixtures.py
python scripts/generate_canonical_direction_boundary_fixtures.py
```

Both regenerated files matched byte for byte. Validate with the existing coverage workflow command, workspace all-target Clippy, rustfmt, and `git diff --check`. Coverage measures execution; it does not establish complete chemical correctness or full RDKit compatibility.
