# chematic and RDKit

chematic is not a drop-in reimplementation of RDKit. It provides a pure-Rust
core, Python wheels, native WASM/Node bindings, and selected RDKit-oriented
APIs. RDKit remains the better choice when broad ecosystem compatibility,
mature 3D workflows, or an unsupported API is required.

See the [migration matrix](rdkit-migration.md) for function-level coverage and
the [compatibility contract](compatibility-scope.md) for stable boundaries.

## At a glance

| Area | chematic | RDKit |
|---|---|---|
| Core implementation | Rust; no C/C++ toolchain on the common path | C++ with Python bindings |
| Browser | Native `wasm32-unknown-unknown` package | Official `@rdkit/rdkit` JavaScript/WASM distribution path |
| Python | Prebuilt `chematic` wheels and a selected RDKit-style subset | Broad, mature reference API |
| Canonical identity | Fail-closed stable-key API; canonical spelling is not a cache-key guarantee | Mature canonicalization with its own conventions |
| Descriptors/fingerprints | Broad native set plus named compatibility modes; parity is metric-specific | Broad reference implementations and ecosystem |
| 3D/MMFF94 | Experimental, typed failure and provenance | Mature ETKDG and force-field workflows |
| Formats/materials | Broad Rust format surface, including several materials/simulation documents | Strong common cheminformatics formats; wider ecosystem integrations |
| Reactions/SMARTS | Bounded implemented subset with typed limits | Broader and more mature chemistry workflow surface |

## Recorded performance

The following are dated macOS arm64 medians, not cross-platform guarantees:

| Operation | chematic | RDKit | Source boundary |
|---|---:|---:|---|
| Canonical SMILES | 24.95 µs/mol | 25.58 µs/mol | v1.0.2 code, 5,000-entry descriptor corpus |
| Canonical SMILES | 18.27 µs/mol | 26.82 µs/mol | v1.0.2 code, independent 5,000-entry corpus |
| SDF graph/property read | 9.48 µs/mol | 99.96 µs/mol | 365 records; chematic graph-only path |
| SDF serialization-only write | 7.62 µs/mol | 79.54 µs/mol | same corpus; 2D layout disabled |
| ECFP4 batch | 54.7 µs/mol | 94.3 µs/mol | historical v0.18.0, 5,000 molecules |
| Parse + compatible Morgan | 3.11–3.67x faster (95% lower bound) | baseline | merged source `7d98dcd3`, official RDKit.js, fixed 10k, GitHub-hosted Chromium/Firefox/WebKit; excludes startup and unsupported Fe(II) row |

The SDF operations are intentionally narrow and do not establish a lead for
layout-enabled writing or every supplier option. Canonical output was produced
for every input, but the two libraries need not choose the same valid spelling.
Exact versions, corpus hashes, commands, and follow-up source A/B results are in
the [benchmark guide](benchmark.md).

## Current compatibility diagnostic

The pinned RDKit 2026.03.6 exposed 10,000-row packet keeps the major surfaces
separate. Public v1.0.19 records 9,770/10,000 correspondence-correct CIP labels,
9,999/10,000 compatible-Morgan rows with one typed Fe(II) refusal, and 14,306
differing SMARTS cells out of 310,000. A post-v1.0.20 source candidate for
issue #632 reduces the 18 SMILES semantic differences to zero while preserving
zero non-isomeric graph differences; the long permutation audit was interrupted,
so it is not a merged or published-package claim. See [validation](validation.md)
for the evidence files and remaining #634/#635 work.

## Accuracy and parity

The opt-in `Mol.rdkit_mw`, `Mol.rdkit_hba`, and
`Mol.rdkit_aromatic_ring_count` profiles use pinned RDKit conventions without
changing the native descriptor defaults. In the v1.0.13 source-built gate,
the current source-built v2 descriptor run matches RDKit for MW and HBA on
5,000/5,000 rows, including the known explicit-isotope row. The current code
uses RDKit-compatible nuclide masses where available and, for an unlisted but
syntactically valid isotope label, RDKit's mass-number fallback. The full per-field
distribution and 12-case structural holdout are recorded in
[`descriptor-rdkit-diagnostics-v1.0.13.json`](https://github.com/kent-tokyo/chematic/blob/main/validation/results/descriptor-rdkit-diagnostics-v1.0.13.json)
and [`descriptor-rdkit-holdout-v1.0.13.json`](https://github.com/kent-tokyo/chematic/blob/main/validation/results/descriptor-rdkit-holdout-v1.0.13.json).

The historical 4,999-molecule descriptor snapshot reports:

- molecular weight: 99.82% within ±0.01 Da;
- HBA, HBD, TPSA, LogP, molar refractivity, Fsp3, and the documented ring
  metrics: 100% at their stated tolerances;
- stereocenter count: 99.96% against the legacy oracle and 98.6% against
  `FindPotentialStereo`;
- opt-in accurate CIP R/S/E/Z labels: 99.64% stable-oracle agreement against
  modern `rdCIPLabeler`; 15 phosphorus rows are excluded from confident output
  as representation-unstable.

These figures cover a subset of the full descriptor surface. Kappa,
Hall-Kier, Bertz, Balaban, BCUT2D, VSA, MQN, SA Score, and other similarly
named values must be treated as chematic-specific unless the validation table
states parity. See [`validation.md`](validation.md).

## Where chematic is a good fit

- Rust-native or browser deployments where a C++ toolchain is undesirable;
- bounded local parsing, 2D analysis, fingerprints, reports, and selected
  materials/simulation formats;
- applications that need typed resource-limit and unsupported outcomes;
- lightweight Python workflows covered by the published binding contract;
- local MCP workflows, with network access limited to the explicit
  PubChem-backed `name_to_smiles` tool.

## Where RDKit is a better fit

- broad Python API and plugin compatibility;
- production-proven ETKDG, conformer, force-field, and advanced stereo use;
- database, workflow-platform, and long-established ecosystem integrations;
- workloads requiring an RDKit behavior that chematic marks partial or
  unsupported;
- exact matching to an RDKit fingerprint, aromaticity, canonicalization, or
  reaction convention outside a named chematic compatibility mode.

## WASM artifact size

The v1.0.15 published-package scorecard measures chematic at **4,005,280 bytes
raw / 1,460,499 bytes gzip** and `@rdkit/rdkit@2026.03.6` at **7,333,095 /
2,379,975 bytes**. The packages have different feature surfaces, so size is a
deployment observation rather than a feature-normalized benchmark. See the
[published-package scorecard](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-16-official-rdkit-js-isolated-browser-v1.0.15.md).

## Interpretation rule

“Faster”, “more accurate”, and “compatible” apply only to the named operation,
version, corpus, configuration, and failure policy. A missing or unsupported
row is not a win, and a microbenchmark is not a replacement claim.
