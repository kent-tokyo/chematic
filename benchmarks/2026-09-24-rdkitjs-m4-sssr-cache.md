# M4 Chrome source-candidate RDKit.js comparison (2026-09-24)

## Scope

This is a source-candidate measurement of the SSSR/cache/VF2 performance
series. It is not evidence for a published package and does not establish
native Python performance across the broader operation matrix.

The candidate is `3650e255e6bb978e9cc3ab2a6e1752ab0d055be1`, built from
`origin/main` `01a86b621f401d22fc9a0e90f85e4c85b6dd49fe` plus the nine
performance commits. The exact candidate diff digest is
`18d38c9841d4626dae9b43cbf0137ec2cda6e78ee5fcb2af92a5ef88f967f8e2`.

## Protocol

- Host: Apple silicon macOS (`arm64`), Rust 1.97.0, wasm-pack 0.13.1,
  Node 24.5.0, Python 3.13.6.
- Browser: Google Chrome 153.0.8010.53, one fresh browser process per arm and
  repetition.
- Comparator: official `@rdkit/rdkit@2026.03.6`.
- Corpus: the first 5,000 rows of `scripts/descriptor_census_corpus.smi`,
  SHA-256 `d6f2ba3f128296f935007f0b0813aa97b6ebcc2457e014ddca2213ddd655276c`.
- Five alternating chematic/RDKit repetitions; each arm warms 20 rows.
- The fingerprint path parses, performs compatible preparation, generates a
  radius-2/2,048-bit fingerprint, and consumes the same packed 256-byte
  LSB-first representation. The prepared lane excludes construction from the
  timed interval for both engines.

The runner rejects a measurement if the per-repetition packed-output checksum
or bit count differs. It completed successfully for all five pairs. No typed
fingerprint refusal occurred in this 5,000-row slice.

## Result

| Operation | chematic mean (ms/mol) | RDKit mean (ms/mol) | Geometric speedup | 95% lower bound |
| --- | ---: | ---: | ---: | ---: |
| SMILES parse | 0.003212 | 0.096348 | 29.996x | — |
| Parse + canonical-SMILES write | 0.054140 | 0.129000 | 2.383x | — |
| Parse + compatible Morgan | 0.027584 | 0.125056 | 4.535x | 4.411x |
| Prepared compatible Morgan | 0.008220 | 0.028036 | 3.411x | 3.321x |

The performance gate passes for the two compatible-Morgan lanes because each
paired 95% lower bound exceeds 1.0. The record is deliberately limited to this
host, browser, corpus, options, and source candidate. It must be repeated from
a registry-installed artifact before being presented as a released-package
claim.

The candidate WebAssembly binary is 4,230,157 bytes (1,548,598 bytes gzip-9;
SHA-256 `ad74c98fe6ee0baab151bd083a02223e5ba1fe4098b6181cecc4a78eb48f06f9`).
The matching machine-readable raw runs, package provenance, and RSS sampling
qualification are in
[`validation/results/rdkitjs-m4-source-candidate-sssr-cache-2026-09-24.json`](../validation/results/rdkitjs-m4-source-candidate-sssr-cache-2026-09-24.json).

## Boundary

This record does **not** verify the reported native-Python 62-operation
matrix: its driver and raw `rdkit_speed_matrix_*` outputs were not present in
the Git worktree. Those claims remain unadopted until their script, corpus
provenance, raw runs, and equivalent-operation definitions are committed and
rerun without concurrent work.
