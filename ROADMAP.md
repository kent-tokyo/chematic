# chematic roadmap

> Release line: **v1.2.5**. The checked Cleanup API reports ambiguous metal transformations.

chematic shares one Rust chemistry core across Rust, Python, Node and WASM. The priorities are correctness, explicit limits and reproducible evidence.

## Current position

- Selected SMILES/SMARTS, descriptor, fingerprint, search and MOL/SDF paths are stable. 3D/MMFF94, pKa/ADMET, IUPAC, broad Markush/polymer support and rich CDXML editing remain bounded or experimental.
- Focused aromatic MCS has a source fix; broad FMCS parity remains open. The earlier frozen v2 holdout is historical candidate evidence.
- Cross-format meaning preservation is checked by format-specific contracts. Unsupported or lossy conversions report refusal or explicit loss.

See [validation](docs/reference/validation.md), the [benchmark index](benchmarks/README.md) and the [open-work ledger](docs/reference/roadmap-open-work.md).

## Priority order

1. **Release integrity:** run Rust, Python, WASM, compatibility, security and package gates on the exact tag; verify crates.io, npm, PyPI and GitHub artifacts independently.
2. **Correctness:** classify Cleanup, MCS, CTfile, stereo, reaction, SMARTS, R-group and interchange residuals. Keep typed refusals visible.
3. **Equivalent-work performance:** separate parsing, first use, prepared reuse, startup, memory and batch work. Use paired blocks and a second host for broad claims.
4. **Workflow and 3D quality:** assess graph and document semantics, 3D typing, geometry, energy, convergence, stereo, clashes, conformer quality, time and memory on unused inputs.

## Product phases

| Phase | Focus | State |
|---|---|---|
| P0 | Release integrity and evidence | Active |
| P1 | Parsing, reactions, SMARTS, descriptors, interchange | Maintenance |
| P2 | Fingerprints, search, stereo, identity | Active |
| P3 | 2D/3D, force fields, conformer quality | Experimental |
| P4 | Bindings, CLI, MCP, browser products | Active |
| P5 | Cross-platform publication and hardening | Active |
| P6 | Advanced chemistry breadth | Selective |

## Accuracy packages

| Package | State |
|---|---|
| A0 Evaluation contract | Complete |
| A1 Perception and descriptors | Active |
| A2 Stereo and identity | Active |
| A3 Fingerprints and retrieval | Maintenance |
| A4 Workflows and interchange | Active |
| A5 Independent adjudication | External |
| A6 3D and force fields | Experimental |

## Evidence rules

Record versions, hashes, corpus, options, host, failure policy and every refusal or skipped row. A source pass is not a package result; coverage is not chemistry accuracy. Once inspected, a sealed dataset becomes development data. See the [benchmark method](docs/reference/benchmark.md) and [compatibility scope](docs/reference/compatibility-scope.md).
