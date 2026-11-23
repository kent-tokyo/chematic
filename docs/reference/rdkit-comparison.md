# chematic and RDKit

chematic is a local-first Rust chemistry library for Rust, Python and
WebAssembly. RDKit has a broader, more mature chemistry and Python ecosystem;
it is the better choice when a workflow needs an unsupported API or production
3D breadth. chematic is not a drop-in reimplementation.

For individual API mappings, see [migration](rdkit-migration.md). For actual
support and refusal boundaries, see [compatibility scope](compatibility-scope.md).

## v1.2.2 results

The published Python comparison uses hash-verified wheels for chematic
1.2.0–1.2.2, RDKit 2026.9.2 and COSMolKit 0.6.3 on one macOS arm64 host.
It covers descriptors, fingerprints and five MCS cases. The
[benchmark report](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-21-v122-published-comparison.md)
gives the timing boundaries and per-case output checks.

The v1.2.2 source gate matches RDKit 2026.09.1 on 4,022/4,069 comparable
Cleanup stress rows and 9,996/10,000 broad rows. Reionize matches all
4,072 stress rows. Later Cleanup candidate results and refusals are in the
[Cleanup audit](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-21-cleanup-v122-regression-audit.md).

## Earlier evidence

| Area | Recorded evidence |
|---|---|
| Chemistry and Python bindings | [v1.1.3 four-wheel packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-v113-published-artifact-packet.md) covers selected descriptors, fingerprints, reaction and 3D rows. |
| Reactions | [reaction follow-up](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-10-v1042-regression-followups.md) record exact outputs and invalid inputs. |
| Browser Morgan | [v1.0.30 three-browser packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-linux-three-browser-paired20.md) covers direct and prepared paths. |
| 3D and force fields | [A6 bounds and conformer follow-up](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-12-a6-rdkit-2026-09-bounds-followup.md) and the [conformer evaluation index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md) retain source and package boundaries. |

Each record describes its artifact, input set, options and failures. The
[validation report](validation.md) links the current release and wider index.

## Practical differences

- **Implementation and deployment:** chematic's Rust core is shared by
  Rust, Python, Node and WASM. RDKit's official JavaScript/WASM package is
  available too, so browser availability alone is not a differentiator.
- **Profiles:** Native ECFP4 and named RDKit-compatible Morgan are different
  definitions. Aromaticity, CIP and canonical identity also have explicit
  modes or fail-closed limits. Matching names do not imply identical defaults.
- **3D:** RDKit's ETKDG and force-field workflows remain broader and more
  mature. chematic now exposes a named fixed-seed RDKit-compatible ETKDGv3
  path with the bounded corpus evidence above. That result does not establish
  universal topology, ensemble-quality, platform-bitwise, or full MMFF94
  parity for the general 3D APIs.
- **Formats:** chematic has selected materials and simulation formats in
  addition to common molecule formats. Read/write and binding coverage vary;
  use the [format matrix](format-capabilities.md) rather than assuming every
  format round-trips.
- **Agent and browser workflows:** chematic emphasizes local execution,
  typed errors, and input accounting. Verify cancellation, memory and
  unsupported behavior against the specific binding you deploy.

For a production migration, pin both versions and options, rebuild persisted
fingerprint indexes under one chosen profile, and compare semantic outputs
plus refusals on representative held-out data. A canonical SMILES string is
not a universal identity key; `canonical_smiles_stable_key()` is intentionally
fail-closed.
