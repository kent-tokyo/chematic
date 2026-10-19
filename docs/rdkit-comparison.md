# chematic and RDKit

chematic is a local-first Rust chemistry library for Rust, Python and
WebAssembly. RDKit has a broader, more mature chemistry and Python ecosystem;
it is the better choice when a workflow needs an unsupported API or production
3D breadth. chematic is not a drop-in reimplementation.

For individual API mappings, see [migration](rdkit-migration.md). For actual
support and refusal boundaries, see [compatibility scope](compatibility-scope.md).

## Current measured boundary

The **v1.1.9 release** retains the named RDKit 2026.03.1-compatible surface,
adds fail-closed canonical-SMILES stability checks, and preserves the v1.1.8
Labute ASA and ring-count behavior. New v1.1.9 performance changes remain
source-only until exact public artifacts are measured. In the recorded v1.0.41 source comparison, chematic matches RDKit at
least as often as COSMolKit 0.5.0rc15 on all 161 operations over ChEMBL 5k and
RDKit.js 10k; the separate stress corpus retains one weaker native canonical
SMILES round-trip result. This is source-build evidence, not a universal or
published-package parity claim. The latest
published-package chemistry and speed comparison packet is pinned to
**v1.0.42** against RDKit 2026.09.1. A separate v1.1.3 Python artifact packet
covers exact public wheels on hosted Linux x86-64, macOS arm64, macOS x86-64,
and Windows x86-64. The v1.1.7 hosted published-wheel chemistry run also
completed all 11 jobs, but its durable per-row packet is not checked in yet.
These are output-validation results, not evidence for npm/WASM, crates.io, or
v1.1.9 package performance.

| Lane | Recorded result | What it does not show |
|---|---|---|
| [v1.1.3 published Python wheels](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-v113-published-artifact-packet.md) | Four exact PyPI wheels reproduce the expected reaction, SMARTS, Morgan, descriptor, CIP, hybridization, MMFF94, MOL stereo, BioTransformer and A6 counts. Morgan is 10,000/10,000 exact on exposed rows and 5,000/5,000 on ChEMBL; A6 is 265/265 independently sound and stereo-clean. | npm/WASM or crates.io chemistry; native per-platform RDKit execution where no matching RDKit wheel exists; universal parity or speed superiority. |
| v1.1.3 unusual-input follow-up | Source candidate `90542885`: MW, Morgan and ring count match 4,072/4,072; every finite comparable Labute ASA value agrees within `1e-9`. CIP has 4,062/4,072 exact rows, with the ten differences classified rather than guessed. | Published artifacts, other corpora, or universal descriptor/CIP parity. |
| Python output audit | 10,000 chemistry rows, 310,000 SMARTS cells and the original 57 reaction fixtures accounted for; 200 SMARTS cells differ and five CIP rows abstain with typed reasons. | Full SMARTS, CIP or reaction compatibility. |
| RDKit-compatible HBA | Published v1.0.30 macOS arm64 wheel: 5,000/5,000 exposed ChEMBL rows match RDKit; v1.0.29 matched 3,641/5,000. | Other descriptors or platforms. |
| Python timing | Of 63 operations, 20 meet exact-output and favorable paired-interval gates in 20 alternating blocks on one host. | A universal speed lead; output-mismatched operations are not wins. |
| Browser Morgan | Published v1.0.30 npm/WASM on Ubuntu 24.04: Chromium, Firefox and WebKit each match 250/250 direct and prepared rows; 20-block speed intervals favor chematic on the measured lane. | Other browser hosts, operations or library-only memory. |
| Reactions | Published v1.0.38: 80/83 exact graph/origin/map rows and three inputs invalid in both engines. | General SMIRKS parity, yield or selectivity prediction. |
| 3D/MMFF94 | Experimental. The hash-pinned published v1.0.42 macOS arm64 wheel and current cross-platform CI pass the 265-row geometry/stereo/clash gate; source-only typing, energy, and conformer records are separate. The current source candidate's equal-count 20-conformer lane is within 1 kcal/mol of RDKit's best on 258/265 rows. | Seven basin/minimization-tail rows, universal convergence or conformer quality, and a public-artifact MMFF94 speed win. |
| Seeded ETKDGv3 | v1.0.40 source lane on Linux x86-64: 4,977/4,977 ChEMBL and 9,947/9,947 RDKit.js successful rows match RDKit 2026.03.1 coordinates exactly. The current source candidate also reproduces RDKit 2026.09.1 raw and triangle-smoothed topology bounds within `1e-12` on all 768,664 A6 cells and includes its forced amide/ester torsions. | Full 2026.09 coordinate parity, other platforms, broad topologies, or published artifacts. |

The [validation report](validation.md) links the exact artifacts, corpus
hashes, failure counts and operation definitions. The [benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks)
keeps historical results versioned; an older speed or size figure is not a
current-package result.

The [2026-10-14 sanitized-model reuse packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-python-sanitized-view-reuse.md)
records the source-wheel descriptor/fingerprint follow-up against RDKit
2026.09.1 and COSMolKit rc22. Its two fixed slices retain first-use, repeated-use
and actual parse-inclusive results, with output agreement required for ranking.
It is separate from the published-package matrix above.

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
