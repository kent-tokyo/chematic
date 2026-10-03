# chematic and RDKit

CheMatic is a local-first Rust chemistry library for Rust, Python and
WebAssembly. RDKit has a broader, more mature chemistry and Python ecosystem;
it is the better choice when a workflow needs an unsupported API or production
3D breadth. CheMatic is not a drop-in reimplementation.

For individual API mappings, see [migration](rdkit-migration.md). For actual
support and refusal boundaries, see [compatibility scope](compatibility-scope.md).

## Current measured boundary

The current release is **v1.0.31**. Release channels are verified, but the
latest published-package chemistry and speed comparison packet is pinned to
**v1.0.30** against RDKit 2026.03.6. The v1.0.31 WASM formula and reaction
JSON fixes still need a package-output rerun.

| Lane | Recorded result | What it does not show |
|---|---|---|
| Python output audit | 10,000 chemistry rows, 310,000 SMARTS cells and the original 57 reaction fixtures accounted for; 200 SMARTS cells differ and five CIP rows abstain with typed reasons. | Full SMARTS, CIP or reaction compatibility. |
| RDKit-compatible HBA | Published v1.0.30 macOS arm64 wheel: 5,000/5,000 exposed ChEMBL rows match RDKit; v1.0.29 matched 3,641/5,000. | Other descriptors or platforms. |
| Python timing | Of 63 operations, 20 meet exact-output and favorable paired-interval gates in 20 alternating blocks on one host. | A universal speed lead; output-mismatched operations are not wins. |
| Browser Morgan | Published v1.0.30 npm/WASM on Ubuntu 24.04: Chromium, Firefox and WebKit each match 250/250 direct and prepared rows; 20-block speed intervals favor chematic on the measured lane. | Other browser hosts, operations or library-only memory. |
| Reactions | Published v1.0.30 Rust graph/origin/template-map identity matches RDKit on 73/83 stratified rows. An unpublished checked-source Python gate has 76 all-axis matches, three typed unsupported, one diagnosed refusal and three invalid in both. | General SMIRKS parity, yield or selectivity prediction; source results are not published-package results. |
| 3D/MMFF94 | Experimental. Published v1.0.31 macOS: 265/265 geometry/stereo/clash, 100/265 converged and 262 comparable same-coordinate energies within 1 kcal/mol. | Cross-platform or independent conformer quality; an MMFF94 speed win. |

The [validation report](validation.md) links the exact artifacts, corpus
hashes, failure counts and operation definitions. The [benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks)
keeps historical results versioned; an older speed or size figure is not a
current-package result.

## Practical differences

- **Implementation and deployment:** CheMatic's Rust core is shared by
  Rust, Python, Node and WASM. RDKit's official JavaScript/WASM package is
  available too, so browser availability alone is not a differentiator.
- **Profiles:** Native ECFP4 and named RDKit-compatible Morgan are different
  definitions. Aromaticity, CIP and canonical identity also have explicit
  modes or fail-closed limits. Matching names do not imply identical defaults.
- **3D:** RDKit's ETKDG and force-field workflows are more mature. CheMatic
  exposes bounded 3D paths with typed failures but does not claim ETKDGv3
  or full MMFF94 parity.
- **Formats:** CheMatic has selected materials and simulation formats in
  addition to common molecule formats. Read/write and binding coverage vary;
  use the [format matrix](format-capabilities.md) rather than assuming every
  format round-trips.
- **Agent and browser workflows:** CheMatic emphasizes local execution,
  typed errors, and input accounting. Verify cancellation, memory and
  unsupported behavior against the specific binding you deploy.

For a production migration, pin both versions and options, rebuild persisted
fingerprint indexes under one chosen profile, and compare semantic outputs
plus refusals on representative held-out data. A canonical SMILES string is
not a universal identity key; `canonical_smiles_stable_key()` is intentionally
fail-closed.
