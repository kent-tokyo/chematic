# chematic-wasm

WebAssembly bindings for [chematic](https://github.com/kent-tokyo/chematic), a pure-Rust cheminformatics library.

Published to npm as [`@kent-tokyo/chematic`](https://www.npmjs.com/package/@kent-tokyo/chematic).

Use the [binding guide](../../docs/reference/language-bindings.md) for
data and error conventions and the [format matrix](../../docs/reference/format-capabilities.md)
for supported readers and writers. See [Cleanup migration](../../docs/reference/cleanup-migration.md)
for refusal codes in the unreleased source.

## Installation

```sh
npm install @kent-tokyo/chematic
```

## Features

- Parse SMILES and common molecule files, calculate descriptors and
  fingerprints, search substructures, and draw SVG.
- Process documented SDF, XYZ and Extended XYZ batches with indexed results and
  explicit partial status.
- Run bounded reaction, R-group and RDKit-compatible workflows.
- Use the opt-in [3D pipeline](#3d-embedding-embed_pipeline_v2_json) with
  per-stage diagnostics.

## Usage

The generated handle declarations include `Symbol.dispose` for explicit
resource release. TypeScript consumers should include `ESNext.Disposable` in
their `tsconfig.json` `lib` list (in addition to their normal browser or Node
libraries):

```json
{
  "compilerOptions": {
    "lib": ["ES2022", "DOM", "ESNext.Disposable"]
  }
}
```

```js
import init, {
  parse_smiles,
  add_hydrogens,
  tanimoto_ecfp4,
  tanimoto_atom_pair,
  tanimoto_torsion,
  brics_fragment_count,
  gasteiger_charges_json,
  slogp_vsa_json,
  smr_vsa_json,
  peoe_vsa_json,
  identify_functional_groups,
  brenk_report_json,
} from '@kent-tokyo/chematic';

await init();

const mol = parse_smiles('CC(=O)Oc1ccccc1C(=O)O'); // aspirin

// Descriptors
console.log(mol.atom_count());          // 13
console.log(mol.molecular_weight());    // ~180.16
console.log(mol.formula());             // "C9H8O4"
console.log(mol.tpsa());               // ~63.6
console.log(mol.logp_crippen());        // ~1.2
console.log(mol.fsp3());               // ~0.111
console.log(mol.qed());                // drug-likeness score [0, 1]
console.log(mol.exact_mass());         // ~180.042
console.log(mol.hbd_count());          // 1
console.log(mol.hba_count());          // 4
console.log(mol.rotatable_bond_count()); // 2 (RDKit Lipinski definition)
console.log(mol.aromatic_ring_count()); // 1
console.log(mol.lipinski_passes());     // true
console.log(mol.canonical_smiles());    // canonical SMILES string

// Brenk alerts, fully local in the browser/Node runtime
console.log(JSON.parse(brenk_report_json('Oc1ccccc1')));
// { passes: false, alerts: ['phenol'] }

// BRICS fragmentation
console.log(brics_fragment_count(mol)); // ≥ 2

// Fingerprint similarity
const caffeine = parse_smiles('Cn1cnc2c1c(=O)n(c(=O)n2C)C');
console.log(tanimoto_ecfp4(mol, caffeine));    // ECFP4 Tanimoto
console.log(tanimoto_atom_pair(mol, caffeine)); // AtomPair Tanimoto
console.log(tanimoto_torsion(mol, caffeine));   // Torsion Tanimoto
```

### RDKit-compatible interoperability

The `rdkit_*` methods are opt-in compatibility paths and do not replace the
native browser APIs. JSON atom indices and coordinate rows are zero-based in
`MolHandle` atom order:

```js
const smarts = mol.rdkit_smarts(true);
const centers = JSON.parse(mol.rdkit_chiral_centers_json(true));
const bitInfo = JSON.parse(mol.rdkit_morgan_bit_info_json(2, 2048, false));

const withH = add_hydrogens(parse_smiles('CCO'));
const coords = JSON.parse(withH.rdkit_embed_json(42, 0));
const bounds = JSON.parse(withH.rdkit_bounds_matrix_json(true, true, false));
```

PDB, XYZ, and MOL2 readers return `RdkitReadResult`; call `molecule()` for a
handle and `coords_json()` for index-aligned coordinates. Browser limits are
explicit: at most 512 atoms for ETKDG/bounds, 4,096 stereoisomers, and
1,000,000 alignment matches. Invalid or unsupported inputs throw a stable JS
error instead of being silently coerced.

### Node.js

The published package is built with wasm-pack's `web` target. In a browser,
`await init()` locates the adjacent WASM asset. Node does not fetch `file:`
URLs, so pass the asset bytes explicitly:

```js
import { readFile } from 'node:fs/promises';
import init, { parse_smiles } from '@kent-tokyo/chematic';

const wasm = await readFile(new URL(
  './node_modules/@kent-tokyo/chematic/chematic_wasm_bg.wasm',
  import.meta.url,
));
await init({ module_or_path: wasm });
const mol = parse_smiles('c1ccccc1');
console.log(mol.formula()); // C6H6
mol.free();
```

### 3D embedding (`embed_pipeline_v2_json`)

Opt-in — does not change behavior of any existing 3D API (`generate_coords`,
`generate_and_minimize_*`, etc.), which remain the defaults.

```js
const response = JSON.parse(embed_pipeline_v2_json(mol, JSON.stringify({
  embedSeed: 7,
  maxAttempts: 8,
  embedTimeoutMs: null,
  useExpTorsions: false,
  useSmallRingTorsions: false,
  useMacrocycleTorsions: false,
  useMacrocycle14Bounds: false,
  includeLegacyTorsionHeuristic: false,
  stereoPolicy: "ignore",
  failOnUnevaluableStereo: false,
  forceFieldPolicy: "none",
  forceFieldMaxIterations: 200,
  gateMmff94TorsionOop: false,
  gateMmff94StretchBend: false,
  ringTorsionPolicy: "fail_closed",
  totalTimeoutMs: null,
  enforceChirality: false,
  expandImplicitHThroughPipeline: false,
})));
// response.ok, response.result / response.error — same shape as
// Mol.embed_pipeline_v2() in the Python binding.
```

For ring-fused declared stereocenters (e.g. testosterone, cholesterol) that
`enforceChirality` alone can't repair, `stereoPolicy: "repair_and_verify"` +
`enforceChirality: true` + `expandImplicitHThroughPipeline: true` are needed
*together* (issue #291/#383) — `pipeline_v2_stereo_safe_config_json` builds
that exact combination so a caller can't set one and forget another:

```js
const configResponse = JSON.parse(pipeline_v2_stereo_safe_config_json(
  "mmff94_with_uff_fallback", "fail_closed"
));
// configResponse.ok, configResponse.config (or configResponse.error, same
// shape as embed_pipeline_v2_json's own error) — modify configResponse.config
// (e.g. a different embedSeed) before passing it to embed_pipeline_v2_json.
const response = JSON.parse(
  embed_pipeline_v2_json(mol, JSON.stringify(configResponse.config))
);
```

Verified working end-to-end under real WASM (both `wasm-pack --target nodejs`
and `--target web`, the latter being what this package's npm build actually
uses) — success, typed-failure, and typed-timeout paths all return real
results. `embedTimeoutMs`/`totalTimeoutMs` use a monotonic clock that's
portable across native and `wasm32-unknown-unknown`
([`web-time`](https://crates.io/crates/web-time), backed by
`Performance.now()` in the browser); this does not claim identical wall-clock
precision across every JS engine, only that the value is finite, non-negative,
and enforced correctly on all of them.

`nearest_neighbors_json` keeps its historical chematic-native ECFP4 profile.
For the separately named RDKit-compatible Morgan profile, use
`rdkit_nearest_neighbors_json(querySmiles, dbSmilesJson, k)`. It returns the
same `{index, tanimoto}` shape and reports preprocessing failures without
silently falling back to native ECFP4.

For repeated queries, construct `new RdkitSearchIndex(dbSmilesJson)` once and
call `index.search_json(querySmiles, k)`. The prepared index is intended for
chunked libraries up to the WASM batch limit and applies the same fail-closed
RDKit-compatible profile without rebuilding database fingerprints per query.

## V3000 SGROUP syntax view

`v3000_sgroups_json(block)` exposes bounded, typed SGROUP syntax without
expanding polymer or Markush semantics. It preserves source order for unknown
attributes and returns `kindToken` for unknown group kinds. Group IDs, parent
references, atom references, and grouped-field counts are validated before
JSON is returned.

```js
const groups = JSON.parse(v3000_sgroups_json(v3000Block));
// [{ id, kind, parentId, atomIds, attributes, kindToken? }]
```

This is a syntax-level API; it does not claim polymer expansion, Markush
interpretation, or cross-engine semantic compatibility.

## Bundle size

Bundle size depends on the release, enabled features, and toolchain. Measure
the installed npm artifact for the version you deploy. The retained
[three-browser packet](../../benchmarks/2026-10-03-v1030-published-linux-three-browser-paired20.md)
documents the bounded browser performance comparison; it is not a current
bundle-size claim.

PNG rasterization (`tiny_skia`) is excluded from the WASM build — use SVG output instead. All SVG depiction APIs remain fully available.

## Versioned document binding boundary

The `*_v1` document APIs provide a stable JSON boundary for downstream editors:

```js
const parsed = JSON.parse(reaction_document_json_v1(JSON.stringify(document)));
const edited = JSON.parse(edit_reaction_document_json_v1(
  JSON.stringify(parsed),
  JSON.stringify({ kind: "set_step_condition", step_id: "step-1", key: "temperature", value: "25 C" }),
));

const cdxmlEnvelope = JSON.parse(cdxml_document_json_v1(cdxml));
const cdxmlAgain = cdxml_document_from_json_v1(JSON.stringify(cdxmlEnvelope));
```

`cdxml_document_json_v1` retains the exact `source` string and returns a
structural `document` summary with opaque objects and `diagnostics`. Use
`edit_cdxml_document_json_v1` for bounded page/object edits; it reparses the
result before returning. Errors are JSON-shaped with stable `code`, `path`, and
`message` fields (`malformed_input`, `resource_limit`,
`unsupported_construct`, `lossy_conversion`, or `serialization_error`).
`reaction_document_to_rxn_v1` and `cdxml_document_projection_json_v1` reject
lossy legacy projections with `lossy_conversion` diagnostics. These APIs do not claim
mechanism correctness, product prediction, complete stoichiometry, or full
ChemDraw/RXN compatibility.

## Building from source

```sh
wasm-pack build --target bundler --release
```
