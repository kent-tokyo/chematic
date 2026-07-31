# chematic-wasm

WebAssembly bindings for [chematic](https://github.com/kent-tokyo/chematic), a pure-Rust cheminformatics library.

Published to npm as [`@kent-tokyo/chematic`](https://www.npmjs.com/package/@kent-tokyo/chematic).

## Installation

```sh
npm install @kent-tokyo/chematic
```

## Features

- Parse SMILES strings into molecule handles
- Molecular descriptors: MW, TPSA, LogP, Fsp3, QED, exact mass, rotatable bonds, HBD/HBA, aromatic ring count, Labute ASA
- Drug-likeness filters: Lipinski, Veber, Egan, REOS, Ghose
- EState indices (Hall & Kier 1991): per-atom values, sum/max/min
- Gasteiger-Marsili PEOE partial charges: per-heavy-atom charges
- VSA descriptors: SlogP_VSA (×12), SMR_VSA (×10), PEOE_VSA (×14)
- SA score: synthetic accessibility estimate [1, 10]
- Functional group identification (Ertl 2017 IFG)
- Canonical SMILES generation
- ECFP4/6, AtomPair, Torsion, and path fingerprints with Tanimoto similarity
- BRICS fragment count
- SDF/MOL block parsing
- Topological descriptors: Wiener index, Hall-Kier κ, χ connectivity indices, Bertz CT
- Shape descriptors (with 3D coordinates): PMI, NPR, radius of gyration, asphericity
- 2D SVG depiction with CPK colors and atom/bond highlighting
- SVG grid layout for multiple molecules
- Reaction SMILES/SMIRKS parsing and transform
- Add/remove explicit hydrogens
- `embed_pipeline_v2_json`: torsion-knowledge-aware 3D embedding + stereo
  verification/repair + policy-gated force field, mirroring the Python
  `Mol.embed_pipeline_v2()` binding ([known blocker](#3d-embedding-pipeline-v2-blocked-see-issue-219))

## Usage

```js
import init, {
  parse_smiles,
  tanimoto_ecfp4,
  tanimoto_atom_pair,
  tanimoto_torsion,
  brics_fragment_count,
  gasteiger_charges_json,
  slogp_vsa_json,
  smr_vsa_json,
  peoe_vsa_json,
  identify_functional_groups,
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
console.log(mol.rotatable_bond_count()); // 3
console.log(mol.aromatic_ring_count()); // 1
console.log(mol.lipinski_passes());     // true
console.log(mol.canonical_smiles());    // canonical SMILES string

// BRICS fragmentation
console.log(brics_fragment_count(mol)); // ≥ 2

// Fingerprint similarity
const caffeine = parse_smiles('Cn1cnc2c1c(=O)n(c(=O)n2C)C');
console.log(tanimoto_ecfp4(mol, caffeine));    // ECFP4 Tanimoto
console.log(tanimoto_atom_pair(mol, caffeine)); // AtomPair Tanimoto
console.log(tanimoto_torsion(mol, caffeine));   // Torsion Tanimoto
```

```js
// Sprint Q: New descriptors (v0.1.15)
console.log(mol.sa_score());                     // synthetic accessibility [1,10]
console.log(mol.labute_asa());                   // Labute approx. surface area (Å²)

// Gasteiger partial charges (per heavy atom)
const charges = JSON.parse(gasteiger_charges_json(mol));
console.log(charges); // [-0.08, 0.12, -0.43, ...]

// VSA descriptor bins
const slogpVsa = JSON.parse(slogp_vsa_json(mol));
const smrVsa   = JSON.parse(smr_vsa_json(mol));
const peoeVsa  = JSON.parse(peoe_vsa_json(mol));

// Functional group identification
const ifg = JSON.parse(identify_functional_groups(mol));
console.log(ifg); // [{"atoms":[1,2,3],"types":"OC=O"}, ...]
```

### 3D embedding (`embed_pipeline_v2_json`) — blocked, see issue #219

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
  ringTorsionPolicy: "fail_closed",
  totalTimeoutMs: null,
})));
// response.ok, response.result / response.error — same shape as
// Mol.embed_pipeline_v2() in the Python binding.
```

**Not functional yet.** The underlying `chematic-3d` pipeline calls
`std::time::Instant::now()` unconditionally for timing/timeout bookkeeping,
which panics under real `wasm32-unknown-unknown` (no host time source) on
every call. Verified working natively (Rust tests, Rust/Python parity
fixtures) but not yet callable from Node or a browser. Tracked in
[issue #219](https://github.com/kent-tokyo/chematic/issues/219); fix is
scoped to `chematic-3d`, not this crate.

## Version History

**v0.1.94** (2026-06-12):
- SA Score corpus expanded: 188 FDA molecules (1415 unique fragments)
- Enhanced fingerprints: True MHFP, True ERG, path FP with bond types
- Full multi-sphere CIP stereochemistry for R/S assignment
- InChI stereo layer round-trip support (tetrahedral and E/Z)

**v0.1.93** (2026-06-12):
- Full multi-sphere CIP priority rules
- Correct stereochemistry assignment for complex chiral centers

**v0.1.92** (2026-06-12):
- InChI stereo layer parsing (tetrahedral `/t` and E/Z `/b`)
- Path fingerprint with bond type interleaving

**v0.1.91** (2026-06-12):
- True MHFP (structural fragment hashing)
- True ERG (Ertl 2017 functional group detection)

## Bundle Size

~500 KB gzip / ~1.3 MB raw (reduced from ~819 KB gzip in v0.4.17, -38.5%).

PNG rasterization (`tiny_skia`) is excluded from the WASM build — use SVG output instead. All SVG depiction APIs remain fully available.

## Building from source

```sh
wasm-pack build --target bundler --release
```
