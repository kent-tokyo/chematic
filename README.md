# chematic

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/chematic.svg)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic.svg)](https://www.npmjs.com/package/@kent-tokyo/chematic)

[日本語](README_ja.md) | [中文](README_zh.md)

A pure-Rust cheminformatics library targeting RDKit feature parity — **with zero C/C++ dependencies**.

> **Why does zero C/C++ matter?**
> RDKit.js, Indigo WASM, and OpenBabel all ship C++ code compiled via Emscripten.
> That means **30–50 MB WASM binaries**, complex build toolchains, and platform-specific build failures.
> chematic compiles to a **~550 KB WASM bundle** with a single `wasm-pack build` — no `cmake`, no `clang`,
> no `-sys` crates, no `build.rs` C compilation anywhere in the dependency tree.

---

## Live Demo

**[https://kent-tokyo.github.io/chematic/](https://kent-tokyo.github.io/chematic/)** — Interactive descriptor calculator, drug-likeness rules, fingerprint similarity, 3D viewer, and reaction schemes running entirely in your browser via WebAssembly.

---

## Design Goals

**Pure Rust, zero C/C++ FFI — guaranteed**
No `rdkit-sys`, no `openbabel-sys`, no `cc` build dependencies, no `bindgen`. Every
algorithm — from SSSR ring perception to ECFP fingerprints to force-field minimization —
is implemented in 100% safe Rust. The entire dependency tree is verified FFI-free.

**WASM-compatible and lightweight**
All crates compile to `wasm32-unknown-unknown` without modification. The npm package
`@kent-tokyo/chematic` is **~550 KB** versus 30–50 MB for C++ FFI alternatives.
No `cmake`, no `emcc`, no Emscripten toolchain required.

**80+ WebAssembly API endpoints**
The WASM layer exposes 80 functions covering descriptors, fingerprints, scaffold analysis,
stereoisomer enumeration, 3D geometry, diversity selection, and more — all callable from
JavaScript/TypeScript with full TypeScript type definitions.

**Domain-specific algorithms**
Rather than wrapping a generic graph library, chematic implements chemistry-specific
algorithms directly: Kekulization, Hückel aromaticity, CIP stereochemistry, SSSR ring
perception, Gasteiger charges, MaxMin/Butina diversity picking.

**Reproducible and deterministic**
Fingerprints use FNV-1a hashing with a fixed invariant ordering. Given the same SMILES
input, the same bits are always produced. No RNG, no platform-specific behavior.

---

## Current Status

All phases complete + Section 4 (WASM, API improvements) + Sprint v0.1.33 (CXSMILES/CXSMARTS + audit) + Sprint v0.1.34 (InChI ring closure + stereo layers) + Sprint v0.1.35 (wasmBridge support) + Sprint v0.1.36 (Issue #1 Audit: BUG-2/3/4 fix) + Sprint v0.1.37 (mol_transforms API + random SMILES) + **Sprint v0.1.69–v0.1.74 (RDKit Gap Analysis: 6 feature implementations)** + **v0.1.88–v0.1.89 (Gap closure 89%: A1–A6, B1–B2 complete)** + **v0.1.91–v0.1.94 (Gap closure 100%: A1–A5, B3 complete)** + **v0.1.95–v0.1.96: MHFP/ERG/Reaction FP true algorithms + MMFF94 BCI charges** — **1,666 library tests, all passing. Zero C/C++ dependencies.**

Latest release: **v0.1.96** (2026-06-13) — MMFF94 BCI partial charges (±0.1e)

| Crate                 | Description                                                                                              | Tests |
|-----------------------|----------------------------------------------------------------------------------------------------------|-------|
| `chematic-core`       | Atom, Bond, Molecule, Element, kekulization (no deps); mutable `add/remove_atom/bond`, `fragments()`, `is_connected()`, `formula_with_isotopes`, `validate_valence`; `StereoGroup`/`StereoGroupKind` | 48    |
| `chematic-smiles`     | OpenSMILES parser, writer, canonical SMILES                                                              | 57    |
| `chematic-perception` | SSSR, Hückel aromaticity + antiaromaticity (4n+2 rule), `apply_aromaticity`, `aromatize`/`kekulize_inplace`, `assign_stereo_from_2d`, `assign_ez_from_2d`, `cip_ez_descriptor` | 34    |
| `chematic-mol`        | MOL/SDF V2000+V3000 (R/W with 2D coords), CML (R/W), CDXML (R); `SdfRecord` with coords+props; MDL RXN R/W; V3000 stereo-group COLLECTION R/W | 63    |
| `chematic-depict`     | 2D SVG (CPK colors, highlighting, grid), DepictData, `detect_crossings`, `render_svg_with_metadata`, reaction SVG; Y-coordinate system documented | 43    |
| `chematic-chem`       | 70+ descriptors, tautomer scoring, scaffold network, BRICS, QED, standardize, mol_hash, stereo (invert/enumerate), CIP, IFG, VSA (EState+Labute), `parse_condensed`, `isotope_distribution`, `num_amide_bonds`, `num_ester_bonds` | 375   |
| `chematic-fp`         | ECFP2/4/6, FCFP4/6, MACCS 166-bit, TopoPF, AtomPair, Torsion — Tanimoto/Dice                           | 50    |
| `chematic-smarts`     | SMARTS, VF2, MCS with chirality matching (`match_chiral_tag`), atom/bond compare modes; Display + Error trait | 87    |
| `chematic-3d`         | 3D coordinate generation, distance geometry constraints, force-field minimization, shape descriptors, ConformerEnsemble with RMSD pruning, PDB/XYZ; WASM RNG seeded | 147   |
| `chematic-rxn`        | Reaction SMILES/SMIRKS, `find_reaction_center` — `run_reactants` with product valence validation        | 30    |
| `chematic-inchi`      | InChI/InChIKey generation; formula/connectivity/hydrogen/stereo/charge/isotope layers; ring closures   | 28    |
| `chematic-wasm`       | **110+ WASM exports** — npm: `@kent-tokyo/chematic` v0.1.95 (~550 KB); InChI API + stereo inversion     | 175   |
| `chematic-iupac`      | Local IUPAC name generation — pure Rust, no network; alkanes/cycloalkanes, alkenes/alkynes, alcohols, amines, halides, aldehydes, ketones, acids, esters, amides, benzene, aromatic heterocycles | 14    |
| `chematic`            | Umbrella crate with feature flags (all sub-crates, incl. `iupac`, `inchi`)                              | 1     |

```
cargo test --workspace --lib   # 1,666 library tests, all passing
```

---

## Quick Start

### Installation

```bash
# Rust
cargo add chematic --git https://github.com/kent-tokyo/chematic --features "smiles,perception,chem,3d,fp"

# JavaScript/TypeScript
npm install @kent-tokyo/chematic@0.1.94
```

### 5-Minute Examples

#### Parse SMILES & check drug-likeness

```rust
use chematic_smiles::parse;
use chematic_chem::*;

let mol = parse("CC(=O)Oc1ccccc1C(=O)O")?;  // aspirin

println!("MW: {:.2}", molecular_weight(&mol));
println!("LogP: {:.2}", logp(&mol));
println!("TPSA: {:.2}", tpsa(&mol));

if lipinski_descriptor_pass(&mol) {
    println!("✓ Passes Lipinski's Rule of Five");
}
```

#### Detect rings & aromaticity

```rust
use chematic_perception::{find_sssr, assign_aromaticity};

let rings = find_sssr(&mol);
let aromatic = assign_aromaticity(&mol);

println!("Rings: {}", rings.ring_count());
// NEW in v0.1.32: Check for antiaromatic systems
if aromatic.has_antiaromaticity(&mol) {
    println!("⚠ Contains antiaromatic rings (unstable)");
}
```

#### Generate 3D coordinates

```rust
use chematic_3d::generate_and_minimize_constrained;

let coords_3d = generate_and_minimize_constrained(&mol);
// NEW in v0.1.32: Constraint satisfaction for better geometry
```

#### Calculate fingerprint similarity

```rust
use chematic_fp::tanimoto_ecfp4;

let benzene = parse("c1ccccc1")?;
let toluene = parse("Cc1ccccc1")?;
let sim = tanimoto_ecfp4(&benzene, &toluene)?;
println!("Similarity: {:.2}", sim);  // ~0.5
```

#### Preserve chemical metadata with CXSMILES

```rust
use chematic_smiles::parse_cxsmiles;

let cx = parse_cxsmiles("CCO |$ethanol$,atomProp:1.role.acceptor,^2:0|")?;
// cx.atom_labels: ["ethanol"]
// cx.atom_props: [(atom: 1, key: "role", value: "acceptor")]
// cx.atom_radicals: [None, 2, None]
```

#### Audit standardization with reports

```rust
use chematic_chem::{StandardizationPipeline, StandardizeOptions};

let opts = StandardizeOptions {
    largest_fragment_only: true,
    neutralize_charges: true,
    ..Default::default()
};
let pipeline = StandardizationPipeline::new(opts);
let (standardized, report) = pipeline.run(&mol);

println!("Status: {:?}", report.status);  // Unchanged | Modified | CompletedWithWarnings
for step in &report.steps {
    println!("  {}: changed={}", step.step.as_str(), step.changed);
}
```

#### Use from WASM/JavaScript

```javascript
import init, { molecule_report_json, parse_cxsmiles_json } from 'chematic-wasm';

await init();

// Parse CXSMILES with metadata
const cx = JSON.parse(parse_cxsmiles_json("CCO |$ethanol$|"));
console.log(cx.atomLabels);  // ["ethanol"]

// Standardize with audit report
const report = JSON.parse(
    molecule_report_json("CC(=O)Oc1ccccc1C(=O)O")
);
console.log(`LogP: ${report.descriptors.logp}`);
console.log(`Lipinski: ${report.filters.lipinski_passes ? '✓' : '✗'}`);
```

### Full Example (Rust)

```rust
use chematic_smiles::parse;
use chematic_perception::{find_sssr, assign_aromaticity};
use chematic_chem::*;
use chematic_3d::generate_and_minimize_dreiding;
use chematic_fp::tanimoto_ecfp4;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse
    let benzene = parse("c1ccccc1")?;
    let toluene = parse("Cc1ccccc1")?;

    // Perception
    let rings = find_sssr(&benzene);
    let arom = assign_aromaticity(&benzene);
    println!("Benzene: {} rings, aromatic: {}", 
        rings.ring_count(), 
        arom.is_aromatic(&benzene));

    // Chemistry
    let mw = molecular_weight(&benzene);
    println!("Benzene MW: {:.2}", mw);

    // 3D
    let coords = generate_and_minimize_dreiding(&benzene);
    println!("3D coordinates generated");

    // Fingerprints
    let sim = tanimoto_ecfp4(&benzene, &toluene)?;
    println!("Benzene-Toluene similarity: {:.2}", sim);

    Ok(())
}
```

---

## SMARTS substructure search

```rust
use chematic_smiles::parse;
use chematic_smarts::{parse_smarts, find_matches};

let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap(); // aspirin
let query = parse_smarts("[$(C(=O)O)]").unwrap();   // carboxylic / ester C
let matches = find_matches(&query, &mol);
println!("C(=O)O groups: {}", matches.len()); // 2
```

---

## Molecular descriptors

```rust
use chematic_smiles::parse;
use chematic_chem::{molecular_weight, tpsa, logp_crippen, fsp3, qed, lipinski_passes};

let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
println!("MW:       {:.2}", molecular_weight(&aspirin)); // ~180.16
println!("TPSA:     {:.2}", tpsa(&aspirin));             // ~63.6
println!("LogP:     {:.2}", logp_crippen(&aspirin));     // ~1.2
println!("Fsp3:     {:.3}", fsp3(&aspirin));             // ~0.111
println!("QED:      {:.3}", qed(&aspirin));              // drug-likeness score
println!("Lipinski: {}", lipinski_passes(&aspirin));     // true
```

---

## BRICS fragmentation

```rust
use chematic_smiles::parse;
use chematic_chem::brics_fragments;

let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
let frags = brics_fragments(&aspirin);
println!("fragments: {}", frags.len()); // ≥ 2
```

---

## Fingerprints

```rust
use chematic_smiles::parse;
use chematic_fp::{ecfp4, atom_pair_fp, torsion_fp};

let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
let caffeine = parse("Cn1cnc2c1c(=O)n(c(=O)n2C)C").unwrap();

let sim_ecfp4    = ecfp4(&aspirin).tanimoto(&ecfp4(&caffeine));
let sim_atompair = atom_pair_fp(&aspirin).tanimoto(&atom_pair_fp(&caffeine));
let sim_torsion  = torsion_fp(&aspirin).tanimoto(&torsion_fp(&caffeine));
```

---

## 2D depiction

```rust
use chematic_smiles::parse;
use chematic_depict::depict_svg;

let caffeine = parse("Cn1cnc2c1c(=O)n(c(=O)n2C)C").unwrap();
let svg = depict_svg(&caffeine);
std::fs::write("caffeine.svg", svg).unwrap();
```

### Highlighted depiction

```rust
use std::collections::HashSet;
use chematic_smiles::parse;
use chematic_depict::depict_svg_highlighted;

let mol = parse("c1ccncc1").unwrap(); // pyridine
let n_idx = mol.atoms().find(|(_, a)| a.element.atomic_number() == 7)
               .map(|(i, _)| i).unwrap();
let svg = depict_svg_highlighted(&mol, &HashSet::from([n_idx]), &HashSet::new());
```

---

## JavaScript / TypeScript (WebAssembly)

> **~550 KB, zero C/C++ dependencies.** Drop-in for browser or Node.js.
> Compare with RDKit.js at ~30 MB built via Emscripten.

```sh
npm install @kent-tokyo/chematic
```

```js
import init, {
  parse_smiles, canonical_tautomer, murcko_scaffold,
  largest_fragment, neutralize_charges,
  tanimoto_ecfp4, tanimoto_ecfp6, tanimoto_maccs,
  brics_fragments_json, mcs_smiles_json,
  get_descriptors_json, sssr_rings_json,
  enumerate_stereo_isomers_json,
  sdf_to_records_json, sdf_from_records_json,
  maxmin_picks_ecfp4_json, butina_cluster_ecfp4_json,
  shape_descriptors_json, generate_3d_minimized_pdb,
} from '@kent-tokyo/chematic';

await init();

// ── Parsing & descriptors ─────────────────────────────────────────
const mol = parse_smiles('CC(=O)Oc1ccccc1C(=O)O'); // aspirin
console.log(mol.molecular_weight()); // ~180.16
console.log(mol.qed());              // drug-likeness [0,1]
console.log(mol.sa_score());         // synthetic accessibility [1,10]
console.log(mol.lipinski_passes());  // true

// All descriptors at once (JSON object)
const desc = JSON.parse(get_descriptors_json(mol));
console.log(desc.mw, desc.tpsa, desc.logP, desc.fsp3);

// ── Molecule processing ───────────────────────────────────────────
const salt = parse_smiles('CC(=O)[O-].[Na+]');
const clean = largest_fragment(salt);        // remove Na+
const neutral = neutralize_charges(clean);   // neutralize [O-]

const tautomer = canonical_tautomer(parse_smiles('Oc1cccc2ccccc12'));
const scaffold = murcko_scaffold(parse_smiles('c1ccc(CC(=O)O)cc1'));

// ── Fingerprints & similarity ─────────────────────────────────────
const caffeine = parse_smiles('Cn1cnc2c1c(=O)n(c(=O)n2C)C');
console.log(tanimoto_ecfp4(mol, caffeine));  // ECFP4 Tanimoto
console.log(tanimoto_ecfp6(mol, caffeine));  // ECFP6 Tanimoto
console.log(tanimoto_maccs(mol, caffeine));  // MACCS Tanimoto

// ── Scaffold / fragmentation / MCS ───────────────────────────────
const frags = JSON.parse(brics_fragments_json(mol));
const mcs = mcs_smiles_json('["CC(=O)O","CC(=O)N"]');

// ── Stereochemistry ───────────────────────────────────────────────
const isomers = JSON.parse(enumerate_stereo_isomers_json(parse_smiles('C(F)(Cl)Br')));
// ["[C@@H](F)(Cl)Br","[C@H](F)(Cl)Br"]

// ── 3D geometry ───────────────────────────────────────────────────
const pdb = generate_3d_minimized_pdb(mol);
const shape = JSON.parse(shape_descriptors_json(mol));
console.log(shape.pmi1, shape.npr1, shape.asphericity);

// ── Diversity selection ───────────────────────────────────────────
const library = '["CC","c1ccccc1","CCO","CCCC","c1ccncc1"]';
const picks = JSON.parse(maxmin_picks_ecfp4_json(library, 3));
const clusters = JSON.parse(butina_cluster_ecfp4_json(library, 0.4));

// ── SDF round-trip with properties ───────────────────────────────
const records = JSON.parse(sdf_to_records_json(sdfString));
// records[0].smiles, records[0].name, records[0].properties.MW

const sdf = sdf_from_records_json(
  '["CC(=O)O"]',
  '["aspirin"]',
  '["MW\t180.16\nSource\tChEMBL"]'
);
```

---

## Comparison with Other Cheminformatics Libraries

| Feature                              | **chematic**             | RDKit (rdkit-sys)   | OpenBabel FFI  | RDKit.js (WASM)   |
|--------------------------------------|--------------------------|---------------------|----------------|-------------------|
| **C/C++ dependencies**               | **None — pure Rust**     | Extensive C++       | Extensive C++  | C++ via Emscripten |
| **WASM binary size**                 | **~550 KB**              | N/A (no WASM)       | N/A (no WASM)  | ~30 MB            |
| **Build requirement**                | `cargo build` only       | cmake + clang       | cmake + clang  | Emscripten SDK    |
| **WASM target support**              | **Full (native)**        | No                  | No             | Yes (Emscripten)  |
| Unsafe Rust                          | **None**                 | Extensive           | Extensive      | N/A               |
| OpenSMILES parser                    | Full                     | Full                | Full           | Full              |
| SMILES writer / canonical            | Yes                      | Yes                 | Yes            | Yes               |
| Kekulization                         | Yes                      | Yes                 | Yes            | Yes               |
| Ring perception (SSSR)               | Yes                      | Yes                 | Yes            | Yes               |
| SDF/MOL V2000+V3000 + SD fields      | Yes                      | Yes                 | Yes            | Yes               |
| 2D depiction (SVG, CPK colors)       | Yes                      | Yes                 | Yes            | Yes               |
| ECFP/FCFP fingerprints (2/4/6)       | **All variants + bitvec**| Yes                 | Yes            | Yes               |
| AtomPair / Torsion / MACCS FP        | Yes                      | Yes                 | Yes            | Yes               |
| Molecular descriptors                | **40+ (MW/LogP/…/SA)**   | ~30                 | ~20            | ~30               |
| BRICS fragmentation                  | Yes (bonds + SMILES)     | Yes                 | No             | Yes               |
| Murcko scaffold                      | Yes                      | Yes                 | No             | Yes               |
| Tautomer normalisation               | Yes                      | Yes                 | No             | Yes               |
| MCS                                  | Yes                      | Yes                 | No             | Yes               |
| Stereoisomer enumeration             | **Yes**                  | Yes                 | No             | Yes               |
| CIP stereo (R/S, E/Z) detail         | **Yes (per-atom JSON)**  | Yes                 | Yes            | Yes               |
| 3D coordinate generation             | Yes (DG + minimization)  | Yes (ETKDG)         | Yes            | Yes               |
| 3D shape descriptors (PMI/NPR/…)     | **Yes**                  | Yes                 | No             | Yes               |
| PDB / XYZ file formats               | Yes                      | Yes                 | Yes            | Yes               |
| MaxMin / Butina diversity picking    | **Yes**                  | Yes                 | No             | No                |
| Reaction SMILES/SMIRKS               | Yes                      | Yes                 | Yes            | Yes               |
| InChI / InChIKey                     | No (C lib required)      | Yes                 | Yes            | Yes               |
| Maintenance (2026)                   | Active                   | Active              | Minimal        | Active            |

Notes:
- chematic WASM binary size measured with `wasm-opt` optimization; RDKit.js is the official WASM build.
- "None" for C/C++ means verified: no `*-sys` crates, no `cc` build dependencies, no `build.rs` C compilation in the entire dependency tree.

---

## Recent Development (v0.1.89–v0.1.95)

**v0.1.96**: MMFF94 BCI partial charges
- `mmff94_charges()` upgraded from electronegativity-only (±0.5e) to BCI table (±0.1e). 25 entries covering C=O, C–O, C–N, N–H, O–H, halogens etc. (Halgren 1996).
- Charge conservation guaranteed: sum(q) == sum(formal charges).
- 1,666 library tests, all passing.

**v0.1.95**: True fingerprint algorithms + IUPAC naming expansion
- **MHFP canonical hash**: Morgan-style circular fragment hashes replace atom-index-dependent byte signatures.
- **ERG pharmacophore node types**: DONOR/ACCEPTOR/POSITIVE/NEGATIVE/HYDROPHOBIC; pyridine N vs pyrrole N-H correctly distinguished.
- **Reaction FP XOR**: confirmed as default.
- **IUPAC naming**: ketones, carboxylic acids, esters, amides, benzene, aromatic heterocycles.

**v0.1.91–v0.1.94: RDKit Gap Closure (A1–A5, B3)**
- **v0.1.91**: True MHFP (structural fragment hashing), True ERG (Ertl 2017 functional groups)
- **v0.1.92**: Path FP with bond type interleaving, InChI stereo layer parsing (`/t`, `/b`)
- **v0.1.93**: Full multi-sphere CIP stereochemistry priority rules (moved to chematic-perception, avoids circular dependency)
- **v0.1.94**: SA Score corpus expanded (145 → 188 FDA molecules, 1034 → 1415 unique fragments)

**v0.1.88–v0.1.90**: InChI stereo layers, Brenk SMARTS, reionization, group normalization

**v0.1.69–v0.1.87**: Initial RDKit gap analysis — SSSR, Kekulization, CIP, 3D geometry, WASM API maturity

For detailed historical roadmap (Phases 1–16, v0.1.14–v0.1.33), see `tasks/todo.md`.

---

## Repository Structure

```
chematic/
├── Cargo.toml               workspace root
├── CHANGELOG.md             version history
├── crates/
│   ├── chematic-core/       Atom, Bond, Molecule, Element, kekulization
│   ├── chematic-smiles/     OpenSMILES parser, writer, canonical SMILES
│   ├── chematic-perception/ SSSR ring perception, Huckel aromaticity
│   ├── chematic-mol/        MOL/SDF V2000+V3000 parser and writer
│   ├── chematic-depict/     2D SVG depiction engine (CPK colors, highlighting)
│   ├── chematic-chem/       Descriptors, BRICS, QED, standardization, scaffold
│   ├── chematic-fp/         ECFP4/6, MACCS, path, AtomPair, Torsion FP
│   ├── chematic-smarts/     SMARTS parser + VF2 subgraph isomorphism, MCS
│   ├── chematic-3d/         3D coordinate generation, PDB/XYZ formats
│   ├── chematic-rxn/        Reaction SMILES parser and writer
│   └── chematic/            Umbrella crate with feature flags
└── tasks/
    ├── todo.md              full roadmap checklist (Japanese)
    └── lessons.md           development lessons learned
```

---

## Development Commands

```bash
cargo build --workspace      # build all crates
cargo test --workspace       # run all tests (736)
cargo check --workspace      # type-check without building
cargo clippy --workspace     # lints
```

---

## License

Licensed under either of Apache License 2.0 or MIT License, at your option.
