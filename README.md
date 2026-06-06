# chematic

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

All phases complete. **869 tests, all passing. Zero C/C++ dependencies.**

| Crate                 | Description                                                                                              | Tests |
|-----------------------|----------------------------------------------------------------------------------------------------------|-------|
| `chematic-core`       | Atom, Bond, Molecule, Element, kekulization (no deps); mutable `with_atom_*` / `with_bond_*` API        | 30    |
| `chematic-smiles`     | OpenSMILES parser, writer, canonical SMILES                                                              | 57    |
| `chematic-perception` | SSSR (Balducci-Pearlman), Huckel aromaticity                                                             | 14    |
| `chematic-mol`        | MOL/SDF V2000+V3000 (R/W with 2D coords), CML (R/W), CDXML multi-fragment + stereo (R), SDF props       | 53    |
| `chematic-depict`     | 2D SVG depiction (CPK colors, highlighting, grid), DepictData for canvas renderers                       | 30    |
| `chematic-chem`       | 40+ descriptors, BRICS, QED, standardization, Murcko scaffold, CIP, IFG, Gasteiger, VSA, SA score, MMP  | 216   |
| `chematic-fp`         | ECFP2/4/6, FCFP4/6, MACCS 166-bit, TopoPF, AtomPair, Torsion — Tanimoto/Dice                           | 50    |
| `chematic-smarts`     | SMARTS parser (recursive, valence, hybridization), VF2 subgraph isomorphism, MCS                        | 82    |
| `chematic-3d`         | 3D coordinate generation, force-field minimization, shape descriptors, ConformerEnsemble, PDB/XYZ       | 68    |
| `chematic-rxn`        | Reaction SMILES parser and writer                                                                        | 26    |
| `chematic-wasm`       | **100+ WASM exports** — npm: `@kent-tokyo/chematic`                                                      | 162   |
| `chematic`            | Umbrella crate with feature flags (all sub-crates)                                                       | 1     |

```
cargo test --workspace   # 869 tests, all passing
```

---

## Quick Start

### Using the umbrella crate

```toml
# Cargo.toml
[dependencies]
chematic = { git = "https://github.com/kent-tokyo/chematic", features = ["smiles", "fp"] }
```

```rust
use chematic::smiles::{parse, canonical_smiles};
use chematic::fp::ecfp4;
```

### Using individual crates

```toml
# Cargo.toml
[dependencies]
chematic-smiles     = { git = "https://github.com/kent-tokyo/chematic" }
chematic-perception = { git = "https://github.com/kent-tokyo/chematic" }
chematic-fp         = { git = "https://github.com/kent-tokyo/chematic" }
```

```rust
use chematic_smiles::{parse, canonical_smiles};
use chematic_perception::{find_sssr, assign_aromaticity};
use chematic_fp::{ecfp4, tanimoto_ecfp4};

fn main() {
    let benzene = parse("c1ccccc1").unwrap();
    let toluene = parse("Cc1ccccc1").unwrap();

    // Ring and aromaticity perception
    let rings = find_sssr(&benzene);
    println!("rings: {}", rings.ring_count()); // 1

    // Fingerprint similarity
    let sim = tanimoto_ecfp4(&benzene, &toluene);
    println!("Tanimoto(benzene, toluene): {sim:.3}"); // ~0.5

    // Canonical SMILES
    println!("{}", canonical_smiles(&benzene)); // c1ccccc1
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

## Roadmap

### Phase 1 — Foundation (complete)
Core types, OpenSMILES parse/write, Kekulization, canonical SMILES.

### Phase 2 — Molecular Perception (complete)
SSSR, Huckel aromaticity, SDF/MOL V2000+V3000, 2D SVG depiction.

### Phase 3 — Chemical Intelligence (complete)
Descriptors (MW, LogP, TPSA, Fsp3, Lipinski), QED, BRICS fragmentation,
ECFP4/6 fingerprints, SMARTS+VF2 (recursive SMARTS, valence, hybridization),
molecular standardization, Murcko scaffold, CIP R/S and E/Z.

### Phase 4 — Similarity and Search (complete)
MACCS 166-bit keys, topological path FP, AtomPair FP, Topological Torsion FP,
MCS, tautomer normalization.

### Phase 5 — 3D Chemistry (complete)
Rule-based 3D coordinate generation, PDB/XYZ formats, UFF-like minimization.

### Phase 6 — RDKit Parity (complete)
Reaction SMILES/SMIRKS ✓, umbrella crate with feature flags ✓,
WASM npm package `@kent-tokyo/chematic` ✓, CPK coloring + highlighted depiction ✓,
ChEMBL 37 full-set validation (2,897,819 molecules, 100.000%) ✓.

### Phase 7 — Extended Descriptors and Diversity (v0.1.14–v0.1.15, complete)
EState indices (Hall & Kier 1991), path fingerprint (DFS path FP, 2048-bit),
SDF/MOL WASM bindings,
functional group identification (Ertl 2017 IFG), Gasteiger-Marsili PEOE partial charges,
VSA descriptors (SlogP_VSA × 12, SMR_VSA × 10, PEOE_VSA × 14),
SA score (complexity-based), MaxMin diversity picking, Butina clustering.

See `tasks/todo.md` for the detailed per-task breakdown.

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
