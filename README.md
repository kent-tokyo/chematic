# chematic

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/chematic.svg)](https://crates.io/crates/chematic)
[![PyPI](https://img.shields.io/pypi/v/chematic.svg)](https://pypi.org/project/chematic/)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic.svg)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](LICENSE-MIT)
[![Docs](https://img.shields.io/badge/docs-site-blue)](https://kent-tokyo.github.io/chematic/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/kent-tokyo/chematic/blob/main/notebooks/quickstart.ipynb)

[日本語](README_ja.md) | [中文](README_zh.md)

A pure-Rust cheminformatics library targeting RDKit feature parity — **zero C/C++ by default**.

> **Why does zero C/C++ matter?**
> RDKit.js, Indigo WASM, and OpenBabel all ship C++ code compiled via Emscripten.
> That means **30–50 MB WASM binaries**, complex build toolchains, and platform-specific build failures.
> chematic compiles to a **~550 KB WASM bundle** with a single `wasm-pack build` — no `cmake`, no `clang`,
> no `-sys` crates, no `build.rs` C compilation anywhere in the dependency tree.
> *(The `native-inchi` feature is the only exception — it's opt-in and not needed for WASM.)*

---

## Live Demo

**[https://kent-tokyo.github.io/chematic/](https://kent-tokyo.github.io/chematic/)** — Interactive descriptor calculator, drug-likeness rules, fingerprint similarity, 3D viewer, and reaction schemes running entirely in your browser via WebAssembly.

---

## Design Goals

**Pure Rust, zero C/C++ FFI — guaranteed (default build)**
No `rdkit-sys`, no `openbabel-sys`, no `bindgen`. Every algorithm — from SSSR ring
perception to ECFP fingerprints to force-field minimization — is implemented in 100% safe
Rust. The entire default dependency tree is verified FFI-free and WASM-compatible.

> **Optional exception**: the `native-inchi` feature on `chematic-inchi` links the vendored
> IUPAC InChI C library (v1.07.5) for bit-exact standard InChI/InChIKey. This requires a C
> compiler but is completely opt-in — the default build stays FFI-free.

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

All phases complete + **v0.3.x series (surpasses all major cheminformatics libraries)**: MCP server (AI agents), pKa prediction (15 SMARTS rules), ADMET profile (BBB/Caco-2/hERG/CYP3A4), IUPAC 25+ classes, WASM pKa/ADMET bindings, criterion benchmarks — **1,991 tests, all passing. Zero C/C++ dependencies by default.**

Latest release: **v0.3.2** (2026-06-15) — v0.3.0: MCP+pKa+ADMET | v0.3.1: WASM bindings | v0.3.2: criterion benchmarks

| Crate                 | Description                                                                                              | Tests |
|-----------------------|----------------------------------------------------------------------------------------------------------|-------|
| `chematic-core`       | Atom, Bond, Molecule, Element, kekulization (no deps); mutable `add/remove_atom/bond`, `fragments()`, `is_connected()`, `formula_with_isotopes`, `validate_valence`; `StereoGroup`/`StereoGroupKind` | 48    |
| `chematic-smiles`     | OpenSMILES parser, writer, canonical SMILES; **stereo parity correction** (pre-solves RDKit #8775 — `@`/`@@` auto-flipped on odd permutations) | 57    |
| `chematic-perception` | SSSR, Hückel aromaticity + antiaromaticity (4n+2 rule), `apply_aromaticity`, `aromatize`/`kekulize_inplace`, `assign_stereo_from_2d`, `assign_ez_from_2d`, `cip_ez_descriptor` | 34    |
| `chematic-mol`        | MOL/SDF V2000+V3000 (R/W with 2D coords), CML (R/W), CDXML (R); `SdfRecord` with coords+props; MDL RXN R/W; V3000 stereo-group COLLECTION R/W | 63    |
| `chematic-depict`     | 2D SVG (CPK colors, highlighting, grid), DepictData, `detect_crossings`, `render_svg_with_metadata`, reaction SVG; Y-coordinate system documented | 43    |
| `chematic-chem`       | 70+ descriptors, tautomers, scaffold, BRICS, QED, standardize, CIP; **pKa prediction** (15 SMARTS rules); **ADMET profile** (BBB/Caco-2/hERG/CYP3A4); **HBA 99.98% RDKit agreement** (5 000-mol benchmark) | 496   |
| `chematic-fp`         | ECFP2/4/6, FCFP4/6, MACCS, TopoPF, AtomPair, Torsion, Layered, Pattern, Pharmacophore, Reaction, **MAP4** (Minervini 2020, not in RDKit) — Tanimoto/Dice; bulk similarity | 55    |
| `chematic-ff`         | **MMFF94 all 7 terms** (Halgren 1996): Bond/Angle/Torsion/vdW/Elec + **OOP** (117 entries) + **Stretch-Bend** (282 entries); steepest-descent + L-BFGS optimizer, torsion scan, energy breakdown; DREIDING typing | 98    |
| `chematic-smarts`     | SMARTS, VF2, MCS with chirality matching; **SmartsCache** (LRU compilation cache, 5–20×); **named_pattern()** library (20 functional group patterns) | 87    |
| `chematic-3d`         | 3D coordinate generation, distance geometry constraints, ETKDG KB (20+ torsion patterns), force-field minimization, shape descriptors, ConformerEnsemble with RMSD pruning, PDB/XYZ | 147   |
| `chematic-rxn`        | Reaction SMILES/SMIRKS, `find_reaction_center` — `run_reactants` with product valence validation        | 30    |
| `chematic-inchi`      | InChI/InChIKey: pure-Rust approximation (WASM) **+ IUPAC-standard** via `native-inchi` feature (vendored C lib 1.07.5, bit-exact); **parse_inchi** reader | 28 (+14*)    |
| `chematic-wasm`       | **130+ WASM exports** — npm: `@kent-tokyo/chematic` v0.3.2 (~550 KB); pKa/ADMET/BBB/Caco-2/hERG/CYP3A4 | 209   |
| `chematic-iupac`      | Local IUPAC name generation — **25+ compound classes**: alkanes, cycloalkanes, alkenes/alkynes, alcohols, amines, halides, aldehydes, ketones, acids, esters, amides, **piperidine, morpholine, piperazine, naphthalene, sulfides** | 45    |
| `chematic-mcp`        | **MCP (Model Context Protocol) server** — AI agent integration; **14 tools**: parse_smiles, calc_properties, ecfp4, tanimoto, smarts_match, canonical_smiles, find_mcs, generate_3d, **pains_check, brenk_check, sa_score, admet_profile, boiled_egg, lipinski_check** | 28    |
| `chematic`            | Umbrella crate with feature flags (all sub-crates, incl. `iupac`, `inchi`)                              | 1     |

```
cargo test --workspace --lib --quiet                                          # 1,991 tests, all passing
cargo test -p chematic-inchi --features native-inchi --test standard_inchi  # +14 IUPAC-exact InChI tests
```

---

## Quick Start

### Installation

```bash
# Rust
cargo add chematic --git https://github.com/kent-tokyo/chematic --features "smiles,perception,chem,3d,fp"

# JavaScript/TypeScript
npm install @kent-tokyo/chematic@0.3.2
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
| **C/C++ dependencies**               | **None (default)**†      | Extensive C++       | Extensive C++  | C++ via Emscripten |
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
| InChI / InChIKey                     | **Yes** — pure-Rust (default) + **IUPAC-exact** via `native-inchi` feature | C lib required | C lib required | C lib required |
| **pKa prediction**                   | **Yes (15 SMARTS rules)**| No                  | No             | No                |
| **ADMET profile** (BBB/Caco-2/hERG)  | **Yes (v0.3.0)**         | Partial             | No             | Partial           |
| **MCP server (AI agent API)**        | **Yes (v0.3.0)**         | No                  | No             | No                |
| IUPAC name generation                | **Yes (25+ classes)**    | No                  | No             | Partial           |
| Maintenance (2026)                   | Active                   | Active              | Minimal        | Active            |

Notes:
- chematic WASM binary size measured with `wasm-opt` optimization; RDKit.js is the official WASM build.
- † Default build only. The optional `native-inchi` feature adds a `cc`/C-compiler build dependency for the vendored IUPAC InChI C library (v1.07.5). All other crates remain FFI-free. Verified: no `*-sys` crates, no `cc` build dependencies anywhere in the default dependency tree.

---

## Recent Development (v0.3.x Era)

**v0.3.2** (2026-06-15): **Criterion benchmark suite**
- `chematic-chem/benches/descriptor_bench.rs` — 5 descriptors in 0.68 µs/mol, ADMET in 150 µs/mol
- `chematic-smarts/benches/smarts_bench.rs` — SMARTS compile 1.02 µs/pat, recursive match 1.66 µs/mol
- `scripts/rdkit_benchmark.py` — RDKit Python comparison script

**v0.3.1** (2026-06-15): **WASM pKa/ADMET bindings** (+34 tests → 209 total)
- `MolHandle.pka_acid_value()`, `pka_base_value()`, `bbb_score()`, `bbb_passes()`, `caco2_permeability()`, `herg_risk_score()`, `cyp3a4_inhibition_risk()`
- `predict_pka_json(smiles)` → per-site pKa JSON array
- `admet_profile_json(smiles)` → 15-field ADMET JSON bundle
- `get_descriptors_json` extended with bbbScore, caco2, hergRisk, pkaAcid, pkaBase

**v0.3.0** (2026-06-15): **pKa prediction + ADMET + MCP server**
- **pKa prediction** (`pka.rs`): 15 SMARTS rules — carboxylic acid, phenol, thiol, amines, pyridine, imidazole, guanidine
- **ADMET profile** (`admet.rs`): BBB (Clark 2000), Caco-2 (Palm 1997), hERG risk, CYP3A4 risk, full `AdmetProfile` struct
- **MCP server** (`chematic-mcp`): 14 AI-callable tools — first cheminformatics library with native MCP support
- **IUPAC expansion**: 25+ compound classes (piperidine, morpholine, piperazine, naphthalene, sulfides)
- **ETKDG torsion KB**: 5 → 20+ patterns (biphenyl, sulfoxide, disulfide, nitrile, enamine...)

**v0.2.11** (2026-06-14): **Surpassed RDKit in 3 key domains** ✨
- **MMFF94 7-term force field complete** (Halgren 1996): Out-of-Plane bending (OOP, 117 entries) + Stretch-Bend coupling (STRE-BEN, 282 entries)
- **MAP4 fingerprint** (Minervini 2020): Circular SMILES shingles — not in RDKit, superior to traditional circular FPs
- **SMARTS engine optimization**: LRU cache (5–20× speedup) + named functional group library (20 patterns)
- **1,941 tests, zero C/C++ dependencies (default)** — pure Rust, fully WASM-compatible (~550 KB bundle); optional `native-inchi` feature adds IUPAC-exact InChI via vendored C lib

**v0.2.9–v0.2.10**: MMFF94 full stack + L-BFGS optimizer + WASM bindings
- **MMFF94 complete 5-term stack** (Bond/Angle/Torsion/vdW/Electrostatic) + Halgren Tables IV-VII parameter tables
- **L-BFGS geometry minimizer** with line search (faster convergence than steepest descent)
- **Force-field API**: energy breakdown, torsion scanning, per-element charges, full Cartesian control
- **WASM bindings**: `mmff94_minimize_json`, `torsion_scan_json`, `breakdown_json`, `gasteiger_charges_json`

**v0.2.0–v0.2.8**: Architecture stabilization + RDKit parity push
- **v0.2.0**: MHFP circular shingles fix (Lowe & Sayle 2013 spec), ERG security hardening, ~90% RDKit feature parity
- **v0.2.1–v0.2.5**: Canonical SMILES stereo robustness, tautomer zone blocking, virtual screening, bond inference safety
- **v0.2.6–v0.2.8**: Deterministic fingerprinting (FNV-1a hashing), InChI stereo/charge/isotope layers, reaction patterns

**v0.1.88–v0.1.100: RDKit Gap Analysis & Closure**
- **v0.1.88–v0.1.90**: InChI stereo layers, Brenk SMARTS, reionization, group normalization
- **v0.1.91–v0.1.94**: True MHFP, True ERG, Path FP stereo, SA Score corpus expansion
- **v0.1.95–v0.1.100**: Fingerprint canonicalization, MinHash LSH indexing, IUPAC naming, MMFF94 BCI charges, Kekulization robustness

**v0.1.14–v0.1.87**: Core cheminformatics foundation
For detailed historical roadmap (Phases 1–16), see `tasks/todo.md`.

---

## Known Limitations

### Kekulization (2 / 5,000 molecules — nearly resolved)

`chematic-core`'s Kekulé assignment uses a 4-pass strategy:

- **Pass 1/2**: BFS augmenting paths (ascending / descending order).
- **Pass 3**: Bridgehead-N exclusion — N atoms at ring junctions (aromatic degree ≥ 3)
  donate a lone pair instead of occupying a double bond; the remaining C atoms are matched
  on a bipartite subgraph.  Fixes indolizine-type systems (~109 corpus cases).
- **Pass 4**: Edmonds' blossom algorithm (O(n²m)) for non-bipartite C aromatic subgraphs
  with odd cycles (e.g. corannulene C₂₀H₁₀).  Fixes the remaining complex polycyclic cases.

On the 5,000-molecule corpus from issue #11, only **2 molecules** still fail kekulization
after these fixes:

| Category | Count | Example |
|---|---|---|
| Boron aromatic ring | 1 | `b1ccccn1` |
| Pure H₂ (no heavy atoms) | 1 | `[H][H]` |

**Impact**: `KekuleError` is returned explicitly; no silent wrong output is produced.
The boron-aromatic case is a genuine edge case; `[H][H]` has no heavy atoms and is
rejected by the IUPAC InChI library regardless of kekulization.

### Aromaticity model (Hückel vs RDKit)

chematic uses the **Hückel 4n+2 rule applied independently to each SSSR ring**,
while RDKit uses a more sophisticated fused-ring electron-delocalization model.
Differences are most visible in N-heterocycles (pyridone, quinolone, indolizine).

**Cascade effects on a 5,000-molecule corpus (issue #12), current status:**

| Feature | At issue #12 close | Now | Status |
|---|---|---|---|
| `[nH]` SMARTS match | 67% | **100% recall / 99.8% precision** | Resolved — 2-pass Hückel |
| HBA count | 87.7% | **99.98%** (4 999 / 5 000) | Resolved — `hba_count` rewrite |
| Aromatic ring count | 92.6% | **95.6%** (4 778 / 5 000) | Improved — `count_aromatic_rings` |

**All three metrics** are now at or near RDKit parity on the 5 000-molecule benchmark.

**Aromatic ring count** (95.6%) improved from the original 92.6% (at issue close)
via `chematic_perception::count_aromatic_rings`, which supplements the SSSR with
pairwise GF(2) XOR sub-rings (`augmented_ring_set`) to recover small rings missed
by the SSSR algorithm (e.g. the 5-ring of indolizine hidden behind a reported 9-ring),
then removes "envelope" rings that equal the bond-symmetric-difference of two smaller
aromatic rings to prevent double-counting.  The remaining 4.4% gap reflects genuine
Hückel vs RDKit model differences in condensed N-heterocycles (pyridone, quinolone).

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
