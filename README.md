# chematic

[日本語](README_ja.md)

A pure-Rust cheminformatics library targeting RDKit feature parity, with no C/C++ FFI.

---

## Design Goals

**Pure Rust, zero C/C++ FFI**
No rdkit-sys, no openbabel bindings. Every algorithm is implemented in safe Rust.

**WASM-compatible and lightweight**
Core crates compile to `wasm32-unknown-unknown` without modification. Binary size is in
the hundreds of KB range, versus tens of MB for C++ FFI wrappers.

**Domain-specific algorithms**
Rather than wrapping a generic graph library, chematic implements chemistry-specific
algorithms directly: Kekulization, Hückel aromaticity, CIP stereochemistry, SSSR ring
perception.

**Reproducible and deterministic**
Fingerprints use FNV-1a hashing with a fixed invariant ordering. Given the same SMILES
input, the same bits are always produced. No RNG, no platform-specific behavior.

---

## Current Status

All phases complete. 544 tests, all passing.

| Crate                 | Description                                                                        | Tests |
|-----------------------|------------------------------------------------------------------------------------|-------|
| `chematic-core`       | Atom, Bond, Molecule, Element, kekulization (no deps)                              | 30    |
| `chematic-smiles`     | OpenSMILES parser, writer, canonical SMILES                                        | 52    |
| `chematic-perception` | SSSR (Balducci-Pearlman), Huckel aromaticity                                       | 14    |
| `chematic-mol`        | MOL/SDF V2000+V3000 parser and writer                                              | 37    |
| `chematic-depict`     | 2D SVG depiction with CPK coloring and atom/bond highlighting                      | 15    |
| `chematic-chem`       | Descriptors, BRICS fragmentation, QED, standardization, Murcko scaffold, CIP      | 216   |
| `chematic-fp`         | ECFP4/6, MACCS 166-bit, topological path, AtomPair, Torsion FP, Tanimoto/Dice     | 44    |
| `chematic-smarts`     | SMARTS parser (recursive, valence, hybridization), VF2 subgraph isomorphism, MCS  | 76    |
| `chematic-3d`         | 3D coordinate generation, PDB/XYZ file formats                                    | 25    |
| `chematic-rxn`        | Reaction SMILES parser and writer                                                  | 15    |
| `chematic-wasm`       | WebAssembly bindings — npm: `@kent-tokyo/chematic`                                 | 18    |
| `chematic`            | Umbrella crate with feature flags (all sub-crates)                                  | 1     |

```
cargo test --workspace   # 544 tests, all passing
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

```sh
npm install @kent-tokyo/chematic
```

```js
import init, { parse_smiles, tanimoto_ecfp4, tanimoto_atom_pair, brics_fragment_count } from '@kent-tokyo/chematic';

await init();

const mol = parse_smiles('CC(=O)Oc1ccccc1C(=O)O'); // aspirin
console.log(mol.molecular_weight()); // ~180.16
console.log(mol.logp_crippen());     // ~1.2
console.log(mol.qed());              // drug-likeness [0,1]
console.log(mol.fsp3());             // fraction sp3 carbons
console.log(brics_fragment_count(mol)); // number of BRICS fragments

const caffeine = parse_smiles('Cn1cnc2c1c(=O)n(c(=O)n2C)C');
console.log(tanimoto_ecfp4(mol, caffeine));    // ECFP4 similarity
console.log(tanimoto_atom_pair(mol, caffeine)); // AtomPair similarity
```

---

## Comparison with Other Cheminformatics Libraries

| Feature                          | chematic                | RDKit (rdkit-sys)  | OpenBabel FFI  | chemcore / purr   |
|----------------------------------|-------------------------|--------------------|----------------|-------------------|
| Language                         | Pure Rust               | Rust + C++ FFI     | Rust + C++ FFI | Pure Rust         |
| WASM target                      | Yes                     | No                 | No             | Partial           |
| Binary size (core)               | ~500 KB                 | ~50 MB             | ~20 MB         | ~200 KB           |
| OpenSMILES parser                | Full                    | Full               | Full           | Partial           |
| SMILES writer / canonical        | Yes                     | Yes                | Yes            | No                |
| Kekulization                     | Yes                     | Yes                | Yes            | No                |
| Aromaticity perception           | Yes (Huckel)            | Yes                | Yes            | Partial           |
| Ring perception (SSSR)           | Yes                     | Yes                | Yes            | No                |
| SDF/MOL V2000+V3000              | Yes                     | Yes                | Yes            | No                |
| 2D depiction (SVG, CPK colors)   | Yes                     | Yes                | Yes            | No                |
| ECFP fingerprints                | Yes (ECFP4/6)           | Yes                | Yes            | No                |
| AtomPair / Torsion fingerprints  | Yes                     | Yes                | Yes            | No                |
| MACCS fingerprints               | Yes (166-bit)           | Yes                | Yes            | No                |
| SMARTS / substructure search     | Yes (VF2 + recursive)   | Yes                | Yes            | No                |
| Molecular descriptors            | Yes (MW/LogP/TPSA/Fsp3/QED/…) | Yes         | Yes            | No                |
| BRICS fragmentation              | Yes                     | Yes                | No             | No                |
| 3D coordinate generation         | Yes (rule-based)        | Yes (ETKDG)        | Yes            | No                |
| PDB/XYZ file formats             | Yes                     | Yes                | Yes            | No                |
| CIP stereochemistry (R/S, E/Z)   | Yes                     | Yes                | Yes            | No                |
| Force field minimization         | Yes (rule-based)        | Yes (UFF/MMFF)     | Yes            | No                |
| Reaction SMILES/SMIRKS           | Yes                     | Yes                | Yes            | No                |
| Unsafe Rust                      | None                    | Extensive          | Extensive      | None              |
| Maintenance (2026)               | Active                  | Active             | Minimal        | Archived          |

Notes:
- Binary sizes are approximate and depend on enabled features.
- chemcore and purr are archived; chematic supersedes their scope.

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
cargo test --workspace       # run all tests (544)
cargo check --workspace      # type-check without building
cargo clippy --workspace     # lints
```

---

## License

Licensed under either of Apache License 2.0 or MIT License, at your option.
