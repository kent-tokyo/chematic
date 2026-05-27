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

Phases 1–3 and Phase 5 (coordinate generation + file I/O) are complete.
Phase 4 (MACCS, topological path, MCS, tautomer normalization) is also done.
332 tests, all passing.

| Crate                 | Description                                                             | Tests |
|-----------------------|-------------------------------------------------------------------------|-------|
| `chematic-core`       | Atom, Bond, Molecule, Element, kekulization (no deps)                   | 30    |
| `chematic-smiles`     | OpenSMILES parser, writer, canonical SMILES                             | 50    |
| `chematic-perception` | SSSR (Balducci-Pearlman), Huckel aromaticity                            | 14    |
| `chematic-mol`        | MOL/SDF V2000+V3000 parser and writer                                   | 36    |
| `chematic-depict`     | 2D SVG depiction (ring+chain templates)                                 | 14    |
| `chematic-chem`       | Descriptors, standardization (salt strip, charge), Murcko scaffold, CIP | 67    |
| `chematic-fp`         | ECFP4/ECFP6, MACCS 166-bit keys, topological path FP, Tanimoto/Dice    | 31    |
| `chematic-smarts`     | SMARTS parser, VF2 subgraph isomorphism, MCS                            | 46    |
| `chematic-3d`         | 3D coordinate generation, PDB/XYZ file formats                          | 15    |
| `chematic-rxn`        | Reaction SMILES parser and writer                                        | 15    |
| `chematic`            | Umbrella crate with feature flags (all sub-crates)                       | 1     |

```
cargo test --workspace   # 332 tests, all passing
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
// Using the umbrella crate
use chematic::smiles::{parse, canonical_smiles};
use chematic::fp::ecfp4;
// chematic = { version = "0.1.0", features = ["smiles", "fp"] }
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
    let arom = assign_aromaticity(&benzene);
    println!("aromatic atoms: {}", arom.aromatic_atom_count()); // 6

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
let query = parse_smarts("C=O").unwrap();
let matches = find_matches(&query, &mol);
println!("C=O groups: {}", matches.len()); // 2
```

---

## Molecular descriptors

```rust
use chematic_smiles::parse;
use chematic_chem::{molecular_weight, tpsa, lipinski_passes};

let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
println!("MW:    {:.2}", molecular_weight(&aspirin)); // ~180.16
println!("TPSA:  {:.2}", tpsa(&aspirin));             // ~63.6
println!("Lipinski: {}", lipinski_passes(&aspirin));  // true
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

---

## Comparison with Other Cheminformatics Libraries

| Feature                       | chematic              | RDKit (rdkit-sys)  | OpenBabel FFI  | chemcore / purr   |
|-------------------------------|-----------------------|--------------------|----------------|-------------------|
| Language                      | Pure Rust             | Rust + C++ FFI     | Rust + C++ FFI | Pure Rust         |
| WASM target                   | Yes                   | No                 | No             | Partial           |
| Binary size (core)            | ~500 KB               | ~50 MB             | ~20 MB         | ~200 KB           |
| OpenSMILES parser             | Full                  | Full               | Full           | Partial           |
| SMILES writer                 | Yes                   | Yes                | Yes            | No                |
| Canonical SMILES              | Yes                   | Yes                | Yes            | No                |
| Kekulization                  | Yes                   | Yes                | Yes            | No                |
| Aromaticity perception        | Yes (Huckel)          | Yes                | Yes            | Partial           |
| Ring perception (SSSR)        | Yes                   | Yes                | Yes            | No                |
| SDF/MOL V2000                 | Yes                   | Yes                | Yes            | No                |
| SDF/MOL V3000                 | Yes                   | Yes                | Yes            | No                |
| 2D depiction (SVG)            | Yes                   | Yes                | Yes            | No                |
| ECFP fingerprints             | Yes (ECFP4/6)         | Yes                | Yes            | No                |
| SMARTS / substructure search  | Yes (VF2)             | Yes                | Yes            | No                |
| Molecular descriptors         | Yes (MW/LogP/TPSA/...) | Yes               | Yes            | No                |
| 3D coordinate generation      | Yes (rule-based)      | Yes (ETKDG)        | Yes            | No                |
| PDB/XYZ file formats          | Yes                   | Yes                | Yes            | No                |
| CIP stereochemistry (R/S)     | Yes (R/S, E/Z)        | Yes                | Yes            | No                |
| MACCS fingerprints            | Yes (166-bit keys)    | Yes                | Yes            | No                |
| Force field minimization      | Yes (rule-based)      | Yes (UFF/MMFF)     | Yes            | No                |
| Reaction SMILES/SMIRKS        | Yes                   | Yes                | Yes            | No                |
| Unsafe Rust                   | None                  | Extensive          | Extensive      | None              |
| Maintenance (2026)            | Active                | Active             | Minimal        | Archived          |

Notes:
- "chematic" column reflects current implementation plus the final planned state.
- Binary sizes are approximate and depend on enabled features.
- chemcore and purr are archived; chematic supersedes their scope.

---

## Roadmap

### Phase 1 — Foundation (complete)
Core types, OpenSMILES parse/write, Kekulization, canonical SMILES. 80 tests.

### Phase 2 — Molecular Perception (complete)
SSSR, Huckel aromaticity, SDF/MOL V2000+V3000, 2D SVG depiction. 63 tests.

### Phase 3 — Chemical Intelligence (complete)
Descriptors (MW, LogP, TPSA, Lipinski), ECFP4/6 fingerprints, SMARTS+VF2,
molecular standardization (salt stripping, charge neutralization), Murcko scaffold,
CIP R/S and E/Z stereochemistry assignment.

### Phase 4 — Similarity and Search (complete)
MACCS 166-bit structural keys ✓, topological path fingerprints ✓, MCS ✓, tautomer normalization ✓.

### Phase 5 — 3D Chemistry (partially complete)
Rule-based 3D coordinate generation, PDB/XYZ formats.
Remaining: UFF force field minimization.

### Phase 6 — RDKit Parity (partially complete)
Reaction SMILES/SMIRKS (chematic-rxn) ✓, umbrella crate with feature flags (chematic) ✓.
Remaining: WASM package (npm: chematic), ChEMBL-scale validation.

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
│   ├── chematic-depict/     2D SVG depiction engine
│   ├── chematic-chem/       Molecular descriptors, standardization, scaffold
│   ├── chematic-fp/         ECFP4/6 fingerprints, Tanimoto/Dice similarity
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
cargo test --workspace       # run all tests (332+)
cargo check --workspace      # type-check without building
cargo clippy --workspace     # lints
```

---

## License

Licensed under either of Apache License 2.0 or MIT License, at your option.
