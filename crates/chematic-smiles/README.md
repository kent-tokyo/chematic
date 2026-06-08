# chematic-smiles

OpenSMILES parser and writer for Rust. Parses arbitrary SMILES strings into a molecular graph and generates canonical SMILES. Pure Rust, no FFI, WASM-compatible.

## Features

- **SMILES parsing**: read OpenSMILES syntax (aromatic, charges, bonds, isotopes, stereo)
- **Canonical SMILES**: deterministic round-trip generation with sorted atom ordering
- **Stereo chemistry**: Up/Down bonds for 3D chirality, E/Z double bonds
- **Aromaticity perception**: automatic aromatic ring detection (Hückel 4n+2)
- **Hydrogen handling**: implicit/explicit hydrogens, valence checks
- **WASM-compatible**: zero C/C++ dependencies

## Quick Start

```rust
use chematic_smiles::parse;

// Parse SMILES
let mol = parse("c1ccccc1O").expect("phenol");
println!("Atoms: {}", mol.atom_count());
println!("Bonds: {}", mol.bond_count());

// Iterate atoms
for i in 0..mol.atom_count() {
    let atom_idx = chematic_core::AtomIdx(i as u32);
    let atom = mol.atom(atom_idx);
    println!("Atom {}: {}", i, atom.element.symbol());
}

// Iterate bonds
for (_bond_idx, bond) in mol.bonds() {
    println!("Bond: {} - {} (order: {:?})", 
        bond.atom1.0, bond.atom2.0, bond.order);
}
```

## API Overview

- `parse(smiles: &str) -> Result<Molecule, ParseError>` — parse SMILES string
- `canonical_smiles(mol: &Molecule) -> String` — generate canonical SMILES
- `Molecule` — molecular graph with atoms, bonds, and bond orders
- `Atom`, `Bond` — low-level types for graph data

## Dependencies

- [`chematic-core`](../chematic-core/README.md) — graph and chemistry types

## References

- OpenSMILES spec: https://opensmiles.org/
- SMILES tutorial: https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html

## See Also

- [`chematic`](../chematic) — main umbrella crate
- [`chematic-inchi`](../chematic-inchi) — InChI generation
- [`chematic-smarts`](../chematic-smarts) — SMARTS matching and reactions
