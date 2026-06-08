# chematic-core

Core types for chematic: `Atom`, `Bond`, `Molecule`, and `Element`. Low-level graph representation for all chematic crates. Pure Rust, WASM-compatible.

## Features

- **Molecule graph**: adjacency-list representation with atoms and bonds
- **Atom types**: element, hybridization, formal charge, aromatic flag, isotope
- **Bond types**: single, double, triple, aromatic bonds with configurable order
- **Element table**: periodic table with atomic number, mass, electronegativity, van der Waals radius
- **Molecular properties**: atom count, bond count, molecular weight, formula
- **WASM-compatible**: zero C/C++ dependencies

## Quick Start

```rust
use chematic_core::{Molecule, AtomIdx, Element, BondOrder};

// Create a molecule manually (or use a parser)
let mut mol = Molecule::new();

// Add atoms (returns their indices)
let c0 = mol.add_atom(Element::C.properties(), false);  // methane C
let h0 = mol.add_atom(Element::H.properties(), false);
let h1 = mol.add_atom(Element::H.properties(), false);
let h2 = mol.add_atom(Element::H.properties(), false);
let h3 = mol.add_atom(Element::H.properties(), false);

// Add bonds
mol.add_bond(c0, h0, BondOrder::Single);
mol.add_bond(c0, h1, BondOrder::Single);
mol.add_bond(c0, h2, BondOrder::Single);
mol.add_bond(c0, h3, BondOrder::Single);

// Query the graph
println!("Atoms: {}", mol.atom_count());      // 5
println!("Bonds: {}", mol.bond_count());      // 4
println!("Mass: {:.2}", mol.molecular_weight());  // 16.04

// Iterate neighbors
for (neighbor_idx, bond_idx) in mol.neighbors(c0) {
    println!("Neighbor: {}", neighbor_idx.0);
}
```

## API Overview

### Types

- `Molecule` — molecular graph: atoms, bonds, adjacency
- `Atom` — atomic data: element, charge, hybridization, aromaticity, isotope
- `Bond` — bond data: atoms, order, stereochemistry (Up/Down, E/Z)
- `Element` — element info: atomic number, mass, vdW radius, electronegativity
- `AtomIdx`, `BondIdx` — newtype wrappers for type safety

### Molecular Properties

```rust
let mw = mol.molecular_weight();
let formula = mol.formula();
let neighbors = mol.neighbors(atom_idx);
let atom = mol.atom(atom_idx);
let bond = mol.bond(bond_idx);
```

## Element Constants

```rust
use chematic_core::Element;

Element::C,     // carbon
Element::H,     // hydrogen
Element::N,     // nitrogen
Element::O,     // oxygen
Element::P,     // phosphorus
Element::S,     // sulfur
Element::Cl,    // chlorine
// ...and all others
```

## Dependencies

None (no external dependencies).

## See Also

- [`chematic-smiles`](../chematic-smiles/README.md) — SMILES parser (builds Molecules)
- [`chematic-smarts`](../chematic-smarts/README.md) — SMARTS matching
- [`chematic`](../chematic) — main umbrella crate
