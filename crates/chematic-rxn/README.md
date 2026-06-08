# chematic-rxn

Reaction SMILES and SMIRKS parser for Rust. Parses chemical transformations (reactants → products) and reaction patterns (transform rules). Pure Rust, RDKit-compatible, WASM-compatible.

## Features

- **Reaction SMILES parsing**: parse reaction equations (e.g., `CC(C)C>>CC(C)[O]`)
- **SMIRKS parsing**: transform patterns for reaction template matching
- **Atom mapping**: track which atoms in reactants map to which atoms in products
- **Reaction properties**: count reactants, products, and agents
- **RDKit compatibility**: parses RDKit reaction SMILES, produces identical results
- **WASM-compatible**: zero C/C++ dependencies

## Quick Start

```rust
use chematic_rxn::parse_reaction;

// Parse a reaction SMILES
let rxn = parse_reaction("CC(C)Br.[Na+].[OH-]>>CC(C)O.[Na+].[Br-]")
    .expect("bimolecular substitution");

println!("Reactants: {}", rxn.reactant_count());
println!("Products: {}", rxn.product_count());

// Access reactants and products as Molecule objects
for mol in rxn.reactants() {
    println!("Atoms: {}", mol.atom_count());
}
```

## API Overview

### Parsing

- `parse_reaction(rxn_smiles: &str) -> Result<Reaction, ParseError>` — parse reaction SMILES
- `parse_smirks(pattern: &str) -> Result<ReactionPattern, ParseError>` — parse SMIRKS pattern

### Reaction Structure

- `Reaction` — contains reactants, agents, products, and optional atom mappings
- `Molecule` — each reactant/product is a standard Molecule
- `ReactionPattern` — SMIRKS pattern for template matching

### Properties

```rust
let rxn = parse_reaction("C.C>>CC")?;
assert_eq!(rxn.reactant_count(), 2);
assert_eq!(rxn.product_count(), 1);
assert_eq!(rxn.agent_count(), 0);

// Get individual molecules
let mol1 = rxn.reactant(0).expect("first reactant");
let prod = rxn.product(0).expect("first product");
```

## Dependencies

- [`chematic-core`](../chematic-core/README.md) — molecular graph
- [`chematic-smiles`](../chematic-smiles/README.md) — SMILES parser

## References

- Reaction SMILES: https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html#reactions
- SMIRKS: https://www.daylight.com/meetings/emug05/Friedman_SMIRKS.pdf

## See Also

- [`chematic`](../chematic) — main umbrella crate
- [`chematic-smiles`](../chematic-smiles) — SMILES parser
- [`chematic-smarts`](../chematic-smarts) — substructure matching and reactions
- [`chematic-depict`](../chematic-depict) — reaction scheme visualization
