# chematic-smarts

SMARTS parser, VF2 subgraph isomorphism, and maximum common subgraph (MCS) for Rust. Pure implementation, RDKit-compatible, WASM-compatible.

## Features

- **SMARTS parsing**: full OpenSMARTS syntax (atom/bond logical expressions, recursive SMARTS)
- **VF2 matching**: fast subgraph isomorphism detection (substructure search)
- **Maximum Common Subgraph (MCS)**: find largest shared structure between two molecules
- **Chirality-aware matching**: optional R/S and E/Z stereo consideration
- **RDKit compatibility**: parses RDKit SMARTS, produces identical matches
- **WASM-compatible**: zero C/C++ dependencies

## Quick Start

```rust
use chematic_smiles::parse;
use chematic_smarts::{parse_smarts, find_matches};

// Parse a query
let query = parse_smarts("[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]1").expect("benzene ring");

// Parse a molecule
let mol = parse("c1ccccc1O").expect("phenol");

// Find all matches
let matches = find_matches(&query, &mol);
println!("Found {} matches", matches.len());
for m in matches {
    println!("Atoms: {:?}", m.atoms);
}
```

## API Overview

### SMARTS Matching

- `find_matches(query: &Molecule, target: &Molecule) -> Vec<AtomMap>` — all substructure matches
- `find_matches_with_config(query: &Molecule, target: &Molecule, config: &MatchConfig) -> Vec<AtomMap>` — with options
- `has_match(query: &Molecule, target: &Molecule) -> bool` — fast boolean check

### Maximum Common Subgraph

- `mcs(mol1: &Molecule, mol2: &Molecule) -> McsResult` — find largest common subgraph
- Result includes atom/bond mappings and fragment count

### Configuration

```rust
use chematic_smarts::MatchConfig;

let config = MatchConfig {
    use_chirality: true,  // consider R/S, E/Z
    max_matches: 1000,    // limit search
};
let matches = chematic_smarts::find_matches_with_config(&query, &mol, &config);
```

## Dependencies

- [`chematic-core`](../chematic-core/README.md) — molecular graph
- [`chematic-perception`](../chematic-perception/README.md) — aromaticity, hybridization

## References

- OpenSMARTS spec: https://opensmarts.org/
- VF2 algorithm: https://en.wikipedia.org/wiki/VF2_(algorithm)
- SMARTS tutorial: https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html

## See Also

- [`chematic`](../chematic) — main umbrella crate
- [`chematic-smiles`](../chematic-smiles) — SMILES parser
- [`chematic-fp`](../chematic-fp) — fingerprints and similarity
