# chematic-perception

Pure Rust molecular structure perception — **ring detection, aromaticity, stereochemistry, charges**.

## Features

### Ring Detection (SSSR)
- **Balducci-Pearlman Algorithm**: Efficiently finds smallest set of smallest rings
- **Gaussian Elimination (GF(2))**: Linear algebra over finite fields
- Usage: `find_sssr(&mol) -> RingSet`

### Aromaticity Model
- **Hückel 4n+2 π-Electron Rule**: Aromatic vs non-aromatic detection
- **Antiaromaticity (4n)** (NEW in v0.1.32): Identifies unstable systems
  - Cyclobutadiene, Cyclooctatetraene, etc.
- **API**: 
  - `assign_aromaticity(&mol) -> AromaticityModel`
  - `ring_classifications(&mol) -> Vec<RingAromaticity>`
  - `antiaromatic_rings(&mol) -> Vec<Vec<AtomIdx>>`
  - `has_antiaromaticity(&mol) -> bool`

### Stereochemistry
- **CIP Priority Rules**: Determine R/S configuration
- **2D → 3D Stereo**: Assign stereoisomers from 2D wedge/dash
- **E/Z Geometry**: Double-bond stereochemistry
- **Chiral Enumeration**: Generate all stereoisomers

### Molecular Properties
- **Implicit Hydrogens**: OpenSMILES valence rules
- **Formal Charges**: Automatic charge assignment
- **Heteroatoms**: N, O, S, P, halogens with correct aromaticity

## Quick Start

### Detect rings

```rust
use chematic_perception::find_sssr;

let rings = find_sssr(&mol);
println!("Rings: {}", rings.ring_count());

for ring in rings.rings() {
    println!("Ring atoms: {:?}", ring.atoms);
}
```

### Classify aromaticity

```rust
use chematic_perception::assign_aromaticity;

let model = assign_aromaticity(&mol);
let classifications = model.ring_classifications(&mol);

for (i, ring_class) in classifications.iter().enumerate() {
    println!("Ring {}: {:?}", i, ring_class);
}
```

### Check for antiaromaticity (NEW)

```rust
if model.has_antiaromaticity(&mol) {
    let antiarom = model.antiaromatic_rings(&mol);
    println!("Antiaromatic rings (unstable): {}", antiarom.len());
}
```

### Assign stereochemistry

```rust
use chematic_perception::assign_stereo_from_2d;

let mol_with_stereo = assign_stereo_from_2d(&mol)?;
// mol_with_stereo: 2D wedge/dash bonds converted to R/S
```

## Algorithms

| Algorithm | Purpose | Complexity |
|-----------|---------|-----------|
| Balducci-Pearlman + GF(2) | SSSR ring detection | O(n³) |
| Hückel 4n+2 Rule | Aromaticity classification | O(n) |
| CIP Priority Rules | R/S stereochemistry | O(n log n) |
| E/Z Geometry Detection | Double-bond stereochemistry | O(n) |

## Crate Dependencies

- `chematic-core` — Atom, Bond, Molecule types
- `chematic-smiles` — SMILES parsing

**Zero FFI**: Pure Rust, WASM-compatible.

## Testing

```bash
cargo test --lib
# 34 tests: SSSR edge cases, aromaticity, stereochemistry
```

## Version History

**v0.1.32** (2026-06-07):
- NEW: Antiaromaticity detection (4n Hückel rule)
- NEW: `ring_classifications()`, `antiaromatic_rings()`, `has_antiaromaticity()` methods
- 16 new aromaticity tests (including antiaromatic systems)

## License

MIT OR Apache-2.0

