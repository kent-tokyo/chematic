# chematic-3d

Pure Rust 3D coordinate generation, geometry optimization, shape analysis — **zero C/C++ dependencies, WASM-compatible**.

## Features

### 3D Coordinate Generation
- **Distance Geometry (DG)**: Rule-based 3D placement from molecular topology
  - Ideal bond lengths (van der Waals radii)
  - Ideal valence angles (tetrahedral, trigonal, etc.)
  - Handles rings, chains, fragments independently
- **Constraint Satisfaction** (NEW in v0.1.32):
  - Iterative constraint projection (O(n²) per iteration)
  - Bond distance enforcement (±0.05 Å tolerance)
  - Angle enforcement (±5° tolerance)
  - Convergence in 5–10 iterations

### Force Field Minimization
- **DREIDING Force Field**: van der Waals + bond/angle energy terms
  - Lorentz-Berthelot combining rules
  - 1-2 and 1-3 exclusions
  - Velocity Verlet integration
- **UFF (Universal Force Field)**: Generic atom types
- **Configurable**: temperature, damping, convergence criteria

### Molecular Dynamics
- **Velocity Verlet Integration**: 2nd-order accuracy
- **Maxwell-Boltzmann Initialization**: Temperature-correct initial velocities
- **WASM RNG Seeding** (NEW in v0.1.32): Cryptographic randomness in browser
  - Previously: fixed seed in WASM → deterministic trajectories ❌
  - Now: js feature enabled → physical randomness ✅
- **Thermostat Support**: Berendsen, Langevin

### 3D File I/O
- **PDB Reading & Writing**: Standard 80-column format
- **XYZ Format**: Simple XYZ coordinate files
- **Coordinate Indexing**: O(1) access via `AtomIdx`

### Shape Descriptors
- **RMSD** (Root Mean Square Deviation): Geometric similarity
- **Alignment**: Kabsch algorithm, rotation matrices
- **Shape Metrics**:
  - Radius of Gyration (Rg)
  - Principal Moments of Inertia (PMI)
  - Asphericity, Eccentricity
  - Normalized PMI (NPR1, NPR2)

### Stereochemistry
- **3D Stereoisomer Assignment**: From 3D coordinates to R/S
- **Chiral Enumeration**: Generate all stereoisomers

### Conformer Ensemble
- **Multiple Conformers**: Manage diverse 3D conformations
- **RMSD Clustering**: Group similar geometries

## Quick Start

### Generate 3D coordinates

```rust
use chematic_3d::generate_coords;
use chematic_smiles::parse;

let mol = parse("c1ccccc1")?;  // benzene
let coords = generate_coords(&mol);

for i in 0..mol.atom_count() {
    let p = coords.get(AtomIdx(i as u32));
    println!("Atom {}: ({:.3}, {:.3}, {:.3})", i, p.x, p.y, p.z);
}
```

### Optimize geometry with force field

```rust
use chematic_3d::generate_and_minimize_dreiding;

let optimized = generate_and_minimize_dreiding(&mol);
// optimized: 3D coordinates minimized with DREIDING
```

### Apply constraint satisfaction (NEW)

```rust
use chematic_3d::{
    generate_coords, build_constraints, satisfy_constraints,
    generate_and_minimize_constrained
};

// One-step: DG → constraints → DREIDING
let coords = generate_and_minimize_constrained(&mol);

// Manual: control constraint iterations
let coords = generate_coords(&mol);
let constraints = build_constraints(&mol);
let projected = satisfy_constraints(&coords, &mol, &constraints, 20);
```

### Run molecular dynamics

```rust
use chematic_3d::run_md;

let config = MDConfig {
    temperature: 300.0,
    timestep: 0.001,  // fs
    steps: 10000,
    ..Default::default()
};

let trajectory = run_md(&mol, &coords, config)?;
println!("MD completed: {} frames", trajectory.frames.len());
```

### Calculate RMSD

```rust
use chematic_3d::rmsd_no_align;

let d = rmsd_no_align(&coords1, &coords2);
println!("RMSD: {:.3} Å", d);
```

## Performance

| Task | Benzene | Naphthalene | Caffeine |
|------|---------|-------------|----------|
| Generate coords | ~20 µs | ~50 µs | ~100 µs |
| DREIDING minimize | ~500 µs | ~2 ms | ~10 ms |
| Constraint projection | ~150 µs | ~400 µs | ~700 µs |
| 100-step MD | ~50 ms | ~100 ms | ~500 ms |

## WASM Compatibility

All functions compile to `wasm32-unknown-unknown` without modification.

```bash
wasm-pack build --target web --release
```

**Bundle Size**: ~550 KB (chematic-wasm including all dependencies)

## Crate Dependencies

- `chematic-core` — Atom, Bond, Molecule types
- `chematic-perception` — Ring perception, aromaticity
- `chematic-ff` — Force field parameters (DREIDING, UFF)
- `fastrand` — RNG (with `js` feature for WASM)

**Zero FFI**: No rdkit-sys, no Boost, no C/C++ at all.

## Documentation

Full API reference:
```bash
cargo doc --open
```

## Examples

See `examples/` directory:
- `3d_generation.rs` — Basic coordinate generation
- `md_simulation.rs` — Molecular dynamics
- `geometry_optimization.rs` — Force field minimization

## Testing

```bash
cargo test --lib
# 80+ tests covering all 3D operations
```

## Version History

**v0.1.32** (2026-06-07):
- NEW: Constraint satisfaction algorithm with bond/angle enforcement
- NEW: WASM RNG seeding fix (MD now uses cryptographic randomness in browser)
- 12 new constraint tests
- All 80 tests passing

**v0.1.30** (2026-06-07):
- DREIDING force field fix (correct VDW radii, combining rules)
- MD velocity initialization fix (correct thermal velocities)
- SPME Ewald summation (experimental)

## License

MIT OR Apache-2.0

## Contributing

Contributions welcome! See [CONTRIBUTING.md](../../CONTRIBUTING.md).
