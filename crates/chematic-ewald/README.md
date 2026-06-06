# chematic-ewald — Ewald Summation for Long-Range Electrostatics

Pure-Rust implementation of Ewald summation and Smooth Particle Mesh Ewald (SPME) for accurate long-range electrostatic energy in molecular dynamics simulations. Full WASM support, zero external FFI.

## Features

### Direct Coulomb Sum (Non-Periodic)
Simple pairwise Coulomb energy for isolated molecules:

```rust
use chematic_ewald::direct_coulomb;

let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
let charges = vec![1.0, -1.0];
let energy = direct_coulomb(&coords, &charges);
// → -332.0637 kcal/mol (K_e * q1 * q2 / r)
```

### Coulomb with Real-Space Cutoff
Efficient short-range Coulomb for periodic systems:

```rust
use chematic_ewald::direct_coulomb_cutoff;

let energy = direct_coulomb_cutoff(&coords, &charges, r_cut);
// Only pairs with r < r_cut contribute
```

### Ewald Damped Coulomb (Complementary Error Function)
Smooth damping for short-range Coulomb in Ewald summation:

```rust
use chematic_ewald::direct_coulomb_damped;

let alpha = 3.5 / r_cut;  // Ewald splitting parameter
let energy = direct_coulomb_damped(&coords, &charges, alpha);
// Uses erfc(α*r) damping
```

### SPME Ewald (Reciprocal Space)
FFT-based long-range electrostatics via Smooth Particle Mesh Ewald:

```rust
use chematic_ewald::{spme_energy, BoxVectors, PmeConfig};

let box_vecs = BoxVectors::cubic(30.0);  // 30 Ų cubic box
let config = PmeConfig {
    alpha: 0.0,           // Auto-compute: 3.5 / r_cut
    kmax: [5, 5, 5],      // Reciprocal space cutoff
    r_cut: 8.0,           // Real-space cutoff (Å)
    mesh: [16, 16, 16],   // PME mesh size
    spline_order: 4,      // B-spline interpolation order
};

let result = spme_energy(&coords, &charges, &box_vecs, &config);
// → EwaldResult { real_energy, reciprocal_energy, self_energy, total_energy }
```

## SPME Algorithm

### Components
1. **Real-space**: Damped Coulomb with erfc cutoff (short-range, computed directly)
2. **Reciprocal-space**: FFT-based long-range via charge interpolation to mesh
3. **Self-energy**: Gaussian charge distribution correction

### Workflow
```
Input: atomic coordinates, charges, periodic box
  ↓
1. Assign Ewald parameter α = 3.5 / r_cut
2. Real-space: ∑_{i<j} K_e * q_i * q_j * erfc(α*r_ij) / r_ij
3. Interpolate charges onto 3D mesh (B-spline order 4)
4. Forward FFT to reciprocal space
5. Reciprocal-space energy via Ewald kernel exp(-k²/4α²) / k²
6. Self-energy correction: -K_e * α / √π * ∑_i q_i²
7. Total: E_real + E_reciprocal + E_self
  ↓
Output: total electrostatic energy
```

## Box Vectors and Periodic Boundary Conditions

```rust
use chematic_ewald::BoxVectors;

// Cubic box
let box = BoxVectors::cubic(30.0);

// Triclinic box (general parallelepiped)
let box = BoxVectors([[
    [30.0,  0.0,  0.0],  // a vector
    [ 0.0, 30.0,  0.0],  // b vector
    [ 0.0,  0.0, 30.0],  // c vector
]]);

let volume = box.volume();  // Box volume in Ų
```

## Configuration and Accuracy

### PmeConfig Parameters

| Parameter | Default | Role |
|-----------|---------|------|
| `alpha` | auto (3.5/r_cut) | Ewald splitting: controls real/reciprocal balance |
| `r_cut` | 8.0 Å | Real-space cutoff |
| `mesh` | [16,16,16] | PME grid size; increase for accuracy |
| `kmax` | [5,5,5] | Reciprocal-space cutoff in k-vectors |
| `spline_order` | 4 | B-spline order (3–5); higher = more accurate |

### Accuracy vs. Performance Trade-off

| Use Case | mesh | kmax | α | Accuracy |
|----------|------|------|---|----------|
| Quick estimate | [8,8,8] | [3,3,3] | auto | ~0.1 kcal/mol |
| Production MD | [24,24,24] | [6,6,6] | auto | ~0.01 kcal/mol |
| High precision | [32,32,32] | [8,8,8] | auto | ~0.001 kcal/mol |

## Integration with chematic-3d

Long-range electrostatics for MD simulations:

```rust
use chematic_ewald::spme_energy;
use chematic_chem::gasteiger_charges;

let charges = gasteiger_charges(&mol);
let result = spme_energy(&coords, &charges, &box_vecs, &config);
let e_elec = result.total_energy;
```

## WASM Usage

```javascript
// via chematic-wasm
const mol = parse_smiles("CCO");
const result = coulomb_energy_json(mol);
// → { "coulomb_energy": -12.34, "unit": "kcal/mol" }
```

## Performance

- **Direct Coulomb**: O(N²), exact, no approximation
- **SPME Real-space**: O(N), grid-independent
- **SPME Reciprocal**: O(N log N) via FFT, scales well to ~10k atoms
- **WASM performance**: N ≤ 50 atoms recommended for interactive use

### Memory

- Direct: ~8 bytes per pair
- SPME: ~8 * mesh³ bytes (e.g., 24³ → ~11 MB) + N * O(1)

## Ewald vs. SPME

| Method | Accuracy | Speed | Memory | Suitable For |
|--------|----------|-------|--------|---|
| Direct Coulomb | Exact | O(N²) | O(N) | Small systems, isolated molecules |
| Ewald (Reciprocal) | ~10⁻⁶ | O(N) | O(1) | Periodic systems |
| SPME | ~10⁻⁶ | O(N log N) | O(mesh³) | Large periodic systems |

## References

- Ewald, P. P. (1921). Die Berechnung optischer und elektrostatischer Gitterpotentiale. *Ann. Phys.*, **64**, 253–287.
- Essmann, U., Perera, L., Berkowitz, M. L., Darden, T., Pedersen, L. G., & Roux, B. (1995). A smooth particle mesh Ewald method. *J. Chem. Phys.*, **103**(19), 8577–8593.

## Limitations & Future Work

- ✅ Real-space Coulomb: Complete
- ✅ Reciprocal-space PME: Implemented (B-spline interpolation basic)
- ⏳ Parallel FFT: Planned
- ⏳ GPU acceleration: Deferred to future phase
- ⏳ Polarizable Ewald (DPME): Phase 4+

## License

MIT OR Apache-2.0
