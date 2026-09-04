# chematic-chem

Pure Rust chemical intelligence library — **descriptors, drug-likeness rules, diversity picking, reaction analysis**.

## Features

### Molecular descriptors

The Python facade exposes 70+ descriptor functions and 190+ individual values;
the Rust crate provides the underlying scalar and vector implementations.
RDKit agreement is metric-specific rather than universal; see
[`docs/validation.md`](../../docs/validation.md).

- **Physicochemical**: MW, LogP (Crippen), TPSA, PSA, MOLAR_REFR, VdW volume
- **Van der Waals volume**: `vabc(mol) -> f64` — Bondi radii + spherical-cap overlaps, no 3D coordinates required
- **Lipophilicity**: LogP (multiple models), MolLogP  
- **Hydrogen Bonding**: HBA (99.98% RDKit agreement on 5 000 molecules), HBD (100% agreement, including S-H thiols), HBA_LIPINSKI, HBD_LIPINSKI
- **Flexibility**: Rotatable bonds, ring count, scaffold RMSD
- **Complexity**: Topological Polar Surface Area, QED, SA Score
- **Topological indices**:
  - `schultz_mti(mol) -> u64` — Schultz Molecular Topological Index: Σ_{i<j}(δᵢ+δⱼ)×dᵢⱼ
  - `gutman_mti(mol) -> u64` — Gutman MTI*: Σ_{i<j}δᵢ×δⱼ×dᵢⱼ
  - `gravitational_index(mol) -> f64` — Σ_{i<j}mᵢmⱼ/dᵢⱼ² (atomic masses, graph distance)
- **Molecular Properties**: Atom/bond counts, isotopes, fragments
- **Fingerprint Similarity**: Tanimoto, Dice (ECFP, MACCS)

### Drug-Likeness Rules
- **Lipinski Rule of Five**: MW < 500, LogP < 5, HBA ≤ 10, HBD ≤ 5
- **Veber's Rules**: TPSA ≤ 140, RotBonds ≤ 10
- **Egan's Rules**: LogP 0.4-5.0, TPSA 20-130
- **Ghose's Rule**: MW 160-480, LogP -0.4-5.6, TPSA 40-130
- **REOS Filters**: Removes likely PAINS, toxic patterns
- **PAINS Filters**: Pan-Assay Interference compounds + alerts
- **Brenk Filters**: Toxicity / instability structural alerts
- **BOILED-Egg** (`boiled_egg()` → `BoiledEggProfile`): GI absorption + BBB penetration prediction (Daina & Zoete 2016)

### Advanced Analysis
- **Murcko Scaffold**: Generic scaffold from molecule
- **Maximum Common Substructure (MCS)**: Molecular similarity
- **Reaction Analysis**: Atom-map detection, reaction center identification
- **Stereochemistry**: Invert/enumerate stereoisomers
- **Isotope Distribution**: Exact mass spectrum simulation

### Diversity Selection
- **MaxMin Picking**: Diverse subset selection
- **Butina Clustering**: Distance-based clustering
- **Tanimoto Similarity**: ECFP-based distances

### Standardization
- **Mol Hash**: Canonical identifier (accounting for salts, tautomers)
- **Standardize**: Remove/identify salts, neutralize charges
- **Parse Condensed**: IUPAC-like condensed formulas

## Quick Start

### Calculate descriptors

```rust
use chematic_chem::molecular_weight;
use chematic_smiles::parse;

let mol = parse("CC(=O)Oc1ccccc1C(=O)O")?; // aspirin
let mw = molecular_weight(&mol);
println!("MW: {:.2}", mw);
```

### Check drug-likeness

```rust
use chematic_chem::{
    lipinski_descriptor_pass,
    veber_descriptor_pass,
    pains_filters
};

if lipinski_descriptor_pass(&mol) {
    println!("✓ Passes Lipinski's Rule of Five");
} else {
    println!("✗ Violates Lipinski");
}

let pains = pains_filters(&mol);
if !pains.is_empty() {
    println!("PAINS alerts: {:?}", pains);
}
```

### Similarity search

```rust
use chematic_chem::compare_molecules;

let results = compare_molecules(&[smiles1, smiles2])?;
for pair in &results.pairwise {
    println!("Similarity: {:.3}", pair.similarities.ecfp4_tanimoto);
}
```

### Diversity picking

```rust
use chematic_chem::maxmin_diversity;

let picks = maxmin_diversity(&molecules, 10)?;
println!("Selected {} diverse molecules", picks.len());
```

## API Reference

| Function | Purpose |
|----------|---------|
| `molecular_weight(mol)` | Average isotope mass |
| `logp(mol)` | XLogP model |
| `tpsa(mol)` | Topological Polar Surface Area |
| `vabc(mol)` | Van der Waals volume (Bondi radii, no 3D needed) |
| `qed(mol)` | Drug-likeness (0-1) |
| `sa_score(mol)` | Synthetic accessibility (1-10) |
| `schultz_mti(mol)` | Schultz Molecular Topological Index |
| `gutman_mti(mol)` | Gutman MTI* |
| `gravitational_index(mol)` | Gravitational index (Σ mᵢmⱼ/dᵢⱼ²) |
| `lipinski_descriptor_pass(mol)` | Lipinski's Rule of Five |
| `compare_molecules(smiles_list)` | Multi-mol similarity |
| `screen_smiles(smiles_list)` | Batch descriptor + filtering |
| `maxmin_diversity(mols, k)` | Diverse subset picking |
| `find_mcs(mols)` | Maximum common substructure |
| `boiled_egg(mol)` | BOILED-Egg GI/BBB prediction (→ `BoiledEggProfile`) |

## Crate Dependencies

- `chematic-core` — Atom, Bond, Molecule types
- `chematic-smiles` — SMILES parsing
- `chematic-perception` — Ring, aromaticity, stereo
- `chematic-fp` — Fingerprint calculation
- `chematic-mol` — File I/O

**Zero FFI**: Pure Rust, WASM-compatible via npm `@kent-tokyo/chematic`.

## Validation

Run `cargo test -p chematic-chem`. Corpus-level agreement and known residuals are documented in [`docs/validation.md`](../../docs/validation.md). Release history is maintained in the repository [`CHANGELOG.md`](../../CHANGELOG.md).

## License

MIT OR Apache-2.0
