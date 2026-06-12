# chematic-chem

Pure Rust chemical intelligence library — **descriptors, drug-likeness rules, diversity picking, reaction analysis**.

## Features

### Molecular Descriptors (40+)
- **Physicochemical**: MW, LogP (XLogP), TPSA, PSA, MOLAR_REFR, VdW volume
- **Lipophilicity**: LogP (multiple models), MolLogP  
- **Hydrogen Bonding**: HBA, HBD, HBA_LIPINSKI, HBD_LIPINSKI
- **Flexibility**: Rotatable bonds, ring count, scaffold RMSD
- **Complexity**: Topological Polar Surface Area, QED, SA Score
- **Molecular Properties**: Atom/bond counts, isotopes, fragments
- **Fingerprint Similarity**: Tanimoto, Dice (ECFP, MACCS)

### Drug-Likeness Rules
- **Lipinski Rule of Five**: MW < 500, LogP < 5, HBA ≤ 10, HBD ≤ 5
- **Veber's Rules**: TPSA ≤ 140, RotBonds ≤ 10
- **Egan's Rules**: LogP 0.4-5.0, TPSA 20-130
- **Ghose's Rule**: MW 160-480, LogP -0.4-5.6, TPSA 40-130
- **REOS Filters**: Removes likely PAINS, toxic patterns
- **PAINS Filters**: Pan-Assay Interference compounds + alerts

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
| `qed(mol)` | Drug-likeness (0-1) |
| `sa_score(mol)` | Synthetic accessibility (1-10) |
| `lipinski_descriptor_pass(mol)` | Lipinski's Rule of Five |
| `compare_molecules(smiles_list)` | Multi-mol similarity |
| `screen_smiles(smiles_list)` | Batch descriptor + filtering |
| `maxmin_diversity(mols, k)` | Diverse subset picking |
| `find_mcs(mols)` | Maximum common substructure |

## Crate Dependencies

- `chematic-core` — Atom, Bond, Molecule types
- `chematic-smiles` — SMILES parsing
- `chematic-perception` — Ring, aromaticity, stereo
- `chematic-fp` — Fingerprint calculation
- `chematic-mol` — File I/O

**Zero FFI**: Pure Rust, WASM-compatible via npm `@kent-tokyo/chematic`.

## Testing

```bash
cargo test --lib
# 248 tests: descriptors, rules, MCS, diversity
```

## Version History

**v0.1.94** (2026-06-12):
- SA Score fragment corpus expanded: 145 → 188 FDA molecules (1034 → 1415 unique fragments)
- Integrated with multi-sphere CIP and enhanced fingerprints

**v0.1.93** (2026-06-12):
- Full multi-sphere CIP priority rules (moved from chematic-chem to chematic-perception)
- Correct R/S stereochemistry assignment (>2 distinct substituents)

**v0.1.32** (2026-06-07):
- 992 total tests (v0.1.30: 948)
- Integrated with v0.1.32 crate updates (3D constraints, aromaticity)

## License

MIT OR Apache-2.0

