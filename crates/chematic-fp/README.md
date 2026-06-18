# chematic-fp

Fast molecular fingerprints for similarity search and clustering. ECFP4/6, RDKit Morgan, MACCS 166-bit, and topological path fingerprints with Tanimoto/Dice similarity. Pure Rust, WASM-compatible.

## Features

- **ECFP fingerprints**: ECFP4 (default), ECFP6, and configurable radius
- **MACCS 166-bit**: public pharmacophore-based bit vector
- **Topological path**: branching-aware path enumeration
- **RDKit-compatible Morgan**: identical to RDKit's default ECFP (FNV-1a 64-bit hash)
- **Similarity metrics**: Tanimoto, Dice, complement similarity
- **MAP4**: MinHashed Atom-Pair fingerprint (Minervini 2020) — `map4(mol, config)` / `tanimoto_map4(a, b)`
- **MHFP / SECFP**: MinHash fingerprint for small molecules — `mhfp(mol, config)`
- **ERG**: Extended Reduced Graph pharmacophore fingerprint — `erg(mol)`
- **Chirality-aware**: optional stereo center consideration (R/S, E/Z)
- **WASM-compatible**: zero C/C++ dependencies

## Quick Start

```rust
use chematic_smiles::parse;
use chematic_fp::ecfp4;
use chematic_fp::similarity::tanimoto;

let mol1 = parse("c1ccccc1O").expect("phenol");
let mol2 = parse("c1ccccc1N").expect("aniline");

let fp1 = ecfp4(&mol1);
let fp2 = ecfp4(&mol2);

let similarity = tanimoto(&fp1, &fp2);
println!("Tanimoto similarity: {:.2}", similarity);
```

## API Overview

### Fingerprint Generators

- `ecfp4(mol: &Molecule) -> Vec<u64>` — 4-bond radius ECFP (most common)
- `ecfp(mol: &Molecule, config: &EcfpConfig) -> Vec<u64>` — configurable ECFP
- `maccs_166(mol: &Molecule) -> BitVec` — MACCS 166 bit pharmacophore vector
- `topo_path(mol: &Molecule) -> Vec<u64>` — topological path fingerprint

### Similarity

- `tanimoto(fp1: &[u64], fp2: &[u64]) -> f64` — Jaccard/Tanimoto coefficient
- `dice(fp1: &[u64], fp2: &[u64]) -> f64` — Dice similarity

## Configuration

```rust
use chematic_fp::EcfpConfig;

let config = EcfpConfig {
    radius: 3,           // 3-bond radius ECFP6
    use_chirality: true, // include R/S, E/Z
    use_bond_types: true,
    min_features: 1,
};
let fp = chematic_fp::ecfp(mol, &config);
```

## Dependencies

- [`chematic-core`](../chematic-core/README.md) — molecular graph
- [`chematic-perception`](../chematic-perception/README.md) — aromaticity

## References

- ECFP/Morgan fingerprints: https://pubs.acs.org/doi/10.1021/ci034667o
- MACCS keys: https://pubs.acs.org/doi/10.1021/c60045a
- Tanimoto similarity: https://en.wikipedia.org/wiki/Jaccard_index

## See Also

- [`chematic`](../chematic) — main umbrella crate
- [`chematic-smiles`](../chematic-smiles) — SMILES parser
- [`chematic-smarts`](../chematic-smarts) — SMARTS substructure matching
