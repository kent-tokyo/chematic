# chematic-wasm

WebAssembly bindings for [chematic](https://github.com/kent-tokyo/chematic), a pure-Rust cheminformatics library.

Published to npm as [`@kent-tokyo/chematic`](https://www.npmjs.com/package/@kent-tokyo/chematic).

## Installation

```sh
npm install @kent-tokyo/chematic
```

## Features

- Parse SMILES strings into molecule handles
- Molecular descriptors: MW, TPSA, LogP, Fsp3, QED, exact mass, rotatable bonds, HBD/HBA, aromatic ring count
- Lipinski Rule-of-Five check
- Canonical SMILES generation
- ECFP4, AtomPair, and Topological Torsion fingerprints with Tanimoto similarity
- BRICS fragment count

## Usage

```js
import init, {
  parse_smiles,
  tanimoto_ecfp4,
  tanimoto_atom_pair,
  tanimoto_torsion,
  brics_fragment_count,
} from '@kent-tokyo/chematic';

await init();

const mol = parse_smiles('CC(=O)Oc1ccccc1C(=O)O'); // aspirin

// Descriptors
console.log(mol.atom_count());          // 13
console.log(mol.molecular_weight());    // ~180.16
console.log(mol.formula());             // "C9H8O4"
console.log(mol.tpsa());               // ~63.6
console.log(mol.logp_crippen());        // ~1.2
console.log(mol.fsp3());               // ~0.111
console.log(mol.qed());                // drug-likeness score [0, 1]
console.log(mol.exact_mass());         // ~180.042
console.log(mol.hbd_count());          // 1
console.log(mol.hba_count());          // 4
console.log(mol.rotatable_bond_count()); // 3
console.log(mol.aromatic_ring_count()); // 1
console.log(mol.lipinski_passes());     // true
console.log(mol.canonical_smiles());    // canonical SMILES string

// BRICS fragmentation
console.log(brics_fragment_count(mol)); // ≥ 2

// Fingerprint similarity
const caffeine = parse_smiles('Cn1cnc2c1c(=O)n(c(=O)n2C)C');
console.log(tanimoto_ecfp4(mol, caffeine));    // ECFP4 Tanimoto
console.log(tanimoto_atom_pair(mol, caffeine)); // AtomPair Tanimoto
console.log(tanimoto_torsion(mol, caffeine));   // Torsion Tanimoto
```

## Building from source

```sh
wasm-pack build --target bundler --release
```
