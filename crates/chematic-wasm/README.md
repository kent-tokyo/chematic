# chematic-wasm

WebAssembly bindings for [chematic](https://github.com/kent-tokyo/chematic), a pure-Rust cheminformatics library.

This crate exposes `#[wasm_bindgen]` bindings so that chematic can be used directly from JavaScript and TypeScript in the browser or Node.js.

## Features

- Parse SMILES strings into molecule handles
- Compute molecular descriptors: molecular weight, TPSA, formula, heavy atom count, H-bond donors/acceptors
- Lipinski Rule-of-Five check
- Canonical SMILES generation
- ECFP4 fingerprints and Tanimoto similarity

## Usage

Build with [wasm-pack](https://rustwasm.github.io/wasm-pack/):

```sh
wasm-pack build --target web
```

Then in JavaScript/TypeScript:

```js
import init, { parse_smiles, tanimoto_ecfp4 } from './pkg/chematic_wasm.js';

await init();

const mol = parse_smiles('c1ccccc1');
console.log(mol.atom_count());       // 6
console.log(mol.molecular_weight()); // ~78.11
console.log(mol.formula());          // "C6H6"
console.log(mol.lipinski_passes());  // true

const aspirin = parse_smiles('CC(=O)Oc1ccccc1C(=O)O');
const sim = tanimoto_ecfp4(mol, aspirin);
console.log(sim); // < 1.0
```
