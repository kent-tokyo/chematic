# Installation

## Requirements

- Python ≥ 3.8
- No C/C++ compiler required
- No conda required

## Install from PyPI

```bash
pip install chematic
```

That's it. `chematic` ships as a pre-built binary wheel for:

| Platform | x86_64 | aarch64 (Apple M-series / AWS Graviton) |
|---|---|---|
| Linux | yes | yes |
| macOS | yes | yes |
| Windows | yes | — |

## Verify the installation

```python
import chematic
print(chematic.__version__)

mol = chematic.from_smiles("c1ccccc1")
print(mol.formula)   # C6H6
```

## Node.js / WASM

The npm package is a browser-oriented WASM build. Browser applications call
`await init()` normally. In Node, pass the adjacent WASM bytes explicitly
because Node does not fetch `file:` URLs:

```js
import { readFile } from 'node:fs/promises';
import init, { parse_smiles } from '@kent-tokyo/chematic';

const wasm = await readFile(new URL(
  './node_modules/@kent-tokyo/chematic/chematic_wasm_bg.wasm',
  import.meta.url,
));
await init({ module_or_path: wasm });
const mol = parse_smiles('c1ccccc1');
console.log(mol.formula()); // C6H6
mol.free();
```

## Optional dependencies

`chematic` itself has no required runtime dependencies beyond Python and numpy.

For the cookbook examples you may also want:

```bash
pip install pandas scikit-learn matplotlib
```

## Build from source

If a pre-built wheel is not available for your platform:

```bash
pip install maturin
git clone https://github.com/kent-tokyo/chematic
cd chematic/crates/chematic-py
maturin develop --release
```

Requires Rust ≥ 1.75 (`rustup install stable`).
