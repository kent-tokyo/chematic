# chematic-inchi

Pure Rust InChI and InChIKey generation for IUPAC standard organic molecules. WASM-compatible, zero C/C++ dependencies.

## Features

- **InChI string generation**: formula, connectivity, hydrogen, stereo (R/S, E/Z), charge, and isotope layers
- **InChIKey generation**: 27-character deterministic identifier
- **Ring closure bonds**: correctly appends back-edge bonds to connectivity layer
- **Stereo chemistry**: tetrahedral (R/S) and double-bond (E/Z) stereo layers via CIP rules
- **WASM-compatible**: no FFI, compiles to WebAssembly

## Quick Start

```rust
use chematic_inchi::inchi;
use chematic_smiles::parse;

let mol = parse("c1ccccc1").expect("benzene");
let inchi_str = inchi(&mol);
assert_eq!(inchi_str, "InChI=1S/C6H6/c1-2-3-4-5-6-1/h1-6H");

// Generate InChIKey
let key = chematic_inchi::inchi_key(&inchi_str);
assert_eq!(key.len(), 27);
```

## Dependencies

- [`chematic-core`](../chematic-core/README.md) — molecular graph types
- [`chematic-smiles`](../chematic-smiles/README.md) — SMILES parser
- [`chematic-chem`](../chematic-chem/README.md) — CIP rules for stereo assignment
- [`sha2`](https://docs.rs/sha2/) — SHA256 for InChIKey hash

## Layers Generated

| Layer | Description | Example |
|---|---|---|
| `/c` | Connectivity (atom adjacency) | `/c1-2-3-4-5-6-1` |
| `/h` | Hydrogen count per atom | `/h1-6H` |
| `/b` | E/Z double bond stereo | `/b2-3-` |
| `/t` | R/S tetrahedral stereo | `/t2-` |
| `/q` | Formal charge (if net ≠ 0) | `/q+1` |
| `/i` | Isotope information (if present) | `/i1D` |

## References

- IUPAC InChI standard: https://iupac.org/what-we-do/inchi/
- InChI layer structure: https://www.inchi-trust.org/

## See Also

- [`chematic`](../chematic) — main umbrella crate
- [`chematic-smiles`](../chematic-smiles) — SMILES I/O
- [`chematic-depict`](../chematic-depict) — SVG 2D structure drawing
