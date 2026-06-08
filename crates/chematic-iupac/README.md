# chematic-iupac

IUPAC systematic name generation for organic molecules. Converts SMILES to systematic nomenclature (e.g., `CC(C)O` → *2-propanol*). Pure Rust, no network required, WASM-compatible.

## Features

- **Systematic IUPAC naming**: alkanes, alkenes, alkynes, alcohols, ethers, aldehydes, ketones, carboxylic acids, amines, amides, and derivatives
- **Locant assignment**: automatic numbering and priority rules
- **Stereochemistry in names**: includes R/S and E/Z nomenclature
- **Common names fallback**: recognizes and outputs common names when appropriate
- **No network dependency**: all rules computed locally
- **WASM-compatible**: zero C/C++ dependencies

## Quick Start

```rust
use chematic_smiles::parse;
use chematic_iupac::iupac_name;

let mol = parse("CC(C)O").expect("2-propanol");
let name = iupac_name(&mol).expect("IUPAC name");
println!("{}", name);  // Output: "2-propanol" or "propan-2-ol"

let mol2 = parse("c1ccccc1O").expect("phenol");
let name2 = iupac_name(&mol2).expect("IUPAC name");
println!("{}", name2);  // Output: "phenol"
```

## API Overview

### Main Function

- `iupac_name(mol: &Molecule) -> Result<String, IupacError>` — generate IUPAC name
- `iupac_name_with_options(mol: &Molecule, opts: &IupacOptions) -> Result<String, IupacError>` — with configuration

### Configuration

```rust
use chematic_iupac::IupacOptions;

let opts = IupacOptions {
    prefer_common_names: true,    // use "phenol" instead of "hydroxybenzene"
    use_systematic_numbering: true,
    include_stereochemistry: true,
};
let name = chematic_iupac::iupac_name_with_options(&mol, &opts)?;
```

## Coverage

Supported functional groups:
- Hydrocarbons: alkanes, alkenes, alkynes, aromatics
- Oxygenated: alcohols, ethers, aldehydes, ketones, carboxylic acids, esters
- Nitrogen: amines, amides, nitriles, imines
- Halogens: alkyl halides, halo-compounds
- Heterocycles: pyridine, furan, thiophene, imidazole (common names)

## Limitations

- No support for complex natural product stereoisomers (tropanes, steroids)
- Limited multi-functional group prioritization (use structure simplification)
- No archaic IUPAC nomenclature (pre-1979 rules)

## Dependencies

- [`chematic-core`](../chematic-core/README.md) — molecular graph
- [`chematic-smiles`](../chematic-smiles/README.md) — SMILES parser

## References

- IUPAC Nomenclature: https://iupac.org/what-we-do/nomenclature/
- Blue Book (Organic Chemistry): https://iupac.org/publications/bluebook/

## See Also

- [`chematic`](../chematic) — main umbrella crate
- [`chematic-chem`](../chematic-chem/README.md) — chemical descriptors and properties
- [`chematic-smiles`](../chematic-smiles/README.md) — SMILES parsing
