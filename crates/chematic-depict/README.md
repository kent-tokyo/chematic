# chematic-depict

2D molecular structure depiction as SVG. Ring templates, wedge/dash stereo bonds, CPK coloring, and automatic layout. Pure Rust, no C/C++ dependencies, WASM-compatible.

## Features

- **SVG output**: vector graphics suitable for web and print
- **PDF output**: `depict_pdf(mol)` / `depict_pdf_opts(mol, opts)` — via `svg2pdf`; enabled by default (`pdf` feature), excluded from WASM
- **EPS output**: `depict_eps(mol)` / `depict_eps_opts(mol, opts)` — pure Rust PostScript, no extra deps, always available
- **Ring templates**: pre-drawn benzene and other common rings for cleaner depiction
- **Stereo bonds**: wedge (Up) and dash (Down) for 3D stereochemistry
- **CPK coloring**: element-specific color scheme (C gray, H white, O red, N blue, etc.)
- **Automatic 2D layout**: uses distance geometry and spring forces for readable coordinates
- **Atom/bond highlighting**: highlight arbitrary atoms or bonds with color overlays
- **Label customization**: hide/show atom symbols, hydrogens, or atom indices
- **PNG rasterization**: optional `png` feature (default=enabled) via `tiny_skia`; excluded from WASM — use SVG there
- **WASM-compatible**: zero C/C++ dependencies (SVG + EPS always available)

## Quick Start

```rust
use chematic_smiles::parse;
use chematic_depict::depict_smiles;

let mol = parse("c1ccccc1O").expect("phenol");

let svg = depict_smiles(
    "c1ccccc1O",
    &Default::default(),
).expect("depiction");

println!("{}", svg);  // Print SVG to stdout or write to file
```

## API Overview

### Main Functions

- `depict_smiles(smiles: &str, options: &DepictOptions) -> Result<String, DepictError>` — SVG from SMILES
- `depict_molecule(mol: &Molecule, options: &DepictOptions) -> String` — SVG from Molecule
- `depict_reaction(rxn_smiles: &str, options: &RenderOptions) -> Result<String, RenderError>` — reaction scheme SVG
- `depict_pdf(mol: &Molecule) -> Vec<u8>` — PDF bytes (`pdf` feature; default on, excluded from WASM)
- `depict_pdf_opts(mol: &Molecule, opts: &DepictOptions) -> Vec<u8>` — PDF with options
- `depict_eps(mol: &Molecule) -> String` — EPS output (always available, pure Rust)
- `depict_eps_opts(mol: &Molecule, opts: &DepictOptions) -> String` — EPS with options

### Configuration

```rust
use chematic_depict::DepictOptions;

let opts = DepictOptions {
    width: 300,
    height: 300,
    atom_color: true,    // CPK colors
    show_atom_labels: true,
    highlight_atoms: vec![0, 1],  // highlight atoms 0 and 1
    highlight_color: "#FF0000",
    ..Default::default()
};
let svg = chematic_depict::depict_molecule(&mol, &opts);
```

## Output Example

For benzene (`c1ccccc1`):
```svg
<svg xmlns="http://www.w3.org/2000/svg" width="300" height="300">
  <!-- regular hexagon with C atoms at vertices -->
  <!-- bonds connecting each vertex -->
</svg>
```

## Dependencies

- [`chematic-core`](../chematic-core/README.md) — molecular graph
- [`chematic-smiles`](../chematic-smiles/README.md) — SMILES parser
- [`chematic-perception`](../chematic-perception/README.md) — 2D coordinate generation

## References

- SVG spec: https://www.w3.org/TR/SVG2/
- CPK coloring: https://en.wikipedia.org/wiki/CPK_coloring

## See Also

- [`chematic`](../chematic) — main umbrella crate
- [`chematic-3d`](../chematic-3d) — 3D structure generation and visualization
- [`chematic-smiles`](../chematic-smiles) — SMILES parser
