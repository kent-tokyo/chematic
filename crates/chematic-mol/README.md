# chematic-mol

Pure Rust molecular file format reader/writer — **SDF, MOL V2000/V3000, CML, CDXML, RXN, MDL V2000/V3000**.

## Features

### MOL Format Support
- **V2000 Parser/Writer**: Full support with `parse_mol` / `write_mol`
- **V2000 with Coords**: 2D coordinates preserved via `parse_mol_with_coords` / `write_mol_with_coords`
- **V3000 Parser/Writer**: Extended blocks (CTAB, ATOM, BOND, COLLECTION)
- **V3000 with Coords** (NEW in v0.1.32): 2D coordinates now recovered via `parse_mol_v3000_with_coords`
  - Previously: coordinates discarded ❌
  - Now: coordinates preserved in Vec<(f64, f64)> ✅

### SDF (Structure-Data Format)
- **Multi-Molecule**: Iterate over molecules with `SdfReader`
- **Coordinate Preservation**: Per-molecule 2D coords via `parse_sdf_with_coords`
- **Property Fields**: Read/write arbitrary data (e.g., activity, source)

### CML (Chemical Markup Language)
- **Parse/Write**: XML-based molecular format
- **Y-Coordinate System** (DOCUMENTED in v0.1.32):
  - CML uses **chemical Y-up** (Y increases upward)
  - SVG rendering requires Y-negation: `svg_y = -cml_y`

### CDXML (ChemDraw XML)
- **Parser**: ChemDraw XML format support
- **Y-Coordinate System** (DOCUMENTED in v0.1.32):
  - CDXML uses **ChemDraw Y-down** (Y increases downward, SVG-compatible)
  - No Y-axis conversion needed for SVG rendering

### MDL RXN (Reaction Format)
- **Parse/Write**: Reaction SMILES/SMIRKS via RXN format
- **Validation**: Product valence checking

### Error Handling (NEW in v0.1.32)
- **Display trait**: All error types implement `std::fmt::Display`
- **Error trait**: CmlError, CdxmlError, Mol2Error, RxnParseError now implement `std::error::Error`

## Quick Start

### Parse a single molecule

```rust
use chematic_mol::parse_mol;

let mol_str = "..."; // MOL V2000 or V3000 string
let (mol, metadata) = parse_mol(mol_str)?;
println!("Atoms: {}", mol.atom_count());
```

### Parse SDF (multi-molecule file)

```rust
use chematic_mol::SdfReader;

let sdf_content = "...";
let mut reader = SdfReader::new(sdf_content.as_bytes())?;

while let Some(record) = reader.next_record()? {
    let (mol, metadata) = record.parse_mol()?;
    println!("Molecule: {}", metadata.name);
    // record.properties: HashMap<String, String>
}
```

### Preserve 2D coordinates

```rust
use chematic_mol::parse_mol_with_coords;

let (mol, metadata, coords) = parse_mol_with_coords(mol_str)?;
// coords: Vec<(f64, f64)> — 2D coordinates for layout/rendering

for (i, &(x, y)) in coords.iter().enumerate() {
    println!("Atom {}: ({:.2}, {:.2})", i, x, y);
}
```

### Parse V3000 with coordinates (NEW)

```rust
use chematic_mol::parse_mol_v3000_with_coords;

let v3000_str = "...";
let (mol, metadata, coords) = parse_mol_v3000_with_coords(v3000_str)?;
// Now V3000 coordinates are available for 2D layout!
```

### Handle Y-coordinate systems

```rust
use chematic_mol::{parse_cml, parse_cdxml};
use chematic_depict::render_svg;

// CML: Y-up → must negate for SVG (Y-down)
let (mol, cml_coords) = parse_cml(cml_str)?;
let svg_coords: Vec<(f64, f64)> = cml_coords.iter()
    .map(|(x, y)| (x, -y))  // Negate Y
    .collect();
render_svg(&mol, &svg_coords)?;

// CDXML: Y-down → compatible with SVG, no conversion needed
let (mol, cdxml_coords) = parse_cdxml(cdxml_str)?;
render_svg(&mol, &cdxml_coords)?;  // Use as-is
```

### Write molecule

```rust
use chematic_mol::write_mol;

let mol_str = write_mol(&mol, &metadata);
println!("{}", mol_str);  // V2000 format
```

## File Format Support Matrix

| Format | Read | Write | Coords | Notes |
|--------|------|-------|--------|-------|
| MOL V2000 | ✅ | ✅ | — | Charged atoms, isotopes, stereo |
| MOL V3000 | ✅ | ✅ | — | Extended atoms, COLLECTION |
| SDF (multi) | ✅ | ✅ | ✅ | Per-record properties |
| CML | ✅ | ✅ | ✅ | Y-up convention (documented) |
| CDXML | ✅ | — | ✅ | Y-down convention (documented) |
| MDL RXN | ✅ | ✅ | — | V2000 reactants/products |

## Coordinate System Handling

**Critical for correct rendering**: Y-axis conventions differ between formats.

| Source | Y-Direction | Required Transform | Documentation |
|--------|-------------|-------------------|---|
| `compute_layout()` | Down (SVG) | None | crates/chematic-depict/src/layout.rs |
| `parse_cml()` | Up (chemical) | Negate: y → -y | crates/chematic-mol/src/cml.rs |
| `parse_cdxml()` | Down (ChemDraw) | None | crates/chematic-mol/src/cdxml.rs |
| `parse_mol_with_coords()` | Depends on source | Check origin | crates/chematic-mol/src/mol2000.rs |
| `parse_mol_v3000_with_coords()` (NEW) | Depends on source | Check origin | crates/chematic-mol/src/mol3000.rs |

## Error Types

All error types implement `std::fmt::Display` and `std::error::Error`:

```rust
use std::error::Error;

fn parse_any(input: &str) -> Result<(), Box<dyn Error>> {
    let (mol, _) = chematic_mol::parse_mol(input)?;
    Ok(())
}
```

## Crate Dependencies

- `chematic-core` — Atom, Bond, Molecule types
- `chematic-smiles` — SMILES parsing (for mol files with SMILES data)

**Zero FFI**: No rdkit-sys, no OpenBabel, no C/C++ libraries.

## Performance

| Format | Size | Parse Time | Notes |
|--------|------|-----------|-------|
| V2000 (100 atoms) | ~5 KB | ~100 µs | Typical molecule |
| V3000 (200 atoms) | ~15 KB | ~200 µs | Extended format |
| SDF (1000 mols) | ~5 MB | ~500 ms | Streaming reader |
| CML (50 atoms) | ~10 KB | ~50 µs | XML parsing |

## Testing

```bash
cargo test --lib
# 63 tests covering all formats, coordinates, round-trips
```

## Version History

**v0.1.32** (2026-06-07):
- NEW: `parse_mol_v3000_with_coords()` recovers 2D coordinates from V3000 atom blocks
- NEW: Y-coordinate system documentation (CML, CDXML, compute_layout)
- NEW: Error trait implementations (CmlError, CdxmlError, Mol2Error, RxnParseError)
- 2 new V3000 coordinate tests

**v0.1.30** (2026-06-07):
- V3000 stereo-group COLLECTION support
- Full V3000 line-continuation handling

## License

MIT OR Apache-2.0

## Contributing

Contributions welcome! File format issues, coordinate system problems, or new format support.

