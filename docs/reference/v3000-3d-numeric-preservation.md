# Native V3000 3D numeric preservation

The native V3000 conformer writers retain supplied finite x/y/z values as
binary64, including signed zero. They keep compact four-decimal fields when
those fields already reparse to the same bits; otherwise they emit round-trip
scientific notation. For example, `1.2345678901234567` no longer becomes
`1.2346` when a 3D conformer is written.

```rust
use chematic_mol::write_mol_v3000_with_conformer;
let output = write_mol_v3000_with_conformer(&molecule, &metadata, &conformer);
```

This applies to `write_mol_v3000_with_conformer` and
`write_mol_v3000_with_conformer_checked`. The existing checked function still
checks square-planar stereo/conformer consistency; it is not a general input,
metadata or finite-coordinate validator. Missing coordinate rows keep the
existing zero-fill behavior. The guarantee concerns supplied finite values.

Numeric tokens are compact, and long atom rows use V3000 continuation lines
between fields. The [MDL CTfile specification](https://mdli.com/downloads/literature/ctfile.pdf)
limits physical V3000 lines to 80 characters and defines dash continuations.
Coordinates and UTF-8/property tokens are not split. This change does not add
validation or a width guarantee for arbitrary oversized opaque property strings.

The native 2D V3000 writer, V2000 writers, immutable source-retention paths,
and Python/WASM topology/2D output APIs retain their existing behavior. This
change does not add a binding/CLI 3D conformer-output route. It does not make
all CTfile metadata, chemistry or stereochemistry lossless.

## Paired diagnostic — 2026-10-24

Nine original records isolate exact four-decimal controls, long mantissas,
tiny displacements, signed zero, large coordinates, one-ULP boundaries,
scientific-token continuation, and a pair of mirrored tetrahedral controls.
Coordinates are authored numeric strings and compared as binary64 bits.
Three additional finite extremes are native-only and excluded from this table.

Each path reads and regenerates the same V3000 source. The native paths use
source x/y/z → conformer → 3D writer. Open Babel 3.2.1 uses
`obabel -imol -omol -x3`, explicitly requesting V3000 through its
[documented output option](https://openbabel.org/docs/FileFormats/MDL_MOL_format.html).
Both outputs are read through RDKit 2026.09.1 with sanitization off, H removal
off and strict parsing on. Every source is first qualified through that same
reader against the authored numeric bits.

| Path | RDKit readable | Coordinates, exact bits | Indexed graph fields | Source CIP agreement | Joint |
|---|---:|---:|---:|---:|---:|
| Native candidate, plain or checked | 9/9 | 9/9 | 9/9 | 9/9 | 9/9 |
| Previous native `aeae3f53`, plain or checked | 9/9 | 2/9 | 9/9 | 9/9 | 2/9 |
| Open Babel 3.2.1, V3000 output | 9/9 | 2/9 | 9/9 | 9/9 | 2/9 |

Indexed fields are atomic number, charge, isotope, atom map, aromatic flag and
bond endpoints/orders. Charges, isotopes and atom maps are all zero in this
packet; nonzero coverage is outside this measurement. Stereo agreement compares atom-indexed RDKit CIP codes
with the source as read by RDKit; two records contain a tetrahedral center,
while seven contain none. This is a coordinate-preservation improvement,
not a measured stereochemistry advantage. A separate raw V3000 observer checks
numeric tokens, including continuations, independently of the external reader.
All candidate V30 lines in the numeric packet are within 80 UTF-8 bytes.

Native tests additionally cover three regenerations, finite extremes, missing
coordinates and deterministic float-bit/rounding-boundary samples. The frozen
3D comparison still checks every non-coordinate logical field while requiring
exact supplied numeric bits. Existing frozen 2D byte comparisons remain in place.

The packet is a numerical regression diagnostic, not a chemical-validity test,
independent holdout, speed measurement, format-wide ranking or published release.
Evidence and reproduction are in `validation/v3000_3d_numeric/README.md` and
`validation/results/v3000-3d-numeric-{baseline,candidate}-2026-10-24.json.gz`.
