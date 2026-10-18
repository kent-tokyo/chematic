# Bounded MOL and SDF titles

[Open Babel #2112](https://github.com/openbabel/openbabel/issues/2112)
reports MOL titles longer than the traditional 80-column header limit.
chematic now applies the same title policy to native V2000 and V3000 writers,
including 2D/3D output, stereo-reporting writers, reused buffers and SDF wrappers.

- The first header line contains at most **80 UTF-8 bytes**. ASCII names up to
  80 bytes remain unchanged. This conservative byte limit also fits an
  80-character reader; it is not a claim that CTfile mandates UTF-8.
- Unicode scalars are never split. The first scalar that cannot fit ends the
  prefix; later characters are not scanned to fill unused bytes.
- Each control character (including CR, LF, tab, NUL, DEL and NEL), and each
  Unicode line/paragraph separator, becomes one ASCII space before measuring
  the byte budget. CRLF therefore becomes two spaces.
- Empty names still produce the required blank first line.
- The caller's `MolMetadata.name` stays unchanged. Output truncation and
  control replacement are deliberate, lossy serialization. Callers that need
  the original name after reading the file must save it as an SD property or
  another sidecar. Stereo-loss reports continue to report stereo only.

```rust
use chematic_mol::{MolMetadata, write_mol};

let mol = chematic_smiles::parse("CO").unwrap();
let metadata = MolMetadata::default().with_name(&"0123456789".repeat(10));
let block = write_mol(&mol, &metadata);
assert_eq!(block.lines().next().unwrap().len(), 80);
assert_eq!(metadata.name.len(), 100);
```

This changes serialization defaults for long/control-containing titles.
There is no unrestricted-title output option. Readers still accept legacy
long titles. The third header line (comment), V3000 logical-line wrapping,
SD property values and other format limits retain their existing behavior;
this is not a general CTfile conformance validator.

## Validation

`crates/chematic-mol/tests/mol_title_width.rs` exercises native output routes,
boundaries, UTF-8, controls, large input, source immutability and buffer reuse.
The independent-reader run uses newly compiled Rust output, not an older
Python extension. Its fixtures, commands and results are recorded in
`validation/results/mol-title-width-2026-10-18/`.
