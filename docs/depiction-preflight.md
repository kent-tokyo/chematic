# SVG depiction publication preflight

`chematic-depict` provides an opt-in preflight for applications that publish
SVG generated from caller-supplied coordinates:

```rust
let report = chematic_depict::preflight_svg(
    &mol,
    &layout,
    &chematic_depict::RenderOptions::default(),
    &chematic_depict::PreflightLimits::default(),
);
if !report.ready {
    // show report.diagnostics and do not publish
}
```

`preflight_svg_json` serializes the same `PreflightReport` for JSON consumers.
Diagnostic order follows molecule atom/bond order and paths such as
`atoms[3].geometry`, `atoms[2].label`, and `bonds[0]/bonds[4]` are stable. The
report also contains the effective depiction style and an FNV-1a fingerprint of
the ordered molecule, coordinates, and style.

The checks are deliberately conservative: non-finite or degenerate geometry,
negative style values, obvious viewBox clipping, label/atom and label/label
overlap, bond crossings, and explicit atom/bond/label limits fail `ready`.
`PreflightLimits` bounds the in-memory graph; the WASM
`preflight_smiles_json(smiles, width, height)` entry point additionally uses the
existing 1 MiB serialized-input and 10,000-atom limits.

The style records `sans-serif`, 12 px text, 1.5 px bonds, and the configured
padding/viewport. Text boxes use fixed conservative estimates. Font selection,
kerning, hinting, and final pixel clipping remain renderer-authoritative, so a
successful preflight is not a browser screenshot or PDF parity claim. The
preflight module has no PDF dependency and does not change the SVG renderer.
