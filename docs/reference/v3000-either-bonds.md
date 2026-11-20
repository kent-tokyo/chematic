# V3000 either bonds and unknown alkene configuration

[Open Babel #2429](https://github.com/openbabel/openbabel/issues/2429) requests
retention of V3000 bond `CFG=2`. The
[BIOVIA CTfile specification](https://www.3ds.com/assets/invest/2026-01/biovia_chemistry_2026_ctfile_formats.pdf)
defines that bond configuration as either (printed page 13). Chematic retains
this marker on single and double bonds in `MolMetadata`. A wavy single bond
also prevents coordinate-based E/Z assignment to double bonds adjacent to
either endpoint. A remote wavy bond does not suppress unrelated alkene stereo.

| Path | Preserved meaning |
|---|---|
| Rust `read_mol_v3000_with_diagnostics` followed by `write_mol_v3000` with its metadata | Original single-bond marker and endpoint order, with normalized output IDs |
| Rust `write_mol_v3000_with_conformer` with source metadata | Original either marker retained in the 3D block |
| Python `from_mol_v3000` / `Mol.to_mol_v3000` | Explicitly unknown alkene remains unknown on reparse; output may use double-bond `CFG=2` |
| WASM `mol_from_v3000_block` / `to_mol_v3000_block` | Same unknown-alkene contract through generated 2D coordinates |

Graph-only bindings retain the semantic fact that an alkene was explicitly
unknown, but not the original marker's carrier-bond position. Their writers may
therefore emit `CFG=2` directly on the double bond. An ordinary graph without
that marker remains free to express caller-supplied E/Z geometry; it is not
silently changed to unknown. Already declared E/Z stereo follows the existing
depiction path. Use Python `roundtrip_mol_v3000_source` or WASM
`roundtrip_mol_v3000_block` for exact unedited source retention; those APIs
predate this change.

This addresses alkene interpretation and V3000 rewriting. It does not add SVG
wavy-bond rendering, general tetrahedral either-bond perception, or changes to
V2000 output.
Atom-line CFG is separate and remains outside this bond configuration profile.

The fixed packet `validation/openbabel-2429-2026-10-24/` contains 16 authored
wavy-alkene variants (left/right carrier, reversed endpoints, mirrored geometry,
2D/3D), a direct double-either control and an unmarked geometry control. Reading,
writing and reparsing preserve the expected defined/unknown alkene state in
18/18 patched Python cases and 18/18 independent RDKit 2026.09.1 cases. The
prior local Chematic runtime and Open Babel 3.2.1 each match 1/18 in their tested
pipelines. Literal SMILES spelling is not compared. These are local fixed-input
regression results, not general stereo coverage or published-release evidence.
