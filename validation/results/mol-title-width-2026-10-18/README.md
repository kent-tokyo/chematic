# Native MOL title width — 2026-10-18

Implements the bounded title behavior requested by Open Babel #2112 in native
MOL/SDF V2000 and V3000 output. The source metadata is never mutated; written
long/control-containing titles deliberately lose content. No unrestricted
output option or full CTfile conformance claim.

## Results

- `before.log`: restored verbatim title emission fails six regression tests;
  four ordinary/buffer cases pass. This is a behavior-restoration control, not
  a build of an entirely unchanged baseline checkout.
- `rust-tests.log`: 864 passing tests (815 library + 11 new + 38 existing
  MOL/SDF stereo integration tests). The frozen comparison writers stay unchanged.
- `after.log`: freshly built native writer output for 7 title fixtures across
  18 public output paths (126 cases), consumed by the independent reader.
- `oracle-summary.json`: RDKit 2026.09.1 strict reads all 126 MOL blocks,
  checks names and the original C–O graph; 35 SD property reads pass. Each
  non-title body is byte-identical to the same writer's ordinary-title case.
  These are synthetic methanol title fixtures, not a general chemical corpus.
- `openbabel-diagnostic.json`: installed Open Babel 3.2.1 emits a synthetic
  100-byte ASCII methane title unchanged at 100 bytes. This separate control
  reproduces the reported defect and is not included in cross-read counts.
- `clippy.log`, `format.log`, `source-gates.log`, `mkdocs.log`: passing local checks.

## Reproduce

From the repository root, with Rust dependencies cached and RDKit installed:

```sh
export CARGO_TARGET_DIR=/private/tmp/chematic-mol-title-target
export CARGO_BUILD_JOBS=1 CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0 CARGO_INCREMENTAL=0
cargo test -p chematic-mol --offline --lib --test mol_title_width --test square_planar_mol_io --test sdf_unspecified_double_stereo --test stereo_reader_integration
cargo test -p chematic-mol --offline --test mol_title_width -- --nocapture > validation/results/mol-title-width-2026-10-18/after.log 2>&1
PYTHONPATH=/private/tmp/chematic-sparse-tanimoto-rdkit python3.13 validation/results/mol-title-width-2026-10-18/oracle.py
cargo clippy -p chematic-mol --offline --lib --test mol_title_width --no-deps -- -D warnings
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-mol-title-site
```

The host-specific RDKit path and `/opt/homebrew/bin/obabel` can be replaced
with matching locally installed tools. No upstream attachments or code were
copied. The new algorithm and all fixtures are authored in this change.
Python/WASM extension binaries were not rebuilt; the Rust code is compiled
and tested directly. Hosted CI remains a separate integration check.
