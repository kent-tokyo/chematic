# Per-input SDF gzip validation

Source base: `d03106b59fdd4d06f40ab267bb540e20508bd988`.
Open Babel issue #3030 was checked live on 2026-10-16 UTC / 2026-10-17 Japan.
The original C++ probe reproduces the sticky gzip flag on Open Babel 3.2.1:
gzip succeeds, replacement plain string/stream reads fail, fresh conversions
read the same plain inputs successfully. Outputs are in `openbabel-probe.tsv`
and `.stderr`. The probe uses an independently authored two-atom SMILES input.

The previous chematic checked SDF API lacked a compression keyword. There is
no corresponding shared conversion session to demonstrate the same defect.
This candidate adds plain/gzip/auto selection per reader and preserves the
plain default. It is a capability addition inspired by the report, not a claim
that the Open Babel session bug already existed in chematic.

`before-change.json` pins the retained older native module and records an
exact Python IO source-blob match between its commit and the current base.
The six initial gzip/path-protocol tests failed with an unexpected compression
keyword. This is a focused API baseline, not a full rebuild of current main.

Local validation:

- `cargo test -p chematic-mol --features crystal --lib --test sdf_compression --test sdf_v3000_recognition --test sdf_selections --offline`: 853 library + 9 new compression + 17 existing SDF integration tests passed.
- Python 3.13.6: 74 passed, including 33 new compression cases; one existing
  OGB reference test skipped because `ogb` was not installed.
- 80 records match the plain checked reader across ten graphs, V2000/V3000
  and four policy/input variants. Names, SD fields, SMILES and stereo diagnostics
  match; multi-member data includes a UTF-8 character split between members.
  These are compression-preservation comparisons, not external chemistry oracles.
- Clippy for chematic-mol/chematic-py, all targets, warnings denied: passed.
- Format, source-stub, refactor-budget, binding inventory/manifest, naming,
  unsafe-surface and strict MkDocs checks: passed.

Rust tests cover short/non-seekable input, interrupted detection, explicit
plain byte preservation, independent reader state after full/partial/failing
inputs, split-member records, CRC/truncation/trailing-junk errors, decoded-byte
limits and malformed-record positions. Python tests cover corresponding file
APIs, path protocols and argument validation. Trailer errors can follow already
emitted records; complete iteration is required for whole-input integrity.
Inherited native SDF limits are after-collection checks, not hard peak-memory
ceilings. Other file APIs, compressed output and generic conversion sessions
are outside this candidate.

The native artifact/source hashes are in `summary.json` and `source-sha256.json`.
Builds used one job, no incremental compilation/debug information, opt-level 1,
strip=none and macOS dynamic-lookup linker flags. Compressed logs preserve actual
results. No package-publication or completed hosted-CI/main-merge claim is made.

Reproduce with the built source package on PYTHONPATH:

```sh
python3 validation/openbabel_issues/sdf-gzip/compare.py
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 -m pytest -q -rs \
  crates/chematic-py/tests/test_sdf_compression.py \
  crates/chematic-py/tests/test_io.py \
  crates/chematic-py/tests/test_ogb_graph.py
```
