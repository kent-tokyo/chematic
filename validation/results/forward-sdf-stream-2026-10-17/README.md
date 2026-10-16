# Binary SDF stream supplier evidence

Pinned base: `1b4603db0a20af91ca598415b48e13c81bee0578`.
This source candidate adds native/compatibility `ForwardSDMolSupplier` and
fixes optional trailing SD-field header information in the shared parser.
The upstream references are RDKit issue #9601 and open proposed PR #9615,
checked live on 2026-10-16 UTC (2026-10-17 Japan).

`before-change.json` records the absence of either supplier on the pinned
base and the intermediate native-module hash that reproduced the SD-field
loss. Its test log shows 251 passes and two failures with the new stream
implementation but the unchanged base field parser. This is not an exact
rebuild of every base API. The intermediate reader also used the existing
lightweight mode; external comparisons led the final new supplier to use
the diagnostic-complete reader, preserving V2000 wedge chirality.

`rdkit-baseline.json` and `rdkit-comparison.json` record RDKit 2025.09.3.
The invalid-input child constructed a supplier for `42` and timed out during
iteration after two seconds; logging was disabled and the subprocess runner
killed/reaped the child. The issue and proposed fix remain upstream work.
The 80 passing comparisons cover ten graphs, V2000/V3000, BytesIO/gzip and
native/compatibility APIs. Both graph outputs are normalized through RDKit;
titles and Activity values are compared separately. Isotopes, charges,
tetrahedral stereo and both E/Z configurations are included.

Local checks:

- Python 3.13.6: 255 passed, zero failures/skips (66 new cases).
- `cargo test -p chematic-mol --lib --offline`: 788 passed.
- `cargo test -p chematic-mol --test sdf_v3000_recognition --test sdf_selections --offline`: 17 passed.
- `cargo clippy -p chematic-py -p chematic-mol --all-targets --offline -- -D warnings`: passed.
- Format, refactor budgets, source stubs, binding inventory/manifest, naming,
  unsafe-surface and strict MkDocs checks: passed.

Builds used one job, no incremental compilation/debug information, opt-level
1 and `strip=none`. Native Python builds used the macOS dynamic-lookup linker
flags; `summary.json` and `source-sha256.json` pin the artifact and source.
Compressed logs retain the actual local results. No published wheel or
completed hosted-CI/main-integration claim is made here.

Reproduce the external comparison with the built source package on PYTHONPATH:

```sh
python3 validation/rdkit_issues/forward-sdf/probe.py --compare
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 -m pytest -q \
  crates/chematic-py/tests/test_forward_sdf.py \
  crates/chematic-py/tests/test_io.py \
  crates/chematic-py/tests/test_rdkit_compat.py
```

The existing path-based supplier's separate lightweight V2000 perception
behavior was observed during this audit and left outside this API addition.
The new supplier's graph comparisons use the complete perception path.
