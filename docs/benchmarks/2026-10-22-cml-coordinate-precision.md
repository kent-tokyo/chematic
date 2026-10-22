# CML supplied 2D coordinate precision

This source change removes fixed four-decimal rounding from `write_cml`.
Finite supplied binary64 values are written with Rust's round-trip decimal
format. This preserves the represented coordinate value; it does not improve
the accuracy of a calculated geometry or promise original decimal spelling.

## Paired development result

The original eight-input packet uses two distinct atoms, C and O, with a single
bond. Layouts isolate fractional positions, signs, translations and small
coordinates near the origin. They are synthetic drawing positions, not
optimized geometries. The packet and acceptance tolerance were authored before
the baseline/candidate measurement. It is development evidence, not an
independent or sealed holdout.

Every output is inspected by Python 3.13.6 ElementTree. The observer matches the
unique C/O elements, checks their neutral/non-isotopic state and the single
bond, then compares their supplied 2D positions. IDs and serialization order
may change. Missing coordinates are not substituted with zero. The acceptance
tolerance is **1e-10 Å absolute**, with no relative tolerance, alignment or
rescaling. Exact binary64 value equality is a separate diagnostic.

| Outcome, eight inputs | Before | Candidate | Open Babel 3.2.1 |
|---|---:|---:|---:|
| Graph and coordinates within tolerance | 2 | 8 | 3 |
| Coordinate error above tolerance | 6 | 0 | 5 |
| Exact represented coordinate values | 2 | 8 | 2 |
| Process/read/write/observer failure | 0 | 0 | 0 |

| Input | Before maximum error (Å) | Candidate | Open Babel maximum error (Å) |
|---|---:|---:|---:|
| Dyadic control | 0 | 0 | 0 |
| Decimal control | 0 | 0 | 0 |
| Fractional layout | 4.3211e-5 | 0 | 2.1099e-7 |
| Opposite signs | 3.2110e-5 | 0 | 1.0988e-7 |
| Near origin | 2.3457e-8 | 0 | 2.3457e-8 |
| Translated layout | 4.3211e-5 | 0 | 2.1106e-7 |
| Near rounding boundary | 5.0000e-5 | 0 | 1.1102e-15 |
| Long decimals | 3.3333e-5 | 0 | 3.3333e-7 |

Open Babel runs with `-icml -ocml`, default options. Its installed `-L cml`
inventory and [official format manual](https://openbabel.org/docs/FileFormats/Chemical_Markup_Language.html)
were inspected; this is a comparison with those defaults, not all possible
Open Babel settings or APIs. No timing or general toolkit superiority follows
from these eight layouts. The small errors are representation losses, not
evidence of a chemically significant geometry change.

## Public boundary and validation

The measured chematic route is the Rust `parse_cml` → `write_cml` adapter in
`file_io_semantic_probe`, not the generic CLI conversion command. Python
`Mol.to_cml(coords)` delegates to the same writer, and its supplied positions
are independently read with ElementTree in the source-built extension tests.
Python `from_cml` and WASM graph handles do not retain source coordinates.
WASM `to_cml` generates a layout; that is not source-coordinate preservation.
No new binding API, input-coordinate storage, 3D or metadata support is added.

Local validation passed 1,478 Rust MOL unit/integration/doc tests, 143 Python
coordinate/named-CML/collection/hydrogen tests (13 new writer cases), and 38
observer/bounded-runner tests. All-target MOL Clippy with warnings denied,
formatting, refactor/script-layout/product-name checks and strict MkDocs also
passed. These are source checks on this host, not published-wheel validation.

Separate Rust tests pass signed zero, subnormals, the normal/subnormal boundary,
maximum finite magnitude and adjacent binary64 values through both CML readers,
compare every coordinate bit, and serialize the read-back positions again.
These are numeric representation checks, not physically meaningful structures.
The writer's existing handling of missing/short coordinate slices and nonfinite
caller values is unchanged; the preservation guarantee concerns supplied finite
values. Output number spelling changes from fixed four-decimal strings.

Baseline source is `d288f6027d71f21130dd119f7cc7df8d58ab8daf`.
This PR is independent of the character-reference reader PR #1126 and does not
include that PR's changes. Both paired reports retain every source, command,
output, stderr, exit/resource status and terminal classification. The validation
index records source, tool, executable and log hashes. The runner bounds its
process streams; the small native output file is checked after completion, not
claimed as a generic hard file-stream bound.

## Reproduction

```sh
CARGO_TARGET_DIR=/tmp/cml-precision-target cargo build --locked \
  -p chematic-mol --example file_io_semantic_probe
python3 -m scripts.validation.check_cml_coordinate_precision \
  --probe /tmp/cml-precision-target/debug/examples/file_io_semantic_probe \
  --obabel /path/to/obabel --output /tmp/cml-precision-new.json
cargo test --locked -p chematic-mol --test cml_coordinate_precision
python3 -m unittest scripts.tests.test_cml_coordinate_precision
```

The runner refuses to overwrite an existing report and requires Open Babel
3.2.1. Rebuild a Python extension from the same source before running
`crates/chematic-py/tests/test_cml_coordinate_precision.py`; an installed wheel
does not attest this source change.

Evidence:
[fixture manifest](https://github.com/kent-tokyo/chematic/blob/fix/cml-coordinate-precision-20261022/validation/cml_coordinate_precision/cases.json),
[baseline report](https://github.com/kent-tokyo/chematic/blob/fix/cml-coordinate-precision-20261022/validation/results/cml-coordinate-precision-baseline-2026-10-22.json),
[candidate report](https://github.com/kent-tokyo/chematic/blob/fix/cml-coordinate-precision-20261022/validation/results/cml-coordinate-precision-candidate-2026-10-22.json),
[validation index](https://github.com/kent-tokyo/chematic/blob/fix/cml-coordinate-precision-20261022/validation/results/cml-coordinate-precision-validation-2026-10-22.json).
