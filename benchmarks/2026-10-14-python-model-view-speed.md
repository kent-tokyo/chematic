# Python core-model stereo preparation follow-up — 2026-10-14

Dropping unused stereo preparation retains all checked outputs, but no incremental lane passes the strict 20/20 gate. First-use TPSA medians improve only about 1% on both slices; this change does not close the substantial first-use deficits. All 27/27 eligible parse-inclusive comparator lanes pass. Full results, including regressions, follow.

This packet follows the [first-use perception candidate](2026-10-14-python-first-use-followup.md).
The exact previous source wheel is the baseline. RDKit 2026.09.1 and COSMolKit
0.5.0-rc.22 remain fixed, as do both unfiltered 1,000-row slices and all seven
operations. First-use, repeated-use and actual parse-inclusive timings remain
separate. These are source wheels on one macOS arm64 host, not published
packages or evidence of toolkit-wide superiority.

## Implementation and chemical boundary

`with_rdkit_model_molecule` converts a sanitized port molecule to a core
`Molecule` for descriptor and fingerprint callbacks. That core representation
does not carry CIP ranks or the port's potential-stereocenter flags. Previously
the conversion requested the complete legacy stereo pass, including ranking
used only to flag potential centers when no stereo was specified.

The conversion now uses the existing `legacy_stereo_perception_unflagged`
path. Specified atom tags and double-bond directions still receive their
established cleanup. Sanitization, hydrogen removal/retention, charges,
aromaticity and conversion remain in place. Public CIP and stereocenter-count
APIs retain complete perception. No descriptor or ring calculation moves into
SMILES parsing, and no result is hard-coded for a benchmark molecule.

Only `crates/chematic-smiles/src/rdkit/mod.rs` changes in production relative
to the previous candidate. Three new tests compare complete and optimized
core views, including atoms, bonds, neighbor/stereo order and directions,
across specified/unspecified/invalid stereo, non-tetrahedral tags, isotopes,
explicit/added hydrogen, cleanup, metal coordination, radicals, errors and
graph edits.

The optimization concerns unused preparation. The descriptor formulas remain
unchanged; RDKit's [MolSurf implementation](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Descriptors/MolSurf.cpp)
also shows that TPSA and Labute ASA consume chemical atom/bond environments.
The installed 2026.09.1 package is the numeric oracle; upstream source is
context, not the measured comparator artifact.

## Corpus, protocol and provenance

Both slices contain all 1,000 nonempty source rows at offsets 0 and 5000 of
`validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`, SHA256
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`.
They are exposed development inputs, not sealed holdouts. No row, timing block
or observed regression is removed.

The unchanged runner rotates four arms across twenty sequential fresh-process
blocks per slice. Native return objects are retained in timed loops;
normalization and comparisons occur afterward. First use excludes parsing and
includes lazy work; repeated use follows one untimed call per fresh molecule;
parse-inclusive loops actually parse and call the operation. All three
boundaries are reported because eager/lazy preparation differs by engine.

A strict win requires equal outputs, the paired 95% bootstrap confidence
interval's lower bound above one, and twenty of twenty faster blocks. Paired
ratios use deterministic 10,000-resample bootstrap intervals. COSMolKit output
differences suppress ranking. Whole-process RSS includes imports, parsing,
preparation and all seven operations; it is not streaming or operation-only
memory evidence.

The wheel was rebuilt in a dedicated initially empty Cargo target with Rust
1.97.0, release opt-level 3, LTO, one codegen unit and `native-inchi`. The exact
old wheel is retained as the comparator baseline. The source manifest records
the pinned base, three prerequisite patches, the new patch, frozen build files,
whole-file workspace integrity markers, wheel hashes and comparator native
library hashes. Other dirty Open Babel, Python API, matrix, chemistry and 3D
changes are excluded from the measured wheel. Source version remains 1.1.3;
wheel hashes distinguish the artifacts.

The wrapper verifies frozen sources, scoped workspace files, tools, patches,
wheels and comparator native libraries before and after each slice. Every
worker checks its installed extension against the supplied wheel. Formal
timing starts after compilation and validation finish. Exploratory workers
run during test compilation are retained as development observations and
never counted as formal blocks.

## Formal measurements

Ratios are paired comparator / candidate medians with 95% bootstrap intervals. Bold cells pass the equal-output, lower-bound and 20/20 gates.

### First use

| Operation | RDKit, slice 0 | COSMolKit, slice 0 | RDKit, slice 5000 | COSMolKit, slice 5000 |
|---|---|---|---|---|
| tpsa | **1.449x [1.434, 1.460]; 20/20** | 0.504x [0.500, 0.513]; 0/20 | 0.319x [0.317, 0.325]; 0/20 | 0.087x [0.084, 0.089]; 0/20 |
| labute_asa | 0.698x [0.691, 0.711]; 0/20 | 0.327x [0.321, 0.336]; 0/20 | 0.156x [0.152, 0.158]; 0/20 | 0.063x [0.062, 0.064]; 0/20 |
| num_rings | 0.019x [0.019, 0.020]; 0/20 | Not ranked: outputs differ | 0.097x [0.095, 0.102]; 0/20 | 0.078x [0.074, 0.082]; 0/20 |
| morgan2_2048 | **3.146x [3.127, 3.167]; 20/20** | **2.044x [2.023, 2.057]; 20/20** | **1.399x [1.390, 1.424]; 20/20** | 0.898x [0.888, 0.903]; 0/20 |
| morgan2_chiral | **1.258x [1.247, 1.268]; 20/20** | 0.824x [0.821, 0.827]; 0/20 | **1.439x [1.425, 1.452]; 20/20** | 0.909x [0.904, 0.916]; 0/20 |
| morgan3_2048 | **3.351x [3.329, 3.393]; 20/20** | **2.254x [2.234, 2.267]; 20/20** | **1.617x [1.604, 1.626]; 20/20** | **1.061x [1.056, 1.068]; 20/20** |
| maccs | **1.882x [1.878, 1.887]; 20/20** | **2.165x [2.161, 2.172]; 20/20** | **1.627x [1.621, 1.631]; 20/20** | **1.742x [1.733, 1.745]; 20/20** |

### Repeated use

| Operation | RDKit, slice 0 | COSMolKit, slice 0 | RDKit, slice 5000 | COSMolKit, slice 5000 |
|---|---|---|---|---|
| tpsa | **9.965x [9.891, 10.112]; 20/20** | **1.725x [1.686, 1.842]; 20/20** | **9.736x [9.541, 10.057]; 20/20** | **1.321x [1.276, 1.330]; 20/20** |
| labute_asa | **2.710x [2.658, 2.734]; 20/20** | **1.412x [1.351, 1.438]; 20/20** | **2.654x [2.572, 2.700]; 20/20** | **1.150x [1.140, 1.177]; 20/20** |
| num_rings | **1.721x [1.682, 1.757]; 20/20** | Not ranked: outputs differ | **1.763x [1.698, 1.783]; 20/20** | 1.086x [1.041, 1.113]; 18/20 |
| morgan2_2048 | **5.229x [5.168, 5.282]; 20/20** | **3.230x [3.200, 3.240]; 20/20** | **4.777x [4.745, 4.815]; 20/20** | **2.904x [2.892, 2.941]; 20/20** |
| morgan2_chiral | **4.976x [4.962, 5.016]; 20/20** | **3.089x [3.085, 3.103]; 20/20** | **4.752x [4.690, 4.799]; 20/20** | **2.859x [2.849, 2.884]; 20/20** |
| morgan3_2048 | **4.870x [4.834, 4.893]; 20/20** | **3.137x [3.109, 3.161]; 20/20** | **4.514x [4.492, 4.564]; 20/20** | **2.899x [2.886, 2.918]; 20/20** |
| maccs | **1.948x [1.944, 1.954]; 20/20** | **2.243x [2.238, 2.249]; 20/20** | **1.725x [1.718, 1.736]; 20/20** | **1.852x [1.840, 1.861]; 20/20** |

### Parse inclusive

| Operation | RDKit, slice 0 | COSMolKit, slice 0 | RDKit, slice 5000 | COSMolKit, slice 5000 |
|---|---|---|---|---|
| tpsa | **12.986x [12.873, 13.129]; 20/20** | **12.091x [12.051, 12.207]; 20/20** | **4.135x [4.104, 4.152]; 20/20** | **4.279x [4.236, 4.292]; 20/20** |
| labute_asa | **12.307x [12.181, 12.368]; 20/20** | **11.615x [11.536, 11.696]; 20/20** | **3.946x [3.929, 3.970]; 20/20** | **4.213x [4.194, 4.220]; 20/20** |
| num_rings | **8.791x [8.707, 8.871]; 20/20** | Not ranked: outputs differ | **11.596x [11.484, 11.618]; 20/20** | **12.679x [12.563, 12.703]; 20/20** |
| morgan2_2048 | **8.642x [8.588, 8.693]; 20/20** | **7.606x [7.589, 7.628]; 20/20** | **4.259x [4.237, 4.277]; 20/20** | **4.115x [4.074, 4.133]; 20/20** |
| morgan2_chiral | **4.440x [4.426, 4.466]; 20/20** | **3.910x [3.900, 3.932]; 20/20** | **4.250x [4.243, 4.274]; 20/20** | **4.124x [4.094, 4.134]; 20/20** |
| morgan3_2048 | **7.999x [7.952, 8.034]; 20/20** | **6.978x [6.948, 7.002]; 20/20** | **4.187x [4.168, 4.213]; 20/20** | **4.011x [3.997, 4.037]; 20/20** |
| maccs | **2.109x [2.102, 2.118]; 20/20** | **2.376x [2.372, 2.380]; 20/20** | **1.787x [1.782, 1.796]; 20/20** | **1.917x [1.914, 1.923]; 20/20** |

### Incremental change from the previous source wheel

| Operation / mode | Previous / candidate, slice 0 | Previous / candidate, slice 5000 |
|---|---|---|
| tpsa/first_use | 1.012x [1.010, 1.017]; 17/20 | 1.012x [1.007, 1.016]; 17/20 |
| tpsa/repeated_use | 0.995x [0.988, 1.006]; 9/20 | 0.985x [0.962, 0.997]; 5/20 |
| tpsa/parse_inclusive | 1.006x [0.995, 1.013]; 12/20 | 1.006x [1.004, 1.011]; 18/20 |
| labute_asa/first_use | 1.007x [0.993, 1.015]; 12/20 | 1.010x [1.004, 1.017]; 18/20 |
| labute_asa/repeated_use | 0.998x [0.988, 1.006]; 8/20 | 0.998x [0.975, 1.013]; 9/20 |
| labute_asa/parse_inclusive | 0.996x [0.989, 1.004]; 7/20 | 1.006x [1.002, 1.015]; 16/20 |
| num_rings/first_use | 1.005x [0.997, 1.014]; 13/20 | 1.018x [1.009, 1.024]; 18/20 |
| num_rings/repeated_use | 1.012x [0.990, 1.024]; 12/20 | 1.001x [0.985, 1.024]; 10/20 |
| num_rings/parse_inclusive | 1.007x [0.999, 1.016]; 13/20 | 0.996x [0.987, 1.009]; 8/20 |
| morgan2_2048/first_use | 1.004x [0.997, 1.008]; 13/20 | 1.004x [0.991, 1.008]; 11/20 |
| morgan2_2048/repeated_use | 1.002x [0.997, 1.006]; 11/20 | 1.001x [0.999, 1.004]; 11/20 |
| morgan2_2048/parse_inclusive | 1.001x [0.998, 1.006]; 10/20 | 0.995x [0.990, 0.998]; 5/20 |
| morgan2_chiral/first_use | 0.990x [0.987, 0.995]; 3/20 | 0.997x [0.993, 1.000]; 5/20 |
| morgan2_chiral/repeated_use | 0.993x [0.986, 0.999]; 5/20 | 0.994x [0.991, 1.001]; 6/20 |
| morgan2_chiral/parse_inclusive | 0.999x [0.994, 1.002]; 10/20 | 0.997x [0.987, 1.006]; 8/20 |
| morgan3_2048/first_use | 1.002x [0.999, 1.006]; 12/20 | 0.997x [0.994, 1.004]; 7/20 |
| morgan3_2048/repeated_use | 1.004x [1.000, 1.011]; 14/20 | 1.005x [1.001, 1.015]; 14/20 |
| morgan3_2048/parse_inclusive | 1.001x [0.998, 1.004]; 13/20 | 1.001x [0.997, 1.011]; 11/20 |
| maccs/first_use | 0.999x [0.997, 1.001]; 9/20 | 0.998x [0.992, 1.001]; 8/20 |
| maccs/repeated_use | 1.000x [0.998, 1.002]; 11/20 | 0.996x [0.991, 1.002]; 7/20 |
| maccs/parse_inclusive | 0.999x [0.996, 1.002]; 7/20 | 0.994x [0.993, 0.998]; 4/20 |

### Whole-process memory

| Slice | Candidate median peak MiB | RDKit MiB | COSMolKit MiB |
|---|---:|---:|---:|
| 0 | 46.54 | 106.56 | 105.06 |
| 5000 | 41.04 | 88.38 | 82.81 |


## Correctness and validation

Normal 10,000 and unusual 4,072 inputs retain 112,576 first/repeated descriptor
outcomes and 56,288 Morgan outcomes byte-for-byte against the prior candidate.
Float bits, errors and nonfinite outcomes remain in the dumps. Identical
nonfinite outcomes are retention evidence, not successful numeric agreement.
The previously recorded RDKit residuals remain; matching the old source does
not erase them. Existing baseline dumps are reused only after their hashes
were verified against the earlier packet.

Isolated release tests for SMILES, chemistry and fingerprints: 1,959 passed,
zero failed and eight existing ignored cases. All three new core-view tests
pass. Release Clippy for the changed SMILES crate and all targets passes with
warnings denied. The selected Python molecule, descriptor, RDKit, fingerprint,
pickle and cross-binding suites pass 590 tests. Five unchanged benchmark-gate
tests pass. Changed Rust files pass rustfmt and the workspace diff has no
whitespace errors.

The isolated copy initially lacked compile-time and runtime reference fixtures.
Those files were restored from the fixed base commit before the successful
full test run; no production source changed to resolve the missing files.
The candidate reconstructs from the base and four patches with all 894 frozen
Rust/workspace-manifest files matching. Paired summaries and retained archive
hashes are independently verified in the final parity record.

## Retained evidence

- [Source, wheel, compiler and comparator manifest](2026-10-14-python-model-view-speed-source.json)
- [Original-slice raw blocks](2026-10-14-python-model-view-speed-slice0-paired20.json)
- [Additional-slice raw blocks](2026-10-14-python-model-view-speed-slice5000-paired20.json)
- [Exact retention, reconstruction and validation](2026-10-14-python-model-view-speed-parity.json)
- [Reconstructable candidate patch](2026-10-14-python-model-view-speed-candidate.patch)
- [Unchanged runner](../scripts/bench_python_speed_frontier.py)

Reconstruct from `ffe682c32cb4c76de0d309649427202a7d8e956d`, the three prerequisite
patches and this candidate patch. Build with the recorded toolchain/settings,
install each exact source wheel into its own directory, and run the unchanged
frontier runner with `--limit 1000 --blocks 20` at each recorded offset. The
retained wrapper adds source and library checks around that runner. This is
local source evidence, not a public release, hosted CI, commit or merge.

Toolkit-wide speed superiority remains unproven. First-use preparation gaps,
broader operations, reactions and 3D, malformed inputs, fresh corpora,
additional hosts and exact published artifacts require separate comparisons.
