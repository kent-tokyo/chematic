# Python sanitized-model aromatic-view reuse — 2026-10-14

First-use TPSA improves 1.032x [1.028, 1.044] and 1.068x [1.057, 1.077]
over the previous source wheel on the two slices, winning all 20 blocks.
Additional-slice first-use Labute ASA improves 1.065x [1.059, 1.070], also
20/20. The corresponding additional-slice read-plus-descriptor workflows
improve 1.045x and 1.047x in 20/20 blocks.

All 27/27 eligible parse-inclusive comparator lanes pass; first-use and
repeated-use results pass 14/27 and 26/27. First-use deficits and small
incremental regressions remain, so toolkit-wide superiority is not established.

This packet follows the [first-use perception candidate](2026-10-14-python-first-use-followup.md).
The exact previous source wheel is the baseline. RDKit 2026.09.1 and COSMolKit
0.5.0-rc.22 remain fixed, as do both unfiltered 1,000-row slices and all seven
operations. First-use, repeated-use and actual parse-inclusive timings remain
separate. These are source wheels on one macOS arm64 host, not published
packages or evidence of toolkit-wide superiority.

## Implementation and chemical boundary

Scalar RDKit-compatible descriptors first sanitize unusual inputs using the
existing pure Rust model port. Their temporary core molecule already carries
that model's aromaticity, bond orders, charges and hydrogen counts. A later
non-identity aromatic-view request previously perceived and sanitized this same
model again.

Scalar model callbacks now seed that temporary molecule's existing
`RdkitParityAromatic` cache with the known sanitized view when its ordinary
identity decision does not apply. The seed holds a clone with empty caches,
so it creates no reference cycle. Identity views need no clone. The descriptor
consumer keeps its existing normalization, including aromatic ring triples.
Generic model callbacks do not receive this seed. Mutation invalidation and
clone isolation use the existing derived-cache contract.

The conversion also keeps the [intermediate unflagged-stereo change](2026-10-14-python-model-view-speed.md):
core molecules cannot carry unused CIP ranks or potential-stereocenter flags,
so the existing unflagged legacy pass keeps specified stereo cleanup without
flag-only ranking. The intermediate packet's complete 20-block observations
are retained: that change alone had no strict incremental win. They are not
replaced by, or counted among, this candidate's new formal blocks.

Only `crates/chematic-smiles/src/rdkit/mod.rs` changes in production relative
to the prior first-use source candidate. Public CIP/stereocenter APIs retain
complete perception. No descriptor or ring calculation moves into parsing,
and no benchmark row receives a special-case result. Sanitization, hydrogen
handling, formulas and numeric accumulation order remain unchanged.

Four regression tests compare complete and optimized core views across
specified/unspecified/invalid stereo, non-tetrahedral tags, isotopes, hydrogen,
cleanup, metal coordination, radicals, errors and graph edits. The new cache
contract checks that generic callbacks do not receive the seed, and that
clones and hydrogen edits discard it.

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

The intermediate wheel used a dedicated initially empty Cargo target. This
wheel reuses that verified checkout and target, rebuilding the edited model
module and its affected dependents. Rust 1.97.0, release opt-level 3, LTO, one
codegen unit and `native-inchi` remain fixed. The exact
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
| tpsa | **1.513x [1.493, 1.518]; 20/20** | 0.521x [0.513, 0.527]; 0/20 | 0.339x [0.332, 0.346]; 0/20 | 0.091x [0.088, 0.095]; 0/20 |
| labute_asa | 0.719x [0.710, 0.730]; 1/20 | 0.335x [0.331, 0.346]; 0/20 | 0.162x [0.160, 0.165]; 0/20 | 0.067x [0.065, 0.068]; 0/20 |
| num_rings | 0.019x [0.019, 0.020]; 0/20 | Not ranked: outputs differ | 0.094x [0.091, 0.095]; 0/20 | 0.074x [0.068, 0.079]; 0/20 |
| morgan2_2048 | **3.137x [3.117, 3.172]; 20/20** | **2.029x [2.009, 2.048]; 20/20** | **1.405x [1.397, 1.417]; 20/20** | 0.895x [0.887, 0.902]; 0/20 |
| morgan2_chiral | **1.268x [1.252, 1.285]; 20/20** | 0.826x [0.820, 0.830]; 0/20 | **1.445x [1.425, 1.458]; 20/20** | 0.908x [0.901, 0.913]; 0/20 |
| morgan3_2048 | **3.374x [3.344, 3.398]; 20/20** | **2.240x [2.206, 2.263]; 20/20** | **1.626x [1.613, 1.639]; 20/20** | **1.066x [1.057, 1.072]; 20/20** |
| maccs | **1.888x [1.884, 1.891]; 20/20** | **2.169x [2.162, 2.180]; 20/20** | **1.632x [1.626, 1.636]; 20/20** | **1.747x [1.739, 1.750]; 20/20** |

### Repeated use

| Operation | RDKit, slice 0 | COSMolKit, slice 0 | RDKit, slice 5000 | COSMolKit, slice 5000 |
|---|---|---|---|---|
| tpsa | **9.952x [9.793, 10.025]; 20/20** | **1.719x [1.629, 1.834]; 20/20** | **9.852x [9.620, 10.065]; 20/20** | **1.346x [1.300, 1.432]; 20/20** |
| labute_asa | **2.751x [2.738, 2.799]; 20/20** | **1.441x [1.359, 1.490]; 20/20** | **2.635x [2.585, 2.670]; 20/20** | **1.264x [1.169, 1.354]; 20/20** |
| num_rings | **1.763x [1.735, 1.806]; 20/20** | Not ranked: outputs differ | **1.733x [1.683, 1.761]; 20/20** | 1.057x [1.032, 1.087]; 18/20 |
| morgan2_2048 | **5.220x [5.163, 5.246]; 20/20** | **3.207x [3.164, 3.271]; 20/20** | **4.796x [4.757, 4.851]; 20/20** | **2.909x [2.880, 2.929]; 20/20** |
| morgan2_chiral | **4.965x [4.935, 4.984]; 20/20** | **3.076x [3.051, 3.082]; 20/20** | **4.740x [4.699, 4.790]; 20/20** | **2.882x [2.858, 2.898]; 20/20** |
| morgan3_2048 | **4.882x [4.859, 4.904]; 20/20** | **3.115x [3.093, 3.132]; 20/20** | **4.583x [4.543, 4.614]; 20/20** | **2.903x [2.869, 2.919]; 20/20** |
| maccs | **1.947x [1.944, 1.955]; 20/20** | **2.243x [2.231, 2.249]; 20/20** | **1.731x [1.727, 1.740]; 20/20** | **1.857x [1.851, 1.863]; 20/20** |

### Parse inclusive

| Operation | RDKit, slice 0 | COSMolKit, slice 0 | RDKit, slice 5000 | COSMolKit, slice 5000 |
|---|---|---|---|---|
| tpsa | **13.086x [12.976, 13.187]; 20/20** | **12.241x [12.183, 12.316]; 20/20** | **4.282x [4.233, 4.291]; 20/20** | **4.440x [4.398, 4.475]; 20/20** |
| labute_asa | **12.316x [12.174, 12.387]; 20/20** | **11.638x [11.558, 11.735]; 20/20** | **4.094x [4.058, 4.129]; 20/20** | **4.389x [4.368, 4.417]; 20/20** |
| num_rings | **8.769x [8.731, 8.825]; 20/20** | Not ranked: outputs differ | **11.714x [11.545, 11.812]; 20/20** | **12.670x [12.555, 12.812]; 20/20** |
| morgan2_2048 | **8.591x [8.528, 8.642]; 20/20** | **7.589x [7.512, 7.622]; 20/20** | **4.265x [4.246, 4.295]; 20/20** | **4.113x [4.093, 4.135]; 20/20** |
| morgan2_chiral | **4.424x [4.390, 4.440]; 20/20** | **3.915x [3.880, 3.921]; 20/20** | **4.279x [4.242, 4.335]; 20/20** | **4.127x [4.114, 4.151]; 20/20** |
| morgan3_2048 | **7.949x [7.924, 7.987]; 20/20** | **6.926x [6.893, 6.968]; 20/20** | **4.199x [4.173, 4.230]; 20/20** | **4.028x [3.992, 4.044]; 20/20** |
| maccs | **2.105x [2.102, 2.112]; 20/20** | **2.373x [2.367, 2.384]; 20/20** | **1.795x [1.791, 1.797]; 20/20** | **1.927x [1.919, 1.931]; 20/20** |

### Incremental change from the previous source wheel

| Operation / mode | Previous / candidate, slice 0 | Previous / candidate, slice 5000 |
|---|---|---|
| tpsa/first_use | **1.032x [1.028, 1.044]; 20/20** | **1.068x [1.057, 1.077]; 20/20** |
| tpsa/repeated_use | 0.986x [0.973, 0.997]; 5/20 | 1.019x [0.989, 1.030]; 12/20 |
| tpsa/parse_inclusive | 1.008x [1.003, 1.016]; 16/20 | **1.045x [1.035, 1.051]; 20/20** |
| labute_asa/first_use | 1.022x [1.020, 1.032]; 19/20 | **1.065x [1.059, 1.070]; 20/20** |
| labute_asa/repeated_use | 0.997x [0.985, 1.005]; 7/20 | 1.005x [0.981, 1.024]; 11/20 |
| labute_asa/parse_inclusive | 1.002x [0.995, 1.011]; 12/20 | **1.047x [1.042, 1.055]; 20/20** |
| num_rings/first_use | 1.002x [0.993, 1.007]; 11/20 | 0.997x [0.991, 1.008]; 9/20 |
| num_rings/repeated_use | 0.995x [0.981, 1.005]; 7/20 | 0.988x [0.975, 1.005]; 8/20 |
| num_rings/parse_inclusive | 0.999x [0.995, 1.008]; 9/20 | 1.002x [0.996, 1.011]; 12/20 |
| morgan2_2048/first_use | 1.002x [0.997, 1.008]; 12/20 | 1.003x [0.996, 1.009]; 11/20 |
| morgan2_2048/repeated_use | 1.002x [0.996, 1.004]; 11/20 | 0.995x [0.992, 1.003]; 7/20 |
| morgan2_2048/parse_inclusive | 0.995x [0.987, 1.006]; 8/20 | 0.994x [0.989, 1.001]; 6/20 |
| morgan2_chiral/first_use | 0.997x [0.990, 1.003]; 8/20 | 0.994x [0.981, 0.999]; 5/20 |
| morgan2_chiral/repeated_use | 0.993x [0.990, 0.995]; 3/20 | 0.990x [0.985, 0.993]; 1/20 |
| morgan2_chiral/parse_inclusive | 0.997x [0.992, 1.001]; 7/20 | 0.998x [0.995, 1.006]; 9/20 |
| morgan3_2048/first_use | 0.996x [0.991, 1.002]; 7/20 | 0.997x [0.995, 1.005]; 9/20 |
| morgan3_2048/repeated_use | 0.998x [0.994, 1.004]; 8/20 | 1.001x [0.996, 1.009]; 10/20 |
| morgan3_2048/parse_inclusive | 0.994x [0.989, 0.999]; 5/20 | 0.998x [0.993, 1.002]; 8/20 |
| maccs/first_use | 1.001x [0.998, 1.004]; 11/20 | 0.997x [0.992, 1.002]; 8/20 |
| maccs/repeated_use | 0.999x [0.996, 1.003]; 9/20 | 0.999x [0.993, 1.001]; 7/20 |
| maccs/parse_inclusive | 0.998x [0.994, 1.002]; 8/20 | 0.994x [0.991, 1.001]; 7/20 |

### Whole-process memory

| Slice | Candidate median peak MiB | RDKit MiB | COSMolKit MiB |
|---|---:|---:|---:|
| 0 | 46.57 | 106.58 | 105.71 |
| 5000 | 40.91 | 88.33 | 83.05 |


The original-slice first-use ASA median improves 1.022x [1.020, 1.032], but
19/20 faster blocks do not satisfy the strict gate. Original-slice repeated
TPSA is about 1.4% slower (ratio 0.986 [0.973, 0.997]). Repeated chiral R2 is
about 0.7% and 1.0% slower on the two slices; original-slice read-plus-R3 and
additional-slice initial chiral R2 are about 0.6% slower. All are retained,
including lanes whose production implementation did not change.

COSMolKit remains faster for initial TPSA/ASA on both slices, initial chiral R2
on both slices and plain R2 on the additional slice. RDKit initial ring count
is substantially faster on both slices; comparable COSMolKit ring count is
faster on the additional slice. Its original-slice ring outputs still differ
on six records and are not ranked. Additional-slice repeated ring count is
1.057x [1.032, 1.087] versus COSMolKit but wins only 18/20 blocks, leaving that
strict gate open. These boundaries remain visible beside the workflow wins.

## Correctness and validation

Normal 10,000 and unusual 4,072 inputs retain 112,576 first/repeated descriptor
outcomes and 56,288 Morgan outcomes byte-for-byte against the prior candidate.
Float bits, errors and nonfinite outcomes remain in the dumps. Identical
nonfinite outcomes are retention evidence, not successful numeric agreement.
The previously recorded RDKit residuals remain; matching the old source does
not erase them. Existing baseline dumps are reused only after their hashes
were verified against the earlier packet.

Final isolated release validation runs all SMILES library tests plus chemistry
scalar-cache and RDKit-reference integration tests: 508 passed, zero failed and
one existing ignored case. All four model-view tests pass. The preceding
unflagged-only candidate also passed the broader 1,959-test SMILES/chemistry/
fingerprint run; that earlier scope is recorded in its separate packet.
Release Clippy for the final changed SMILES crate and all targets passes with
warnings denied. Final selected Python molecule, descriptor, RDKit, fingerprint,
pickle and cross-binding suites pass 590 tests. Five unchanged benchmark-gate
tests pass; changed Rust files pass rustfmt and the diff has no whitespace
errors.

The final candidate reconstructs from the base, three prerequisites and its
alternative patch with all 894 frozen Rust/workspace-manifest files matching.
Reference fixtures are pinned to the base commit. Paired summaries and retained
archive hashes are independently verified in the parity record.

## Retained evidence

- [Source, wheel, compiler and comparator manifest](2026-10-14-python-sanitized-view-reuse-source.json)
- [Original-slice raw blocks](2026-10-14-python-sanitized-view-reuse-slice0-paired20.json)
- [Additional-slice raw blocks](2026-10-14-python-sanitized-view-reuse-slice5000-paired20.json)
- [Exact retention, reconstruction and validation](2026-10-14-python-sanitized-view-reuse-parity.json)
- [Reconstructable candidate patch](https://github.com/kent-tokyo/chematic/blob/0bc130f24a7d163535e0f4ad59eb5a155847cf07/benchmarks/2026-10-14-python-sanitized-view-reuse-candidate.patch)
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
