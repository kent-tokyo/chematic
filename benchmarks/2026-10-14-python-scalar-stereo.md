# Python scalar and stereo reuse follow-up — 2026-10-14

The source candidate improves prepared TPSA and Labute ASA calls by about
10–18x over the previous source candidate on two fixed 1,000-row slices.
Prepared chiral Morgan improves 3.335x on slice 0; on the low-stereo slice 5000
its baseline interval crosses parity. Prepared chiral Morgan beats both comparators in
all 20 blocks on both slices. Long scalar batches have favorable paired 95%
intervals against both comparators, but still lose individual COSMolKit blocks.
The strict all-block scalar gate and the complete speed target remain open.

## Source and implementation

The baseline is the prior [Python speed-frontier candidate](2026-10-14-python-speed-frontier.md),
not the clean v1.1.3 baseline: base `ffe682c32cb4c76de0d309649427202a7d8e956d`
plus that packet's six-file production patch. Its exact wheel SHA256 is
`7e77ff905a393315fa4036063e349e44e9691a3596486404fad62ae356890709`.
The new wheel SHA256 is
`b3a7b1a050f6517420a4f3e627a0f19ab407b310914c2848ad2461af87cd1408`.
Both say 1.1.3; these hashes and patches distinguish source artifacts. Neither
result establishes the behavior of a new public release.

- Three typed inline f64 cache cells reuse default TPSA, RDKit-default TPSA
  and whole-molecule Labute ASA. No f64 bit pattern is reserved as a sentinel.
- The RDKit-model scalar helper memoizes the outer calculation on its input
  molecule. A nested model calculation bypasses this cache, preserving the
  existing rule that nested callbacks calculate on the supplied graph.
- Morgan caches prepared tetrahedral/CIP and double-bond stereo. Radius and
  fold width reuse it; the fingerprint itself is still recomputed.
- Typical TPSA atom environments keep up to four neighbors on the stack;
  higher degrees retain a spilling fallback.

All caches belong to the existing molecule derived cache. Graph and stereo
edits clear it and cloning starts empty. Concurrent cold readers may compute
the same deterministic value, then publish it safely. No calculation moves
into SMILES parsing to improve a prepared-call benchmark. The tradeoff is
additional per-molecule cache cells and retained stereo vectors until edit,
clone or drop; whole-process RSS is not a per-molecule memory claim.

The [source manifest](2026-10-14-python-scalar-stereo-source.json) and
[follow-up patch](https://github.com/kent-tokyo/chematic/blob/0bc130f24a7d163535e0f4ad59eb5a155847cf07/benchmarks/2026-10-14-python-scalar-stereo-candidate.patch) reconstruct
these changes after the prior patch. Every reconstructed file matches its
recorded hash. The isolated build excludes unrelated dirty workspace changes,
including the Open Babel, 3D and compatibility work. Release wheels use Rust
1.97.0, opt-level 3, LTO, one codegen unit and native-inchi, matching the previous
wheel's build lane. A concurrent distance-matrix limit change entered the shared
SMILES module after this snapshot; it is preserved outside the measured patch.
The new scalar helper is byte-identical, while the whole shared module is not
the measured artifact. The comparator artifacts remain RDKit 2026.09.1 and
COSMolKit 0.5.0rc22 on macOS arm64 / CPython 3.13.6 / 10 logical CPUs.

## Original three-mode comparison

Each slice contains exactly 1,000 nonempty rows from
`validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`, at offsets 0 and
5000, SHA256 `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`.
Offsets were fixed before either timing run. Both slices belong to an exposed
existing corpus; this is not a sealed holdout or an additional host.

Twenty sequential four-arm blocks rotate every arm through every position
five times. Native return objects are retained during timing; normalization,
output agreement and hashing occur afterward. Fresh molecules are used for
each operation and preparation mode. First use excludes parsing but includes
lazy work not performed by chematic's parser. Repeated use gives each molecule
one untimed operation. Parse inclusive performs actual parse and operation in
the same timed loop. Output failures or disagreements suppress speed ranking.
The COSMolKit ring count differs on six slice-0 rows; none is filtered out.

All table entries are paired comparator / candidate median ratios,
bootstrap 95% intervals, and winning blocks. Above 1 favors the candidate.
A strict win needs equal outputs, lower confidence bound above 1, and 20/20
winning blocks. A favorable median or interval alone is not that strict win.

### First use

| Operation | RDKit, slice 0 | COSMolKit, slice 0 | RDKit, slice 5000 | COSMolKit, slice 5000 |
|---|---|---|---|---|
| tpsa | 1.431x [1.410, 1.493]; 15/20 | 0.489x [0.390, 0.500]; 0/20 | 0.327x [0.294, 0.375]; 0/20 | 0.088x [0.073, 0.108]; 0/20 |
| labute_asa | 0.704x [0.683, 0.717]; 1/20 | 0.323x [0.274, 0.329]; 0/20 | 0.149x [0.140, 0.172]; 0/20 | 0.059x [0.054, 0.067]; 0/20 |
| num_rings | 0.019x [0.018, 0.021]; 0/20 | Not ranked: outputs differ | 0.101x [0.086, 0.109]; 0/20 | 0.095x [0.084, 0.130]; 0/20 |
| morgan2_2048 | 3.104x [3.039, 3.146]; 20/20 | 1.991x [1.954, 2.027]; 19/20 | 1.477x [1.332, 1.549]; 20/20 | 0.932x [0.863, 1.096]; 8/20 |
| morgan2_chiral | 1.276x [1.234, 1.304]; 16/20 | 0.822x [0.808, 0.836]; 2/20 | 1.415x [1.326, 1.568]; 20/20 | 0.923x [0.843, 1.082]; 9/20 |
| morgan3_2048 | 3.341x [3.303, 3.385]; 20/20 | 2.231x [2.211, 2.287]; 20/20 | 1.614x [1.546, 1.845]; 20/20 | 1.061x [1.017, 1.223]; 15/20 |
| maccs | 1.889x [1.882, 1.895]; 20/20 | 2.174x [2.165, 2.202]; 20/20 | 1.622x [1.412, 1.702]; 20/20 | 1.721x [1.531, 1.957]; 20/20 |

### Repeated use

| Operation | RDKit, slice 0 | COSMolKit, slice 0 | RDKit, slice 5000 | COSMolKit, slice 5000 |
|---|---|---|---|---|
| tpsa | 10.076x [9.552, 10.351]; 20/20 | 1.703x [1.596, 1.828]; 19/20 | 8.477x [7.829, 9.471]; 20/20 | 1.250x [1.117, 1.795]; 17/20 |
| labute_asa | 2.735x [2.668, 2.798]; 19/20 | 1.381x [1.221, 1.412]; 16/20 | 2.269x [1.765, 2.606]; 20/20 | 1.041x [0.839, 1.385]; 12/20 |
| num_rings | 1.725x [1.558, 1.817]; 15/20 | Not ranked: outputs differ | 1.761x [1.450, 1.964]; 19/20 | 1.003x [0.855, 1.309]; 10/20 |
| morgan2_2048 | 5.248x [5.137, 5.300]; 20/20 | 3.169x [3.111, 3.242]; 20/20 | 4.713x [4.527, 5.370]; 20/20 | 3.074x [2.758, 3.668]; 20/20 |
| morgan2_chiral | 4.960x [4.853, 5.209]; 20/20 | 3.093x [3.042, 3.196]; 20/20 | 4.679x [4.426, 5.014]; 20/20 | 3.210x [2.785, 3.663]; 20/20 |
| morgan3_2048 | 4.840x [4.703, 4.899]; 20/20 | 3.121x [3.063, 3.183]; 20/20 | 4.437x [4.204, 4.812]; 20/20 | 3.045x [2.767, 3.620]; 20/20 |
| maccs | 1.941x [1.589, 1.952]; 20/20 | 2.230x [2.012, 2.243]; 20/20 | 1.702x [1.481, 1.790]; 20/20 | 1.883x [1.547, 2.244]; 20/20 |

### Parse inclusive

| Operation | RDKit, slice 0 | COSMolKit, slice 0 | RDKit, slice 5000 | COSMolKit, slice 5000 |
|---|---|---|---|---|
| tpsa | 13.086x [12.909, 13.286]; 20/20 | 11.964x [11.337, 12.164]; 20/20 | 3.918x [3.638, 4.501]; 20/20 | 4.421x [3.682, 5.157]; 20/20 |
| labute_asa | 12.385x [12.268, 12.828]; 20/20 | 11.470x [10.522, 11.665]; 20/20 | 3.885x [3.516, 4.316]; 20/20 | 4.069x [3.756, 5.057]; 20/20 |
| num_rings | 8.567x [8.449, 8.690]; 20/20 | Not ranked: outputs differ | 11.020x [10.179, 12.679]; 20/20 | 11.833x [11.033, 16.324]; 20/20 |
| morgan2_2048 | 8.689x [8.587, 8.799]; 20/20 | 7.537x [7.416, 7.606]; 20/20 | 4.212x [3.904, 4.628]; 20/20 | 4.253x [3.699, 5.015]; 20/20 |
| morgan2_chiral | 4.440x [4.414, 4.549]; 20/20 | 3.935x [3.880, 4.012]; 20/20 | 4.268x [3.980, 4.653]; 20/20 | 4.003x [3.636, 4.727]; 20/20 |
| morgan3_2048 | 7.952x [7.807, 7.986]; 20/20 | 6.949x [6.824, 6.990]; 20/20 | 4.129x [3.940, 4.540]; 20/20 | 4.097x [3.567, 4.763]; 20/20 |
| maccs | 2.109x [1.641, 2.124]; 20/20 | 2.348x [2.156, 2.387]; 20/20 | 1.820x [1.717, 1.933]; 20/20 | 1.971x [1.692, 2.071]; 20/20 |

Raw records: [slice 0](2026-10-14-python-scalar-stereo-slice0-paired20.json),
[slice 5000](2026-10-14-python-scalar-stereo-slice5000-paired20.json).
These records are retained unchanged by the additional protocol below.
The first-use ASA and plain R2 slice-0 baseline intervals show about 1–2%
extra cost. Other first-use baseline intervals are mixed or cross parity.
This is a reuse optimization, not evidence that every first-use lane improved.
The independent second slice also exposes weaker first-use plain-Morgan and
scalar ratios against COSMolKit; the earlier prefix-only results cannot be
generalized to every row or preparation boundary.

## Additional long scalar-call batches

The original scalar timers were around 35 microseconds per 1,000 prepared
calls and included substantial block outliers. After inspecting those results,
a separate, fixed protocol was declared: 1,000 calls on each of the same 1,000
prepared molecules, one million actual property calls per operation and arm.
Both fixed slices and 20 rotating blocks are retained. Every returned native
float remains live during timing, and every one is checked outside timing
against the preparation result. Both engines' output agreement remains a gate.
This accounts for 320 million timed scalar calls across both operations,
four arms, 20 blocks and two slices. It measures repeated-call throughput;
it does not replace first-use, parse-inclusive or short-batch results.

The [additional source/protocol manifest](2026-10-14-python-scalar-reuse-source.json)
records that this protocol was chosen after the short-batch result. Production
code and wheel artifacts were unchanged. The additional
[runner](../scripts/bench_python_scalar_reuse.py) and its tests are also retained
in a [tools patch](https://github.com/kent-tokyo/chematic/blob/0bc130f24a7d163535e0f4ad59eb5a155847cf07/benchmarks/2026-10-14-python-scalar-reuse-tools.patch), with its hash in
the parity record. No block is discarded or rerun to obtain an all-block win.

| Operation | RDKit, slice 0 | COSMolKit, slice 0 | RDKit, slice 5000 | COSMolKit, slice 5000 |
|---|---|---|---|---|
| tpsa | 10.054x [9.903, 10.201]; 20/20 | 1.237x [1.209, 1.268]; 19/20 | 9.851x [9.783, 10.103]; 20/20 | 1.239x [1.220, 1.258]; 19/20 |
| labute_asa | 2.598x [2.529, 2.623]; 20/20 | 1.090x [1.076, 1.107]; 18/20 | 2.559x [2.461, 2.606]; 20/20 | 1.092x [1.069, 1.110]; 19/20 |

Raw long-batch records: [slice 0](2026-10-14-python-scalar-reuse-slice0-paired20.json),
[slice 5000](2026-10-14-python-scalar-reuse-slice5000-paired20.json).
Both scalar lanes beat RDKit in every long-batch block. COSMolKit intervals
favor the candidate, but TPSA wins 19/20 blocks on each slice and ASA wins
18/20 and 19/20. They still fail this packet's stricter all-block gate.

The source-candidate improvements over the previous candidate are:

| Prepared operation | Previous candidate / new, slice 0 | Previous candidate / new, slice 5000 |
|---|---|---|
| tpsa (long batches) | 9.977x [9.692, 10.258]; 20/20 | 16.276x [15.921, 16.483]; 20/20 |
| labute_asa (long batches) | 13.388x [13.113, 13.685]; 20/20 | 17.779x [17.035, 18.031]; 20/20 |
| morgan2_chiral (original protocol) | 3.335x [3.262, 3.364]; 20/20 | 0.994x [0.935, 1.070]; 8/20 |

## Correctness and validation

The [parity record](2026-10-14-python-scalar-stereo-parity.json) preserves exact
float bits and every error/nonfinite outcome. Four descriptor calls, first and
repeated, on normal 10,000 plus unusual 4,072 inputs give 112,576 results that
are bit-identical to the previous candidate. First and repeated values also
match exactly. Morgan R2/R3 with chirality off/on gives 56,288 cells identical
to both the previous candidate and the pinned RDKit reference. The new Morgan
dump is byte-identical to the already retained prior dump, so it is not duplicated.

Independent RDKit comparison must be kept separate from baseline retention:

- RDKit-default TPSA: 14,072/14,072 matches within absolute 1e-9.
- Labute ASA: 14,069 finite matches, plus three nonfinite reference/source
  outcomes; these three are not counted as successful agreement or speed wins.
- Ring count: normal 10,000/10,000; unusual 4,068/4,072. Four pre-existing
  differences remain, and this change does not modify ring perception.
- Default S/P-inclusive TPSA: normal 10,000/10,000; unusual 4,041/4,072.
  Thirty-one pre-existing model differences remain. The speed tables use the
  separate RDKit-default TPSA API, not this default extended model.

Final isolated Rust core/chemistry/fingerprint/SMILES tests: 2,117 passed,
zero failed, eight existing ignored cases. Seven new contract tests cover
variant/radius/width separation, atom/bond/topology/stereo edits, empty clone
caches, explicit-H and cleanup models, nested model contexts and concurrent
readers. Release Clippy passes for all four crates and all targets, warnings
denied. Python descriptor, RDKit, stable fingerprint, pickle and cross-binding
suites: 526 passed. Five frontier-gate tests and two long-batch validation tests
pass; the latter reject missing calls, nonfinite/wrong types and a differing
value in any retained repetition. Changed Rust formatting and diff checks pass.
Logs and the shared baseline/candidate descriptor dump are retained as gzip
artifacts with hashes in the parity record.

All timing records include whole-process peak RSS. For the long protocol it
includes imports, parsing, preparation and one million retained floats; for the
original protocol it includes all seven operations. These are distinct memory
boundaries and must not be interpreted as streaming or operation-only RSS.

## Reproduce

Apply the prior production patch, then this follow-up patch to the pinned base.
Preserve the exact prior candidate wheel as the new baseline. Build and install
the new candidate in a separate directory with the recorded settings; do not
reuse a stale target across archived source snapshots. Source/patch and loaded
extension hashes are checked against the recorded artifacts.

```sh
python scripts/bench_python_speed_frontier.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --offset 0 --limit 1000 --blocks 20 \
  --baseline-wheel /path/to/previous-candidate.whl --baseline-path /path/to/baseline \
  --candidate-wheel /path/to/new-candidate.whl --candidate-path /path/to/candidate \
  --base-revision ffe682c32cb4c76de0d309649427202a7d8e956d \
  --source-manifest benchmarks/2026-10-14-python-scalar-stereo-source.json \
  --output /path/to/slice0.json
```

Repeat with offset 5000. For long scalar batches, apply the tools patch and use
`bench_python_scalar_reuse.py` with the same corpus, offsets, rows, wheels and
installed paths, `--blocks 20 --repeats 1000`, and the separate scalar-reuse
source manifest. Run `verify_python_descriptor_frontier.py` for baseline,
candidate and RDKit in their respective environments, and
`verify_python_morgan_frontier.py` for the candidate; compare the retained exact
bit dumps and declared reference tolerances.

## What remains

Toolkit-wide superiority is unproven. First-use TPSA/ASA, initial rings and
first-use chiral Morgan still have measured deficits. COSMolKit's scalar
all-block gates remain open even after longer timing. First-use plain Morgan
on the second slice also lacks a COSMolKit win. The comparison needs broader
common functions, reactions/3D, memory-specific work, sealed corpora, other
hosts and exact published artifacts before a broader claim. No release,
commit, hosted CI or merge is implied by these source-candidate results.
