# Compact Morgan linear-scan candidate: rejected

The candidate replaces a binary search for every compact Morgan environment
with one monotone scan over the sorted previous-round environments. It preserves
outputs, but **fails the predeclared source performance gate and must not be
merged on this evidence**. Both primary intervals cross one and six of 42 source
paired medians fall below the 0.99 retention floor. No timing epoch was rerun or pooled.

This is a completed, one-host source experiment. It does not establish a new
main, release, public-wheel, whole-toolkit, or cross-host speed improvement.

## Fixed identities and boundaries

- Baseline: `c84c68f297f95562d1d29c75234f2190cd5f9d75`.
- Measured candidate: `f331e3a9353b22027671570c66c2558f736aa821`; later commits add evidence only.
- Python 3.13.6 on macOS arm64; RDKit 2026.09.1, COSMolKit 0.5.0, NumPy 2.5.3.
- Rust release settings: opt-level3, LTO=true, codegen-units1, no stripping,
  native InChI, one build job, no incremental compilation, locked offline dependencies.
  Baseline target began empty. Candidate reused dependency cache and actually
  recompiled fingerprint/Python production libraries. No benchmark-specific production
  feature, dependency, cache, preparation relocation, API change, or unsafe code was added.
- Corpus: unfiltered 1,000-row windows at offsets0 and5000 of
  `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`.
- Timed boundary: each public API's native return (chematic packed bytes versus
  RDKit/COSMolKit fingerprint objects). Canonical on-bit conversion and checking are
  outside the timers. These numbers describe those APIs, including their different
  return-object costs; they do not isolate the native fingerprint algorithm.
- Global warmup uses disjoint objects. First-use reparses outside **every** timed
  batch and calls the operation once per fresh molecule. Repeated-use prepares a
  fresh batch, calls that operation once untimed, then reuses those objects.
  Parse-inclusive performs fresh parsing and the operation inside each timer.
- Every full batch returns all1,000 outputs. Every batch is checked against the
  same engine's canonical outputs, and all modes are checked against each other.
  Primary operations accumulate at least50ms of timed full-batch work per mode;
  source controls at least20ms. Durations, counts, and output hashes are retained.
  The average is sum(raw batch ns)/actual complete batches, floored by less than1ns.
  Preparation and checking gaps are excluded; this is not a contiguous steady-state
  or cold-interpreter benchmark. Whole-process RSS includes imports/preparation/checks;
  it is recorded but not ranked.
- Each profile runs all24 permutations of the four arms once, shuffled with
  seeds20261017+offset. Each arm occupies each position six times. Samples stay
  paired by block. Statistics use10,000 fixed-seed bootstrap resamples; rounded
  display values never decide the median floor. There are no filtered rows,
  missing observations, interrupted profiles, or performance-informed reruns.

The [plan](2026-10-17-morgan-linear-artifacts/measurement-plan.json) was saved
before production edits and candidate measurements. The [source manifest](2026-10-17-morgan-linear-artifacts/source-manifest.json)
seals all production sources, loaded adapters, driver, corpus-linked plans,
wheel hashes, module paths/native hashes, and1,883 non-bytecode runtime files.
These were checked before and after both profiles and again during packaging.
Package versions identify the installed comparator artifacts; this was not a
registry update check or a claim about another version.

- baseline wheel SHA-256: `7a3fc614b887fdd3fc93801a8bce1e3df31100cddeacaa8444bc268190d0c91b`.
- candidate wheel SHA-256: `99aaf790b8e124dd714fd0e634557feb245bd5815b06850f1d5a029ba99d99d8`.
- offset0 corpus SHA-256: `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`; exact row-window SHA-256: `844703091ec45b9c0833f705c9b2e532425c79edd71d874cf74b8d2f1e018c82`.
- offset5000 corpus SHA-256: `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`; exact row-window SHA-256: `a84109c44220cdbc44e011b2745e2e3d6253e0d6a448191c977678e543f8186b`.

## Source decision

Adoption required both predeclared primary source CI lower bounds>1, every one
of42 unrounded source paired medians>=0.99, complete data, and exact source outputs.
Neither primary CI is positive. The six floor misses below are binding even where
a CI also includes one. Unchanged descriptor control bodies do not establish that
the algorithm itself caused those fluctuations; the preregistered retention test
still fails, so no exception or favorable resampling was used.

| Offset | Operation | Mode | Baseline/candidate paired median | 95% CI | Faster blocks |
|---:|---|---|---:|---|---:|
| 0 | morgan2_chiral | repeated_use | 1.007344 | [0.996787, 1.024552] | 15/24 |
| 5000 | morgan2_chiral | repeated_use | 1.004094 | [0.987829, 1.017851] | 14/24 |
| 0 | labute_asa | first_use | 0.986769 | [0.974611, 1.007378] | 9/24 |
| 0 | num_rings | first_use | 0.986465 | [0.963830, 1.003296] | 8/24 |
| 0 | num_rings | parse_inclusive | 0.989029 | [0.973683, 1.003387] | 10/24 |
| 5000 | morgan2_2048 | parse_inclusive | 0.989949 | [0.968006, 1.005395] | 8/24 |
| 5000 | tpsa | first_use | 0.989667 | [0.981885, 1.000728] | 8/24 |
| 5000 | labute_asa | parse_inclusive | 0.985738 | [0.976510, 1.015684] | 9/24 |

Additionally offset0 `num_rings/repeated_use` has an entirely negative interval:
0.991299 [0.983534,0.998288],7/24 faster. Its median clears the floor but its
negative interval is retained. The0.989949 Morgan parse-inclusive median above
fails despite rounding to0.990 at three decimals.

All42 source lanes, including every negative or inconclusive interval:

| Offset | Operation | Mode | Paired median | 95% CI | Faster blocks |
|---:|---|---|---:|---|---:|
| 0 | morgan2_2048 | first_use | 1.005750 | [0.995007, 1.011663] | 14/24 |
| 0 | morgan2_2048 | repeated_use | 1.003322 | [0.995936, 1.026661] | 14/24 |
| 0 | morgan2_2048 | parse_inclusive | 1.005044 | [0.995472, 1.009780] | 14/24 |
| 0 | morgan2_chiral | first_use | 1.001527 | [0.994846, 1.012580] | 12/24 |
| 0 | morgan2_chiral | repeated_use | 1.007344 | [0.996787, 1.024552] | 15/24 |
| 0 | morgan2_chiral | parse_inclusive | 0.998128 | [0.995098, 1.003906] | 10/24 |
| 0 | morgan3_2048 | first_use | 1.000614 | [0.993524, 1.005782] | 13/24 |
| 0 | morgan3_2048 | repeated_use | 1.007848 | [0.985381, 1.013360] | 16/24 |
| 0 | morgan3_2048 | parse_inclusive | 0.997671 | [0.993734, 1.007562] | 11/24 |
| 0 | tpsa | first_use | 1.003967 | [0.987630, 1.016790] | 14/24 |
| 0 | tpsa | repeated_use | 1.007884 | [0.999351, 1.011712] | 16/24 |
| 0 | tpsa | parse_inclusive | 0.991451 | [0.981486, 1.005511] | 10/24 |
| 0 | labute_asa | first_use | 0.986769 | [0.974611, 1.007378] | 9/24 |
| 0 | labute_asa | repeated_use | 1.000268 | [0.990808, 1.006259] | 13/24 |
| 0 | labute_asa | parse_inclusive | 0.991353 | [0.982067, 1.007641] | 9/24 |
| 0 | num_rings | first_use | 0.986465 | [0.963830, 1.003296] | 8/24 |
| 0 | num_rings | repeated_use | 0.991299 | [0.983534, 0.998288] | 7/24 |
| 0 | num_rings | parse_inclusive | 0.989029 | [0.973683, 1.003387] | 10/24 |
| 0 | maccs | first_use | 0.995937 | [0.994019, 1.001180] | 9/24 |
| 0 | maccs | repeated_use | 1.000284 | [0.991464, 1.004654] | 12/24 |
| 0 | maccs | parse_inclusive | 0.998391 | [0.992826, 1.005464] | 10/24 |
| 5000 | morgan2_2048 | first_use | 0.998963 | [0.981546, 1.016721] | 12/24 |
| 5000 | morgan2_2048 | repeated_use | 1.006601 | [0.993939, 1.025963] | 14/24 |
| 5000 | morgan2_2048 | parse_inclusive | 0.989949 | [0.968006, 1.005395] | 8/24 |
| 5000 | morgan2_chiral | first_use | 0.991165 | [0.981147, 1.009782] | 9/24 |
| 5000 | morgan2_chiral | repeated_use | 1.004094 | [0.987829, 1.017851] | 14/24 |
| 5000 | morgan2_chiral | parse_inclusive | 1.004194 | [0.983186, 1.018420] | 14/24 |
| 5000 | morgan3_2048 | first_use | 0.997228 | [0.981108, 1.002740] | 9/24 |
| 5000 | morgan3_2048 | repeated_use | 0.990703 | [0.970633, 1.007600] | 10/24 |
| 5000 | morgan3_2048 | parse_inclusive | 0.996952 | [0.980645, 1.003835] | 10/24 |
| 5000 | tpsa | first_use | 0.989667 | [0.981885, 1.000728] | 8/24 |
| 5000 | tpsa | repeated_use | 1.005459 | [0.990212, 1.022493] | 14/24 |
| 5000 | tpsa | parse_inclusive | 1.008730 | [0.998926, 1.020401] | 16/24 |
| 5000 | labute_asa | first_use | 1.000980 | [0.989186, 1.008564] | 12/24 |
| 5000 | labute_asa | repeated_use | 1.012385 | [0.995365, 1.025790] | 15/24 |
| 5000 | labute_asa | parse_inclusive | 0.985738 | [0.976510, 1.015684] | 9/24 |
| 5000 | num_rings | first_use | 0.993941 | [0.970735, 1.014036] | 10/24 |
| 5000 | num_rings | repeated_use | 1.008171 | [0.998693, 1.016094] | 16/24 |
| 5000 | num_rings | parse_inclusive | 0.992116 | [0.977315, 1.016799] | 10/24 |
| 5000 | maccs | first_use | 1.003042 | [0.992050, 1.008150] | 13/24 |
| 5000 | maccs | repeated_use | 1.004731 | [1.001017, 1.026560] | 17/24 |
| 5000 | maccs | parse_inclusive | 0.997379 | [0.991496, 1.007841] | 11/24 |

## Fresh external observations on the rejected candidate

For each external lane the gate requires equivalent outputs, CI lower>1, and24/24
faster paired blocks. All ordinary timed outputs match.32 of36 lanes pass that
stricter speed gate. These measurements do not rescue the failed source decision.
The offset0 first-use chiral lane is slower than COSMolKit:0.820467
[0.813832,0.831359],0/24 faster (native API batches roughly12.321ms versus10.106ms).
The other three non-wins have positive medians but21/24,21/24,or23/24 faster
blocks; they are not labelled consistent wins.

Each cell is comparator/candidate paired median [95% CI]; faster blocks.

| Offset | Operation | Mode | RDKit | COSMolKit |
|---:|---|---|---|---|
| 0 | morgan2_2048 | first_use | 3.349268 [3.301601, 3.384595];24/24 | 2.107447 [2.085815, 2.119830];24/24 |
| 0 | morgan2_2048 | repeated_use | 5.459820 [5.420918, 5.516055];24/24 | 3.285710 [3.269921, 3.324281];24/24 |
| 0 | morgan2_2048 | parse_inclusive | 8.718187 [8.654798, 8.768805];24/24 | 7.653291 [7.609552, 7.747052];24/24 |
| 0 | morgan2_chiral | first_use | 1.297235 [1.289657, 1.305535];24/24 | 0.820467 [0.813832, 0.831359];0/24 |
| 0 | morgan2_chiral | repeated_use | 5.396659 [5.321039, 5.450128];24/24 | 3.265829 [3.248260, 3.301860];24/24 |
| 0 | morgan2_chiral | parse_inclusive | 4.485394 [4.460156, 4.508210];24/24 | 3.938083 [3.901311, 3.966230];24/24 |
| 0 | morgan3_2048 | first_use | 3.466386 [3.448772, 3.480091];24/24 | 2.238358 [2.231400, 2.255877];24/24 |
| 0 | morgan3_2048 | repeated_use | 5.054921 [4.949983, 5.096114];24/24 | 3.163830 [3.147236, 3.189922];24/24 |
| 0 | morgan3_2048 | parse_inclusive | 8.092720 [8.014045, 8.149680];24/24 | 7.049885 [6.956573, 7.123678];24/24 |
| 5000 | morgan2_2048 | first_use | 1.601333 [1.586472, 1.637141];24/24 | 1.011954 [1.007096, 1.033980];21/24 |
| 5000 | morgan2_2048 | repeated_use | 5.033747 [4.891531, 5.085057];24/24 | 3.009291 [2.880286, 3.091650];24/24 |
| 5000 | morgan2_2048 | parse_inclusive | 4.563935 [4.396937, 4.612959];24/24 | 4.368171 [4.272584, 4.449053];24/24 |
| 5000 | morgan2_chiral | first_use | 1.605616 [1.587573, 1.645895];24/24 | 1.027897 [1.014800, 1.040781];21/24 |
| 5000 | morgan2_chiral | repeated_use | 4.988879 [4.823370, 5.128287];24/24 | 2.983106 [2.926063, 3.059350];24/24 |
| 5000 | morgan2_chiral | parse_inclusive | 4.577449 [4.477129, 4.632291];24/24 | 4.384786 [4.350446, 4.477999];24/24 |
| 5000 | morgan3_2048 | first_use | 1.805440 [1.773463, 1.838818];24/24 | 1.181307 [1.167015, 1.201038];23/24 |
| 5000 | morgan3_2048 | repeated_use | 4.619374 [4.442131, 4.771084];24/24 | 2.883112 [2.807106, 2.989921];24/24 |
| 5000 | morgan3_2048 | parse_inclusive | 4.528959 [4.455486, 4.593018];24/24 | 4.204612 [4.180258, 4.308648];24/24 |

Only the three Morgan operations have fresh external timings. TPSA, Labute ASA,
ring count and MACCS are **source-only controls** here. No other workflow is ranked.

## Correctness and native-code validation

- Existing fingerprint suite:387 tests passed, zero failures/ignored. Includes
 24,576 exhaustive comparisons (all1,024 five-vertex simple undirected carbon
 graphs ×3 ring-flag choices ×2 chirality choices ×4 radii), with ordered
 emissions and fallback results compared against the complete frozen old function.
 Directed tests cover empty/isolated/symmetric/ring/stereo cases,64/65 atom inline
 spill and128/129 bond compact/fallback boundaries. Existing general-path comparison,
 mutation and radius/width cache tests remain in the full suite.
- The frozen complete reference matches baseline exactly after only the test
 function name, visibility and module indentation change. Its hash/check is saved.
- On10,000 ordinary plus4,072 stress inputs, baseline and candidate rows (including
 error/refusal fields) are identical. All56,288 Morgan values (R2/R3 ×chirality
 off/on) also match RDKit2026.09.1 exactly; there are no source refusals.
- COSMolKit matches all40,000 ordinary cells. On stress inputs, each plain radius
 matches4,066/4,072 with6 parse refusals. Each chiral radius matches4,059/4,072 with
 the same6 refusals and7 different fingerprints. All exceptional rows and their
 values/errors remain in the differences file. They are not counted as agreement
 or silently omitted from a speed claim. The timing windows contain ordinary rows.
- All192 compressed per-arm journals are present. Independent audit checks2,880
 accounting cells,177,202 complete timed batches and177,202,000 processed records,
 all canonical outputs/digests, exact permutations/positions, and the raw statistics.
 Every source and external timed output matches over all modes and batches.
- Native production libraries actually recompiled. Distinct extension hashes and
 the exact Morgan closure symbol's disassembly (3,095 versus3,071 lines) show the
 candidate exists in the measured code, rather than only in a reference path.
 The [code-realization record](2026-10-17-morgan-linear-artifacts/code-realization.json)
 identifies both symbols, commands, code and extension hashes.
- Local Clippy all-targets, rustfmt check and refactor budget passed.12 benchmark
 driver/frontier tests pass. New tests verify first-use freshness, repeated-use
 priming/reuse, parse-inclusive timer placement, fractional batch accounting,
 refusal/mismatch/nonpositive-time rejection, unfiltered corpus windows, and
 preservation of the prior journal when atomic replacement fails.

## Diagnostic follow-up, excluded from formal timing

After both formal profiles finished, a separate5-second macOS sample profiled
baseline parsing plus first-use chiral Morgan on the first1,000 ordinary rows.
The sampled fingerprint method has2,693 inclusive stacks;1,596 are inside stereo
preparation (~59.3%). Within that subtree, observed work includes sanitization,
RDKit ring perception and CIP ranking. These are inclusive sampling counts, not
independent additive times, a causal speedup estimate, or another formal trial.
The sample supports investigating that preparation next; it does not justify
bypassing sanitization, changing stereo semantics, or adopting this candidate.
The diagnostic worker, sample and output log are retained.

## Evidence and reproduction

- [Independent audit and all78 comparisons](2026-10-17-morgan-linear-artifacts/formal-audit.json),
  [unrounded CSV](2026-10-17-morgan-linear-artifacts/formal-lanes.csv),
  [broad audit](2026-10-17-morgan-linear-artifacts/broad-audit.json).
- Both `profile-*/` directories contain the original plan, result/progress records
  and96 atomic journals apiece. `probe-*.json.gz` are labelled preflight-only and
  never enter the formal statistics. `broad-*.jsonl.gz` retain every validation row.
- [Hash inventory](2026-10-17-morgan-linear-artifacts/manifest.json) covers every
  archived file. Baseline/candidate wheels and installed native extensions remain
  locally under `/private/tmp/chematic-morgan-linear-evidence-20261017`; this Git
  packet records their hashes, full production source manifests, flags, build logs,
  and changed function disassembly. It does not vendor comparator wheels or assert
  that source rebuilds in different paths/toolchains reproduce the identical hash.
- The original checkout's unrelated dirty documentation and other tasks were
  preserved. Only this task's inactive, regenerated test/build caches were removed
  after process/open-file checks; cleanup records distinguish disk readings from
  the bytes measured in those caches. Historical source/evidence remained intact.

To independently audit the archived formal journals in place:

```sh
python3 benchmarks/2026-10-17-morgan-linear-artifacts/audit_formal.py
```

This recomputes the audit outputs; use a copy to preserve the sealed originals.
For broad verification first decompress the four `broad-*.jsonl.gz` in a copy of
the artifact directory, then run `audit_broad.py` there. Source build commands and
all settings are in `baseline-source.json` and `candidate-build-source.json`.
The timing command below documents the recorded local host paths. For a new
host/epoch, rebuild both exact sources, run all-arm preflight, generate fresh
source/runtime/artifact seals, and create new empty output directories; do not
reuse this experiment's samples, modify its completed directories or silently
reuse its host-specific expected module paths.

```sh
PYTHONDONTWRITEBYTECODE=1 /private/tmp/chematic-cosmolkit-stereo-venv/bin/python \
  scripts/bench_python_morgan_linear.py \
  --root /Users/k_tanabe/.codex/worktrees/rdkit-cosmolkit-speed-20261017/chematic \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --offset 0 --limit 1000 \
  --source-manifest /private/tmp/chematic-morgan-linear-evidence-20261017/source-manifest.json \
  --baseline-wheel /private/tmp/chematic-morgan-linear-evidence-20261017/baseline-wheel/chematic-1.1.6-cp313-cp313-macosx_11_0_arm64.whl \
  --candidate-wheel /private/tmp/chematic-morgan-linear-evidence-20261017/candidate-wheel/chematic-1.1.6-cp313-cp313-macosx_11_0_arm64.whl \
  --baseline-path /private/tmp/chematic-morgan-linear-evidence-20261017/baseline-site \
  --candidate-path /private/tmp/chematic-morgan-linear-evidence-20261017/candidate-site \
  --output /private/tmp/chematic-morgan-linear-evidence-20261017/profile-0
```

The second original invocation used `--offset5000` (separate arguments
`--offset 5000`) and `profile-5000`. Existing plans prevent rerunning either epoch.
Main integration and a universal speed claim remain unjustified by this packet.
