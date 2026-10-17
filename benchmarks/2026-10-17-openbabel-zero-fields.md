# V2000 zero-attribute tail: source performance gate passed

The V2000 atom-row writer appends one identical 33-byte constant when mass
difference, encoded charge, valence and atom map are all zero. It still computes
every attribute normally and keeps the complete old path for nonzero fields.
Numeric formatting, stereo perception, public APIs, dependencies and caches are
unchanged. The fast return exits only the atom-row helper, preserving all later
bond/property/group/terminator output.

Both predeclared SDF-write primary confidence lower bounds exceed one and all
30 unrounded paired source medians clear the 0.99 retention floor. The completed
local source gate therefore **passes**. There is a disclosed control retreat of
about 0.6%; this does not mean every operation became faster. Integration remains
subject to the exact final commit's hosted CI and protected-main checks.

This is one-host native source evidence, not a published-package, fresh holdout,
cross-host or whole-Open-Babel superiority claim. Earlier rejected experiments
in [PR #984](https://github.com/kent-tokyo/chematic/pull/984) and
[PR #993](https://github.com/kent-tokyo/chematic/pull/993) remain rejected;
none of their elapsed observations are reused or multiplied into these results.

## Fixed source, artifacts and host

- Baseline `c84c68f297f95562d1d29c75234f2190cd5f9d75` (including the main 3D stereo
  correctness follow-up). Measured candidate
  `a0b47272beb93af19bed1f2fe303237affbc0970`. Later commits add evidence only.
- Both initially empty independent targets use opt-level3, full LTO, one CGU,
  one build job, incremental=false, strip=none, locked offline dependencies and
  Rust 1.97.0/LLVM 22.1.6. Compiler invocation logs prove both MOL production
  compilations. Seven probe-only files are byte-identical between arms; the
  Python measurement adapter is also fixed. No benchmark-specific product path
  is introduced.
- The actual SDF writer owner changes from 2,893 to 2,925 instructions. Both
  disassemblies and distinct native binary hashes are retained. This proves
  realization of the change, not a speed prediction or causal explanation of
  unrelated control fluctuations.
- macOS 27.0.1 arm64, Apple M4, Mac16,12, 10 CPUs, 16 GiB RAM. CPU details were
  collected after measurement. Host load averages were recorded at every block;
  external activity/CPU frequency were not controlled.
- Open Babel 3.2.1 (`Jul 11 2026 -- 19:21:27`) native C++ -O3 probes, installed
  prefix `/opt/homebrew/Cellar/open-babel/3.2.1`. Exact prebuilt probes and their
  identical source/header hashes are reused; no earlier times are used.
  Package-local SDK/runtime/data files are sealed before/after profiles.
  System and package-external runtime files are outside that seal.
- Independent semantic observer: RDKit 2025.09.3, CPython 3.13.6, NumPy 2.5.3.
  These are pinned installed comparator/observer identities, not a claim about
  the latest upstream release.

The [prospective plan](2026-10-17-openbabel-zero-fields-artifacts/measurement-plan.json),
[source manifest](2026-10-17-openbabel-zero-fields-artifacts/source-manifest.json),
[artifact/runtime seal](2026-10-17-openbabel-zero-fields-artifacts/seal.json) and
[native realization](2026-10-17-openbabel-zero-fields-artifacts/code-realization.json)
retain exact sources, inputs, flags, paths, hashes and commands. All 2,110 source
files and 427 artifact/runtime files matched their seals again after all timings.

## Workload and timing boundary

The same ten **unfiltered** profiles total 4,384 records: two diverse 1,000-row
SDF slices, corresponding V3000 slices, 64 linear chains and continued versions,
plus two arbitrary-coordinate 3D chain profiles and two F-branched variants.
Synthetic sizes span 16..1024 carbons; their deterministic nonplanar coordinates
are not physical conformers or a 3D-quality evaluation. Every original corpus
and SHA-256 is in the plan and archived compressed input.

Each profile has 24 fresh-process blocks. SDF uses all six three-arm
permutations four times. The other eight profiles use 12 AB and 12 BA source
orders. Every profile's exact fixed order is retained and independently checked.
Fresh external speed coverage is only the two SDF slices. V3000/3D are source
retention controls here; this packet adds no fresh external ranking for them.

Each native mode warms for at least 500 ms, then times at least 100 ms and at least
three full batches. Raw duration, actual repetitions, actual processed records,
output units, warmup counts/budgets and whole-process peak RSS stay in each
journal. Per-batch ns is raw elapsed divided by actual repetitions, floored by
less than 1 ns. Faster implementations can process more full batches. Ratios are
paired by block; marginal medians need not have the same ratio as paired medians.
Bootstrap uses 10,000 resamples with seed 0 source/seed 1 Open Babel. Unrounded
source medians decide retention; all raw negative and inconclusive observations
are retained. No timings were rerun, filtered or pooled across epochs.

Input loading, native startup and plugin discovery are outside timers. Parse and
roundtrip create fresh result objects each batch; prepared write reuses its
parsed objects. SDF framing is timed, V3000 framing is outside timing. Rust SDF
reuses a caller-owned String; V3000 and Open Babel allocate. This describes these
native API boundaries, not an isolated universal serialization algorithm. RSS
includes preparation/warmup/loading and is not isolated-operation or streaming
memory evidence; it is not ranked. These are warmed throughput measurements,
not cold-start latency.

Every native call was journalled atomically before the next call, with at least
1 GiB free checked before each block. Full source snapshots and full rewritten
streams were checked in preflight after the same write warmup/timing budget.
Formal calls check counts and output-unit totals against those preflight modes;
this protocol does **not** retain a separate content digest for every inner
native repetition. That boundary is distinct from byte equality of the saved
full preflight outcomes.

## Complete source matrix

Adoption requires both `sdf-0/write` and `sdf-5000/write` CI lower bounds>1,
every one of 30 unrounded source median ratios>=0.99, exact source outcomes,
and all planned profiles complete. All conditions pass. Both primary intervals
are positive, with 22/24 and 23/24 faster blocks; source adoption did not require
24/24 wins. No threshold exception was applied.

Milliseconds are per complete batch from actual repetition counts.

| Profile | Mode | Base ms | Candidate ms | Paired base/candidate [95% CI] | Faster blocks |
|---|---|---:|---:|---|---:|
| sdf-0 | parse | 12.180 | 12.204 | 0.998706 [0.994223, 1.002089] | 9/24 |
| sdf-0 | write | 3.729 | 3.588 | 1.038636 [1.026344, 1.059909] | 22/24 |
| sdf-0 | roundtrip | 16.297 | 16.113 | 1.010892 [1.008310, 1.014802] | 22/24 |
| sdf-5000 | parse | 4.444 | 4.456 | 0.998004 [0.986758, 1.004428] | 9/24 |
| sdf-5000 | write | 1.529 | 1.469 | 1.040083 [1.036150, 1.048394] | 23/24 |
| sdf-5000 | roundtrip | 6.117 | 6.076 | 1.006501 [0.999908, 1.015400] | 16/24 |
| v3000-0 | parse | 12.905 | 12.889 | 1.000418 [0.998776, 1.005269] | 12/24 |
| v3000-0 | write | 2.692 | 2.686 | 1.003536 [0.999769, 1.007853] | 16/24 |
| v3000-0 | roundtrip | 15.917 | 15.931 | 0.999025 [0.994913, 1.003102] | 11/24 |
| v3000-5000 | parse | 4.914 | 4.903 | 1.004461 [0.995772, 1.012532] | 13/24 |
| v3000-5000 | write | 0.696 | 0.695 | 1.002243 [0.999932, 1.006619] | 16/24 |
| v3000-5000 | roundtrip | 5.651 | 5.635 | 1.003230 [0.998038, 1.009104] | 16/24 |
| v3000-chains | parse | 5.539 | 5.518 | 1.006068 [0.995643, 1.008412] | 14/24 |
| v3000-chains | write | 1.231 | 1.236 | 0.998669 [0.991071, 1.004904] | 12/24 |
| v3000-chains | roundtrip | 6.820 | 6.820 | 1.002466 [0.993018, 1.008929] | 13/24 |
| v3000-continued | parse | 7.734 | 7.748 | 0.996771 [0.991504, 1.000194] | 8/24 |
| v3000-continued | write | 1.229 | 1.231 | 0.997794 [0.993333, 1.005812] | 10/24 |
| v3000-continued | roundtrip | 9.033 | 9.027 | 1.000997 [0.993814, 1.005387] | 13/24 |
| v3000_3d-spiral | parse | 7.454 | 7.440 | 1.000244 [0.998618, 1.008507] | 12/24 |
| v3000_3d-spiral | write | 2.425 | 2.441 | 1.001627 [0.985351, 1.011833] | 12/24 |
| v3000_3d-spiral | roundtrip | 9.912 | 9.913 | 1.001803 [0.995617, 1.007115] | 13/24 |
| v3000_3d-branched-spiral | parse | 84.144 | 83.603 | 1.009221 [0.999397, 1.016105] | 16/24 |
| v3000_3d-branched-spiral | write | 2.422 | 2.434 | 0.993975 [0.987098, 1.004820] | 9/24 |
| v3000_3d-branched-spiral | roundtrip | 86.543 | 86.226 | 1.001979 [0.998994, 1.006314] | 15/24 |
| v3000_3d-offset | parse | 7.555 | 7.554 | 1.000507 [0.996422, 1.003209] | 12/24 |
| v3000_3d-offset | write | 2.414 | 2.413 | 0.996118 [0.991244, 1.007359] | 10/24 |
| v3000_3d-offset | roundtrip | 9.979 | 10.018 | 0.996140 [0.993971, 1.001074] | 9/24 |
| v3000_3d-branched-offset | parse | 84.650 | 84.014 | 1.005696 [0.998732, 1.015912] | 15/24 |
| v3000_3d-branched-offset | write | 2.403 | 2.416 | 1.004074 [0.987525, 1.013936] | 13/24 |
| v3000_3d-branched-offset | roundtrip | 87.009 | 86.496 | 1.004514 [0.997031, 1.013309] | 15/24 |

The lowest source median is branched-spiral V3000 3D **write**:
0.993975 [0.987098,1.004820],9/24 faster. Its ~0.6% observed median retreat is
inside the predeclared 1% allowance; its interval is inconclusive. There are no
source intervals with an unrounded upper endpoint below one. Every source median
below one and every unrounded interval is retained in the audit/CSV above.

## Fresh Open Babel SDF observations

External ranking requires matching accepted canonical graphs for parse;
write/roundtrip also require indexed atom/bond fields, coordinates to 1e-4,
names/properties and preservation of explicit unknown stereo markers. Additional
unknown markers are allowed by the declared normalization. Source-to-source
snapshot and rewrite checks remain byte-exact.

All 1,000 rows in each SDF slice pass source parse and rewritten-field checks.
Both Open Babel rewritten-field checks also pass. Its first-slice canonical-output
parse check differs at rows 90 and 745, so that lane remains unranked. No row is
removed. Full original/rewritten observations and every mismatch index are saved;
this boundary does not attribute the difference exclusively to Open Babel's
reader rather than its canonical writer/observer path.

Strict external wins additionally require CI lower>1 and 24/24 faster blocks.
All five eligible lanes pass. These are direct new observations, not products
of earlier source/external ratios.

| Profile | Mode | Open Babel/candidate paired ratio [95% CI] | Faster blocks |
|---|---|---|---:|
| sdf-0 | parse | unranked: 2 canonical-output disagreements | — |
| sdf-0 | write | 6.508426 [6.451382, 6.519333] | 24/24 |
| sdf-0 | roundtrip | 7.240532 [7.224030, 7.250846] | 24/24 |
| sdf-5000 | parse | 9.217602 [9.135989, 9.248472] | 24/24 |
| sdf-5000 | write | 9.201869 [9.125039, 9.305987] | 24/24 |
| sdf-5000 | roundtrip | 9.026827 [8.995872, 9.106416] | 24/24 |

## Correctness and audit

- 1,388 MOL/SDF tests passed in 67 suites, zero failed/ignored. A complete frozen
  baseline atom-row function differs only in its declaration name and matches
  every atom of 14,072 normal/unusual graphs, 15,120 metadata combinations with
  four coordinate choices, and R-group-symbol boundaries. The fixture covers
  natural/explicit isotopes, charges, explicit H/valence, zero and overflow maps,
  wildcards, negative zero, rounding boundaries and nonfinite formatting fallback.
  This proves the tested row bytes are retained; it is not a new broad chemistry
  parity claim for unsupported V2000 fields.
- All 4,384 full native parse snapshots and warmed rewrites match between arms,
  including coordinates, atom/bond state, diagnostics, metadata and property
  content. No molecular row or output field was dropped to obtain speed.
- All-target release Clippy, rustfmt, unchanged refactor ratchet and 24 existing
  measurement/semantic runner tests passed. The normal commit-hook Clippy also
  passed. Test-only API spellings were corrected before formal timing; initial
  diagnostics and the correction record remain. Identical release commands after
  those cfg(test)-only edits returned the same three copied native binaries.
- Independent audit verifies every one of 1,584 atomic native timing journals,
  all 24 orders/profile and both-arm/full-batch counts, warmup/timing budgets,
  elapsed normalization and constant per-batch output units. It checks 40,777
  complete timed batches and 20,536,000 processed records and recomputes all 30
  source plus six external statistical vectors. The unranked vector is retained
  as raw data with `ranked=false`; it cannot support a speed claim.
- A second portable audit was run on a separate archive copy, without original
  corpus files or SDK paths. It checks the archive SHA inventory, decompresses
  inputs and source outcomes, and recomputes the same statistical decision.
- The original dirty checkout's four documentation edits and unrelated work are
  untouched. Only owned inactive regenerable caches were removed after process,
  open-file and copied-artifact checks. Source, source archives, binaries,
  diagnostics, logs and old experiments remain intact. The cleanup logs separate
  measured cache sizes from free-space changes caused by concurrent activity.

## Evidence and reproduction

- [Independent audit](2026-10-17-openbabel-zero-fields-artifacts/audit.json),
  [unrounded 36-vector CSV](2026-10-17-openbabel-zero-fields-artifacts/lanes.csv),
  [decision](2026-10-17-openbabel-zero-fields-artifacts/decision.json),
  [validation](2026-10-17-openbabel-zero-fields-artifacts/validation.json), and
  [semantic counts/indices](2026-10-17-openbabel-zero-fields-artifacts/semantics-summary.json).
- All ten `paired-*.json` records and 1,584 individual `*-journals/*.json` files
  retain the original completed experiment. Every negative/inconclusive source
  result, all loads/RSS and all unranked external data stay available.
- Compressed inputs, native full preflight streams and observer outputs preserve
  their original bytes. Observer compression provenance gives uncompressed SHA
  values. The [manifest](2026-10-17-openbabel-zero-fields-artifacts/manifest.json)
  hashes every archived file.
- Exact six source binaries, two comparator binaries and the baseline source
  archive remain locally under `/private/tmp/chematic-v2000-zero-fields-evidence-20261017`;
  their hashes/paths are recorded, not vendored into this Git packet. Build
  invocations, complete source manifests, probe patch, native disassembly,
  compiler identities and the candidate diff are retained. Different paths or
  toolchains are not asserted to reproduce the same binary hash.

To recheck statistics and archived outputs, copy the artifact directory and run:

```sh
python3 audit.py --archive-only
```

The portable mode checks the archived SHA inventory and writes
`archive-audit.json`/`archive-lanes.csv`, preserving the original results. It
cannot re-run native code or certify the SDK currently installed on another
host. Native reproduction needs the exact baseline/candidate sources, identical
common probes, pinned SDK, build flags and fresh source/runtime seals documented
in the manifests. The original local sequence was:

```sh
python3 build.py baseline
python3 build.py candidate
python3 formal.py --preflight-only
# Validate source/tests/code realization and create the recorded seal.
python3 formal.py
python3 audit.py
```

Those helper paths refer to the retained local evidence directory and source
checkout paths recorded in their files. New experiments require new empty build
and result directories and fresh seals; existing completed epochs reject reruns.
The complete local gate and hosted CI/main integration are separate decisions.
