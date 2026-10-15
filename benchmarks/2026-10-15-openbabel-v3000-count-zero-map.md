# V3000 bounded capacities and zero atom-map writing — 2026-10-15

The production optimization is adopted in the attached source candidate. No toolkit-wide superiority, current-main or release claim follows from this run.

Decision: Adopt the combined bounded rounded-capacity reservation and zero atom-map writer shortcut as the tenth unmerged source candidate. NCI parse 1.057x [1.050,1.067], write 1.096x [1.082,1.105] and roundtrip 1.064x [1.057,1.068] against the original nine-candidate baseline all win 21/21 blocks. Every original and continuation source-relative CI lower is above 1; no profile selects a regression follow-up. This does not make every source-relative lane a strict win: original-slice write/roundtrip, chain parse and continuation parse/roundtrip have fewer than 21 faster blocks. All 2,064 original snapshots/rewrites and 64 continuation variants agree exactly, with 848 native tests and 10 runner gates passing. Controlled chains beat pinned native Open Babel by 126.171x parse, 6.709x write and 89.918x roundtrip, each with semantic agreement and all21 faster. Diverse Open Babel profiles remain unranked for graph/stereo/rewrite differences. The two preceding capacity-only trials are rejected, separately archived and never pooled. This incremental source comparison does not resolve older unrelated Python/native negatives or establish toolkit-wide superiority, current main, hosted CI or publication.

The V3000 reader reserves known atom/bond capacities in its builder, atom/bond source-ID maps, XY coordinates and Z-coordinate storage. Reservation occurs once at BEGIN ATOM after the final COUNTS row, so repeated COUNTS does not allocate or discard stores repeatedly. Requested capacities approach the original geometric Vec growth size and remain bounded by the actual number of logical V30 rows; a huge declared count in a short input does not trigger a huge allocation. Counts remain validation values, never a limit on parsed records. Parsing order, state transitions, typed-error precedence, coordinates, stereo/geometry and metadata are unchanged. The 2D writer also bypasses integer formatting for the common zero atom-map field, appending the same exact bytes; nonzero maps keep the existing formatter. Existing frozen-writer comparisons exercise absent, zero, nonzero and maximum atom maps. No public API, dependency, tolerance or unsafe code is added.

The immediate baseline contains all nine preceding accepted candidates, including the six-decimal positive-zero parser and descriptor normalization storage reuse. All rejected direct-line, flat-geometry and two-spelling zero-z trials remain excluded. The baseline native release/LTO executable was built before both capacity trials; its exact binary and all baseline source hashes are verified and reused. The candidate is built freshly after this revision; source, binaries, corpus and Open Babel libraries are hashed. All original 2,064 records retain byte-identical parser snapshots and rewrites between the two source arms, and all 64 physical-continuation variants match the original chains in both arms.

Scoped Rust validation: 848 passed, no failures or ignored tests. This includes the full MOL library and four integration suites for 2D E/Z, wedge stereo, 3D conformers and square-planar I/O. Runner semantic gates: 10 passed. Release all-target Clippy, formatting and refactor checks passed. The frozen original COUNTS transition is compared on 27,000 count-lexeme, row-budget and source-line combinations, with a prior COUNTS row to check repeated declarations and partial mutations. Incomplete rows and block markers are also checked. A separate test repeats maximum-platform-integer declarations 16 times without allocation, then checks exact bounded capacities at BEGIN ATOM. Full typed-error Debug representations, expected counts and state discriminants agree. The original synthetic recipe, license, minimization, pinned source oracle and normalized payload hash are recorded in the new count-capacity fixture manifest. The exact-capacity trial is archived separately and rejected for reproduced prepared-writer regressions. This revision rounds reservations toward the original Vec growth capacities, still bounded by actual logical rows, and reruns every validation and fixed timing profile. A new test compares capacity hints with original Vec growth for counts 0 through 1,024. The exact same verified nine-candidate baseline executable is reused, with its build log and origin hash; new timing observations are not pooled. No new molecular records were added.

## Protocol

All fixed profiles are retained without filtering: 1,000 records each at offsets 0 and 5,000, and 64 linear chains of lengths 16 through 1,024 in steps of 16. This is an exposed development corpus on one macOS arm64 host, without a sealed holdout. Open Babel 3.2.1 native library/plugins and RDKit 2026.09.1 observer are pinned. Each original profile has 21 rotating three-arm fresh-process blocks and five full-batch repeats per mode. All builds and validation precede timing; profiles and diagnostics run sequentially. Only this run owned, completed build target was removed after retaining both measured binaries.

Ratios are comparator time / candidate time; above 1 means faster. Brackets are paired median-ratio 95% intervals from 10,000 deterministic bootstrap resamples. Strict wins require semantic agreement, CI lower > 1 and all 21 blocks faster. Every record, block, mismatch, failure and initial/adaptive negative is retained. The 2D writer directly appends the exact constant zero-map suffix for absent or zero atom maps; nonzero maps retain the original formatter. The parser and writer effects are combined in this measured candidate, without an isolated attribution of each effect.

Timers exclude startup, input loading/framing and plugin discovery. The native harness prepares a parsed batch and performs an untimed warm write before timing. Parse retains the complete diagnostic batch, prepared writing creates output strings, and roundtrip parses and writes fresh records. RSS is whole-process high-water usage including preparation and plugins; it is not per-operation allocation evidence. The avoided geometric vector growth follows from the implementation and capacity assertions; no operation-level allocation count is measured.

## Original full matrix

| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |
|---|---|---:|---:|---|---|
| 0 | parse | 162.315 | 158.413 | 1.025 [1.022, 1.025]; 21/21 **strict win** | unranked (semantic mismatch) |
| 0 | write | 21.051 | 20.149 | 1.052 [1.043, 1.058]; 20/21 | unranked (semantic mismatch) |
| 0 | roundtrip | 184.667 | 179.527 | 1.029 [1.024, 1.031]; 20/21 | unranked (semantic mismatch) |
| 5000 | parse | 46.739 | 44.286 | 1.057 [1.050, 1.067]; 21/21 **strict win** | unranked (semantic mismatch) |
| 5000 | write | 7.751 | 7.068 | 1.096 [1.082, 1.105]; 21/21 **strict win** | unranked (semantic mismatch) |
| 5000 | roundtrip | 54.711 | 51.400 | 1.064 [1.057, 1.068]; 21/21 **strict win** | unranked (semantic mismatch) |
| chains | parse | 33.270 | 32.596 | 1.023 [1.018, 1.028]; 20/21 | 126.171 [125.685, 126.476]; 21/21 **strict win** |
| chains | write | 15.582 | 14.185 | 1.096 [1.091, 1.104]; 21/21 **strict win** | 6.709 [6.643, 6.747]; 21/21 **strict win** |
| chains | roundtrip | 48.494 | 46.654 | 1.040 [1.035, 1.043]; 21/21 **strict win** | 89.918 [89.434, 90.099]; 21/21 **strict win** |

## Semantic exclusions

Open Babel disagreements remain unranked. No stereo marker is removed and no row is dropped to obtain a speed result. Diverse inputs include graph mismatches, explicit unknown-double-stereo loss and rejected rewrites; those lanes cannot support equivalent-output speed wins. Controlled-chain agreement does not establish broad chemical correctness or superiority.

| Profile | Open Babel parse graph mismatches | Parse rejected | Rewrite rejected | Rewrite graph mismatches | Unknown-stereo mismatches |
|---|---:|---:|---:|---:|---:|
| 0 | 3 | 0 | 153 | 1 | 1 |
| 5000 | 137 | 0 | 0 | 137 | 137 |
| chains | 0 | 0 | 0 | 0 | 0 |

## Separate continuation diagnostic

The predeclared source-only diagnostic splits every numeric V30 atom and bond row of all controlled chains after its id into continued physical lines. All original and continued snapshots and rewrites agree between both source arms. It uses 21 alternating two-arm fresh-process blocks, all modes and five repeats; it has no Open Babel ranking or RSS measurement and never replaces an original observation.

| Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks |
|---|---:|---:|---|
| parse | 51.279 | 50.693 | 1.010 [1.006, 1.018]; 17/21 |
| write | 15.325 | 14.149 | 1.083 [1.068, 1.087]; 21/21 **strict win** |
| roundtrip | 66.475 | 64.773 | 1.031 [1.023, 1.036]; 20/21 |

## Adaptive regression checks

Any initial baseline-relative interval wholly below 1 selects that profile for one separate 21-block source-only follow-up in all three modes. All initial results remain above; follow-ups are not pooled or substituted, and new follow-up negatives do not cause recursive retries.

```json
{
  "policy": "One separate 21-block two-source-arm follow-up per selected profile, all three modes; retain initial results; no pooling or recursive retries.",
  "profiles": {}
}
```

## Reproduction and status

[The evidence packet](2026-10-15-openbabel-v3000-count-zero-map-artifacts/README.md) contains both measured source executables, every initial/follow-up result and semantic output, full build/test logs, source manifests, fixture recipe/provenance, the exact patch and a complete build/input archive. Every 10057 candidate archive hash and 10054 reconstructed baseline hash was verified. The patch also reconstructs the changed files in an independent tiny Git repository. The [exact-capacity trial](2026-10-15-openbabel-v3000-count-capacity.md) and [rounded-capacity trial](2026-10-15-openbabel-v3000-count-rounded.md) are rejected and retained separately. No earlier negative is removed or pooled with this candidate. The immediate baseline excludes both capacity trials; this candidate revises the rounded-capacity implementation by adding the writer shortcut. The earlier direct-line, flat-geometry and two-spelling zero-z implementations remain excluded. The fixed inputs and benchmark policy are unchanged; no observation from an earlier run is pooled with this run.

Last verified shared Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared refs were not refreshed; hashes identify version 1.1.4 source candidates, not current main. This run has no hosted CI, merge or publication. Previous source-candidate negatives in unrelated Python or native lanes are not resolved by an incremental comparison. The earlier rejected direct-line candidate remains archived and excluded.
