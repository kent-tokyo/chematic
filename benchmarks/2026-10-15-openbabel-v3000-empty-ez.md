# Skip unused E/Z work without double bonds — 2026-10-15

The shared 2D E/Z reader now returns immediately when its literal-double-bond list is empty. The old function also returned no directions or diagnostics in that case, but first computed aromaticity and ring information that no later phase used. Existing atom/bond annotations remain untouched. For nonempty lists, the same aromaticity query stays in the caller with its original lifetime; the unchanged per-bond classification loop moves into a small private helper. Carrier reconciliation, diagnostic ordering and commits are unchanged.

This source-version 1.1.4 candidate excludes the rejected direct-line experiment. The baseline still contains exactly the five preceding candidates: cycle-root pruning, V3000 writer buffers, shared Kekule descriptor views, V3000 parser token reuse and ring-component reuse. Both release/LTO source executables were freshly built in a separate snapshot before/after this change. Last verified Git base is `cac480fb23cd8adeec6e5e9bed1f18677f87b392`; shared refs were not refreshed after earlier stalls. Hashes identify source candidates, not current main or a registry release.

All 2,064 original fixed records have byte-identical baseline/candidate parser snapshots and rewrites. All 64 continued-chain variants also match the original-chain snapshots and rewrites in both arms. Source framing and parser token storage are unchanged in this candidate.

Scoped Rust validation: 1,063 passed, 0 failed, 1 existing ignored test. The complete perception/MOL library suites and both 2D-stereo reader integration suites ran. Runner semantic-gate tests: 10 passed. Release all-target Clippy for both crates, formatting and refactor checks passed. A frozen copy of the old complete E/Z function matches all E/Z diagnostics, directions, graph fields and annotations on 10,000 fixed inputs and targeted cases: no-double graphs, existing directions, aromatic/query/dative orders, explicitly unspecified double bonds, missing/nonfinite coordinates and graph edits.

The earlier [direct-line candidate](2026-10-15-openbabel-v3000-direct-lines.md) was rejected: it gave no clear original-profile parse improvement and reproduced writer regressions in separate follow-ups. All of its source, tests, timings and negatives remain archived; none of its production changes is included here. Results from separate candidates are not pooled.

## Protocol

Unchanged and unfiltered fixed profiles: 1,000 records each at offsets 0 and 5,000, and all 64 linear chains of lengths 16, 32, ..., 1,024. Exposed development corpora on one macOS arm64 host; no sealed holdout or toolkit-wide claim. Open Babel 3.2.1 native C++ library/plugins and RDKit 2026.09.1 observer are pinned. Each profile has 21 rotating three-arm fresh-process blocks, five full-batch repeats per mode. Builds/tests precede all timing; profiles, diagnostics and any adaptive follow-ups run sequentially. Source, binaries, libraries and inputs are verified before/after measurement.

Native API timers exclude startup, input loading, delimiter splitting and plugin discovery. The unchanged harness prepares a parsed batch and performs an untimed warm write even before parse timing. Parse retains the full diagnostic batch; prepared writing allocates strings after one warm write; roundtrip actually parses and writes a fresh batch. RSS is the whole-process high-water mark including preparation/plugins/allocator state, not per-operation allocation.

Ratios are comparator elapsed time / candidate elapsed time; greater than 1 means faster. Brackets show paired median-ratio 95% intervals from 10,000 deterministic bootstrap resamples. A strict speed win requires semantic agreement, CI lower > 1 and all 21 blocks faster. All input records, blocks, errors, mismatches and output markers are retained. Intervals crossing 1 do not establish an improvement. The writer is unchanged; writer movement is recorded without attributing it to a new writer algorithm.

## Original full matrix

| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |
|---|---|---:|---:|---|---|
| 0 | parse | 163.624 | 163.535 | 1.001 [0.998, 1.006]; 12/21 | unranked (semantic mismatch) |
| 0 | write | 21.259 | 21.354 | 1.001 [0.991, 1.007]; 11/21 | unranked (semantic mismatch) |
| 0 | roundtrip | 186.210 | 186.291 | 1.001 [0.996, 1.005]; 12/21 | unranked (semantic mismatch) |
| 5000 | parse | 48.191 | 47.326 | 1.014 [1.013, 1.021]; 20/21 | unranked (semantic mismatch) |
| 5000 | write | 7.692 | 7.759 | 0.992 [0.980, 1.003]; 7/21 | unranked (semantic mismatch) |
| 5000 | roundtrip | 55.838 | 55.227 | 1.013 [1.005, 1.014]; 18/21 | unranked (semantic mismatch) |
| chains | parse | 36.115 | 34.460 | 1.046 [1.040, 1.054]; 21/21 **strict win** | 119.877 [119.126, 120.666]; 21/21 **strict win** |
| chains | write | 15.348 | 15.248 | 0.998 [0.993, 1.011]; 9/21 | 6.232 [6.179, 6.312]; 21/21 **strict win** |
| chains | roundtrip | 51.279 | 49.593 | 1.038 [1.029, 1.042]; 19/21 | 85.171 [84.250, 85.858]; 21/21 **strict win** |

## Semantic exclusions

Open Babel disagreements remain unranked. No marker was removed or input filtered to make its results pass. The diverse-profile original observations include unknown-double-bond stereo loss and rejected rewrites; they are neither equivalent-output speed wins nor failures of the new source candidate.

| Profile | Open Babel parse graph mismatches | Parse rejected | Rewrite rejected | Rewrite graph mismatches | Unknown-stereo mismatches |
|---|---:|---:|---:|---:|---:|
| 0 | 3 | 0 | 153 | 1 | 1 |
| 5000 | 137 | 0 | 0 | 137 | 137 |
| chains | 0 | 0 | 0 | 0 | 0 |

## Continuation fallback diagnostic

The separately declared source-only diagnostic splits every numeric V30 atom and bond record of all 64 controlled chains after its id into two physical lines. Original and continued snapshots and rewrites match byte for byte in both source arms. It runs 21 alternating two-source-arm fresh-process blocks with the same five-repeat native timers and all three modes, without RSS or an Open Babel arm. It checks the same E/Z shortcut under continued physical framing; all legacy parser behavior is identical in both arms. It never replaces an original result.

| Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks |
|---|---:|---:|---|
| parse | 55.067 | 53.006 | 1.035 [1.032, 1.047]; 20/21 |
| write | 15.586 | 15.405 | 1.010 [0.996, 1.022]; 11/21 |
| roundtrip | 70.628 | 68.383 | 1.028 [1.022, 1.038]; 21/21 **strict win** |

## Adaptive negative observations

Every original baseline-relative interval wholly below 1 selects that profile for one separate source-only 21-block follow-up, retaining all three modes. The same rule applies to the initial continuation diagnostic. No pooling, substitution of original results or recursive retries. Newly appearing follow-up negatives remain unresolved observations.

```json
{
  "policy": "One separate 21-block two-source-arm follow-up per selected profile, all three modes; retain initial results; no pooling or recursive retries.",
  "profiles": {}
}
```

## Limits and evidence

Complete Open Babel-wide speed superiority is not established. Controlled-chain wins apply to that profile; diverse semantic exclusions, other formats, rich query/SGROUP chemistry, cold startup, filesystem/CLI/Python overhead, general 3D, independent corpora, other hosts and published artifacts remain separate evaluations. Earlier candidates’ regressions against older checkpoints are not resolved by this incremental V3000 comparison. No main merge, hosted CI, release or publication was performed.

[Evidence packet](2026-10-15-openbabel-v3000-empty-ez-artifacts/README.md): every original/adaptive block, load observation, mismatch, source and binary/library/input hash, exact output streams, logs, measured source executables, reconstructable patch and complete frozen build/input archive. Archive bytes and both reconstructed source trees are verified. The [existing corpus provenance](../validation/benchmark_corpora/chemistry-speed-followups-v1.manifest.json) and [third-party notices](../THIRD_PARTY_NOTICES.md) apply.
