# V3000 ASCII token splitting — 2026-10-15

The production optimization is rejected; only its report and evidence are delivered. No toolkit-wide superiority, current-main or release claim follows.

Decision: Reject the ASCII-tokenization prototype for ordinary-read acceleration. Source-relative parse medians are 1.004x [0.995,1.008] for slice0, 0.992x [0.936,1.013] for slice5000, and 0.995x [0.983,1.043] for controlled chains; no ordinary parse improvement is established. The separate continuation diagnostic improves parse 1.024x [1.016,1.031] (20/21), but this small non-strict diagnostic gain does not justify the added branch/scans for the main intended workload. Writer movement is not attributed to a new algorithm because writer code is unchanged. All original/adaptive outcomes, 872 native tests, ten runner gates and 894 frozen-reader cases pass and remain archived. No initial CI selects a negative follow-up. The chain collection interruption, lost unpersisted block15, all fourteen unchanged original saved blocks, timestamp-only Open Babel revalidation difference and exact resumption are disclosed. No production source is delivered. The next independent candidate targets direct integer output; no measurements are pooled.

The V3000 reader uses split_ascii_whitespace only when the logical payload is entirely ASCII and has no vertical tab. Rust ASCII whitespace excludes vertical tab, whereas Unicode whitespace includes it; those payloads and all non-ASCII payloads retain split_whitespace. Token-vector reuse, physical-line validation, continuation joining, parser state/error precedence, coordinates, stereo, diagnostics, chemistry and metadata are unchanged. The writer is unchanged; prepared-write timing movement is not attributed to a new writer algorithm. No public API, dependency, unsafe code or tolerance is added.

The immediate baseline contains all eleven previously accepted source candidates, including packed private aromatic-component facts. All rejected experiments remain excluded. Both native release/LTO executables were built freshly: baseline before editing with all source hashes checked, candidate in a separate fresh target. All builds, validation and cleanup of owned completed targets preceded timing. All 2,064 original records retain byte-identical snapshots and rewritten bytes between source arms; all 64 continued chains match their original records in both arms.

Scoped Rust tests: 872 passed, zero failed/ignored. The full MOL library and six integration suites cover E/Z and wedge stereo, 3D, square planar, V3000 semantic preservation and source-ID lookup. The frozen original reader is compared on 894 recipe cases: 149 separators (all ASCII characters, every Unicode whitespace, two Unicode non-whitespace boundaries), direct/continued framing and valid/invalid-index/invalid-coordinate variants of the existing ethanol fixture. Exact typed errors or indexed atoms/bonds/neighbor order, coordinate bits, metadata, geometry/stereo diagnostics and rewrite bytes agree. The recipe/manifest record source hashes, provenance, license, minimization and the pinned source oracle. No new molecular records are added. Ten benchmark runner gates, release all-target Clippy, formatting and the unchanged refactor budget pass.

## Collection interruption

The controlled-chain initial collection stopped after 14 saved blocks when opening the block-15 checkpoint failed with transient ENOSPC. The completed in-memory observation for block 15 could not be recovered and is explicitly unavailable. The 14-block incomplete matrix and failure logs are retained. resume-native.py imported the unchanged timing, observer, eligibility and bootstrap functions, revalidated all source/input/library/binary hashes and complete semantic outcomes/source output hashes, then recollected block 15 and continued the identical rotation through block 21. All fourteen previously saved blocks are unchanged. Both diverse profiles completed without interruption and were not repeated. These boundaries and host load limit small-difference interpretation. A temporary owned checkpoint reserve is released on later ENOSPC, and computed observations remain in memory until saved. Open Babel rewrites changed only their program-line timestamps on all 64 records; every other byte and all observer results were identical on revalidation, as recorded in resume-semantic-verification.json. The initial timestamp-hash assertion failure and repair are retained. No timed native operation changed.

## Protocol

Fixed unfiltered development profiles: 1,000 records at offsets 0 and 5,000, plus 64 controlled linear chains of lengths 16 through 1,024. One macOS arm64 host; Open Babel 3.2.1 native libraries/plugins and RDKit 2026.09.1 observer pinned by hashes. Each profile has 21 rotating three-arm fresh-process blocks, five batch repeats, and all parse/prepared-write/roundtrip modes. No sealed holdout or toolkit-wide coverage.

Ratios are comparator elapsed time / candidate elapsed time; above 1 is faster. Brackets are paired median-ratio 95% intervals from 10,000 deterministic bootstrap resamples. Strict wins require semantic agreement, CI lower > 1 and all 21 blocks faster. All successfully saved original and follow-up observations are retained without pooling; the unavailable observation and recollection are disclosed above.

Timers exclude startup, input loading/framing and plugin discovery. Preparation retains a full parsed batch and includes one untimed warm write; parse evaluates new batches with default diagnostics, prepared-write creates output strings, roundtrip parses and writes fresh records. RSS is whole-process high-water usage including preparation/plugins; it is not per-operation allocation evidence. Source, input, native executable and library hashes are checked.

## Original full matrix

| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |
|---|---|---:|---:|---|---|
| 0 | parse | 162.459 | 161.247 | 1.004 [0.995, 1.008]; 12/21 | unranked (semantic mismatch) |
| 0 | write | 20.412 | 20.185 | 1.019 [1.012, 1.031]; 18/21 | unranked (semantic mismatch) |
| 0 | roundtrip | 181.471 | 182.600 | 0.996 [0.993, 1.003]; 9/21 | unranked (semantic mismatch) |
| 5000 | parse | 83.112 | 82.736 | 0.992 [0.936, 1.013]; 9/21 | unranked (semantic mismatch) |
| 5000 | write | 12.300 | 11.743 | 1.024 [0.940, 1.066]; 13/21 | unranked (semantic mismatch) |
| 5000 | roundtrip | 96.746 | 95.612 | 1.036 [1.009, 1.042]; 15/21 | unranked (semantic mismatch) |
| chains | parse | 39.525 | 33.807 | 0.995 [0.983, 1.043]; 9/21 | 124.186 [122.632, 125.312]; 21/21 **strict win** |
| chains | write | 16.708 | 14.668 | 1.048 [1.018, 1.076]; 18/21 | 7.031 [6.804, 7.175]; 21/21 **strict win** |
| chains | roundtrip | 57.164 | 48.789 | 0.998 [0.995, 1.007]; 9/21 | 90.317 [87.708, 91.293]; 21/21 **strict win** |

## Semantic exclusions

Diverse Open Babel lanes remain unranked when internal graphs, stereo or rewrite fields differ. No row or stereo marker is removed. Controlled-chain agreement establishes only that checked workload.

| Profile | Open Babel parse graph mismatches | Parse rejected | Rewrite rejected | Rewrite graph mismatches | Unknown-stereo mismatches |
|---|---:|---:|---:|---:|---:|
| 0 | 3 | 0 | 153 | 1 | 1 |
| 5000 | 137 | 0 | 0 | 137 | 137 |
| chains | 0 | 0 | 0 | 0 | 0 |

## Separate continuation diagnostic

Every numeric atom/bond row of all 64 controlled chains is split after its id into continued physical lines. All snapshots and rewrites match the original chains in both arms. The source-only diagnostic uses 21 alternating two-arm fresh-process blocks, five repeats and all modes, without Open Babel ranking or RSS measurement. It never replaces an original profile.

| Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks |
|---|---:|---:|---|
| parse | 56.682 | 55.074 | 1.024 [1.016, 1.031]; 20/21 |
| write | 15.297 | 15.021 | 1.019 [1.016, 1.027]; 18/21 |
| roundtrip | 71.885 | 70.199 | 1.028 [1.013, 1.030]; 18/21 |

## Adaptive regression checks

The predeclared rule selects an entire profile for one separate 21-block source-only follow-up in all modes if any initial source-relative CI is wholly below 1. Initial and follow-up outcomes remain separate. Newly appearing negatives do not trigger recursive retries.

```json
{
  "policy": "One separate 21-block two-source-arm follow-up per selected profile, all three modes; retain initial results; no pooling or recursive retries.",
  "profiles": {}
}
```

## Reproduction and status

[The evidence packet](2026-10-15-openbabel-v3000-ascii-artifacts/README.md) retains exact measured binaries, all matrices/snapshots/rewrites, failures if any, test/build logs, full frozen source/input archive and exact patch. Archive readback verifies all 10061 candidate hashes and reconstructs all 10058 baseline hashes. A separate tiny Git repository verifies patch reconstruction. Source reconstruction is not bit-reproducibility of rebuilt binaries.

Last verified shared Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`; refs were not refreshed. These are version 1.1.4 unmerged source candidates, with no hosted CI or publication. Prior unrelated source-candidate regressions against older baselines remain unresolved by this incremental comparison.
