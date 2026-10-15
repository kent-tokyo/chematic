# V3000 direct payload consumption — 2026-10-15

Rejected for adoption: the direct-line candidate did not improve original-profile parsing and repeated source-only follow-ups reproduced writer regressions of roughly 3–4%. All results and the source patch are retained. No production source or changelog from this candidate is applied to the attached worktree.

Ordinary V3000 payloads now go directly from validated physical lines into the parser and reused token buffer. The logical-line vector is deferred until the first continuation; that entire suffix uses the unchanged continuation folder and a reused suffix token buffer. This removes the intermediate logical-line array on inputs without continuations. Physical-line validation still inspects the whole input before parsing, preserving error precedence. Unicode token splitting, source line numbers, metadata, diagnostics and chemical interpretation are unchanged. This is not a constant-memory file-streaming parser.

Unmerged source version 1.1.4 candidate. The immediate baseline includes all five previous candidates: cycle-root pruning, V3000 writer buffers, shared Kekulé descriptor views, V3000 parser token reuse and ring-component reuse. Both native release/LTO binaries were freshly built; the baseline source was verified before editing. Last verified Git base is `cac480fb23cd8adeec6e5e9bed1f18677f87b392`; shared Git refs were not refreshed after earlier stalls. Source and binary hashes identify this experiment, not current main, a release or a registry artifact.

All 2,064 original fixed records retain byte-identical baseline/candidate internal snapshots and rewrites. All 64 continued-chain variants also retain the original-chain snapshots and rewrite bytes in both source arms. The continuation diagnostic transforms physical framing only; it is not an additional diverse-chemistry or Open Babel ranking claim.

Full MOL tests: 1,018 passed, 0 failed, 0 ignored. Runner semantic-gate tests: 10 passed. Release all-target MOL Clippy, workspace formatting and the refactor budget passed. After Clippy identified one redundant borrow in an added test, that style-only fix passed all five added tests again; the full 1,018-test run preceded the fix and its log is retained. Frozen-reader differential tests compare complete graphs, coordinate bits, metadata, stereo/group fields, 2D/3D diagnostics, rewrites and displayed errors. They cover every continuation split position, malformed/dangling continuations, first-physical-line error precedence, CRLF/no final newline, Unicode whitespace, malformed numbers and duplicate/opaque attributes.

## Protocol

Unchanged and unfiltered fixed profiles: 1,000 records each at offsets 0 and 5,000, and all 64 linear chains of lengths 16, 32, ..., 1,024. Exposed development corpora on one macOS arm64 host; no sealed holdout or toolkit-wide claim. Open Babel 3.2.1 native C++ library/plugins and RDKit 2026.09.1 observer are pinned. Each profile has 21 rotating three-arm fresh-process blocks, five full-batch repeats per mode. Builds/tests precede all timing; profiles, diagnostics and any adaptive follow-ups run sequentially. Source, binaries, libraries and inputs are verified before/after measurement.

Native API timers exclude startup, input loading, delimiter splitting and plugin discovery. The unchanged harness prepares a parsed batch and performs an untimed warm write even before parse timing. Parse retains the full diagnostic batch; prepared writing allocates strings after one warm write; roundtrip actually parses and writes a fresh batch. RSS is the whole-process high-water mark including preparation/plugins/allocator state, not per-operation allocation.

Ratios are comparator elapsed time / candidate elapsed time; greater than 1 means faster. Brackets show paired median-ratio 95% intervals from 10,000 deterministic bootstrap resamples. A strict speed win requires semantic agreement, CI lower > 1 and all 21 blocks faster. All input records, blocks, errors, mismatches and output markers are retained. Intervals crossing 1 do not establish an improvement. The writer is unchanged; writer movement is recorded without attributing it to a new writer algorithm.

## Original full matrix

| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |
|---|---|---:|---:|---|---|
| 0 | parse | 289.790 | 292.817 | 0.995 [0.970, 1.029]; 9/21 | unranked (semantic mismatch) |
| 0 | write | 37.241 | 39.003 | 0.967 [0.915, 0.977]; 6/21 | unranked (semantic mismatch) |
| 0 | roundtrip | 322.314 | 337.592 | 0.992 [0.937, 1.007]; 10/21 | unranked (semantic mismatch) |
| 5000 | parse | 62.919 | 63.045 | 0.997 [0.986, 1.018]; 9/21 | unranked (semantic mismatch) |
| 5000 | write | 9.952 | 10.351 | 0.965 [0.951, 0.975]; 3/21 | unranked (semantic mismatch) |
| 5000 | roundtrip | 73.106 | 73.487 | 0.994 [0.989, 1.002]; 8/21 | unranked (semantic mismatch) |
| chains | parse | 36.153 | 36.106 | 0.999 [0.997, 1.013]; 10/21 | 113.197 [111.970, 114.186]; 21/21 **strict win** |
| chains | write | 15.387 | 15.899 | 0.961 [0.946, 0.972]; 1/21 | 5.997 [5.952, 6.028]; 21/21 **strict win** |
| chains | roundtrip | 51.761 | 51.998 | 0.995 [0.989, 0.998]; 5/21 | 80.714 [79.943, 81.447]; 21/21 **strict win** |

## Semantic exclusions

Open Babel disagreements remain unranked. No marker was removed or input filtered to make its results pass. The diverse-profile original observations include unknown-double-bond stereo loss and rejected rewrites; they are neither equivalent-output speed wins nor failures of the new source candidate.

| Profile | Open Babel parse graph mismatches | Parse rejected | Rewrite rejected | Rewrite graph mismatches | Unknown-stereo mismatches |
|---|---:|---:|---:|---:|---:|
| 0 | 3 | 0 | 153 | 1 | 1 |
| 5000 | 137 | 0 | 0 | 137 | 137 |
| chains | 0 | 0 | 0 | 0 | 0 |

## Continuation fallback diagnostic

The separately declared source-only diagnostic splits every numeric V30 atom and bond record of all 64 controlled chains after its id into two physical lines. Original and continued snapshots and rewrites match byte for byte in both source arms. It runs 21 alternating two-source-arm fresh-process blocks with the same five-repeat native timers and all three modes, without RSS or an Open Babel arm. It diagnoses the legacy suffix path rather than replacing any original result.

| Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks |
|---|---:|---:|---|
| parse | 54.210 | 54.723 | 0.990 [0.987, 0.995]; 4/21 |
| write | 15.247 | 15.906 | 0.957 [0.950, 0.962]; 1/21 |
| roundtrip | 69.778 | 70.761 | 0.983 [0.977, 0.993]; 2/21 |

## Adaptive negative observations

Every original baseline-relative interval wholly below 1 selects that profile for one separate source-only 21-block follow-up, retaining all three modes. The same rule applies to the initial continuation diagnostic. No pooling, substitution of original results or recursive retries. Newly appearing follow-up negatives remain unresolved observations.

```json
{
  "policy": "One separate 21-block two-source-arm follow-up per selected profile, all three modes; retain initial results; no pooling or recursive retries.",
  "profiles": {
    "0": {
      "negative_initial_lanes": {
        "write": {
          "median": 0.967,
          "ci95": [
            0.915,
            0.977
          ],
          "winning_blocks": 6,
          "blocks": 21
        }
      },
      "corpus": "/private/tmp/chematic-openbabel-v3000-20261014/corpus-0.molblocks",
      "records": 1000,
      "gate": {
        "source_snapshot_sha256": {
          "baseline": "e2fce68d8c33616317f46397cf79ee877ef856460a35b84a622ddff1d8eced74",
          "candidate": "e2fce68d8c33616317f46397cf79ee877ef856460a35b84a622ddff1d8eced74"
        },
        "byte_identical_rewrites": true,
        "source_parse_snapshots_equal": true
      }
    },
    "5000": {
      "negative_initial_lanes": {
        "write": {
          "median": 0.965,
          "ci95": [
            0.951,
            0.975
          ],
          "winning_blocks": 3,
          "blocks": 21
        }
      },
      "corpus": "/private/tmp/chematic-openbabel-v3000-20261014/corpus-5000.molblocks",
      "records": 1000,
      "gate": {
        "source_snapshot_sha256": {
          "baseline": "7fa169eebb0d112936fd4e135d08702eb1150faa6791111f52792c51d504d67d",
          "candidate": "7fa169eebb0d112936fd4e135d08702eb1150faa6791111f52792c51d504d67d"
        },
        "byte_identical_rewrites": true,
        "source_parse_snapshots_equal": true
      }
    },
    "chains": {
      "negative_initial_lanes": {
        "write": {
          "median": 0.961,
          "ci95": [
            0.946,
            0.972
          ],
          "winning_blocks": 1,
          "blocks": 21
        },
        "roundtrip": {
          "median": 0.995,
          "ci95": [
            0.989,
            0.998
          ],
          "winning_blocks": 5,
          "blocks": 21
        }
      },
      "corpus": "/private/tmp/chematic-openbabel-v3000-20261014/corpus-chains.molblocks",
      "records": 64,
      "gate": {
        "source_snapshot_sha256": {
          "baseline": "5add029fe89ad301067f53ee404c87ebb54221733bf4c5280de463eabe5e1550",
          "candidate": "5add029fe89ad301067f53ee404c87ebb54221733bf4c5280de463eabe5e1550"
        },
        "byte_identical_rewrites": true,
        "source_parse_snapshots_equal": true
      }
    },
    "continued": {
      "negative_initial_lanes": {
        "parse": {
          "median": 0.99,
          "ci95": [
            0.987,
            0.995
          ],
          "winning_blocks": 4,
          "blocks": 21
        },
        "write": {
          "median": 0.957,
          "ci95": [
            0.95,
            0.962
          ],
          "winning_blocks": 1,
          "blocks": 21
        },
        "roundtrip": {
          "median": 0.983,
          "ci95": [
            0.977,
            0.993
          ],
          "winning_blocks": 2,
          "blocks": 21
        }
      },
      "corpus": "/private/tmp/chematic-openbabel-lines-20261015/corpus-continued-chains.molblocks",
      "records": 64,
      "gate": {
        "all_original_and_continued_source_snapshots_and_rewrites_byte_identical": true,
        "sha256": {
          "continued-baseline.snapshot": "5add029fe89ad301067f53ee404c87ebb54221733bf4c5280de463eabe5e1550",
          "original-candidate.molblocks": "b67b1fb9cdb4c2e8977a34dd309a13ebc035050ba7e41473114c62f2c7d5ca95",
          "continued-candidate.snapshot": "5add029fe89ad301067f53ee404c87ebb54221733bf4c5280de463eabe5e1550",
          "original-baseline.molblocks": "b67b1fb9cdb4c2e8977a34dd309a13ebc035050ba7e41473114c62f2c7d5ca95",
          "continued-candidate.molblocks": "b67b1fb9cdb4c2e8977a34dd309a13ebc035050ba7e41473114c62f2c7d5ca95",
          "original-baseline.snapshot": "5add029fe89ad301067f53ee404c87ebb54221733bf4c5280de463eabe5e1550",
          "continued-baseline.molblocks": "b67b1fb9cdb4c2e8977a34dd309a13ebc035050ba7e41473114c62f2c7d5ca95",
          "original-candidate.snapshot": "5add029fe89ad301067f53ee404c87ebb54221733bf4c5280de463eabe5e1550"
        },
        "records": 64
      }
    }
  }
}
```

### focused21-chains

| Lane | Baseline / candidate [95% CI]; faster blocks |
|---|---|
| parse | 1.002 [0.997, 1.010]; 12/21 |
| write | 0.968 [0.956, 0.971]; 0/21 |
| roundtrip | 0.992 [0.982, 0.997]; 3/21 |

### focused21-0

| Lane | Baseline / candidate [95% CI]; faster blocks |
|---|---|
| parse | 0.946 [0.942, 0.951]; 0/21 |
| write | 0.968 [0.964, 0.979]; 1/21 |
| roundtrip | 0.994 [0.986, 1.002]; 7/21 |

### focused21-5000

| Lane | Baseline / candidate [95% CI]; faster blocks |
|---|---|
| parse | 1.013 [1.002, 1.018]; 15/21 |
| write | 0.961 [0.957, 0.969]; 0/21 |
| roundtrip | 1.006 [1.002, 1.007]; 15/21 |

### focused21-continued

| Lane | Baseline / candidate [95% CI]; faster blocks |
|---|---|
| parse | 0.999 [0.991, 1.002]; 10/21 |
| write | 0.962 [0.954, 0.967]; 0/21 |
| roundtrip | 0.987 [0.981, 0.993]; 2/21 |

## Limits and evidence

Complete Open Babel-wide speed superiority is not established. Controlled-chain wins apply to that profile; diverse semantic exclusions, other formats, rich query/SGROUP chemistry, cold startup, filesystem/CLI/Python overhead, general 3D, independent corpora, other hosts and published artifacts remain separate evaluations. Earlier candidates’ regressions against older checkpoints are not resolved by this incremental V3000 comparison. No main merge, hosted CI, release or publication was performed.

[Evidence packet](2026-10-15-openbabel-v3000-direct-lines-artifacts/README.md): every original/adaptive block, load observation, mismatch, source and binary/library/input hash, exact output streams, logs, measured source executables, reconstructable patch and complete frozen build/input archive. Archive bytes and both reconstructed source trees are verified. The [existing corpus provenance](../validation/benchmark_corpora/chemistry-speed-followups-v1.manifest.json) and [third-party notices](../THIRD_PARTY_NOTICES.md) apply.
