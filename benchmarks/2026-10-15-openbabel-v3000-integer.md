# V3000 direct integer writing — 2026-10-15

The production optimization is adopted in the attached source candidate. No toolkit-wide superiority, current-main or release claim follows.

Decision: Adopt direct V3000 2D integer output as the twelfth unmerged source candidate, with an explicit read/write tradeoff. Initial source-relative prepared-write speedups are 1.448x [1.441,1.455] (slice0), 1.880x [1.875,1.886] (slice5000), and 2.064x [2.040,2.092] (chains), all 21 faster. Initial roundtrip speedups 1.028x,1.062x,1.184x also win all 21. Separate predeclared follow-ups retain write gains of 1.441x and 1.875x (both all 21), while read regressions reproduce at 0.985x [0.983,0.991] and 0.980x [0.976,0.987]. These read negatives are not removed or explained away; the writer gain and positive whole parse/write workflow outweigh them for this candidate. The reader source is unchanged, so no new parse algorithm claim is made and read-only throughput remains a limitation. Source-only continuation write 2.033x and roundtrip 1.133x both win all 21. Controlled-chain native Open Babel speedups are parse 126.559x, write 14.138x and roundtrip 108.602x, all with semantic agreement and all 21 faster. Diverse Open Babel lanes remain unranked for graph/stereo/rewrite differences. All 2,064 original and 64 continuation snapshots/rewrites match byte-for-byte between source arms; 872 native tests, ten runner gates, 143,702 integer-oracle fields, Clippy, formatting and unchanged refactor budget pass. Baseline11 executable/source identity is verified; candidate built freshly; no timing samples pooled. This candidate collection completed without a checkpoint interruption. The earlier rejected ASCII experiment, including its interruption and lost observation, remains separately archived and excluded. No complete toolkit-wide superiority, current-main, hosted-CI or publication claim follows.

The V3000 2D writer appends atom/bond ids using the existing private decimal formatter with zero padding width; the existing three-column field table now also supplies its unpadded suffix. Bond-order codes, all currently 0 through 9, are appended as one decimal digit. Existing whitespace, element symbols, atom maps, properties, wedge orientation, CFG fields and floating-point formatting remain unchanged. Values above 999 retain the original general integer path. Three-column V2000 output remains unchanged. The formatter is exposed only within the crate. No external public API, dependency, unsafe code or tolerance is added. The rejected ASCII lexer is excluded and the reader is unchanged; parse timing movement is not attributed to a new reader algorithm.

The immediate baseline contains all eleven previously accepted source candidates, including packed private aromatic-component facts. All rejected experiments remain excluded. The exact baseline11 native release/LTO executable is reused from the preceding rejected ASCII trial, where it was built before either edit. All 10058 baseline source hashes and its binary hash match; baseline-artifact-origin.json records the origin. Candidate is built freshly in a separate target. No timing observation from the ASCII trial is pooled. All builds, validation and cleanup of owned completed targets preceded timing. All 2,064 original records retain byte-identical snapshots and rewritten bytes between source arms; all 64 continued chains match their original records in both arms.

Scoped Rust tests: 872 passed, zero failed/ignored. The full MOL library and six integration suites cover E/Z and wedge stereo, 3D, square planar, V3000 semantic preservation and source-ID lookup. The decimal formatter is checked against the standard integer formatter for 143,702 fields: widths 0 through 12, all values 0 through 1,024, decimal-power boundaries, maximum u32 values, and 10,000 deterministic u32 samples. Existing frozen V3000 writer tests and full raw corpus snapshots/rewrites check unchanged output, including maps, properties, stereo and coordinates. No new molecular input records are added. Ten benchmark runner gates, release all-target Clippy, formatting and the unchanged refactor budget pass.

## Protocol

Fixed unfiltered development profiles: 1,000 records at offsets 0 and 5,000, plus 64 controlled linear chains of lengths 16 through 1,024. One macOS arm64 host; Open Babel 3.2.1 native libraries/plugins and RDKit 2026.09.1 observer pinned by hashes. Each profile has 21 rotating three-arm fresh-process blocks, five batch repeats, and all parse/prepared-write/roundtrip modes. No sealed holdout or toolkit-wide coverage.

Ratios are comparator elapsed time / candidate elapsed time; above 1 is faster. Brackets are paired median-ratio 95% intervals from 10,000 deterministic bootstrap resamples. Strict wins require semantic agreement, CI lower > 1 and all 21 blocks faster. Every original and follow-up observation is retained without pooling or substitution. An external checkpoint/artifact-write wrapper retries ENOSPC using an owned 16MiB reserve, preserving computed reports in memory; source benchmark timing/observer/bootstrap functions and all native timers are unchanged. Driver source and any retry messages are retained.

Timers exclude startup, input loading/framing and plugin discovery. Preparation retains a full parsed batch and includes one untimed warm write; parse evaluates new batches with default diagnostics, prepared-write creates output strings, roundtrip parses and writes fresh records. RSS is whole-process high-water usage including preparation/plugins; it is not per-operation allocation evidence. Source, input, native executable and library hashes are checked.

## Original full matrix

| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |
|---|---|---:|---:|---|---|
| 0 | parse | 162.235 | 163.254 | 0.989 [0.986, 0.996]; 4/21 | unranked (semantic mismatch) |
| 0 | write | 20.490 | 14.109 | 1.448 [1.441, 1.455]; 21/21 **strict win** | unranked (semantic mismatch) |
| 0 | roundtrip | 182.701 | 178.265 | 1.028 [1.025, 1.031]; 21/21 **strict win** | unranked (semantic mismatch) |
| 5000 | parse | 45.296 | 45.451 | 0.997 [0.989, 0.999]; 3/21 | unranked (semantic mismatch) |
| 5000 | write | 7.037 | 3.740 | 1.880 [1.875, 1.886]; 21/21 **strict win** | unranked (semantic mismatch) |
| 5000 | roundtrip | 52.569 | 49.383 | 1.062 [1.055, 1.070]; 21/21 **strict win** | unranked (semantic mismatch) |
| chains | parse | 33.068 | 32.874 | 1.007 [1.000, 1.015]; 15/21 | 126.559 [125.339, 127.807]; 21/21 **strict win** |
| chains | write | 14.226 | 6.862 | 2.064 [2.040, 2.092]; 21/21 **strict win** | 14.138 [14.055, 14.358]; 21/21 **strict win** |
| chains | roundtrip | 47.347 | 39.830 | 1.184 [1.175, 1.199]; 21/21 **strict win** | 108.602 [106.333, 109.167]; 21/21 **strict win** |

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
| parse | 58.630 | 58.206 | 1.006 [1.000, 1.012]; 14/21 |
| write | 15.908 | 7.784 | 2.033 [2.017, 2.062]; 21/21 **strict win** |
| roundtrip | 74.570 | 65.895 | 1.133 [1.125, 1.137]; 21/21 **strict win** |

## Adaptive regression checks

The predeclared rule selects an entire profile for one separate 21-block source-only follow-up in all modes if any initial source-relative CI is wholly below 1. Initial and follow-up outcomes remain separate. Newly appearing negatives do not trigger recursive retries.

```json
{
  "policy": "One separate 21-block two-source-arm follow-up per selected profile, all three modes; retain initial results; no pooling or recursive retries.",
  "profiles": {
    "0": {
      "negative_initial_lanes": {
        "parse": {
          "median": 0.989,
          "ci95": [
            0.986,
            0.996
          ],
          "winning_blocks": 4,
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
        "parse": {
          "median": 0.997,
          "ci95": [
            0.989,
            0.999
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
    }
  }
}
```

### focused21-0

| Lane | Baseline / candidate [95% CI]; faster blocks |
|---|---:|
| parse | 0.985 [0.983, 0.991]; 2/21 |
| write | 1.441 [1.420, 1.448]; 21/21 **strict win** |
| roundtrip | 1.041 [1.037, 1.044]; 20/21 |

### focused21-5000

| Lane | Baseline / candidate [95% CI]; faster blocks |
|---|---:|
| parse | 0.980 [0.976, 0.987]; 1/21 |
| write | 1.875 [1.866, 1.895]; 21/21 **strict win** |
| roundtrip | 1.062 [1.059, 1.073]; 21/21 **strict win** |

## Reproduction and status

[The evidence packet](2026-10-15-openbabel-v3000-integer-artifacts/README.md) retains exact measured binaries, all matrices/snapshots/rewrites, failures if any, test/build logs, full frozen source/input archive and exact patch. Archive readback verifies all 10059 candidate hashes and reconstructs all 10058 baseline hashes. A separate tiny Git repository verifies patch reconstruction. Source reconstruction is not bit-reproducibility of rebuilt binaries.

Last verified shared Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`; refs were not refreshed. These are version 1.1.4 unmerged source candidates, with no hosted CI or publication. Prior unrelated source-candidate regressions against older baselines remain unresolved by this incremental comparison.
