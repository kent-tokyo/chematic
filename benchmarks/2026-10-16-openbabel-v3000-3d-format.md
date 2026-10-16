# V3000 3D coordinate formatting — 2026-10-16/17 (JST)

Local adoption gate: **PASS**. Open Babel-wide superiority remains unproved. These are native source results on one shared macOS arm64 host.

## Implementation and preservation

Base `faa8092e6050ad6250c68498970245cfbc6579e9`, source version 1.1.6. Only `append_atoms_3d` changes in production: it appends atom/map integers and x/y/z using the existing integer and `push_fixed4` helpers already used by V2000 and 2D V3000. Ties, large magnitudes and nonfinite values still use the exact Rust formatter fallback. Atom properties, missing-coordinate zeros, checked-writer validation, bond writing, stereo, collections, headers, allocation policy, parsing and public APIs are unchanged. No unsafe code or new dependency. Earlier unpromoted candidates are excluded.

The complete original 3D atom function is frozen in a private regression test and independently checked against base tokens after its name is restored. Tests compare whole atom sections, retained prefixes, signed zeros, numeric/map/id boundaries, missing/extra points, pseudoatoms and Unicode/duplicate opaque properties. An additional 20,000 deterministic f64 bit patterns and 60,000 rounding-tie/adjacent-value cases match exactly. Existing full-writer 2D/3D/checked-writer differential, rich metadata and stereo tests remain active. All 4,256 profile-input snapshots and rewrite byte streams match the baseline; chemical structures overlap across profiles.

MOL validation: 1,345 passed, 0 failed, 0 ignored across 63 suites (`--lib --tests`). All-target release Clippy, formatting, refactor budget and ten runner semantic gates passed. Correctness tests use opt3/LTO=false/16 codegen units; timing binaries use full LTO/one codegen unit. Source/build snapshots contain Git-tracked inputs plus the named 3D probe, excluding ignored Python caches.

## Frozen workload and limits

The six existing retention profiles are unchanged: two unfiltered 1,000-record SDF slices, two unfiltered 1,000-record V3000 slices (offsets 0 and 5,000), 64 linear alkanes with 16–1,024 atoms, and continuation-only copies of those chains. The two added primary profiles each reuse the 64 achiral chains with explicit deterministic nonplanar finite coordinates. `controlled-3d-recipe.json` records source/output hashes and formulas. These are synthetic coordinate-serialization controls, not physically optimized conformers, diverse 3D chemistry, embedding quality or an independent holdout. No rows were filtered after measuring outputs.

The new 3D probe calls the real public conformer writer and refuses an absent parsed conformer before timing during warmup. The same probe is built against both source trees. Raw source differs only through Rustfmt normalization (closure blocks/trailing commas); `frozen-writer-proof.json` verifies equality after applying that formatter and retains both raw hashes. Two-dimensional probes and the Open Babel harnesses are unchanged.

Preparation and the plan were frozen on 2026-10-16 JST; formal collection crossed midnight into 2026-10-17 without changing the workload or binaries.

## Protocol and adoption rule

Twenty-one rotating three-arm fresh-process blocks per profile; three native batch repetitions for parse, prepared write and parse-plus-write. Startup/input loading/plugin discovery are excluded; V3000 framing is outside the timer, SDF framing inside. Parse retains the whole diagnostic batch. Prepared write follows one warm write. SDF reuses a Rust String; V3000 and Open Babel return allocating strings. The 3D probe retains all xyz fields. RSS is whole-process high-water memory including preparation and warmup, not isolated operation or streaming memory. External host activity is uncontrolled; load averages are retained.

Rust 1.97.0, opt3/full LTO/one codegen unit/one job, incremental=false, strip=none. Open Babel 3.2.1 C++ -O3 executables are reused by exact hashes; historical timing samples are not reused. RDKit 2026.09.1 independently observes graph, indexed atoms/bonds, xyz within 1e-4, names, properties and explicit-unknown stereo. Internal canonical graphs gate parse rankings; accepted equivalent rewrites gate write/roundtrip. Rejected/mismatching rows stay in place and unranked. The 281 package-local native-file entries are rehashed before/after profiles; aliases may repeat contents. System/transitive libraries and package-external data are outside this seal.

Ratios are reference time / candidate time; intervals use the unchanged deterministic 10,000-resample paired bootstrap, rounded to three decimals. Adoption was fixed before timing: exact source bytes/full errors, both controlled 3D write CI lower>1, every one of 24 source-relative medians >=0.99, and disclosure of all negative intervals. A strict Open Babel win additionally needs matching declared outputs and 21/21 faster blocks. Single-host source results do not establish published-package or toolkit-wide superiority.

An initial startup attempt reached one baseline SDF parse child but `/usr/bin/time -l` failed to obtain `sysctl kern.clockrate` inside the sandbox; no block or valid timing sample was accepted. Its incomplete report and failure logs are retained separately. After obtaining process-statistics permission, exactly one complete 21-block epoch per profile ran with the unchanged frozen sources, binaries, corpus and plan. No performance-informed rerun, adaptive filtering, sample pooling or substitution occurred. A prior integrity-only launch also stopped before workers because Rustfmt completed after a test-file hash was sealed; that test-only manifest was re-sealed before measurement.

## Complete matrix

Median milliseconds per full batch (native elapsed divided by three repetitions). Ratios are paired, not ratios of displayed medians.

| Profile | Mode | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; wins | Open Babel / candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
| sdf-0 | parse | 12.359 | 12.347 | 1.003 [0.997, 1.009]; 13/21 | unranked: semantic prerequisite fails |
| sdf-0 | write | 3.817 | 3.804 | 1.005 [0.993, 1.013]; 11/21 | 6.206 [6.095, 6.256]; 21/21 **strict win** |
| sdf-0 | roundtrip | 16.392 | 16.369 | 0.999 [0.995, 1.007]; 10/21 | 7.176 [7.152, 7.205]; 21/21 **strict win** |
| sdf-5000 | parse | 4.689 | 4.682 | 0.999 [0.993, 1.003]; 10/21 | 8.885 [8.847, 8.910]; 21/21 **strict win** |
| sdf-5000 | write | 1.591 | 1.598 | 0.994 [0.985, 1.004]; 9/21 | 8.547 [8.495, 8.644]; 21/21 **strict win** |
| sdf-5000 | roundtrip | 6.359 | 6.357 | 0.997 [0.995, 1.003]; 10/21 | 8.770 [8.727, 8.825]; 21/21 **strict win** |
| v3000-0 | parse | 13.095 | 13.145 | 0.998 [0.993, 1.004]; 9/21 | unranked: semantic prerequisite fails |
| v3000-0 | write | 2.721 | 2.730 | 1.000 [0.990, 1.006]; 10/21 | unranked: semantic prerequisite fails |
| v3000-0 | roundtrip | 16.037 | 16.036 | 1.002 [0.995, 1.004]; 12/21 | unranked: semantic prerequisite fails |
| v3000-5000 | parse | 5.176 | 5.179 | 1.001 [0.991, 1.012]; 13/21 | unranked: semantic prerequisite fails |
| v3000-5000 | write | 0.702 | 0.704 | 0.999 [0.992, 1.010]; 9/21 | unranked: semantic prerequisite fails |
| v3000-5000 | roundtrip | 5.857 | 5.889 | 0.996 [0.993, 1.001]; 9/21 | unranked: semantic prerequisite fails |
| v3000-chains | parse | 5.688 | 5.656 | 1.004 [0.999, 1.013]; 14/21 | 144.702 [144.320, 145.196]; 21/21 **strict win** |
| v3000-chains | write | 1.217 | 1.207 | 1.013 [0.998, 1.017]; 14/21 | 15.759 [15.425, 15.874]; 21/21 **strict win** |
| v3000-chains | roundtrip | 6.924 | 6.889 | 0.999 [0.996, 1.015]; 10/21 | 121.364 [120.893, 121.851]; 21/21 **strict win** |
| v3000-continued | parse | 7.788 | 7.771 | 1.001 [0.995, 1.008]; 11/21 | 106.241 [105.256, 107.111]; 21/21 **strict win** |
| v3000-continued | write | 1.213 | 1.214 | 1.001 [0.994, 1.007]; 11/21 | unranked: semantic prerequisite fails |
| v3000-continued | roundtrip | 9.059 | 9.018 | 1.004 [1.002, 1.005]; 15/21 | unranked: semantic prerequisite fails |
| v3000-3d-spiral | parse | 50.783 | 49.946 | 1.000 [0.970, 1.007]; 11/21 | 18.113 [17.231, 18.291]; 21/21 **strict win** |
| v3000-3d-spiral | write | 14.913 | 2.767 | 5.330 [5.289, 5.510]; 21/21 **strict win** | 8.556 [8.358, 8.752]; 21/21 **strict win** |
| v3000-3d-spiral | roundtrip | 65.987 | 53.101 | 1.229 [1.197, 1.253]; 19/21 | 17.459 [17.038, 17.873]; 21/21 **strict win** |
| v3000-3d-offset | parse | 60.979 | 59.202 | 1.030 [0.993, 1.114]; 14/21 | 18.506 [17.536, 20.333]; 21/21 **strict win** |
| v3000-3d-offset | write | 18.391 | 3.337 | 5.543 [5.434, 5.915]; 21/21 **strict win** | 8.560 [8.173, 9.167]; 21/21 **strict win** |
| v3000-3d-offset | roundtrip | 80.012 | 62.350 | 1.280 [1.241, 1.361]; 20/21 | 17.595 [16.999, 19.834]; 21/21 **strict win** |

## Decision, every negative interval and every failed floor

```json
{
  "local_adoption_gate_passed": true,
  "all_source_snapshots_and_rewrites_exact": true,
  "primary": {
    "v3000-3d-spiral": {
      "speedup": {
        "median": 5.33,
        "ci95": [
          5.289,
          5.51
        ],
        "winning_blocks": 21,
        "blocks": 21
      },
      "semantic_prerequisite": true,
      "strict_speed_win": true
    },
    "v3000-3d-offset": {
      "speedup": {
        "median": 5.543,
        "ci95": [
          5.434,
          5.915
        ],
        "winning_blocks": 21,
        "blocks": 21
      },
      "semantic_prerequisite": true,
      "strict_speed_win": true
    }
  },
  "primary_ci_pass": true,
  "all_24_source_median_floor": 0.99,
  "below_floor": {
    "sdf-0": {},
    "sdf-5000": {},
    "v3000-0": {},
    "v3000-5000": {},
    "v3000-chains": {},
    "v3000-continued": {},
    "v3000-3d-spiral": {},
    "v3000-3d-offset": {}
  },
  "unranked_source_lanes": [],
  "all_negative_source_intervals": {
    "sdf-0": {},
    "sdf-5000": {},
    "v3000-0": {},
    "v3000-5000": {},
    "v3000-chains": {},
    "v3000-continued": {},
    "v3000-3d-spiral": {},
    "v3000-3d-offset": {}
  },
  "openbabel_strict_wins": {
    "sdf-0": {
      "parse": false,
      "write": true,
      "roundtrip": true
    },
    "sdf-5000": {
      "parse": true,
      "write": true,
      "roundtrip": true
    },
    "v3000-0": {
      "parse": false,
      "write": false,
      "roundtrip": false
    },
    "v3000-5000": {
      "parse": false,
      "write": false,
      "roundtrip": false
    },
    "v3000-chains": {
      "parse": true,
      "write": true,
      "roundtrip": true
    },
    "v3000-continued": {
      "parse": true,
      "write": false,
      "roundtrip": false
    },
    "v3000-3d-spiral": {
      "parse": true,
      "write": true,
      "roundtrip": true
    },
    "v3000-3d-offset": {
      "parse": true,
      "write": true,
      "roundtrip": true
    }
  },
  "eight_complete_single_epochs_no_sample_pooling": true,
  "initial_failed_rss_launch_completed_blocks": 0,
  "universal_superiority": false
}
```

## Whole-process RSS

| Profile | Mode | Baseline MiB | Candidate MiB | Open Babel MiB |
|---|---|---:|---:|---:|
| sdf-0 | parse | 15.14 | 15.17 | 52.73 |
| sdf-0 | write | 10.66 | 10.64 | 39.67 |
| sdf-0 | roundtrip | 15.20 | 15.17 | 52.84 |
| sdf-5000 | parse | 10.91 | 10.89 | 41.33 |
| sdf-5000 | write | 7.73 | 7.73 | 33.33 |
| sdf-5000 | roundtrip | 10.92 | 10.94 | 41.39 |
| v3000-0 | parse | 16.50 | 16.52 | 51.92 |
| v3000-0 | write | 10.91 | 10.92 | 38.98 |
| v3000-0 | roundtrip | 16.59 | 16.59 | 52.03 |
| v3000-5000 | parse | 11.06 | 11.06 | 40.88 |
| v3000-5000 | write | 7.72 | 7.72 | 32.94 |
| v3000-5000 | roundtrip | 11.08 | 11.08 | 40.92 |
| v3000-chains | parse | 13.70 | 13.70 | 50.75 |
| v3000-chains | write | 9.34 | 9.36 | 39.34 |
| v3000-chains | roundtrip | 13.72 | 13.75 | 51.48 |
| v3000-continued | parse | 14.08 | 14.08 | 51.03 |
| v3000-continued | write | 9.72 | 9.70 | 39.70 |
| v3000-continued | roundtrip | 14.11 | 14.11 | 51.80 |
| v3000-3d-spiral | parse | 16.42 | 16.48 | 50.52 |
| v3000-3d-spiral | write | 10.58 | 10.77 | 39.16 |
| v3000-3d-spiral | roundtrip | 16.47 | 16.44 | 51.28 |
| v3000-3d-offset | parse | 16.59 | 16.56 | 50.61 |
| v3000-3d-offset | write | 10.67 | 10.84 | 39.27 |
| v3000-3d-offset | roundtrip | 16.64 | 16.56 | 51.34 |

## Retained proof

The [2026-10-16-openbabel-v3000-3d-format-artifacts](2026-10-16-openbabel-v3000-3d-format-artifacts/README.md) packet preserves every raw matrix, decision, mismatch and negative interval, frozen plan/recipe, source/library/corpus/binary hashes, reconstructable patch and complete logs/drivers. Full binaries, generated profile inputs, parse snapshots and rewrites remain locally at `/private/tmp/chematic-v3000-coordinate-format-evidence-20261016` with hashes. The audit reopens all saved outputs, checks fresh post-timing native dump/snapshot/SMILES observers against the same declared gates, and recomputes all 48 paired comparisons; source reconstruction does not imply bit-reproducible rebuilds.

Pre-timing test setup corrections are retained: the empty graph constructor was corrected and an unnecessary trait import removed. Tiny fixtures absent from sparse checkout were restored unchanged from Git. Only completed dedicated build targets were removed after process/open-file checks; the eight measured executables, all source and evidence were preserved and rehashed. The checkpoint wrapper retries only the same in-memory write after ENOSPC, never a worker or timer.

Open items include diverse/realistic 3D conformers, other formats, larger/streaming memory, additional hosts, sealed holdout chemistry, cold CLI startup, Python/WASM and published artifacts. Unchanged lanes moving in the matrix do not isolate the cause, but still count under the frozen retention rule.

