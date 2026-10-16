# 3D tetrahedral candidate guard — 2026-10-17

Local adoption gate: **PASS**. Open Babel-wide superiority remains unproved. One shared macOS arm64 host, native source APIs.

The measured revision is `a093ce86d8c18eae2c17db3d33405489dc309ace`. Later integration includes main `4a97a599` (writer PR #960 and MOL2 stereo/bond-preservation PR #963). The guard function is byte-identical to the measured version. Combined native validation has 1,356 passes, 0 failures and 0 ignored tests; hosted CI is repeated for the final head. The frozen packet and all timing ratios remain assigned to the measured revision, not the combined head.

## Implementation and preservation

Measured baseline `7d6b092f448913a4ef87f1a3d6353789dad69cb5` and candidate `a093ce86d8c18eae2c17db3d33405489dc309ace`, source version 1.1.6. The sole production change is an early return in `perceive_tetrahedral_from_3d`: it requires at least one untagged, non-aromatic atom with exactly three or four graph neighbors before invoking global stereo/CIP/ring perception. The original output filter can emit a tetrahedral tag only on that shape. If any atom satisfies this necessary condition, the complete original implementation runs unchanged. Both V2000 and V3000 diagnostic readers use this helper.

The guard is a necessary output-shape condition, not a substitute stereocenter classifier. Ring bridgeheads, paired ring centers, sulfur/selenium lone pairs, explicit graph H, phosphorus, existing tetrahedral and square-planar tags remain on their previous paths. Coordinate parsing, dimensional diagnostics, square-planar perception, wedge-conflict correction, ring-alkene correction, and writers retain their implementations. No public API or dependency change. The skipped no-output branch can avoid filling derived graph caches; all prepared-write lanes and whole-process RSS remain measured.

The complete original tetrahedral function is frozen and verified against baseline Git after its name is restored. The candidate suffix after the guard is byte-identical to the original. Differential tests compare atoms, bonds, directional anchors, stereo-neighbor order, tags, R-group labels, and enhanced stereo membership on all 14,072 pinned graphs, every five-vertex simple graph with reflected coordinates, and directed lone-pair/ring/metal/explicit-H cases, existing tags, empty or absent coordinates. These synthetic coordinates test preservation; they are not physical conformers.

All 4,256 profile-input snapshots and rewrite streams match the baseline; chemistry overlaps between profiles. MOL validation: 1,348 passed, 0 failed, 0 ignored across 63 suites. Release all-target MOL Clippy, format, refactor budget and ten runner semantic gates passed. Native correctness uses opt3/LTO=false/16 codegen units; timing uses full LTO/one codegen unit. The 8,333 checked-in fuzz/benchmark fixtures have a separate validation manifest.

## Frozen workload and protocol

Eight unchanged profiles: two unfiltered 1,000-record SDF slices, two unfiltered 1,000-record V3000 slices (offsets 0/5,000), 64 linear alkanes of 16–1,024 atoms, continuation-only copies, and two deterministic nonplanar-coordinate versions of those chains. The last two are the primary absent-candidate controls. They are synthetic, achiral chains, not diverse physical 3D conformers, embedding-quality tests or an independent holdout. Inputs were reused by exact hashes; previous timing samples were not reused. No rows were removed.

Twenty-one rotating three-arm fresh-process blocks per profile, three native batch repetitions for parse, prepared write and parse-plus-write. Loading, startup and plugin discovery are excluded. V3000 framing is outside timing; SDF framing is inside. Parse retains complete diagnostic batches. Prepared write has one untimed warm write. SDF reuses a Rust String; V3000 and Open Babel allocate strings. The identical 3D probe, verified after Rustfmt normalization, calls the real conformer writer and rejects an absent conformer before timing.

Rust 1.97.0, opt3/full LTO/one codegen unit/one job, incremental=false, strip=none. Open Babel 3.2.1 C++ -O3 probes and pinned corpora are reused by hashes. RDKit 2026.09.1 independently observes declared fields. Parse ranking requires equivalent canonical graphs; write/roundtrip ranking requires accepted equivalent rewrites, with indexed atoms/bonds, xyz tolerance 1e-4, names/properties and unknown stereo. Rejected/mismatching rows remain unranked. Source/input, executable and 281 package-local native-file entries are checked before/after profiles. System libraries and package-external runtime data are outside that seal.

Ratios are reference time / candidate time; intervals use the unchanged deterministic paired bootstrap, 10,000 resamples. Adoption requires exact source outputs, both primary parse interval lower bounds >1 and all 24 source-relative medians >=0.99. Every negative interval is disclosed. Strict Open Babel wins additionally require the declared semantic prerequisite and 21/21 faster blocks. Exactly one source-valid epoch per profile; no pooled epochs or sample substitution. An earlier eight-profile epoch is invalid because its candidate retained baseline machine code; all original observations are retained separately and excluded from every source ratio here. The unchanged candidate and unchanged gates were rerun only after binary-source proof corrected that build defect. External host activity is uncontrolled; load averages and whole-process high-water RSS are retained. RSS includes preparation and warmup, not isolated operation or streaming memory.

The first baseline compilation failed with ENOSPC and supplied zero timing samples. The first complete epoch then exposed a shared-target Cargo cache defect: baseline and candidate tetrahedral machine-code bodies were identical despite different source. Source ratios from that epoch are invalid. The candidate-only target was isolated, the entire MOL package was explicitly cleaned, and verbose compilation plus disassembly verified the corrected library and three executables before the one valid epoch per profile. No Rust source, workload, primary, block count or acceptance gate changed. Full invalid observations, code/profile proofs and build logs are retained in `invalid-first-epoch/` and at `/private/tmp/chematic-3d-stereo-guard-evidence-20261017`. Baseline executables are reused by exact hashes, with original-source machine-code proof; old timing samples are not reused. After process/open-file checks, 500,468 KiB of the prior completed experiment's dedicated compiler cache was removed; its exact wheels/extensions, source and evidence were retained and rehashed. Compilation then succeeded with the same source/flags. Native validation needed restoration of checked-in fuzz fixtures in the sparse worktree. Those failures remain separate from completed checks and timings. The checkpoint wrapper retries only an identical in-memory write after ENOSPC, never a worker or timer.

The successful writer optimization in PR #960 and all earlier unpromoted candidates are excluded from this measured base. These timings must not be multiplied by earlier writer ratios or assigned to a later combined revision. Python/WASM and published packages are separate evidence lanes.

## Complete matrix

Median milliseconds per full batch; native elapsed divided by three repetitions. Ratios are paired, not ratios of displayed medians.

| Profile | Mode | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; wins | Open Babel / candidate [95% CI]; wins |
|---|---|---:|---:|---|---|
| sdf-0 | parse | 12.415 | 12.431 | 1.002 [0.998, 1.005]; 13/21 | unranked: semantic prerequisite fails |
| sdf-0 | write | 3.811 | 3.845 | 0.995 [0.985, 1.001]; 7/21 | 6.187 [6.128, 6.231]; 21/21 **strict win** |
| sdf-0 | roundtrip | 16.480 | 16.379 | 1.000 [0.996, 1.011]; 11/21 | 7.150 [7.096, 7.215]; 21/21 **strict win** |
| sdf-5000 | parse | 4.792 | 4.815 | 0.992 [0.975, 1.001]; 7/21 | 8.991 [8.875, 9.098]; 21/21 **strict win** |
| sdf-5000 | write | 1.640 | 1.636 | 1.001 [0.986, 1.010]; 11/21 | 8.724 [8.494, 8.826]; 21/21 **strict win** |
| sdf-5000 | roundtrip | 6.640 | 6.608 | 0.999 [0.987, 1.005]; 10/21 | 8.851 [8.756, 8.995]; 21/21 **strict win** |
| v3000-0 | parse | 13.386 | 13.434 | 0.998 [0.991, 1.007]; 9/21 | unranked: semantic prerequisite fails |
| v3000-0 | write | 2.818 | 2.843 | 0.991 [0.981, 0.996]; 6/21 | unranked: semantic prerequisite fails |
| v3000-0 | roundtrip | 16.528 | 16.549 | 0.998 [0.989, 1.004]; 8/21 | unranked: semantic prerequisite fails |
| v3000-5000 | parse | 5.437 | 5.487 | 1.000 [0.993, 1.008]; 11/21 | unranked: semantic prerequisite fails |
| v3000-5000 | write | 0.743 | 0.738 | 0.999 [0.996, 1.010]; 9/21 | unranked: semantic prerequisite fails |
| v3000-5000 | roundtrip | 6.130 | 6.199 | 1.000 [0.995, 1.006]; 10/21 | unranked: semantic prerequisite fails |
| v3000-chains | parse | 5.687 | 5.677 | 1.005 [0.995, 1.012]; 13/21 | 142.742 [142.490, 144.466]; 21/21 **strict win** |
| v3000-chains | write | 1.257 | 1.254 | 1.007 [1.002, 1.019]; 15/21 | 15.254 [15.124, 15.316]; 21/21 **strict win** |
| v3000-chains | roundtrip | 6.986 | 6.965 | 1.004 [1.000, 1.012]; 14/21 | 119.203 [118.665, 119.611]; 21/21 **strict win** |
| v3000-continued | parse | 7.789 | 7.801 | 0.999 [0.995, 1.001]; 9/21 | 105.630 [105.112, 106.086]; 21/21 **strict win** |
| v3000-continued | write | 1.269 | 1.245 | 1.010 [0.994, 1.038]; 13/21 | unranked: semantic prerequisite fails |
| v3000-continued | roundtrip | 9.060 | 9.086 | 1.000 [0.991, 1.003]; 10/21 | unranked: semantic prerequisite fails |
| v3000-3d-spiral | parse | 45.060 | 6.762 | 6.641 [6.601, 6.695]; 21/21 **strict win** | 119.730 [119.153, 120.103]; 21/21 **strict win** |
| v3000-3d-spiral | write | 13.416 | 13.190 | 1.019 [1.013, 1.024]; 20/21 | 1.604 [1.584, 1.611]; 21/21 **strict win** |
| v3000-3d-spiral | roundtrip | 58.301 | 19.862 | 2.935 [2.915, 2.947]; 21/21 **strict win** | 41.859 [41.645, 41.983]; 21/21 **strict win** |
| v3000-3d-offset | parse | 44.939 | 6.858 | 6.547 [6.518, 6.569]; 21/21 **strict win** | 118.775 [118.013, 118.975]; 21/21 **strict win** |
| v3000-3d-offset | write | 13.577 | 13.339 | 1.019 [1.015, 1.022]; 20/21 | 1.627 [1.615, 1.633]; 21/21 **strict win** |
| v3000-3d-offset | roundtrip | 58.192 | 20.173 | 2.882 [2.872, 2.904]; 21/21 **strict win** | 41.375 [41.256, 41.688]; 21/21 **strict win** |

## Decision and every negative interval/floor failure

```json
{
  "local_adoption_gate_passed": true,
  "all_source_snapshots_and_rewrites_exact": true,
  "primary": {
    "v3000-3d-spiral": {
      "speedup": {
        "median": 6.641,
        "ci95": [
          6.601,
          6.695
        ],
        "winning_blocks": 21,
        "blocks": 21
      },
      "semantic_prerequisite": true,
      "strict_speed_win": true
    },
    "v3000-3d-offset": {
      "speedup": {
        "median": 6.547,
        "ci95": [
          6.518,
          6.569
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
    "v3000-0": {
      "write": {
        "median": 0.991,
        "ci95": [
          0.981,
          0.996
        ],
        "winning_blocks": 6,
        "blocks": 21
      }
    },
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
  "eight_complete_source_valid_epochs_no_sample_pooling": true,
  "invalid_source_contrast_epochs_retained": 8,
  "invalid_epoch_samples_used": 0,
  "binary_source_provenance_valid": true,
  "failed_initial_build_timing_samples": 0,
  "universal_superiority": false
}
```

## Whole-process RSS

| Profile | Mode | Baseline MiB | Candidate MiB | Open Babel MiB |
|---|---|---:|---:|---:|
| sdf-0 | parse | 15.28 | 15.14 | 52.75 |
| sdf-0 | write | 10.73 | 10.61 | 39.67 |
| sdf-0 | roundtrip | 15.28 | 15.22 | 52.88 |
| sdf-5000 | parse | 11.02 | 10.94 | 41.33 |
| sdf-5000 | write | 7.83 | 7.88 | 33.36 |
| sdf-5000 | roundtrip | 11.08 | 10.95 | 41.42 |
| v3000-0 | parse | 16.31 | 16.59 | 51.95 |
| v3000-0 | write | 10.80 | 10.97 | 38.98 |
| v3000-0 | roundtrip | 16.36 | 16.64 | 52.05 |
| v3000-5000 | parse | 10.95 | 11.08 | 40.89 |
| v3000-5000 | write | 7.67 | 7.80 | 32.94 |
| v3000-5000 | roundtrip | 10.97 | 11.12 | 40.97 |
| v3000-chains | parse | 13.86 | 13.97 | 50.75 |
| v3000-chains | write | 9.27 | 9.38 | 39.36 |
| v3000-chains | roundtrip | 13.91 | 14.02 | 51.47 |
| v3000-continued | parse | 14.25 | 14.34 | 51.00 |
| v3000-continued | write | 9.66 | 9.77 | 39.64 |
| v3000-continued | roundtrip | 14.28 | 14.36 | 51.75 |
| v3000-3d-spiral | parse | 16.45 | 16.23 | 50.48 |
| v3000-3d-spiral | write | 10.66 | 10.31 | 39.16 |
| v3000-3d-spiral | roundtrip | 16.50 | 16.28 | 51.27 |
| v3000-3d-offset | parse | 16.52 | 16.33 | 50.59 |
| v3000-3d-offset | write | 10.72 | 10.44 | 39.20 |
| v3000-3d-offset | roundtrip | 16.58 | 16.36 | 51.36 |

## Retained evidence and limits

The [2026-10-17-openbabel-3d-stereo-guard-artifacts](2026-10-17-openbabel-3d-stereo-guard-artifacts/README.md) packet retains every raw matrix, mismatch and negative interval, frozen plan/recipe, source/native/corpus/binary hashes, reconstructable patch and validation/build drivers. Full executables, source snapshots and saved output streams remain at `/private/tmp/chematic-3d-stereo-guard-valid-evidence-20261017` with hashes. The independent audit reconstructs the patch, verifies measured source against Git, reopens saved outputs, checks fresh post-timing native observers, and recomputes all 48 comparisons. Source reconstruction does not imply bit-reproducible rebuilds.

Broader physical 3D conformers, other formats/workloads, large streaming memory, sealed holdouts, additional hosts, cold CLI startup, Python/WASM and exact published packages remain open. Movement in unchanged lanes does not isolate a causal regression but still counts under the frozen retention rule. A finite one-host comparison does not prove toolkit-wide superiority.


See the [main integration validation](2026-10-17-openbabel-3d-stereo-guard-integration/README.md): 1,356 tests and all 4,256 untimed snapshots/rewrites pass.
