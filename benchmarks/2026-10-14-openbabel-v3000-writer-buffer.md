# V3000 writer buffer reuse (2026-10-14)

V3000 2D/3D writing now appends atom rows, bond rows and ordinary properties directly to the final reserved String. It removes temporary per-row format Strings, heap-allocated CFG fragments and the CFG-filter Vec/join. Coordinate formatting, field order, duplicate-property first-match handling, whitespace normalization, stereo depiction, checked-write errors and unknown metadata behavior are unchanged.

## Source and output contract

This is a frozen local source comparison, version 1.1.4. The attached source worktree was last verified at `cac480fb23cd8adeec6e5e9bed1f18677f87b392` and contains the preceding unmerged relevant-cycle candidate in both arms. Original checkout/common Git reads stalled, so an explicit source copy and file hashes replace fresh Git/main verification. This is not a new merge, release or exact published-package comparison. Only `crates/chematic-mol/src/mol3000/writer.rs` changes among the 1,662 pre-existing snapshot files; two test files are added.

The test-only frozen previous writer verifies exact bytes for rounding ties, negative zero, large/nonfinite coordinates, missing coordinate rows, 2D/3D/checked entry points, long chains, raw metadata, duplicate IDs and historical CFG/whitespace filtering. Full MOL validation: **1,010 passed, zero failed/ignored**, plus **10 runner-gate tests**. Release Clippy (all MOL targets, warnings denied), workspace format and refactor checks pass. Missing copied fixtures were restored unchanged before the complete successful retry.

Every one of the 2,064 measured records has byte-identical baseline/candidate rewrites and internal parser snapshots. Independent RDKit 2026.09.1 checks indexed atoms/bonds, canonical isomeric graph, finite coordinates within the declared tolerance, names/properties and explicit unknown-double-bond markers. Retained error outputs never become passing rows.

## Protocol

The unchanged fixed, unfiltered profiles are 1,000 records at each SMILES offset 0 and 5,000, plus all 64 linear chains of lengths 16, 32, ..., 1,024. Existing chemistry-speed-followups corpus provenance and licenses apply. These are exposed development data; the chains are a controlled scalability profile, not diverse chemistry or a sealed holdout. Exact corpus hashes and generation settings are in the source manifest.

Native Rust release/LTO builds and installed Open Babel 3.2.1 C++ library/plugin binaries run on one macOS arm64 host. Every profile uses 21 fresh-process, sequential, rotating three-arm blocks and five full-batch repeats per process. Startup, file loading, delimiter splitting and initial plugin discovery are excluded. Parsing retains the full diagnostic batch; prepared writing follows a warm write; roundtrip actually parses a new batch and writes it within one timer. The unchanged harness warms a prepared write even before the parse timer, so parser comparisons describe that allocator/warmup context. Whole-process RSS includes all preparation, imports/plugins and allocator state.

All sources, binaries, corpora and native library hashes are checked before/after timing. Darwin time -l required OS-statistics access: the initial restricted attempt is retained as an incomplete report with zero completed blocks, and the full authorized run is reported below. No row or block is filtered.

Ratios are paired comparator/candidate medians, with deterministic 10,000-resample 95% bootstrap intervals. Above 1 means candidate faster. A strict win requires output agreement, CI lower endpoint above 1, and 21/21 faster blocks.

## Complete timing matrix

| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |
|---|---|---:|---:|---|---|
| 0 | parse | 35.590 | 36.097 | 0.988 [0.986, 0.993]; 2/21 | unranked |
| 0 | write | 7.571 | 4.413 | 1.708 [1.694, 1.726]; 21/21 | unranked |
| 0 | roundtrip | 43.804 | 40.730 | 1.077 [1.069, 1.082]; 21/21 | unranked |
| 5000 | parse | 11.272 | 11.450 | 0.983 [0.977, 0.995]; 4/21 | unranked |
| 5000 | write | 3.472 | 1.627 | 2.123 [2.110, 2.130]; 21/21 | unranked |
| 5000 | roundtrip | 14.848 | 13.046 | 1.141 [1.133, 1.147]; 21/21 | unranked |
| chains | parse | 12.008 | 11.852 | 1.011 [0.994, 1.037]; 13/21 | 86.188 [84.952, 89.027]; 21/21 |
| chains | write | 9.619 | 3.935 | 2.426 [2.324, 2.464]; 21/21 | 6.063 [5.900, 6.343]; 21/21 |
| chains | roundtrip | 21.799 | 16.014 | 1.349 [1.328, 1.372]; 21/21 | 66.668 [64.927, 67.383]; 21/21 |

Open Babel retains its prior diverse-profile parser/rewrite mismatches. All six diverse Open Babel lanes remain unranked; their raw speed is not an equivalent-output win. No output markers are repaired and no failing inputs removed. Only the separately declared controlled-chain profile can support direct Open Babel speed comparisons.

## Regressions and limits

- Profile 0, `parse`: baseline/candidate 0.988 [0.986, 0.993]; 2/21.
- Profile 5000, `parse`: baseline/candidate 0.983 [0.977, 0.995]; 4/21.

The full matrix retains these costs even though only writer source changes. It does not identify whether code layout, allocator preparation or host conditions caused them. This candidate is not a uniformly faster implementation or evidence of complete Open Babel superiority. Other formats, rich query/SGROUP chemistry, general 3D, cold startup, filesystem/CLI/Python overhead, streaming memory, fresh diverse corpora, additional hosts and published artifacts remain separate comparisons.

## Whole-process peak RSS

| Profile | Lane | Baseline MiB | Candidate MiB | Open Babel MiB |
|---|---|---:|---:|---:|
| 0 | parse | 17.14 | 17.30 | 51.95 |
| 0 | write | 11.14 | 11.31 | 38.98 |
| 0 | roundtrip | 17.22 | 17.34 | 52.05 |
| 5000 | parse | 11.52 | 11.59 | 40.92 |
| 5000 | write | 7.94 | 8.06 | 32.94 |
| 5000 | roundtrip | 11.50 | 11.58 | 40.95 |
| chains | parse | 14.64 | 14.80 | 50.81 |
| chains | write | 10.08 | 10.17 | 39.39 |
| chains | roundtrip | 14.66 | 14.78 | 51.55 |

## Evidence and reproduction

- [Frozen source/compiler/library/corpus manifest](2026-10-14-openbabel-v3000-writer-buffer-source.json)
- [Reconstructable candidate patch including tests](https://github.com/kent-tokyo/chematic/blob/0bc130f24a7d163535e0f4ad59eb5a155847cf07/benchmarks/2026-10-14-openbabel-v3000-writer-buffer-candidate.patch)
- [Exact rewrite/snapshot archives and validation hashes](2026-10-14-openbabel-v3000-writer-buffer-verification.json)
- [Restricted-environment attempt, zero completed blocks](2026-10-14-openbabel-v3000-writer-buffer-sandbox-attempt.json)
- [Profile 0: all 21 raw blocks](2026-10-14-openbabel-v3000-writer-buffer-0-paired21.json)
- [Profile 5000: all 21 raw blocks](2026-10-14-openbabel-v3000-writer-buffer-5000-paired21.json)
- [Profile chains: all 21 raw blocks](2026-10-14-openbabel-v3000-writer-buffer-chains-paired21.json)

For baseline reconstruction, apply the retained [prior relevant-cycle candidate patch](https://github.com/kent-tokyo/chematic/blob/0bc130f24a7d163535e0f4ad59eb5a155847cf07/benchmarks/2026-10-14-python-cycle-root-pruning-candidate.patch) to the last verified base; its changed-file hashes match this frozen baseline. Then reconstruct the three changed/added writer files from that baseline plus the new candidate patch (reconstruction was checked byte-for-byte). Build `openbabel_v3000_throughput` with `cargo build --release --locked -p chematic-mol --example openbabel_v3000_throughput` in a dedicated target directory; keep each executable before rebuilding. Use the unchanged native C++ harness and pinned Open Babel library/plugin hashes. Run:

```sh
python scripts/bench_openbabel_v3000_throughput.py \
  --baseline <baseline-binary> --candidate <candidate-binary> \
  --openbabel <openbabel-binary> --corpus <exact-profile.molblocks> \
  --records 1000 --blocks 21 --repeats 5 \
  --source-manifest benchmarks/2026-10-14-openbabel-v3000-writer-buffer-source.json \
  --output <profile.json>
```

Use `--records 64` for chains. OS-statistics access is required for Darwin `time -l`. The source manifest names all exact input hashes; existing corpus generation and license records are unchanged.


## Adaptive parser follow-up

The small CI-below-1 parse ratios in the original two diverse profiles prompted a separate two-source-arm follow-up: 21 fresh-process alternating blocks per profile, five repeats, identical binaries, corpora and native warmup. Open Babel is omitted only in this additional incremental comparison. Original three-arm results remain unchanged above; no samples are pooled or discarded. This is the same exposed data, not an independent corpus.

| Profile | Baseline / candidate [95% CI] | Faster blocks |
|---|---|---:|
| 0 | 1.006 [1.002, 1.051] | 17/21 |
| 5000 | 1.021 [0.998, 1.052] | 14/21 |

The small parser slowdowns were not independently reproduced here. The different arm context and host conditions prevent a causal claim; this does not erase the original results or establish uniform improvement. [Every focused-repeat sample](2026-10-14-openbabel-v3000-writer-buffer-parse-repeat21.json) is retained.

Existing [corpus provenance and licensing](../validation/benchmark_corpora/chemistry-speed-followups-v1.manifest.json) and [third-party notices](../THIRD_PARTY_NOTICES.md) apply.
