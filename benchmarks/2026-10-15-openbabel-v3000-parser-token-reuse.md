# V3000 parser token buffer reuse (2026-10-15)

V3000 reading now reuses one token vector within a molecule, clearing it for each logical line instead of allocating it for every atom, bond and block marker. The physical-line validation and logical continuation prepass remain unchanged, and their payloads outlive the borrowed tokens. Token contents and Unicode-whitespace splitting are unchanged. No chemical interpretation, metadata retention, diagnostic handling or error policy changes.


On the two fixed diverse profiles, baseline-relative parser median ratios are 1.072 and 1.149 (16/21 and 19/21 faster blocks). Controlled-chain parsing is 1.303 (20/21). These baseline-relative results do not meet the stricter 21/21 demonstrated-win criterion. In the controlled-chain Open Babel comparison, parse/write/roundtrip ratios are 113.158 / 6.234 / 81.891, with semantic agreement and 21/21 faster blocks for all three lanes. No original baseline-relative interval is wholly below 1, so no adaptive repeat was selected.

## Source and correctness

Fresh local native release builds, source version 1.1.4. Baseline is the frozen attached worktree including the prior unmerged [cycle](2026-10-14-python-cycle-root-pruning.md), [V3000 writer](2026-10-14-openbabel-v3000-writer-buffer.md), and [Kekulé/descriptor-view](2026-10-14-python-kekule-view-reuse.md) candidates in both arms. Last verified Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Git refs were not refreshed after the preceding common-metadata stalls. Source hashes identify this run; it is not a new main merge, release or registry-artifact comparison. Only `crates/chematic-mol/src/mol3000.rs` and the changelog change among the existing frozen files; one test file is added.

- Frozen per-line-allocation reference versus the new complete reader: long/short/empty/continued rows, Unicode whitespace, malformed physical and logical lines, duplicate/unknown properties and invalid/nonfinite coordinates. Full graphs, coordinate bits, metadata, stereo/group fields, 2D/3D diagnostics, rewrite bytes and displayed errors are compared.
- Full MOL tests: 1,013 passed, zero failed/ignored; runner-gate tests: 10 passed. All-target release Clippy for MOL, workspace formatting and refactor checks pass. The first full test run lacked existing fuzz inputs in the isolated copy; that failure is retained, and the complete rerun passed after copying the unchanged inputs.
- Every one of the 2,064 measured records has byte-identical baseline/candidate rewrites and internal parser snapshots. Independent RDKit 2026.09.1 observes indexed atom/bond identity, canonical isomeric graph, finite coordinates within the declared tolerance, names/properties and explicit unknown-double-bond markers. Rejected rows remain errors.

## Protocol

Unchanged, unfiltered fixed profiles: 1,000 records at SMILES offsets 0 and 5,000, and all 64 linear chains of lengths 16, 32, ..., 1,024. The profiles are exposed development data, not sealed holdout data; chains are a separately declared controlled scalability profile, not diverse chemistry. Existing [provenance/licenses](../validation/benchmark_corpora/chemistry-speed-followups-v1.manifest.json) and [third-party notices](../THIRD_PARTY_NOTICES.md) apply. Exact generated MOL blocks are retained with hashes and generation settings.

Rust release/LTO native MOL APIs and installed Open Babel 3.2.1 C++ library/plugin binaries on one macOS arm64 host. Each profile has 21 sequential fresh-process blocks with rotating three-arm order and five full-batch repeats. Startup, loading, delimiter splitting and initial plugin discovery are excluded. Parse retains the diagnostic batch; prepared writing follows an untimed warm write; roundtrip actually parses and writes a new batch within one timer. The unchanged harness warms writing even before the parse timer, so parser results include that allocator/warmup context. Whole-process peak RSS includes preparation, imports/plugins and allocator state. Builds/tests precede timing; profiles are sequential without concurrent compilation or compression. Darwin time -l uses authorized OS-statistics access.

All source, binary, native library and corpus hashes are checked before/after measurement. Paired comparator/candidate medians use 10,000 deterministic bootstrap resamples for 95% intervals. Above 1 means candidate faster. A strict win requires the semantic prerequisite, CI lower endpoint above 1 and 21/21 faster blocks. No input, block, output marker or mismatch is removed or repaired.

## Complete timing matrix

| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |
|---|---|---:|---:|---|---|
| 0 | parse | 333.790 | 304.403 | 1.072 [1.044, 1.113]; 16/21 | unranked |
| 0 | write | 39.928 | 39.106 | 1.038 [1.003, 1.066]; 15/21 | unranked |
| 0 | roundtrip | 370.724 | 333.348 | 1.075 [1.055, 1.096]; 17/21 | unranked |
| 5000 | parse | 87.201 | 75.254 | 1.149 [1.114, 1.170]; 19/21 | unranked |
| 5000 | write | 12.343 | 11.824 | 1.042 [0.995, 1.073]; 13/21 | unranked |
| 5000 | roundtrip | 99.124 | 87.442 | 1.137 [1.118, 1.143]; 19/21 | unranked |
| chains | parse | 47.384 | 36.337 | 1.303 [1.294, 1.313]; 20/21 | 113.158 [112.322, 113.782]; 21/21 |
| chains | write | 16.173 | 15.277 | 1.048 [1.026, 1.066]; 18/21 | 6.234 [6.166, 6.282]; 21/21 |
| chains | roundtrip | 63.834 | 51.542 | 1.233 [1.228, 1.240]; 20/21 | 81.891 [80.560, 82.273]; 21/21 |

Diverse Open Babel lanes with parser/rewrite disagreements remain unranked. Their raw timings do not establish an equivalent-output speed win. Controlled-chain results support only their declared profile. The writer algorithm is unchanged in this candidate; its timings are nevertheless retained and can vary with code layout, preparation and allocator state. This experiment does not isolate those causes.

## Regressions and limits

All original baseline-relative intervals wholly below 1:
- None in this full matrix. Intervals crossing 1 are not improvement claims.

Other formats, rich query/SGROUP chemistry, general 3D, cold startup, filesystem/CLI/Python overhead, streaming memory, independent diverse corpora, other hosts and published artifacts remain separate comparisons. Earlier candidates' unresolved regressions are not settled by this incremental native V3000 test. This is not evidence of complete Open Babel-wide superiority.

## Whole-process peak RSS

| Profile | Lane | Baseline MiB | Candidate MiB | Open Babel MiB |
|---|---|---:|---:|---:|
| 0 | parse | 17.31 | 17.84 | 52.02 |
| 0 | write | 11.30 | 11.28 | 39.00 |
| 0 | roundtrip | 17.95 | 17.36 | 52.12 |
| 5000 | parse | 11.66 | 11.59 | 41.03 |
| 5000 | write | 8.00 | 8.02 | 32.95 |
| 5000 | roundtrip | 11.56 | 11.62 | 41.06 |
| chains | parse | 14.84 | 14.84 | 50.75 |
| chains | write | 10.03 | 10.05 | 39.41 |
| chains | roundtrip | 14.88 | 14.88 | 51.53 |

## Evidence and reproduction

- [Frozen source/compiler/library/corpus manifest](2026-10-15-openbabel-v3000-parser-token-reuse-source.json), [reconstructable patch](https://github.com/kent-tokyo/chematic/blob/0bc130f24a7d163535e0f4ad59eb5a155847cf07/benchmarks/2026-10-15-openbabel-v3000-parser-token-reuse-candidate.patch), [verification/hashes](2026-10-15-openbabel-v3000-parser-token-reuse-verification.json).
- [Profile 0](2026-10-15-openbabel-v3000-parser-token-reuse-0-paired21.json), [profile 5000](2026-10-15-openbabel-v3000-parser-token-reuse-5000-paired21.json), [controlled chains](2026-10-15-openbabel-v3000-parser-token-reuse-chains-paired21.json): every raw block, mismatch, native timing, output hash and RSS.
- [Self-contained tested source and exact input blocks](../validation/results/openbabel-v3000-parser-token-reuse-source-2026-10-15.tar.gz). Extract `candidate/`; reconstruct baseline by restoring the two `baseline-overrides/` files and removing the new token-buffer test file. Verify all baseline hashes. Preceding candidate prerequisites are included.
- Exact rewrite/snapshot streams and logs: `validation/results/openbabel-v3000-parser-token-reuse-*`. Measured native executable paths and hashes are retained; Open Babel libraries/plugins are pinned separately.

Build each source with `cargo build --release --locked -p chematic-mol --example openbabel_v3000_throughput` in a dedicated target, retaining each executable. Use the unchanged C++ harness and pinned Open Babel binaries/plugins. In the pinned observer environment:

```sh
python scripts/bench_openbabel_v3000_throughput.py \
  --baseline <baseline-binary> --candidate <candidate-binary> \
  --openbabel <openbabel-harness> --corpus <exact-profile.molblocks> \
  --records 1000 --blocks 21 --repeats 5 \
  --source-manifest benchmarks/2026-10-15-openbabel-v3000-parser-token-reuse-source.json \
  --output <profile.json>
```

Use `--records 64` for controlled chains. Darwin time -l requires OS-statistics access.
