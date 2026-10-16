# Python scalar property cache — 2026-10-16

The candidate fails its predeclared local adoption gate and remains unpromoted. First-use baseline/candidate median ratios: offset 0: TPSA 1.082x, ASA 1.074x; offset 5000: TPSA 1.038x, ASA 1.041x. Prepared TPSA baseline/candidate median ratios: offset 0 1.022x; offset 5,000 1.019x. Complete superiority over RDKit or COSMolKit is not established.

## Change and preservation

This combined candidate keeps the prior conservative aromatic-edge identity proof and improves the normalization path when a new aromatic view is required. Final bond orders are stored in a dense index-aligned vector instead of a temporary hash map. Before adding each atom to the one final builder, the existing local implicit-H helper evaluates the exact target aromatic flag and incident bond orders. When that count differs from the original, the original H value is stored immediately. This removes the temporary pre-H vector and the second complete graph rebuild, while preserving the same integer valence calculation, bond order, promotion/authoritative semantics, adjacency order and copied metadata.

The fourth delta makes Python TPSA/ASA properties read the existing graph-specific scalar slot directly. On a miss they compute through the unchanged chemistry API; its result populates that same slot. The core scalar initializer computes outside its OnceLock initialization lock, so this nested cold initialization does not deadlock. Python molecule graphs are immutable, and these property reads occur outside the native model callbacks; metadata updates do not change graph-derived values. The library's thread-local nested-model bypass remains unchanged for native consumers. No new cache field, Python return-object reuse, graph mutation or parser precomputation is introduced. The same f64 return still creates a Python float on each read.

Two Python regressions cover cold calls in a fresh process with a timeout, repeated values, metadata changes, pickle restoration and new graphs produced by salt removal. The boundary audit records the immutable Arc storage and absence of Python-side native-model callback entry points.

The combined source also keeps the scalar cache-miss sanitizer callback out of line and inlines the TPSA/ASA scalar entry points. The previous compiled TPSA wrapper reserved a 0x6d0-byte working frame plus saved registers before its cached-value check. The new callback boundary keeps that model-building frame on the cache-miss path. The thread-local nested-model guard remains before scalar caching; callback evaluation, outside-lock initialization, cache invalidation and fallback behavior are unchanged. Assembly is diagnostic evidence, not a wall-time prediction.

The previous H-correction rebuild omitted R-group sidecars. The new single builder retains that historical conditional omission, avoiding an unrelated behavior change. Stereo groups with bond memberships, neighbor order, atom tags, ring-closure metadata and direction stashing/anchors are retained exactly. No public API, dependency or unsafe code is added. The existing aromaticity algorithms and model membership sets remain unchanged.

- A frozen rebuild and independent old bond-order map compare 112,576 complete graph/metadata results across all 14,072 pinned inputs, parsed and attempted Kekule forms, both aromaticity algorithms and promote-only/authoritative application. The oracle exercises 2,766 corrected atom counts. An enriched pyrrole case explicitly checks tags, R-group sidecars, enhanced stereo, directions, closure metadata and empty derived caches.
- The edge-proof oracle also matches all 14,072 inputs cold/warm/cloned, 2,048 small graphs, mutations and valid externally seeded flags. All 5,166 positive proofs imply the original component identity verdict (666 and 302 positives in the timed slices); some already satisfy the earlier saturated/forest shortcut, so these are coverage counts rather than newly accelerated call counts.
- Source wheels retain exact bytes for 112,576 descriptor outcomes, 56,288 Morgan outputs and 28,144 MACCS outcomes, with no parse/evaluation exceptions. Six pre-existing NaNs remain in first/repeated Labute ASA on three stress rows (698, 707, 721). These are retained explicitly, not chemical successes or timed wins. The fixed timed slices contain no nonfinite outputs.
- Unchanged Rust library/test sources retain the third candidate's 1,627 passing native unit checks (perception 276, SMILES 393, chemistry 958; zero failed, three existing ignored), reused only after source hashes are verified. These native tests were not rerun for the Python-only fourth delta. Python: 2,304 passed, two skipped. Runner gates: nine passed. The third candidate's release all-target native Clippy checks remain bound to the unchanged Rust library sources; this candidate separately passes release all-target Python-binding Clippy, workspace formatting and refactor checks. These local checks are distinct from hosted CI and package publication.

The standalone [edge proof](2026-10-16-python-aromatic-edge-proof.md) did not establish gains on offset 5,000. The [normalization follow-up](2026-10-16-python-aromatic-normalization.md) established cold gains on both profiles but missed the source median floor for offset-0 prepared TPSA: 0.985 [0.966, 1.023]. That interval is inconclusive about a slowdown. Its separate MACCS first-use negative interval remains recorded. These complete matrices are retained; no observations are pooled or substituted. The [outlined scalar candidate](2026-10-16-python-scalar-outline.md) passed the original cold-gain/median-floor conditions but not its additional prepared-TPSA condition on offset 0: 1.005 [1.000, 1.021]. The fourth source candidate adds the Python property cache boundary before any new measurement, keeping all adoption conditions unchanged.

## Frozen conditions

Baseline: clean `19f07cd198d55c31afc99c2fae437a4b891b109c`, source version 1.1.5. The previous unpromoted ring/MACCS and Open Babel candidates are excluded. The baseline wheel is reused by exact wheel/installed-extension hash from the first experiment, where it was freshly built from the named clean main. The fourth candidate wheel is freshly built. Both use full LTO, one codegen unit, one build job, `native-inchi` and `CARGO_PROFILE_RELEASE_STRIP=none`. CPython 3.13.6, Rust 1.97.0, macOS arm64; RDKit 2026.09.1 and COSMolKit 0.5.0. All source, corpus, wheel and installed-extension identities are checked; the independent comparator manifest additionally pins all 149 RDKit package native files (63 .so plus 86 bundled .dylib) and one COSMolKit .so before/after each profile; system libraries outside the package directories are excluded.

Offsets 0 and 5,000 of the existing pinned 10k corpus each contribute 1,000 unfiltered inputs. These are exposed development inputs on one shared host, not a sealed holdout. Each profile uses 20 sequential rotating four-arm fresh-process blocks. First use excludes parsing on fresh objects and can compare different eager/lazy preparation. Repeated use warms each operation once per object. Parse-inclusive times real parsing and evaluation together. Native return objects are retained; normalization/output checks occur outside timing. Whole-process high-water RSS includes imports, preparation and every operation. Builds/tests and owned-cache cleanup finish before timing. Other workloads on this shared host are not controlled; raw load averages remain recorded.

Ratios are reference time / candidate time, above 1 favors the candidate. Intervals use 10,000 deterministic paired bootstrap resamples. A strict comparator win requires matching finite outputs, CI lower >1 and 20/20 faster blocks. All observations are retained. COSMolKit ring-count disagreements at offset 0 suppress that comparator's rank in all modes.

The adoption policy was frozen before candidate timing: exact full source outputs; TPSA and ASA first-use CI lower >1 on both profiles; prepared TPSA CI lower >1 on both profiles; every source-relative lane median >=0.99. Every negative interval remains disclosed, including unchanged operations. There are no adaptive timing reruns, substitutions or relaxed gates. See [decision](2026-10-16-python-scalar-property-cache-artifacts/decision.json).

## Complete matrix: offset 0

| Operation / preparation | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.082 [1.077, 1.092]; 19/20 | 1.606 [1.586, 1.642]; 20/20 **strict win** | 0.582 [0.577, 0.591]; 0/20 |
| tpsa/repeated_use | 1.022 [1.002, 1.033]; 14/20 | 9.885 [9.771, 10.205]; 20/20 **strict win** | 1.759 [1.699, 1.831]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.026 [1.019, 1.030]; 20/20 **strict win** | 13.295 [13.203, 13.411]; 20/20 **strict win** | 12.523 [12.413, 12.632]; 20/20 **strict win** |
| labute_asa/first_use | 1.074 [1.073, 1.082]; 19/20 | 0.786 [0.772, 0.798]; 0/20 | 0.378 [0.368, 0.389]; 0/20 |
| labute_asa/repeated_use | 1.000 [0.977, 1.027]; 10/20 | 2.740 [2.681, 2.815]; 20/20 **strict win** | 1.471 [1.436, 1.511]; 19/20 |
| labute_asa/parse_inclusive | 1.014 [1.004, 1.027]; 14/20 | 12.516 [12.289, 12.584]; 20/20 **strict win** | 11.814 [11.689, 11.937]; 20/20 **strict win** |
| num_rings/first_use | 1.005 [0.993, 1.014]; 12/20 | 0.038 [0.037, 0.039]; 0/20 | unranked: output mismatch |
| num_rings/repeated_use | 0.984 [0.972, 0.994]; 5/20 | 1.746 [1.715, 1.854]; 20/20 **strict win** | unranked: output mismatch |
| num_rings/parse_inclusive | 1.000 [0.991, 1.007]; 10/20 | 11.557 [11.366, 11.658]; 20/20 **strict win** | unranked: output mismatch |
| morgan2_2048/first_use | 1.021 [1.017, 1.027]; 17/20 | 3.183 [3.156, 3.216]; 20/20 **strict win** | 2.071 [2.054, 2.081]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 0.998 [0.988, 1.007]; 9/20 | 5.162 [5.137, 5.195]; 20/20 **strict win** | 3.206 [3.166, 3.218]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.007 [1.000, 1.018]; 14/20 | 8.681 [8.623, 8.805]; 20/20 **strict win** | 7.660 [7.613, 7.756]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.007 [1.003, 1.011]; 15/20 | 1.273 [1.270, 1.285]; 20/20 **strict win** | 0.836 [0.833, 0.841]; 0/20 |
| morgan2_chiral/repeated_use | 0.997 [0.991, 1.000]; 5/20 | 5.013 [4.978, 5.050]; 20/20 **strict win** | 3.108 [3.085, 3.144]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.001 [0.996, 1.011]; 11/20 | 4.468 [4.432, 4.524]; 20/20 **strict win** | 3.956 [3.917, 3.981]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.014 [1.010, 1.021]; 18/20 | 3.410 [3.383, 3.425]; 20/20 **strict win** | 2.268 [2.252, 2.299]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.999 [0.998, 1.003]; 8/20 | 4.833 [4.818, 4.896]; 20/20 **strict win** | 3.114 [3.099, 3.149]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.007 [1.001, 1.017]; 14/20 | 8.036 [7.997, 8.098]; 20/20 **strict win** | 7.058 [7.017, 7.079]; 20/20 **strict win** |
| maccs/first_use | 1.001 [0.996, 1.004]; 11/20 | 1.802 [1.797, 1.808]; 20/20 **strict win** | 2.068 [2.059, 2.077]; 20/20 **strict win** |
| maccs/repeated_use | 1.000 [0.996, 1.003]; 12/20 | 1.854 [1.848, 1.860]; 20/20 **strict win** | 2.131 [2.126, 2.138]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.999 [0.998, 1.005]; 8/20 | 2.019 [2.016, 2.024]; 20/20 **strict win** | 2.274 [2.271, 2.279]; 20/20 **strict win** |

## Complete matrix: offset 5000

| Operation / preparation | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.038 [1.030, 1.044]; 19/20 | 0.460 [0.451, 0.468]; 0/20 | 0.124 [0.122, 0.126]; 0/20 |
| tpsa/repeated_use | 1.019 [1.008, 1.035]; 15/20 | 10.155 [9.850, 10.231]; 20/20 **strict win** | 1.398 [1.331, 1.464]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.028 [1.021, 1.035]; 19/20 | 5.354 [5.301, 5.376]; 20/20 **strict win** | 5.513 [5.461, 5.528]; 20/20 **strict win** |
| labute_asa/first_use | 1.041 [1.035, 1.052]; 19/20 | 0.218 [0.214, 0.224]; 0/20 | 0.093 [0.088, 0.094]; 0/20 |
| labute_asa/repeated_use | 1.009 [0.994, 1.037]; 12/20 | 2.719 [2.644, 2.838]; 20/20 **strict win** | 1.218 [1.190, 1.231]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 1.033 [1.021, 1.041]; 18/20 | 5.036 [4.983, 5.063]; 20/20 **strict win** | 5.324 [5.265, 5.376]; 20/20 **strict win** |
| num_rings/first_use | 1.011 [1.001, 1.018]; 14/20 | 0.103 [0.102, 0.106]; 0/20 | 0.088 [0.086, 0.092]; 0/20 |
| num_rings/repeated_use | 1.005 [0.986, 1.020]; 12/20 | 1.736 [1.704, 1.804]; 20/20 **strict win** | 1.139 [1.107, 1.163]; 20/20 **strict win** |
| num_rings/parse_inclusive | 1.015 [1.006, 1.023]; 15/20 | 11.832 [11.723, 12.009]; 20/20 **strict win** | 12.862 [12.756, 12.893]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.034 [1.029, 1.040]; 20/20 **strict win** | 1.615 [1.598, 1.625]; 20/20 **strict win** | 1.030 [1.023, 1.042]; 18/20 |
| morgan2_2048/repeated_use | 1.014 [1.005, 1.019]; 19/20 | 4.735 [4.701, 4.756]; 20/20 **strict win** | 2.877 [2.862, 2.913]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.031 [1.027, 1.035]; 19/20 | 4.753 [4.716, 4.778]; 20/20 **strict win** | 4.557 [4.528, 4.584]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.035 [1.027, 1.043]; 20/20 **strict win** | 1.646 [1.635, 1.657]; 20/20 **strict win** | 1.070 [1.058, 1.093]; 20/20 **strict win** |
| morgan2_chiral/repeated_use | 1.004 [0.996, 1.007]; 11/20 | 4.729 [4.711, 4.756]; 20/20 **strict win** | 2.917 [2.896, 2.946]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.027 [1.015, 1.032]; 18/20 | 4.717 [4.686, 4.742]; 20/20 **strict win** | 4.561 [4.530, 4.612]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.028 [1.020, 1.035]; 19/20 | 1.838 [1.824, 1.851]; 20/20 **strict win** | 1.217 [1.200, 1.224]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.005 [0.999, 1.011]; 14/20 | 4.554 [4.523, 4.602]; 20/20 **strict win** | 2.888 [2.868, 2.915]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.020 [1.008, 1.032]; 18/20 | 4.621 [4.581, 4.640]; 20/20 **strict win** | 4.390 [4.345, 4.414]; 20/20 **strict win** |
| maccs/first_use | 1.000 [0.996, 1.004]; 10/20 | 1.577 [1.574, 1.580]; 20/20 **strict win** | 1.683 [1.681, 1.693]; 20/20 **strict win** |
| maccs/repeated_use | 1.001 [0.998, 1.005]; 12/20 | 1.654 [1.649, 1.659]; 20/20 **strict win** | 1.771 [1.764, 1.779]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.001 [0.999, 1.004]; 10/20 | 1.738 [1.729, 1.744]; 20/20 **strict win** | 1.863 [1.860, 1.868]; 20/20 **strict win** |

## All source-relative intervals wholly below 1

```json
{
  "0": {
    "num_rings/repeated_use": {
      "median": 0.984,
      "ci95": [
        0.972,
        0.994
      ],
      "winning_blocks": 5,
      "blocks": 20
    }
  },
  "5000": {}
}
```

Movement in an unchanged operation does not establish a new algorithm or its cause. Source improvement, semantic equivalence, equivalent-output competitor wins and universal toolkit superiority are separate conclusions.

## Process memory and parser context

| Offset | Baseline MiB | Candidate MiB | RDKit MiB | COSMolKit MiB |
|---|---:|---:|---:|---:|
| 0 | 46.70 | 46.66 | 106.62 | 106.23 |
| 5000 | 40.02 | 39.93 | 88.36 | 83.73 |

Parser-only medians from the same workers (context only; paired operation gates above remain the adoption criteria):

| Offset | Baseline ms | Candidate ms | RDKit ms | COSMolKit ms |
|---|---:|---:|---:|---:|
| 0 | 3.150 | 3.165 | 52.933 | 53.368 |
| 5000 | 1.548 | 1.536 | 23.293 | 27.057 |

## Evidence and reproduction

- Raw complete [offset 0](2026-10-16-python-scalar-property-cache-artifacts/frontier-0.json.gz) and [offset 5000](2026-10-16-python-scalar-property-cache-artifacts/frontier-5000.json.gz) matrices include every arm, operation, preparation mode, mismatch, extension/output hash, load average and RSS.
- [Source manifest](2026-10-16-python-scalar-property-cache-artifacts/source-manifest.json), [predeclared plan](2026-10-16-python-scalar-property-cache-artifacts/measurement-plan.json), [comparator libraries](2026-10-16-python-scalar-property-cache-artifacts/comparator-manifest.json), [reconstructable patch](2026-10-16-python-scalar-property-cache-artifacts/candidate.patch.gz), [verification](2026-10-16-python-scalar-property-cache-artifacts/verification.json), [exact output hashes](2026-10-16-python-scalar-property-cache-artifacts/parity.json), [retained nonfinite outcomes](2026-10-16-python-scalar-property-cache-artifacts/preserved-nonfinite-outcomes.json). The patch was independently reconstructed from the pinned base and every changed-file hash verified.
- Build, test, Clippy, formatting, refactor and timing logs and the exact validation/launch drivers are compressed or saved in this artifact directory. Full source wheels, installed extensions and generated output dumps remain locally in `/private/tmp/chematic-property-cache-evidence-20261016` and are not duplicated into Git. Their hashes and generation scripts are retained.
- Earlier fixture/import failures and untimed intermediate tests belong to the separately retained first two candidates. The unchanged native source checks are reused by hash; the new Python binding is rebuilt, checked and tested before timing. The exact final source is reconstructable from the pinned patch. Completed owned build caches are removed only after preserving wheel/extension identities. No failed-build or failed-test observation enters a performance matrix.

Reproduce the baseline from the named Git commit in an isolated checkout. Decompress/apply `candidate.patch.gz` for the candidate. Build both with `CARGO_BUILD_JOBS=1 CARGO_INCREMENTAL=0 CARGO_PROFILE_RELEASE_STRIP=none maturin build --release --locked --features native-inchi --manifest-path crates/chematic-py/Cargo.toml --interpreter <python3.13> --out <wheel-dir>` using separate targets or preserving each wheel before the next build. Install to separate `pip --target` directories with `--no-deps --no-compile`. Run the three `scripts/verify_python_*_frontier.py` generators separately for each wheel through `PYTHONPATH`; descriptors/Morgan require `--engine chematic`. All use the checked-in pinned 10k and stress corpora.

With the pinned comparator environment, run the unchanged `scripts/bench_python_speed_frontier.py --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi --offset <0-or-5000> --limit 1000 --blocks 20 --baseline-wheel <wheel> --candidate-wheel <wheel> --baseline-path <target> --candidate-path <target> --base-revision 19f07cd198d55c31afc99c2fae437a4b891b109c --source-manifest <manifest.json> --output <result.json>`. The retained launch wrapper only reserves 64 MiB and retries an owned checkpoint write with the same in-memory observation on ENOSPC, leaving workers, timers, output gates and statistics unchanged. Adjust its local paths when reproducing.

The timed TPSA operation is RDKit-default N/O-only TPSA; the broad preservation dump additionally checks S/P-inclusive TPSA. Prepared descriptor/ring deficits, the COSMolKit ring-count mismatch, independent diverse holdout data, additional hosts, 3D, reactions, multithreading, WASM and exact published artifacts remain separate work. This frozen packet does not validate later main commits or published packages.
