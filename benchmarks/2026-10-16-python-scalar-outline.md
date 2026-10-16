# Python scalar cache-path outlining — 2026-10-16

The candidate fails its predeclared local adoption gate and remains unpromoted. First-use baseline/candidate median ratios: offset 0: TPSA 1.079x, ASA 1.072x; offset 5000: TPSA 1.043x, ASA 1.047x. Prepared TPSA baseline/candidate median ratios: offset 0 1.005x; offset 5,000 1.013x. Complete superiority over RDKit or COSMolKit is not established.

## Retained source after the fourth experiment

This review branch restores this exact third candidate after the fourth
[Python property cache experiment](2026-10-16-python-scalar-property-cache.md)
fails its prepared ring-count regression gate. The fourth source remains
reconstructable from its separate patch; no timing is pooled, replaced or
repeated after this selection. This third candidate still fails its additional
prepared-TPSA condition and is not approved for main adoption. The
[final source selection](2026-10-16-python-scalar-outline-artifacts/final-source-selection.json)
verifies all retained source hashes and the restoration of the two Python
binding files from the pinned base.

The [packet audit](2026-10-16-python-scalar-outline-artifacts/packet-audit.json) independently reconstructs all four patches, verifies retained source hashes and checks all eight complete matrices against their frozen decisions. The [audit script](2026-10-16-python-scalar-outline-artifacts/verify_packets.py) creates no timing samples.

Later frozen manifests carry the first experiment's legacy source_scope prose.
The file hashes and reconstructable cumulative patches establish their actual
scope; additional_primary records the prepared-TPSA condition for candidates
three and four. Those frozen manifests have not been rewritten.

## Change and preservation

This combined candidate keeps the prior conservative aromatic-edge identity proof and improves the normalization path when a new aromatic view is required. Final bond orders are stored in a dense index-aligned vector instead of a temporary hash map. Before adding each atom to the one final builder, the existing local implicit-H helper evaluates the exact target aromatic flag and incident bond orders. When that count differs from the original, the original H value is stored immediately. This removes the temporary pre-H vector and the second complete graph rebuild, while preserving the same integer valence calculation, bond order, promotion/authoritative semantics, adjacency order and copied metadata.

The third source candidate also keeps the scalar cache-miss sanitizer callback out of line and inlines the TPSA/ASA scalar entry points. The previous compiled TPSA wrapper reserved a 0x6d0-byte working frame plus saved registers before its cached-value check. The new callback boundary keeps that model-building frame on the cache-miss path. The thread-local nested-model guard remains before scalar caching; callback evaluation, outside-lock initialization, cache invalidation and fallback behavior are unchanged. Assembly is diagnostic evidence, not a wall-time prediction.

The previous H-correction rebuild omitted R-group sidecars. The new single builder retains that historical conditional omission, avoiding an unrelated behavior change. Stereo groups with bond memberships, neighbor order, atom tags, ring-closure metadata and direction stashing/anchors are retained exactly. No public API, dependency or unsafe code is added. The existing aromaticity algorithms and model membership sets remain unchanged.

- A frozen rebuild and independent old bond-order map compare 112,576 complete graph/metadata results across all 14,072 pinned inputs, parsed and attempted Kekule forms, both aromaticity algorithms and promote-only/authoritative application. The oracle exercises 2,766 corrected atom counts. An enriched pyrrole case explicitly checks tags, R-group sidecars, enhanced stereo, directions, closure metadata and empty derived caches.
- The edge-proof oracle also matches all 14,072 inputs cold/warm/cloned, 2,048 small graphs, mutations and valid externally seeded flags. All 5,166 positive proofs imply the original component identity verdict (666 and 302 positives in the timed slices); some already satisfy the earlier saturated/forest shortcut, so these are coverage counts rather than newly accelerated call counts.
- Source wheels retain exact bytes for 112,576 descriptor outcomes, 56,288 Morgan outputs and 28,144 MACCS outcomes, with no parse/evaluation exceptions. Six pre-existing NaNs remain in first/repeated Labute ASA on three stress rows (698, 707, 721). These are retained explicitly, not chemical successes or timed wins. The fixed timed slices contain no nonfinite outputs.
- Native unit tests in all three changed crates: 1,627 passed (perception 276, SMILES 393, chemistry 958), zero failed, three existing ignored. Python: 2,302 passed, two skipped. Runner gates: nine passed. Release all-target perception Clippy, workspace formatting and refactor checks pass. These local checks are distinct from hosted CI and package publication.

The standalone [edge proof](2026-10-16-python-aromatic-edge-proof.md) did not establish gains on offset 5,000. The [normalization follow-up](2026-10-16-python-aromatic-normalization.md) established cold gains on both profiles but missed the source median floor for offset-0 prepared TPSA: 0.985 [0.966, 1.023]. That interval is inconclusive about a slowdown. Its separate MACCS first-use negative interval remains recorded. These complete matrices are retained; the third source candidate specifically changes scalar cache-path code generation before new timing. No observations are pooled or substituted.

## Frozen conditions

Baseline: clean `19f07cd198d55c31afc99c2fae437a4b891b109c`, source version 1.1.5. The previous unpromoted ring/MACCS and Open Babel candidates are excluded. The baseline wheel is reused by exact wheel/installed-extension hash from the first experiment, where it was freshly built from the named clean main. The third candidate wheel is freshly built. Both use full LTO, one codegen unit, one build job, `native-inchi` and `CARGO_PROFILE_RELEASE_STRIP=none`. CPython 3.13.6, Rust 1.97.0, macOS arm64; RDKit 2026.09.1 and COSMolKit 0.5.0. All source, corpus, wheel and installed-extension identities are checked; the independent comparator manifest additionally pins all 149 RDKit package native files (63 .so plus 86 bundled .dylib) and one COSMolKit .so before/after each profile; system libraries outside the package directories are excluded.

Offsets 0 and 5,000 of the existing pinned 10k corpus each contribute 1,000 unfiltered inputs. These are exposed development inputs on one shared host, not a sealed holdout. Each profile uses 20 sequential rotating four-arm fresh-process blocks. First use excludes parsing on fresh objects and can compare different eager/lazy preparation. Repeated use warms each operation once per object. Parse-inclusive times real parsing and evaluation together. Native return objects are retained; normalization/output checks occur outside timing. Whole-process high-water RSS includes imports, preparation and every operation. Builds/tests and owned-cache cleanup finish before timing. Other workloads on this shared host are not controlled; raw load averages remain recorded.

Ratios are reference time / candidate time, above 1 favors the candidate. Intervals use 10,000 deterministic paired bootstrap resamples. A strict comparator win requires matching finite outputs, CI lower >1 and 20/20 faster blocks. All observations are retained. COSMolKit ring-count disagreements at offset 0 suppress that comparator's rank in all modes.

The adoption policy was frozen before candidate timing: exact full source outputs; TPSA and ASA first-use CI lower >1 on both profiles; prepared TPSA CI lower >1 on both profiles; every source-relative lane median >=0.99. Every negative interval remains disclosed, including unchanged operations. There are no adaptive timing reruns, substitutions or relaxed gates. See [decision](2026-10-16-python-scalar-outline-artifacts/decision.json).

## Complete matrix: offset 0

| Operation / preparation | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.079 [1.076, 1.085]; 19/20 | 1.622 [1.606, 1.641]; 20/20 **strict win** | 0.577 [0.564, 0.588]; 0/20 |
| tpsa/repeated_use | 1.005 [1.000, 1.021]; 13/20 | 10.113 [10.040, 10.200]; 20/20 **strict win** | 1.776 [1.737, 1.838]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.024 [1.010, 1.029]; 17/20 | 13.314 [13.236, 13.428]; 20/20 **strict win** | 12.454 [12.356, 12.557]; 20/20 **strict win** |
| labute_asa/first_use | 1.072 [1.066, 1.081]; 19/20 | 0.785 [0.766, 0.793]; 0/20 | 0.370 [0.361, 0.379]; 0/20 |
| labute_asa/repeated_use | 1.025 [1.006, 1.034]; 15/20 | 2.746 [2.682, 2.802]; 20/20 **strict win** | 1.487 [1.436, 1.555]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 1.022 [1.014, 1.026]; 17/20 | 12.534 [12.462, 12.634]; 20/20 **strict win** | 11.881 [11.789, 11.930]; 20/20 **strict win** |
| num_rings/first_use | 1.000 [0.996, 1.004]; 10/20 | 0.037 [0.037, 0.040]; 0/20 | unranked: output mismatch |
| num_rings/repeated_use | 1.002 [0.991, 1.008]; 11/20 | 1.743 [1.711, 1.780]; 20/20 **strict win** | unranked: output mismatch |
| num_rings/parse_inclusive | 1.003 [0.989, 1.009]; 12/20 | 11.598 [11.531, 11.672]; 20/20 **strict win** | unranked: output mismatch |
| morgan2_2048/first_use | 1.021 [1.016, 1.025]; 20/20 **strict win** | 3.207 [3.173, 3.233]; 20/20 **strict win** | 2.079 [2.065, 2.087]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 0.996 [0.994, 1.002]; 7/20 | 5.202 [5.160, 5.245]; 20/20 **strict win** | 3.215 [3.201, 3.224]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.010 [1.000, 1.014]; 13/20 | 8.714 [8.657, 8.771]; 20/20 **strict win** | 7.686 [7.643, 7.718]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.010 [1.005, 1.014]; 20/20 **strict win** | 1.279 [1.272, 1.288]; 20/20 **strict win** | 0.839 [0.836, 0.844]; 0/20 |
| morgan2_chiral/repeated_use | 1.000 [0.996, 1.006]; 9/20 | 5.031 [5.014, 5.070]; 20/20 **strict win** | 3.148 [3.131, 3.168]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.005 [0.999, 1.009]; 13/20 | 4.495 [4.471, 4.504]; 20/20 **strict win** | 3.948 [3.938, 3.973]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.023 [1.010, 1.026]; 18/20 | 3.419 [3.403, 3.449]; 20/20 **strict win** | 2.285 [2.266, 2.305]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.998 [0.997, 1.002]; 7/20 | 4.867 [4.845, 4.907]; 20/20 **strict win** | 3.124 [3.110, 3.135]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.008 [1.003, 1.016]; 16/20 | 8.056 [8.009, 8.100]; 20/20 **strict win** | 6.994 [6.961, 7.062]; 20/20 **strict win** |
| maccs/first_use | 1.003 [1.001, 1.004]; 15/20 | 1.804 [1.799, 1.806]; 20/20 **strict win** | 2.068 [2.064, 2.072]; 20/20 **strict win** |
| maccs/repeated_use | 0.998 [0.996, 1.004]; 9/20 | 1.855 [1.853, 1.857]; 20/20 **strict win** | 2.136 [2.132, 2.143]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.002 [0.998, 1.004]; 12/20 | 2.013 [2.008, 2.018]; 20/20 **strict win** | 2.269 [2.264, 2.274]; 20/20 **strict win** |

## Complete matrix: offset 5000

| Operation / preparation | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.043 [1.035, 1.047]; 20/20 **strict win** | 0.463 [0.458, 0.465]; 0/20 | 0.124 [0.122, 0.126]; 0/20 |
| tpsa/repeated_use | 1.013 [1.002, 1.023]; 14/20 | 9.919 [9.840, 9.990]; 20/20 **strict win** | 1.360 [1.312, 1.418]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.035 [1.028, 1.038]; 20/20 **strict win** | 5.350 [5.306, 5.387]; 20/20 **strict win** | 5.502 [5.449, 5.528]; 20/20 **strict win** |
| labute_asa/first_use | 1.047 [1.040, 1.055]; 20/20 **strict win** | 0.214 [0.210, 0.218]; 0/20 | 0.089 [0.086, 0.091]; 0/20 |
| labute_asa/repeated_use | 1.024 [1.005, 1.044]; 15/20 | 2.682 [2.623, 2.754]; 20/20 **strict win** | 1.222 [1.153, 1.266]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 1.039 [1.034, 1.043]; 19/20 | 4.999 [4.952, 5.056]; 20/20 **strict win** | 5.340 [5.290, 5.372]; 20/20 **strict win** |
| num_rings/first_use | 1.018 [1.009, 1.031]; 15/20 | 0.102 [0.098, 0.106]; 0/20 | 0.085 [0.082, 0.090]; 0/20 |
| num_rings/repeated_use | 1.007 [0.999, 1.017]; 13/20 | 1.707 [1.675, 1.742]; 20/20 **strict win** | 1.095 [1.064, 1.119]; 20/20 **strict win** |
| num_rings/parse_inclusive | 1.005 [0.995, 1.016]; 12/20 | 11.806 [11.684, 11.937]; 20/20 **strict win** | 12.855 [12.623, 12.948]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.036 [1.031, 1.041]; 19/20 | 1.610 [1.597, 1.618]; 20/20 **strict win** | 1.029 [1.019, 1.042]; 19/20 |
| morgan2_2048/repeated_use | 1.011 [1.006, 1.016]; 17/20 | 4.785 [4.750, 4.817]; 20/20 **strict win** | 2.890 [2.864, 2.910]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.020 [1.016, 1.027]; 18/20 | 4.738 [4.697, 4.773]; 20/20 **strict win** | 4.530 [4.516, 4.585]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.037 [1.025, 1.043]; 20/20 **strict win** | 1.660 [1.634, 1.667]; 20/20 **strict win** | 1.058 [1.053, 1.075]; 20/20 **strict win** |
| morgan2_chiral/repeated_use | 0.998 [0.994, 1.004]; 9/20 | 4.728 [4.715, 4.799]; 20/20 **strict win** | 2.907 [2.889, 2.928]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.026 [1.013, 1.034]; 18/20 | 4.739 [4.712, 4.785]; 20/20 **strict win** | 4.550 [4.522, 4.582]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.032 [1.026, 1.039]; 19/20 | 1.838 [1.821, 1.857]; 20/20 **strict win** | 1.214 [1.203, 1.219]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.006 [0.998, 1.009]; 13/20 | 4.546 [4.506, 4.616]; 20/20 **strict win** | 2.914 [2.889, 2.947]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.029 [1.019, 1.034]; 18/20 | 4.629 [4.595, 4.647]; 20/20 **strict win** | 4.406 [4.370, 4.434]; 20/20 **strict win** |
| maccs/first_use | 1.001 [0.994, 1.008]; 11/20 | 1.576 [1.568, 1.581]; 20/20 **strict win** | 1.685 [1.675, 1.690]; 20/20 **strict win** |
| maccs/repeated_use | 1.001 [0.999, 1.004]; 13/20 | 1.660 [1.654, 1.666]; 20/20 **strict win** | 1.768 [1.763, 1.776]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.004 [0.996, 1.007]; 11/20 | 1.734 [1.731, 1.742]; 20/20 **strict win** | 1.855 [1.852, 1.861]; 20/20 **strict win** |

## All source-relative intervals wholly below 1

```json
{
  "0": {},
  "5000": {}
}
```

Movement in an unchanged operation does not establish a new algorithm or its cause. Source improvement, semantic equivalence, equivalent-output competitor wins and universal toolkit superiority are separate conclusions.

## Process memory and parser context

| Offset | Baseline MiB | Candidate MiB | RDKit MiB | COSMolKit MiB |
|---|---:|---:|---:|---:|
| 0 | 46.67 | 46.62 | 106.58 | 106.26 |
| 5000 | 40.09 | 39.88 | 88.38 | 83.78 |

Parser-only medians from the same workers (context only; paired operation gates above remain the adoption criteria):

| Offset | Baseline ms | Candidate ms | RDKit ms | COSMolKit ms |
|---|---:|---:|---:|---:|
| 0 | 3.131 | 3.160 | 52.773 | 53.236 |
| 5000 | 1.544 | 1.538 | 23.328 | 27.020 |

## Evidence and reproduction

- Raw complete [offset 0](2026-10-16-python-scalar-outline-artifacts/frontier-0.json.gz) and [offset 5000](2026-10-16-python-scalar-outline-artifacts/frontier-5000.json.gz) matrices include every arm, operation, preparation mode, mismatch, extension/output hash, load average and RSS.
- [Source manifest](2026-10-16-python-scalar-outline-artifacts/source-manifest.json), [predeclared plan](2026-10-16-python-scalar-outline-artifacts/measurement-plan.json), [comparator libraries](2026-10-16-python-scalar-outline-artifacts/comparator-manifest.json), [reconstructable patch](2026-10-16-python-scalar-outline-artifacts/candidate.patch.gz), [verification](2026-10-16-python-scalar-outline-artifacts/verification.json), [exact output hashes](2026-10-16-python-scalar-outline-artifacts/parity.json), [retained nonfinite outcomes](2026-10-16-python-scalar-outline-artifacts/preserved-nonfinite-outcomes.json). The patch was independently reconstructed from the pinned base and every changed-file hash verified.
- Build, test, Clippy, formatting, refactor and timing logs and the exact validation/launch drivers are compressed or saved in this artifact directory. Full source wheels, installed extensions and generated output dumps remain locally in `/private/tmp/chematic-scalar-outline-evidence-20261016` and are not duplicated into Git. Their hashes and generation scripts are retained.
- Earlier fixture/import failures and untimed intermediate tests belong to the separately retained first two candidates. This final source passes the native/Python checks before timing; its exact source is reconstructable from the pinned patch. Completed owned build caches are removed only after preserving wheel/extension identities. No failed-build or failed-test observation enters a performance matrix.

Reproduce the baseline from the named Git commit in an isolated checkout. Decompress/apply `candidate.patch.gz` for the candidate. Build both with `CARGO_BUILD_JOBS=1 CARGO_INCREMENTAL=0 CARGO_PROFILE_RELEASE_STRIP=none maturin build --release --locked --features native-inchi --manifest-path crates/chematic-py/Cargo.toml --interpreter <python3.13> --out <wheel-dir>` using separate targets or preserving each wheel before the next build. Install to separate `pip --target` directories with `--no-deps --no-compile`. Run the three `scripts/verify_python_*_frontier.py` generators separately for each wheel through `PYTHONPATH`; descriptors/Morgan require `--engine chematic`. All use the checked-in pinned 10k and stress corpora.

With the pinned comparator environment, run the unchanged `scripts/bench_python_speed_frontier.py --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi --offset <0-or-5000> --limit 1000 --blocks 20 --baseline-wheel <wheel> --candidate-wheel <wheel> --baseline-path <target> --candidate-path <target> --base-revision 19f07cd198d55c31afc99c2fae437a4b891b109c --source-manifest <manifest.json> --output <result.json>`. The retained launch wrapper only reserves 64 MiB and retries an owned checkpoint write with the same in-memory observation on ENOSPC, leaving workers, timers, output gates and statistics unchanged. Adjust its local paths when reproducing.

The timed TPSA operation is RDKit-default N/O-only TPSA; the broad preservation dump additionally checks S/P-inclusive TPSA. Prepared descriptor/ring deficits, the COSMolKit ring-count mismatch, independent diverse holdout data, additional hosts, 3D, reactions, multithreading, WASM and exact published artifacts remain separate work. This frozen packet does not validate later main commits or published packages.
