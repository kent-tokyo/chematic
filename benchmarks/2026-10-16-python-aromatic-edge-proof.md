# Python aromatic edge proof — 2026-10-16

The candidate fails its predeclared local adoption gate and remains unpromoted. First-use baseline/candidate median ratios: offset 0: TPSA 1.095x, ASA 1.076x; offset 5000: TPSA 1.003x, ASA 0.998x. Complete superiority over RDKit or COSMolKit is not established.

This standalone source is followed by a distinct [normalization candidate](2026-10-16-python-aromatic-normalization.md). Its measurements are retained separately; no samples are pooled or substituted.

## Change and preservation

Before constructing cyclic-component scratch tables, a read-only edge scan proves a conservative subset of the existing identity verdict. Every aromatic edge must be cyclic and join flagged atoms. Every other cyclic edge must be single (up/down included), join unflagged atoms and have no double/triple/quadruple bond incident on either endpoint. Those components cannot mix aromatic and non-aromatic cyclic edges or require the existing candidate-cycle check. The earlier model, sanitization, saturated/forest and macrocycle-donor guards remain. When the proof fails, the original component and fallback algorithms run. No new aromaticity rule, public API, dependency or unsafe code is added, and no work is deliberately moved into parsing to improve a getter benchmark.

The ordinary model-disagreement predicate was already memoized; this candidate reduces first-use scratch construction, not repeated model-gate computation. The production change is confined to `crates/chematic-perception/src/rdkit_parity/view.rs`. The patch also adds a frozen-reference test and a MACCS output generator, and applies rustfmt-only changes to the existing Kekulé-view test; test logic and assertions are unchanged.

- The frozen pre-change shortcut matches on all 14,072 pinned inputs, cold/warm/cloned states, 2,048 five-atom graphs, mixed/exocyclic systems, macrocycle donors, valid externally seeded ring flags and mutations. All 5,166 positive edge proofs also imply the original component-based identity verdict (666 and 302 positives in the timed slices). Some positives already satisfy the earlier saturated/forest shortcut; these counts are diagnostic coverage, not counts of newly accelerated calls.
- Source wheels retain exact bytes for 112,576 descriptor outcomes (first/repeated float bits), 56,288 Morgan outputs and 28,144 MACCS outcomes, without parse/evaluation exceptions. Six descriptor outcomes are pre-existing NaNs: first/repeated Labute ASA on three stress rows (698, 707, 721). They are retained explicitly, not counted as chemical successes or timed speed wins. The fixed timed slices contain no nonfinite outputs.
- Native perception tests: 274 passed, zero failed, one existing ignored. Python: 2,302 passed, two skipped. Benchmark-runner gates: nine passed. Release all-target perception Clippy, workspace formatting and refactor-budget checks pass. These are local checks, separate from hosted CI and publication.

## Frozen conditions

Baseline: clean `19f07cd198d55c31afc99c2fae437a4b891b109c`, source version 1.1.5. The previous unpromoted ring/MACCS and Open Babel candidates are excluded. Both Python wheels are fresh release builds with full LTO, one codegen unit, one build job, `native-inchi` and `CARGO_PROFILE_RELEASE_STRIP=none`. CPython 3.13.6, Rust 1.97.0, macOS arm64; RDKit 2026.09.1 and COSMolKit 0.5.0. All source, corpus, wheel and installed-extension identities are checked; the independent comparator manifest additionally pins all 63 RDKit and one COSMolKit package .so files before/after each profile. Bundled .dylib dependencies were not frozen in this first matrix.

Offsets 0 and 5,000 of the existing pinned 10k corpus each contribute 1,000 unfiltered inputs. These are exposed development inputs on one shared host, not a sealed holdout. Each profile uses 20 sequential rotating four-arm fresh-process blocks. First use excludes parsing on fresh objects and can compare different eager/lazy preparation. Repeated use warms each operation once per object. Parse-inclusive times real parsing and evaluation together. Native return objects are retained; normalization/output checks occur outside timing. Whole-process high-water RSS includes imports, preparation and every operation. Builds/tests and owned-cache cleanup finish before timing. Other workloads on this shared host are not controlled; raw load averages remain recorded.

Ratios are reference time / candidate time, above 1 favors the candidate. Intervals use 10,000 deterministic paired bootstrap resamples. A strict comparator win requires matching finite outputs, CI lower >1 and 20/20 faster blocks. All observations are retained. COSMolKit ring-count disagreements at offset 0 suppress that comparator's rank in all modes.

The adoption policy was frozen before candidate timing: exact full source outputs; TPSA and ASA first-use CI lower >1 on both profiles; every source-relative lane median >=0.99. Every negative interval remains disclosed, including unchanged operations. There are no adaptive timing reruns, substitutions or relaxed gates. See [decision](2026-10-16-python-aromatic-edge-proof-artifacts/decision.json).

## Complete matrix: offset 0

| Operation / preparation | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.095 [1.080, 1.121]; 19/20 | 1.646 [1.616, 1.666]; 19/20 | 0.578 [0.569, 0.594]; 0/20 |
| tpsa/repeated_use | 0.999 [0.982, 1.024]; 10/20 | 9.823 [9.674, 10.042]; 20/20 **strict win** | 1.770 [1.635, 1.845]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.024 [1.014, 1.035]; 17/20 | 13.241 [13.180, 13.393]; 20/20 **strict win** | 12.427 [12.362, 12.548]; 20/20 **strict win** |
| labute_asa/first_use | 1.076 [1.063, 1.101]; 19/20 | 0.793 [0.779, 0.804]; 0/20 | 0.385 [0.363, 0.395]; 0/20 |
| labute_asa/repeated_use | 1.004 [0.981, 1.031]; 10/20 | 2.707 [2.639, 2.824]; 20/20 **strict win** | 1.484 [1.408, 1.520]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 1.020 [1.012, 1.037]; 19/20 | 12.538 [12.417, 12.626]; 20/20 **strict win** | 11.939 [11.831, 12.012]; 20/20 **strict win** |
| num_rings/first_use | 1.004 [0.995, 1.011]; 12/20 | 0.039 [0.037, 0.041]; 0/20 | unranked: output mismatch |
| num_rings/repeated_use | 0.999 [0.990, 1.020]; 9/20 | 1.755 [1.723, 1.816]; 20/20 **strict win** | unranked: output mismatch |
| num_rings/parse_inclusive | 1.005 [0.994, 1.014]; 12/20 | 11.593 [11.547, 11.680]; 20/20 **strict win** | unranked: output mismatch |
| morgan2_2048/first_use | 1.026 [1.019, 1.036]; 19/20 | 3.206 [3.178, 3.238]; 20/20 **strict win** | 2.069 [2.042, 2.085]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 0.998 [0.990, 1.006]; 9/20 | 5.201 [5.170, 5.285]; 20/20 **strict win** | 3.186 [3.168, 3.198]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.012 [0.997, 1.020]; 13/20 | 8.684 [8.618, 8.759]; 20/20 **strict win** | 7.678 [7.624, 7.722]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.014 [1.010, 1.020]; 17/20 | 1.281 [1.267, 1.294]; 20/20 **strict win** | 0.842 [0.836, 0.859]; 0/20 |
| morgan2_chiral/repeated_use | 1.000 [0.997, 1.004]; 10/20 | 5.058 [5.021, 5.106]; 20/20 **strict win** | 3.109 [3.088, 3.134]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.012 [1.009, 1.014]; 18/20 | 4.517 [4.502, 4.534]; 20/20 **strict win** | 3.971 [3.957, 4.000]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.019 [1.011, 1.030]; 19/20 | 3.448 [3.424, 3.459]; 20/20 **strict win** | 2.281 [2.276, 2.302]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.000 [0.996, 1.004]; 10/20 | 4.873 [4.848, 4.891]; 20/20 **strict win** | 3.111 [3.083, 3.119]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.017 [1.003, 1.019]; 15/20 | 8.066 [7.985, 8.102]; 20/20 **strict win** | 7.059 [6.978, 7.084]; 20/20 **strict win** |
| maccs/first_use | 1.000 [0.995, 1.004]; 10/20 | 1.804 [1.796, 1.809]; 20/20 **strict win** | 2.066 [2.061, 2.077]; 20/20 **strict win** |
| maccs/repeated_use | 0.998 [0.994, 1.005]; 9/20 | 1.862 [1.849, 1.868]; 20/20 **strict win** | 2.141 [2.124, 2.147]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.001 [0.999, 1.004]; 12/20 | 2.018 [2.009, 2.030]; 20/20 **strict win** | 2.275 [2.263, 2.278]; 20/20 **strict win** |

## Complete matrix: offset 5000

| Operation / preparation | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.003 [0.996, 1.013]; 12/20 | 0.446 [0.436, 0.451]; 0/20 | 0.119 [0.117, 0.122]; 0/20 |
| tpsa/repeated_use | 1.005 [0.995, 1.021]; 12/20 | 9.792 [9.663, 9.899]; 20/20 **strict win** | 1.344 [1.325, 1.364]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.001 [0.996, 1.006]; 11/20 | 5.196 [5.172, 5.235]; 20/20 **strict win** | 5.345 [5.313, 5.377]; 20/20 **strict win** |
| labute_asa/first_use | 0.998 [0.995, 1.005]; 9/20 | 0.206 [0.201, 0.211]; 0/20 | 0.085 [0.083, 0.086]; 0/20 |
| labute_asa/repeated_use | 1.001 [0.988, 1.028]; 10/20 | 2.652 [2.601, 2.695]; 20/20 **strict win** | 1.194 [1.169, 1.207]; 19/20 |
| labute_asa/parse_inclusive | 0.999 [0.991, 1.005]; 9/20 | 4.861 [4.830, 4.875]; 20/20 **strict win** | 5.156 [5.126, 5.226]; 20/20 **strict win** |
| num_rings/first_use | 1.005 [0.996, 1.012]; 11/20 | 0.100 [0.099, 0.105]; 0/20 | 0.084 [0.081, 0.086]; 0/20 |
| num_rings/repeated_use | 1.011 [1.002, 1.018]; 16/20 | 1.719 [1.694, 1.781]; 20/20 **strict win** | 1.135 [1.098, 1.159]; 20/20 **strict win** |
| num_rings/parse_inclusive | 0.996 [0.993, 1.007]; 9/20 | 11.725 [11.619, 11.796]; 20/20 **strict win** | 12.727 [12.646, 12.893]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.003 [1.001, 1.006]; 15/20 | 1.557 [1.546, 1.561]; 20/20 **strict win** | 0.992 [0.988, 1.006]; 8/20 |
| morgan2_2048/repeated_use | 1.005 [0.998, 1.010]; 13/20 | 4.744 [4.680, 4.759]; 20/20 **strict win** | 2.866 [2.854, 2.904]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.002 [0.993, 1.006]; 11/20 | 4.602 [4.570, 4.650]; 20/20 **strict win** | 4.450 [4.410, 4.483]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.007 [1.001, 1.013]; 14/20 | 1.596 [1.585, 1.616]; 20/20 **strict win** | 1.026 [1.021, 1.038]; 20/20 **strict win** |
| morgan2_chiral/repeated_use | 1.004 [0.999, 1.011]; 14/20 | 4.766 [4.730, 4.810]; 20/20 **strict win** | 2.908 [2.889, 2.940]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 0.998 [0.992, 1.000]; 6/20 | 4.626 [4.594, 4.664]; 20/20 **strict win** | 4.439 [4.403, 4.463]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.004 [1.002, 1.006]; 15/20 | 1.798 [1.783, 1.808]; 20/20 **strict win** | 1.174 [1.168, 1.179]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.008 [1.003, 1.014]; 16/20 | 4.571 [4.538, 4.589]; 20/20 **strict win** | 2.901 [2.886, 2.942]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.001 [0.997, 1.006]; 11/20 | 4.546 [4.515, 4.587]; 20/20 **strict win** | 4.287 [4.276, 4.321]; 20/20 **strict win** |
| maccs/first_use | 1.000 [0.997, 1.002]; 9/20 | 1.576 [1.573, 1.579]; 20/20 **strict win** | 1.684 [1.678, 1.687]; 20/20 **strict win** |
| maccs/repeated_use | 0.998 [0.994, 1.002]; 8/20 | 1.656 [1.650, 1.663]; 20/20 **strict win** | 1.768 [1.762, 1.776]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.001 [0.996, 1.005]; 11/20 | 1.736 [1.730, 1.744]; 20/20 **strict win** | 1.861 [1.853, 1.870]; 20/20 **strict win** |

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
| 0 | 46.67 | 46.70 | 106.58 | 106.20 |
| 5000 | 40.09 | 39.99 | 88.43 | 83.75 |

Parser-only medians from the same workers (context only; paired operation gates above remain the adoption criteria):

| Offset | Baseline ms | Candidate ms | RDKit ms | COSMolKit ms |
|---|---:|---:|---:|---:|
| 0 | 3.178 | 3.165 | 53.081 | 53.511 |
| 5000 | 1.548 | 1.541 | 23.269 | 26.861 |

## Evidence and reproduction

- Raw complete [offset 0](2026-10-16-python-aromatic-edge-proof-artifacts/frontier-0.json.gz) and [offset 5000](2026-10-16-python-aromatic-edge-proof-artifacts/frontier-5000.json.gz) matrices include every arm, operation, preparation mode, mismatch, extension/output hash, load average and RSS.
- [Source manifest](2026-10-16-python-aromatic-edge-proof-artifacts/source-manifest.json), [predeclared plan](2026-10-16-python-aromatic-edge-proof-artifacts/measurement-plan.json), [comparator libraries](2026-10-16-python-aromatic-edge-proof-artifacts/comparator-manifest.json), [reconstructable patch](2026-10-16-python-aromatic-edge-proof-artifacts/candidate.patch.gz), [verification](2026-10-16-python-aromatic-edge-proof-artifacts/verification.json), [exact output hashes](2026-10-16-python-aromatic-edge-proof-artifacts/parity.json), [retained nonfinite outcomes](2026-10-16-python-aromatic-edge-proof-artifacts/preserved-nonfinite-outcomes.json). The patch was independently reconstructed from the pinned base and every changed-file hash verified.
- Build, test, Clippy, formatting, refactor and timing logs and the exact validation/launch drivers are compressed or saved in this artifact directory. Full source wheels, installed extensions and generated output dumps remain locally in `/private/tmp/chematic-aromatic-proof-evidence-20261016` and are not duplicated into Git. Their hashes and generation scripts are retained.
- Initial native test compilation and Python collection lacked existing fixtures in the sparse checkout. Their logs remain; unchanged tracked fixtures were unfolded and the complete successful reruns above passed. No performance observations came from those failed launches. Only unused incremental/owned completed build caches were removed, preserving sources, binaries and evidence.

Reproduce the baseline from the named Git commit in an isolated checkout. Decompress/apply `candidate.patch.gz` for the candidate. Build both with `CARGO_BUILD_JOBS=1 CARGO_INCREMENTAL=0 CARGO_PROFILE_RELEASE_STRIP=none maturin build --release --locked --features native-inchi --manifest-path crates/chematic-py/Cargo.toml --interpreter <python3.13> --out <wheel-dir>` using separate targets or preserving each wheel before the next build. Install to separate `pip --target` directories with `--no-deps --no-compile`. Run the three `scripts/verify_python_*_frontier.py` generators separately for each wheel through `PYTHONPATH`; descriptors/Morgan require `--engine chematic`. All use the checked-in pinned 10k and stress corpora.

With the pinned comparator environment, run the unchanged `scripts/bench_python_speed_frontier.py --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi --offset <0-or-5000> --limit 1000 --blocks 20 --baseline-wheel <wheel> --candidate-wheel <wheel> --baseline-path <target> --candidate-path <target> --base-revision 19f07cd198d55c31afc99c2fae437a4b891b109c --source-manifest <manifest.json> --output <result.json>`. The retained launch wrapper only reserves 64 MiB and retries an owned checkpoint write with the same in-memory observation on ENOSPC, leaving workers, timers, output gates and statistics unchanged. Adjust its local paths when reproducing.

Prepared descriptor/ring deficits, the COSMolKit ring-count mismatch, independent diverse holdout data, additional hosts, 3D, reactions, multithreading, WASM and exact published artifacts remain separate work. This frozen packet does not validate later main commits or published packages.
