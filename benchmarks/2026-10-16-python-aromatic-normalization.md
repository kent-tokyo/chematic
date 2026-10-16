# Python aromatic normalization — 2026-10-16

The candidate fails its predeclared local adoption gate and remains unpromoted. First-use baseline/candidate median ratios: offset 0: TPSA 1.087x, ASA 1.072x; offset 5000: TPSA 1.041x, ASA 1.046x. Complete superiority over RDKit or COSMolKit is not established.

## Change and preservation

This combined candidate keeps the prior conservative aromatic-edge identity proof and improves the normalization path when a new aromatic view is required. Final bond orders are stored in a dense index-aligned vector instead of a temporary hash map. Before adding each atom to the one final builder, the existing local implicit-H helper evaluates the exact target aromatic flag and incident bond orders. When that count differs from the original, the original H value is stored immediately. This removes the temporary pre-H vector and the second complete graph rebuild, while preserving the same integer valence calculation, bond order, promotion/authoritative semantics, adjacency order and copied metadata.

The previous H-correction rebuild omitted R-group sidecars. The new single builder retains that historical conditional omission, avoiding an unrelated behavior change. Stereo groups with bond memberships, neighbor order, atom tags, ring-closure metadata and direction stashing/anchors are retained exactly. No public API, dependency or unsafe code is added. The existing aromaticity algorithms and model membership sets remain unchanged.

- A frozen rebuild and independent old bond-order map compare 112,576 complete graph/metadata results across all 14,072 pinned inputs, parsed and attempted Kekule forms, both aromaticity algorithms and promote-only/authoritative application. The oracle exercises 2,766 corrected atom counts. An enriched pyrrole case explicitly checks tags, R-group sidecars, enhanced stereo, directions, closure metadata and empty derived caches.
- The edge-proof oracle also matches all 14,072 inputs cold/warm/cloned, 2,048 small graphs, mutations and valid externally seeded flags. All 5,166 positive proofs imply the original component identity verdict (666 and 302 positives in the timed slices); some already satisfy the earlier saturated/forest shortcut, so these are coverage counts rather than newly accelerated call counts.
- Source wheels retain exact bytes for 112,576 descriptor outcomes, 56,288 Morgan outputs and 28,144 MACCS outcomes, with no parse/evaluation exceptions. Six pre-existing NaNs remain in first/repeated Labute ASA on three stress rows (698, 707, 721). These are retained explicitly, not chemical successes or timed wins. The fixed timed slices contain no nonfinite outputs.
- Native perception tests: 276 passed, zero failed, one existing ignored. Python: 2,302 passed, two skipped. Runner gates: nine passed. Release all-target perception Clippy, workspace formatting and refactor checks pass. These local checks are distinct from hosted CI and package publication.

The standalone [edge-proof candidate](2026-10-16-python-aromatic-edge-proof.md) did not establish first-use gains on offset 5,000. Its complete matrices and failed adoption decision remain separate. The second candidate adds normalization changes before any new timing; no observations are pooled or substituted.

## Frozen conditions

Baseline: clean `19f07cd198d55c31afc99c2fae437a4b891b109c`, source version 1.1.5. The previous unpromoted ring/MACCS and Open Babel candidates are excluded. The baseline wheel is reused by exact wheel/installed-extension hash from the first experiment, where it was freshly built from the named clean main. The second candidate wheel is freshly built. Both use full LTO, one codegen unit, one build job, `native-inchi` and `CARGO_PROFILE_RELEASE_STRIP=none`. CPython 3.13.6, Rust 1.97.0, macOS arm64; RDKit 2026.09.1 and COSMolKit 0.5.0. All source, corpus, wheel and installed-extension identities are checked; the independent comparator manifest additionally pins all 149 RDKit package native files (63 .so plus 86 bundled .dylib) and one COSMolKit .so before/after each profile; system libraries outside the package directories are excluded.

Offsets 0 and 5,000 of the existing pinned 10k corpus each contribute 1,000 unfiltered inputs. These are exposed development inputs on one shared host, not a sealed holdout. Each profile uses 20 sequential rotating four-arm fresh-process blocks. First use excludes parsing on fresh objects and can compare different eager/lazy preparation. Repeated use warms each operation once per object. Parse-inclusive times real parsing and evaluation together. Native return objects are retained; normalization/output checks occur outside timing. Whole-process high-water RSS includes imports, preparation and every operation. Builds/tests and owned-cache cleanup finish before timing. Other workloads on this shared host are not controlled; raw load averages remain recorded.

Ratios are reference time / candidate time, above 1 favors the candidate. Intervals use 10,000 deterministic paired bootstrap resamples. A strict comparator win requires matching finite outputs, CI lower >1 and 20/20 faster blocks. All observations are retained. COSMolKit ring-count disagreements at offset 0 suppress that comparator's rank in all modes.

The adoption policy was frozen before candidate timing: exact full source outputs; TPSA and ASA first-use CI lower >1 on both profiles; every source-relative lane median >=0.99. Every negative interval remains disclosed, including unchanged operations. There are no adaptive timing reruns, substitutions or relaxed gates. See [decision](2026-10-16-python-aromatic-normalization-artifacts/decision.json).

## Complete matrix: offset 0

| Operation / preparation | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.087 [1.077, 1.093]; 20/20 **strict win** | 1.642 [1.609, 1.661]; 20/20 **strict win** | 0.583 [0.576, 0.600]; 0/20 |
| tpsa/repeated_use | 0.985 [0.966, 1.023]; 8/20 | 9.728 [9.575, 9.917]; 20/20 **strict win** | 1.674 [1.570, 1.741]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.023 [1.019, 1.031]; 17/20 | 13.332 [13.190, 13.439]; 20/20 **strict win** | 12.441 [12.353, 12.501]; 20/20 **strict win** |
| labute_asa/first_use | 1.072 [1.069, 1.078]; 19/20 | 0.773 [0.762, 0.788]; 0/20 | 0.372 [0.362, 0.385]; 0/20 |
| labute_asa/repeated_use | 0.999 [0.980, 1.009]; 10/20 | 2.752 [2.686, 2.799]; 20/20 **strict win** | 1.444 [1.383, 1.523]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 1.024 [1.015, 1.033]; 19/20 | 12.490 [12.441, 12.644]; 20/20 **strict win** | 11.889 [11.807, 11.946]; 20/20 **strict win** |
| num_rings/first_use | 0.993 [0.989, 1.005]; 7/20 | 0.037 [0.036, 0.038]; 0/20 | unranked: output mismatch |
| num_rings/repeated_use | 0.998 [0.975, 1.009]; 9/20 | 1.735 [1.709, 1.761]; 20/20 **strict win** | unranked: output mismatch |
| num_rings/parse_inclusive | 1.006 [0.989, 1.010]; 12/20 | 11.581 [11.470, 11.640]; 20/20 **strict win** | unranked: output mismatch |
| morgan2_2048/first_use | 1.027 [1.019, 1.032]; 19/20 | 3.171 [3.140, 3.200]; 20/20 **strict win** | 2.064 [2.056, 2.080]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 0.998 [0.994, 1.002]; 8/20 | 5.176 [5.146, 5.219]; 20/20 **strict win** | 3.169 [3.154, 3.197]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.007 [0.999, 1.012]; 13/20 | 8.672 [8.636, 8.708]; 20/20 **strict win** | 7.633 [7.566, 7.661]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.011 [1.006, 1.014]; 15/20 | 1.282 [1.273, 1.292]; 20/20 **strict win** | 0.835 [0.827, 0.840]; 0/20 |
| morgan2_chiral/repeated_use | 1.000 [0.995, 1.002]; 11/20 | 5.036 [5.000, 5.066]; 20/20 **strict win** | 3.110 [3.087, 3.127]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.003 [0.997, 1.013]; 13/20 | 4.473 [4.455, 4.498]; 20/20 **strict win** | 3.953 [3.923, 3.956]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.023 [1.011, 1.032]; 18/20 | 3.418 [3.392, 3.440]; 20/20 **strict win** | 2.264 [2.250, 2.277]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.000 [0.995, 1.004]; 10/20 | 4.854 [4.821, 4.906]; 20/20 **strict win** | 3.092 [3.077, 3.105]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.009 [1.005, 1.015]; 20/20 **strict win** | 8.058 [8.024, 8.100]; 20/20 **strict win** | 7.023 [7.000, 7.058]; 20/20 **strict win** |
| maccs/first_use | 0.998 [0.994, 0.999]; 4/20 | 1.798 [1.794, 1.803]; 20/20 **strict win** | 2.064 [2.058, 2.069]; 20/20 **strict win** |
| maccs/repeated_use | 0.996 [0.990, 1.000]; 5/20 | 1.850 [1.845, 1.855]; 20/20 **strict win** | 2.124 [2.120, 2.129]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.996 [0.993, 1.000]; 5/20 | 2.010 [2.003, 2.020]; 20/20 **strict win** | 2.266 [2.259, 2.274]; 20/20 **strict win** |

## Complete matrix: offset 5000

| Operation / preparation | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.041 [1.036, 1.050]; 19/20 | 0.462 [0.454, 0.469]; 0/20 | 0.123 [0.120, 0.128]; 0/20 |
| tpsa/repeated_use | 1.010 [0.978, 1.020]; 12/20 | 9.754 [9.618, 10.021]; 20/20 **strict win** | 1.326 [1.282, 1.351]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.034 [1.026, 1.044]; 19/20 | 5.337 [5.293, 5.389]; 20/20 **strict win** | 5.500 [5.450, 5.553]; 20/20 **strict win** |
| labute_asa/first_use | 1.046 [1.042, 1.054]; 19/20 | 0.216 [0.211, 0.218]; 0/20 | 0.091 [0.088, 0.095]; 0/20 |
| labute_asa/repeated_use | 1.022 [0.996, 1.078]; 13/20 | 2.617 [2.500, 2.652]; 20/20 **strict win** | 1.207 [1.186, 1.240]; 19/20 |
| labute_asa/parse_inclusive | 1.033 [1.024, 1.044]; 18/20 | 4.999 [4.959, 5.028]; 20/20 **strict win** | 5.370 [5.320, 5.414]; 20/20 **strict win** |
| num_rings/first_use | 0.994 [0.986, 1.008]; 8/20 | 0.101 [0.096, 0.105]; 0/20 | 0.085 [0.082, 0.091]; 0/20 |
| num_rings/repeated_use | 1.005 [0.991, 1.016]; 12/20 | 1.697 [1.662, 1.756]; 20/20 **strict win** | 1.123 [1.077, 1.151]; 19/20 |
| num_rings/parse_inclusive | 1.009 [0.990, 1.029]; 12/20 | 11.684 [11.637, 11.817]; 20/20 **strict win** | 12.769 [12.666, 12.876]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.030 [1.027, 1.037]; 18/20 | 1.600 [1.582, 1.610]; 20/20 **strict win** | 1.020 [1.010, 1.027]; 15/20 |
| morgan2_2048/repeated_use | 1.011 [1.004, 1.019]; 17/20 | 4.721 [4.675, 4.783]; 20/20 **strict win** | 2.843 [2.818, 2.888]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.028 [1.020, 1.032]; 18/20 | 4.748 [4.704, 4.779]; 20/20 **strict win** | 4.533 [4.511, 4.562]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.028 [1.021, 1.036]; 17/20 | 1.650 [1.630, 1.665]; 20/20 **strict win** | 1.048 [1.034, 1.064]; 18/20 |
| morgan2_chiral/repeated_use | 0.999 [0.988, 1.006]; 9/20 | 4.742 [4.699, 4.784]; 20/20 **strict win** | 2.861 [2.839, 2.894]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.021 [1.010, 1.035]; 16/20 | 4.729 [4.650, 4.760]; 20/20 **strict win** | 4.531 [4.508, 4.557]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.022 [1.016, 1.028]; 17/20 | 1.837 [1.807, 1.851]; 20/20 **strict win** | 1.200 [1.187, 1.212]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.003 [0.998, 1.010]; 14/20 | 4.508 [4.478, 4.558]; 20/20 **strict win** | 2.877 [2.864, 2.906]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.027 [1.004, 1.035]; 15/20 | 4.599 [4.543, 4.657]; 20/20 **strict win** | 4.365 [4.300, 4.402]; 20/20 **strict win** |
| maccs/first_use | 0.996 [0.989, 1.000]; 7/20 | 1.563 [1.560, 1.569]; 20/20 **strict win** | 1.669 [1.662, 1.681]; 20/20 **strict win** |
| maccs/repeated_use | 0.996 [0.994, 1.001]; 7/20 | 1.649 [1.643, 1.660]; 20/20 **strict win** | 1.766 [1.760, 1.768]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.999 [0.995, 1.006]; 9/20 | 1.733 [1.726, 1.742]; 20/20 **strict win** | 1.856 [1.847, 1.865]; 20/20 **strict win** |

## All source-relative intervals wholly below 1

```json
{
  "0": {
    "maccs/first_use": {
      "median": 0.998,
      "ci95": [
        0.994,
        0.999
      ],
      "winning_blocks": 4,
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
| 0 | 46.69 | 46.63 | 106.57 | 106.24 |
| 5000 | 40.00 | 39.91 | 88.41 | 83.77 |

Parser-only medians from the same workers (context only; paired operation gates above remain the adoption criteria):

| Offset | Baseline ms | Candidate ms | RDKit ms | COSMolKit ms |
|---|---:|---:|---:|---:|
| 0 | 3.146 | 3.165 | 53.057 | 53.409 |
| 5000 | 1.560 | 1.557 | 23.302 | 27.067 |

## Evidence and reproduction

- Raw complete [offset 0](2026-10-16-python-aromatic-normalization-artifacts/frontier-0.json.gz) and [offset 5000](2026-10-16-python-aromatic-normalization-artifacts/frontier-5000.json.gz) matrices include every arm, operation, preparation mode, mismatch, extension/output hash, load average and RSS.
- [Source manifest](2026-10-16-python-aromatic-normalization-artifacts/source-manifest.json), [predeclared plan](2026-10-16-python-aromatic-normalization-artifacts/measurement-plan.json), [comparator libraries](2026-10-16-python-aromatic-normalization-artifacts/comparator-manifest.json), [reconstructable patch](2026-10-16-python-aromatic-normalization-artifacts/candidate.patch.gz), [verification](2026-10-16-python-aromatic-normalization-artifacts/verification.json), [exact output hashes](2026-10-16-python-aromatic-normalization-artifacts/parity.json), [retained nonfinite outcomes](2026-10-16-python-aromatic-normalization-artifacts/preserved-nonfinite-outcomes.json). The patch was independently reconstructed from the pinned base and every changed-file hash verified.
- Build, test, Clippy, formatting, refactor and timing logs and the exact validation/launch drivers are compressed or saved in this artifact directory. Full source wheels, installed extensions and generated output dumps remain locally in `/private/tmp/chematic-hydrogen-patch-evidence-20261016` and are not duplicated into Git. Their hashes and generation scripts are retained.
- An intermediate shared-include import edit failed compilation and was corrected before the final native tests, wheel build and timing. Untimed intermediate tests are retained locally; the final tested source is pinned by the reconstructable patch. Sparse-fixture failures belong to the separately retained first candidate. Completed owned build caches were removed only after preserving wheel/extension identities. No performance observations came from failed builds or tests.

Reproduce the baseline from the named Git commit in an isolated checkout. Decompress/apply `candidate.patch.gz` for the candidate. Build both with `CARGO_BUILD_JOBS=1 CARGO_INCREMENTAL=0 CARGO_PROFILE_RELEASE_STRIP=none maturin build --release --locked --features native-inchi --manifest-path crates/chematic-py/Cargo.toml --interpreter <python3.13> --out <wheel-dir>` using separate targets or preserving each wheel before the next build. Install to separate `pip --target` directories with `--no-deps --no-compile`. Run the three `scripts/verify_python_*_frontier.py` generators separately for each wheel through `PYTHONPATH`; descriptors/Morgan require `--engine chematic`. All use the checked-in pinned 10k and stress corpora.

With the pinned comparator environment, run the unchanged `scripts/bench_python_speed_frontier.py --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi --offset <0-or-5000> --limit 1000 --blocks 20 --baseline-wheel <wheel> --candidate-wheel <wheel> --baseline-path <target> --candidate-path <target> --base-revision 19f07cd198d55c31afc99c2fae437a4b891b109c --source-manifest <manifest.json> --output <result.json>`. The retained launch wrapper only reserves 64 MiB and retries an owned checkpoint write with the same in-memory observation on ENOSPC, leaving workers, timers, output gates and statistics unchanged. Adjust its local paths when reproducing.

Prepared descriptor/ring deficits, the COSMolKit ring-count mismatch, independent diverse holdout data, additional hosts, 3D, reactions, multithreading, WASM and exact published artifacts remain separate work. This frozen packet does not validate later main commits or published packages.
