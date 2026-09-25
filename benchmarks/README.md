# chematic benchmark records

This directory contains dated, reproducible measurement records. Numbers are
scoped to the source revision, package versions, corpus, hardware, runtime,
and operation boundary written in each record. They are not universal speed,
accuracy, or compatibility claims.

## Start here

| Need | Start with |
|---|---|
| Review the post-v1.0.23 output-identical speed branch | [`2026-09-24-perf-speed3-output-identical.md`](2026-09-24-perf-speed3-output-identical.md) |
| Review the M4 source-candidate SSSR/cache/VF2 browser rerun | [`2026-09-24-rdkitjs-m4-sssr-cache.md`](2026-09-24-rdkitjs-m4-sssr-cache.md) |
| Understand the rules and how to report a result | [`docs/benchmark.md`](../docs/benchmark.md) |
| Review current public-package fingerprint and 3D evidence | [`2026-09-23-public-package-fingerprint-3d-v1.0.20.md`](2026-09-23-public-package-fingerprint-3d-v1.0.20.md) |
| Review the post-freeze MMFF94 stereo-safe source speed gate | [`2026-09-23-mmff94-stereo-safe-performance.md`](2026-09-23-mmff94-stereo-safe-performance.md) |
| Review current-source MMFF94 same-coordinate energy evidence | [`2026-09-23-mmff94-current-source-energy.md`](2026-09-23-mmff94-current-source-energy.md) |
| Review the analytic MMFF94 source speed/convergence candidate | [`2026-09-23-a6-analytic-mmff94-source-candidate.md`](2026-09-23-a6-analytic-mmff94-source-candidate.md) |
| Review the MMFF94 stereo-safe source quality candidate | [`2026-09-23-a6-mmff94-stereo-safe-quality.md`](2026-09-23-a6-mmff94-stereo-safe-quality.md) |
| Review current Parse + compatible Morgan evidence | [`2026-09-20-parse-morgan-rdkitjs.md`](2026-09-20-parse-morgan-rdkitjs.md) |
| Compare current similarity search with RDKit | [`2026-09-12-similarity-search-a3-v1.0.13.md`](2026-09-12-similarity-search-a3-v1.0.13.md) |
| Reproduce the 1.10x hot-path gate | [`2026-09-05-hotpath-110.md`](2026-09-05-hotpath-110.md) |
| Check file-streaming contracts | [v1.0.12 validation matrix](../validation/results/cross-engine-matrix-v1.0.12.json) |
| Check current streaming safety gate | [`2026-09-11-streaming-safety-v1.0.13.json`](2026-09-11-streaming-safety-v1.0.13.json) |
| Check current isolated official RDKit.js browser comparison | [`2026-09-16-official-rdkit-js-isolated-browser-v1.0.15.md`](2026-09-16-official-rdkit-js-isolated-browser-v1.0.15.md) |
| Historical v1.0.12 official RDKit.js Node rerun | [`../validation/results/competitive-benchmark-rdkitjs-2026-09-11-v1.0.12.json`](../validation/results/competitive-benchmark-rdkitjs-2026-09-11-v1.0.12.json) |
| Historical v1.0.12 official RDKit.js Node report | [`2026-09-11-official-rdkit-js-v1.0.12.md`](2026-09-11-official-rdkit-js-v1.0.12.md) |
| Check official RDKit.js browser gate | [`2026-09-10-official-rdkit-js-browser-gate-v1.0.11.md`](2026-09-10-official-rdkit-js-browser-gate-v1.0.11.md) |
| Check WASM artifact size | [`2026-09-09-wasm-size-v1.0.10.md`](2026-09-09-wasm-size-v1.0.10.md) |
| Find older measurements | [Historical snapshots](#historical-snapshots) |

The current release line is v1.0.20. Older records remain versioned historical
measurements where their headers say so; a release does not imply that an older
measurement was rerun.

## Performance and scaling

### Search and fingerprint paths

| Record | Scope |
|---|---|
| [`2026-09-24-rdkitjs-m4-sssr-cache.md`](2026-09-24-rdkitjs-m4-sssr-cache.md) | Apple-silicon Chrome source-candidate rerun: exact packed-output contract across 5,000 rows; parse-inclusive compatible Morgan 4.535x geometric speedup (95% lower bound 4.411x) and prepared lane 3.411x (lower bound 3.321x); not a registry-package or native-Python-matrix claim |
| [`2026-09-23-mmff94-stereo-safe-performance.md`](2026-09-23-mmff94-stereo-safe-performance.md) | Source `e9f178fa`: fixed-265 and post-freeze 100-row MMFF94 stereo-safe quality/speed gates; both pass, but publication and remaining A6 numerical exits stay separate |
| [`2026-09-23-public-package-fingerprint-3d-v1.0.20.md`](2026-09-23-public-package-fingerprint-3d-v1.0.20.md) | Registry-installed v1.0.20: exact supported-domain compatible Morgan speed wins, UFF best-of-10 speed/coverage win, and 265/265 clash-free MMFF94 stereo-safe quality with a remaining speed deficit |
| [`2026-09-22-public-package-fingerprint-3d.md`](2026-09-22-public-package-fingerprint-3d.md) | Historical published v1.0.19 baseline and source-candidate evidence |
| [`2026-09-23-a6-analytic-mmff94-source-candidate.md`](2026-09-23-a6-analytic-mmff94-source-candidate.md) | Source `af7c0c44` analytic MMFF94/300-iteration packet against the published RDKit 2026.3.6 wheel; speed and quality boundaries remain separate |
| [`2026-09-23-a6-mmff94-stereo-safe-quality.md`](2026-09-23-a6-mmff94-stereo-safe-quality.md) | Historical source `dd7fe3e9` result that first closed the 12-failure/four-clash MMFF94 cohort; the public-package confirmation is linked above |
| [`2026-09-20-parse-morgan-rdkitjs.md`](2026-09-20-parse-morgan-rdkitjs.md) | Merged source candidate, equivalent-output Parse + compatible Morgan comparison with official RDKit.js: exact supported-domain bits, local four-lane records, and three-browser GitHub-hosted paired-CI gate; not a published-package claim |
| [`2026-09-09-similarity-search-v1.0.9.md`](2026-09-09-similarity-search-v1.0.9.md) | 4,500-entry library / 500-query exact top-k comparison with RDKit; latency and ranking overlap are separate axes |
| [`2026-09-09-similarity-search-v1.0.9.json`](2026-09-09-similarity-search-v1.0.9.json) | Machine-readable similarity-search measurements and ranking checks |
| [`2026-09-11-similarity-search-v1.0.12.md`](2026-09-11-similarity-search-v1.0.12.md) | Three-lane native/native, RDKit-compatible/RDKit, and cross-profile top-10 gate with separate failure accounting |
| [`2026-09-11-similarity-search-v1.0.12.json`](2026-09-11-similarity-search-v1.0.12.json) | Machine-readable v1.0.12 similarity-search gate |
| [`2026-09-12-similarity-search-a3-v1.0.13.md`](2026-09-12-similarity-search-a3-v1.0.13.md) | A3 rerun after explicit-aromatic recovery; 4,500/500 coverage and exact compatible top-10 gate |
| [`2026-09-12-similarity-search-a3-v1.0.13.json`](2026-09-12-similarity-search-a3-v1.0.13.json) | Machine-readable A3 rerun |
| [`rdkit-ecfp4-cross-binding-parity-5000-v1.0.13.json`](../validation/results/rdkit-ecfp4-cross-binding-parity-5000-v1.0.13.json) | Current Rust/Python/Node/WASM ECFP4 cross-binding gate: 5,000/5,000 |
| [`rdkit-ecfp4-cross-binding-parity-5000-parse-morgan-candidate-2026-09-20.json`](../validation/results/rdkit-ecfp4-cross-binding-parity-5000-parse-morgan-candidate-2026-09-20.json) | Frozen Parse + Morgan candidate rebuilt in isolated Python and Node-WASM bindings: every Rust/Python/Node-WASM pair is 5,000/5,000; this is binding parity, not a separate RDKit oracle |
| [`rdkit-ecfp4-sparse-cross-binding-parity-5000-v1.0.13.json`](../validation/results/rdkit-ecfp4-sparse-cross-binding-parity-5000-v1.0.13.json) | Raw sparse identifier/count cross-binding gate: 5,000/5,000 |
| [`rdkit-ecfp4-raw-bitinfo-cross-binding-parity-5000-v1.0.13.json`](../validation/results/rdkit-ecfp4-raw-bitinfo-cross-binding-parity-5000-v1.0.13.json) | Raw identifier bitInfo cross-binding gate: 5,000/5,000 |
| [`rdkit-ecfp4-bitinfo-cross-binding-parity-5000-v1.0.13.json`](../validation/results/rdkit-ecfp4-bitinfo-cross-binding-parity-5000-v1.0.13.json) | Folded bitInfo cross-binding gate: 5,000/5,000 |
| [`native-descriptor-regression-v1.0.13.json`](../validation/results/native-descriptor-regression-v1.0.13.json) | Native descriptor default regression gate: 4 fixed fixtures, 4/4 |
| [`descriptor-rotatable-rdkit-parity-v1.0.13.json`](../validation/results/descriptor-rotatable-rdkit-parity-v1.0.13.json) | A1 additional-family gate: RDKit Strict rotatable bonds, 5,000/5,000, zero failures/unsupported rows |
| [`descriptor-stereocenter-rdkit-parity-v1.0.13.json`](../validation/results/descriptor-stereocenter-rdkit-parity-v1.0.13.json) | A1 diagnostic potential-stereocenter gate: 5,000/5,000 compared, 4,995 exact; five bridged/ring-tied residuals remain unadopted |
| [`rdkit-search-cross-binding-parity-v1.0.13.json`](../validation/results/rdkit-search-cross-binding-parity-v1.0.13.json) | Full chunked 500-query/4,500-entry RDKit-compatible top-k search contract across Rust/Python/Node/WASM; independent RDKit oracle is separate |
| [`descriptor-cross-binding-parity-5000-v1.0.13.json`](../validation/results/descriptor-cross-binding-parity-5000-v1.0.13.json) | Five-field descriptor binding contract: 5,000/5,000 |
| [`2026-09-22-canonical-orbit-perf-v1.0.19.md`](2026-09-22-canonical-orbit-perf-v1.0.19.md) | Exact RENKIN target 2 canonical-search differential: 24 leaves to 1, 4.86x paired-median local measurement, zero output mismatches |
| [`2026-09-11-canonical-orbit-perf-v1.0.13.md`](2026-09-11-canonical-orbit-perf-v1.0.13.md) | Historical exact canonical-search/orbit-pruning differential and instrumentation for issue #372 |
| [`2026-09-09-wasm-rdkit-gate.md`](2026-09-09-wasm-rdkit-gate.md) | Same-corpus Node/WASM comparison with the installed official RDKit.js package |
| [`2026-09-09-wasm-rdkit-gate.json`](2026-09-09-wasm-rdkit-gate.json) | Machine-readable WASM comparison output and exact fingerprint parity count |
| [`2026-09-09-wasm-rdkit-paired.md`](2026-09-09-wasm-rdkit-paired.md) | Same-process paired Node/WASM timing follow-up |
| [`2026-09-09-wasm-rdkit-paired.json`](2026-09-09-wasm-rdkit-paired.json) | Machine-readable paired timing and fingerprint parity output |
| [`2026-09-11-official-rdkit-js-v1.0.12.md`](2026-09-11-official-rdkit-js-v1.0.12.md) | v1.0.12 same-process Node comparison against `@rdkit/rdkit@2025.3.4-1.0.0`, including artifact digests |
| [`../validation/results/competitive-benchmark-rdkitjs-2026-09-11-v1.0.12.json`](../validation/results/competitive-benchmark-rdkitjs-2026-09-11-v1.0.12.json) | Machine-readable official RDKit.js comparison result |
| [`2026-09-10-official-rdkit-js-browser-gate-v1.0.11.md`](2026-09-10-official-rdkit-js-browser-gate-v1.0.11.md) | Real Playwright Chromium comparison against the official package |
| [`../validation/results/competitive-browser-rdkitjs-2026-09-10-v1.0.11.json`](../validation/results/competitive-browser-rdkitjs-2026-09-10-v1.0.11.json) | Machine-readable browser comparison result |
| [`2026-09-16-official-rdkit-js-isolated-browser-v1.0.15.md`](2026-09-16-official-rdkit-js-isolated-browser-v1.0.15.md) | Fresh-browser-process, fixed exposed 10k comparison against `@rdkit/rdkit@2026.03.6`; it records operation-specific results and 9,999/9,999 configured-bit agreement on the declared supported domain across three browsers, plus one typed Fe(II) coordination boundary |
| [`../validation/results/competitive-browser-rdkitjs-isolated-v1.0.15-2026-09-16.json`](../validation/results/competitive-browser-rdkitjs-isolated-v1.0.15-2026-09-16.json) | Machine-readable isolated browser comparison result |
| [`../validation/results/competitive-browser-rdkitjs-ecfp4-parity-5k-v1.0.15-2026-09-16.json`](../validation/results/competitive-browser-rdkitjs-ecfp4-parity-5k-v1.0.15-2026-09-16.json) | Direct browser ECFP4/Morgan comparison against `@rdkit/rdkit@2026.03.6`: 5,000/5,000 exact after the large polycyclic-aromatic explicit-aromatic regression fix; this does not generalize to unmeasured corpus or option configurations |
| [`../validation/results/competitive-browser-rdkitjs-isolated-10k-rss-v1.0.15-2026-09-16.json`](../validation/results/competitive-browser-rdkitjs-isolated-10k-rss-v1.0.15-2026-09-16.json) | Chrome-owned-process 10k memory lane: browser JS heap, chematic linear-memory allocation, and sampled process-tree RSS; summed RSS is not unique physical memory and is not cross-browser evidence |
| [`../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-chromium-10k.json`](../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-chromium-10k.json) | Merged source candidate, Chrome 10k/20-repetition equivalent-output Parse + compatible Morgan comparison: 9,999 supported rows and paired median 2.63x; see the 2026-09-20 record for the completed second-host gate |
| [`../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-firefox-10k.json`](../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-firefox-10k.json) | Merged source candidate, Firefox 10k/10-repetition lane: same configured output and paired median 2.06x; browser timer granularity is retained in the raw record |
| [`../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-webkit-10k.json`](../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-webkit-10k.json) | Merged source candidate, WebKit 10k/10-repetition lane: same configured output and paired median 2.96x; the separately linked GitHub-hosted gate is the second-host measurement |
| [`../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-chromium-chembl-5k.json`](../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-chromium-chembl-5k.json) | Merged source candidate, independent ChEMBL 5k Chrome/20-repetition lane: 5,000/5,000 exact and paired median 6.17x; it is independent-corpus evidence, not a second host |
| [`../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-parity-10k-chromium.json`](../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-parity-10k-chromium.json), [`…-firefox.json`](../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-parity-10k-firefox.json), [`…-webkit.json`](../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-parity-10k-webkit.json) | Direct configured-bit gates for the three local browser lanes: each is 9,999/9,999 supported with one unchanged typed Fe(II) refusal |
| [`../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-parity-chembl-5k-chromium.json`](../validation/results/parse-morgan-rdkitjs-local-candidate-2026-09-20-parity-chembl-5k-chromium.json) | Direct independent ChEMBL configured-bit gate: 5,000/5,000 exact; it proves output equality but not cross-host performance |
| [`2026-09-10-generated-wasm-v3000-gate-v1.0.11.md`](2026-09-10-generated-wasm-v3000-gate-v1.0.11.md) | Generated Web-WASM V3000 metadata preservation gate |
| [`2026-09-11-v3000-cross-engine-v1.0.12.md`](2026-09-11-v3000-cross-engine-v1.0.12.md) | Bounded V3000 canonical/semantic probe across schematic, RDKit, and Open Babel |
| [`2026-09-11-v3000-cross-engine-v1.0.12.json`](2026-09-11-v3000-cross-engine-v1.0.12.json) | Machine-readable V3000 cross-engine probe and explicit Indigo availability boundary |
| [`2026-09-05-prepared-index.md`](2026-09-05-prepared-index.md) | Reusable prepared fingerprint index on the pinned ten-molecule fixture |
| [`2026-09-05-tanimoto-parallel.md`](2026-09-05-tanimoto-parallel.md) | Serial/parallel dense Tanimoto matrix parity and scaling |
| [`2026-09-05-hotpath-110.md`](2026-09-05-hotpath-110.md) | Alternating source A/B gate for canonical SMILES, SDF, and parsing |
| [`2026-09-05-hot-path-follow-up.md`](2026-09-05-hot-path-follow-up.md) | Historical follow-up gate with rejected experiments and exact-output checks |

### Parsing, descriptors, and chemistry workloads

| Record | Scope |
|---|---|
| [`2026-09-24-perf-speed3-output-identical.md`](2026-09-24-perf-speed3-output-identical.md) | v1.0.23 `6f8645ef` vs perf branch `6bd3fae0`: byte-identical outputs (29 operations x 46,736 molecules, 0 differing rows), paired Rust cold timing, and the Python matrix vs RDKit 2026.03.6 re-run on both revisions; RDKit agreement unchanged. 2-vCPU VM, not a package claim |
| [`2026-09-24-perf-digest-diff-v1.0.23-20b2cbc3.json`](2026-09-24-perf-digest-diff-v1.0.23-20b2cbc3.json) | Output-identity digest summary (all operations) for the record above |
| [`2026-09-24-perf-digest-diff-v1.0.23-6bd3fae0-smarts.json`](2026-09-24-perf-digest-diff-v1.0.23-6bd3fae0-smarts.json) | Output-identity digest summary (SMARTS-affected operations) after the warm-up fix |
| [`2026-09-24-perf-digest-time-pair-v1.0.23-6bd3fae0.json`](2026-09-24-perf-digest-time-pair-v1.0.23-6bd3fae0.json) | Paired Rust cold timings for the record above |
| [`2026-09-24-python-op-matrix-vs-rdkit-v1.0.23-6f8645ef.json`](2026-09-24-python-op-matrix-vs-rdkit-v1.0.23-6f8645ef.json) | Python matrix, v1.0.23 wheel |
| [`2026-09-24-python-op-matrix-vs-rdkit-branch-20b2cbc3.json`](2026-09-24-python-op-matrix-vs-rdkit-branch-20b2cbc3.json) | Python matrix, superseded branch wheel `20b2cbc3` (single-query `has_substructure` regression, fixed in `6bd3fae0`) |
| [`2026-09-24-python-op-matrix-vs-rdkit-branch-6bd3fae0.json`](2026-09-24-python-op-matrix-vs-rdkit-branch-6bd3fae0.json) | Python matrix, branch wheel `6bd3fae0` |
| [`2026-09-24-rdkit-agreement-branch-6bd3fae0.json`](2026-09-24-rdkit-agreement-branch-6bd3fae0.json) | RDKit agreement counts for `6bd3fae0` (identical to the accuracy-branch record) |
| [`2026-09-25-perf-digest-diff-v1.0.23-6cb3adce-atom-order.json`](2026-09-25-perf-digest-diff-v1.0.23-6cb3adce-atom-order.json) | Issue #650 atom-order branch `6cb3adce` vs v1.0.23: canonical SMILES, largest fragment, standardize, InChI, CIP and stereocenter outputs byte-identical on the seven benchmark corpora (280,416 rows, 0 differing); the atom-order contract itself is checked by `crates/chematic-smiles/tests/atom_order.rs` (`CHEMATIC_ATOM_ORDER_CORPUS`) |
| [`2026-09-24-rdkit-agreement-accuracy-branch.md`](2026-09-24-rdkit-agreement-accuracy-branch.md) | v1.0.22 `297fec4d` vs accuracy branch `ea967d47`: row-level RDKit 2026.03.6 agreement for 15 RDKit-defined operations on three 5k corpora (atom pair/torsion/MACCS/QED/Murcko/Kekulé descriptors), plus a branch timing re-run; 2-vCPU VM, not a package claim |
| [`2026-09-24-rdkit-agreement-v1.0.22-297fec4d.json`](2026-09-24-rdkit-agreement-v1.0.22-297fec4d.json) | Raw agreement counts for v1.0.22 in the record above |
| [`2026-09-24-rdkit-agreement-branch-ea967d47.json`](2026-09-24-rdkit-agreement-branch-ea967d47.json) | Raw agreement counts for the branch in the record above |
| [`2026-09-24-python-op-matrix-vs-rdkit-branch-ea967d47.json`](2026-09-24-python-op-matrix-vs-rdkit-branch-ea967d47.json) | Raw per-repeat timing and agreement for the branch timing re-run |
| [`2026-09-24-python-op-matrix-vs-rdkit-perf-branch.md`](2026-09-24-python-op-matrix-vs-rdkit-perf-branch.md) | Source branch `8cc01365` vs base `01a86b62` vs RDKit 2026.03.6: 62 Python operations with equivalence classes and row-level output agreement; 21 checked operations faster with full agreement (base: 6); base-vs-branch outputs byte-identical. 2-vCPU cloud VM, not a package claim |
| [`2026-09-24-python-op-matrix-vs-rdkit-base-01a86b62.json`](2026-09-24-python-op-matrix-vs-rdkit-base-01a86b62.json) | Raw per-repeat timings and agreement for the base run of the record above |
| [`2026-09-24-python-op-matrix-vs-rdkit-branch-8cc01365.json`](2026-09-24-python-op-matrix-vs-rdkit-branch-8cc01365.json) | Raw per-repeat timings and agreement for the branch run of the record above |
| [`2026-09-04-canonical-fast-path.md`](2026-09-04-canonical-fast-path.md) | Canonical SMILES on two 5,000-molecule corpora |
| [`2026-09-04-sdf-fast-path.md`](2026-09-04-sdf-fast-path.md) | SDF graph/property read and serialization-only write |
| [`2026-09-04-rdkit-openbabel.md`](2026-09-04-rdkit-openbabel.md) | Earlier RDKit/Open Babel comparison with explicit operation boundaries |
| [`2026-09-04-wasm-size.md`](2026-09-04-wasm-size.md) | Earlier WASM artifact-size measurement |
| [`2026-09-05-descriptor-scaling.md`](2026-09-05-descriptor-scaling.md) | `descriptors_array` column selection, digest, and allocation contract |
| [`2026-09-05-descriptor-streaming.md`](2026-09-05-descriptor-streaming.md) | Descriptor provenance and streaming fixture contract |
| [`2026-09-05-descriptor-topology.md`](2026-09-05-descriptor-topology.md) | Wiener/Kappa/Chi topology descriptors |
| [`2026-09-05-distance-descriptors.md`](2026-09-05-distance-descriptors.md) | AutoCorr2D, Moran, and Geary distance descriptors |
| [`2026-09-05-rdkit-openbabel-speed.md`](2026-09-05-rdkit-openbabel-speed.md) | In-process chematic/RDKit and separately scoped Open Babel CLI timing |
| [`2026-09-05-hotpath-110.json`](2026-09-05-hotpath-110.json) | Machine-readable hot-path gate output |
| [`2026-09-05-mmff94-prepared-nonbonded.md`](2026-09-05-mmff94-prepared-nonbonded.md) | Prepared MMFF94 nonbonded terms and energy parity |
| [`2026-09-05-mmff94-gradient-parallel.md`](2026-09-05-mmff94-gradient-parallel.md) | Bounded parallel finite-difference gradient probes |
| [`2026-09-09-mmff94-nonbonded-gradient-v1.0.10.md`](2026-09-09-mmff94-nonbonded-gradient-v1.0.10.md) | Prepared MMFF94 vdW/electrostatic analytic-gradient parity |
| [`2026-09-09-uff-prepared-topology-v1.0.10.md`](2026-09-09-uff-prepared-topology-v1.0.10.md) | UFF prepared topology, all 45 declared type-parameter soundness, and analytic-gradient parity |
| [v1.0.10 RDKit MMFF94 availability oracle](../validation/results/mmff94-rdkit-availability-oracle-v1.0.10.json) | Independent 265-molecule parse, embed, force-field construction, and finite-energy availability boundary |
| [`2026-09-11-streaming-safety-v1.0.12.md`](2026-09-11-streaming-safety-v1.0.12.md) | Historical v1.0.12 ten-format malformed, oversized, gzip, and generated parser-entry safety gate |
| [`2026-09-11-streaming-safety-v1.0.12.json`](2026-09-11-streaming-safety-v1.0.12.json) | Machine-readable historical v1.0.12 streaming safety gate evidence |
| [`2026-09-11-streaming-safety-v1.0.13.json`](2026-09-11-streaming-safety-v1.0.13.json) | Latest completed v1.0.13 streaming safety gate evidence; not a v1.0.20 rerun |
| [v1.0.12 streaming failure taxonomy](../validation/results/streaming-failure-taxonomy-v1.0.12.json) | Machine-readable typed failure taxonomy for the 120-case base plus twenty Extended XYZ, mmCIF, MOL2, CML, and CDXML supplemental cases |
| [v1.0.10 reaction SMARTS contract](../validation/results/reaction-smarts-bounded-contract-v1.0.10.json) | Bounded 20-case aromatic, bond-order, mapped-agent, agent-OR, disconnected-component assignment/rejection, hydrogen-count, pipe-alternative, and embedding-selection presence contract |

## Streaming and cross-engine contracts

These records primarily measure record accounting, failure behavior, and
boundary semantics. They must not be read as like-for-like throughput claims
unless the record explicitly says that the APIs and process boundaries match.

### Aggregate matrices

| Record | Scope |
|---|---|
| [`2026-09-08-streaming-matrix-v1.0.9.json`](2026-09-08-streaming-matrix-v1.0.9.json) | Rust-only ten-format matrix, plain/gzip stages, limits, digests, and repetitions |
| [`2026-09-08-streaming-cross-engine-matrix-v1.0.9.json`](2026-09-08-streaming-cross-engine-matrix-v1.0.9.json) | Same-input ten-format agreement across chematic and installed RDKit/Open Babel lanes |
| [`2026-09-09-streaming-cross-engine-matrix-v1.0.10.md`](2026-09-09-streaming-cross-engine-matrix-v1.0.10.md) | v1.0.10 same-input contract refresh with explicitly non-ranking throughput context |
| [`2026-09-09-streaming-cross-engine-matrix-v1.0.10.json`](2026-09-09-streaming-cross-engine-matrix-v1.0.10.json) | Historical machine-readable ten-format contract matrix |
| [`2026-09-11-streaming-cross-engine-matrix-v1.0.12.md`](2026-09-11-streaming-cross-engine-matrix-v1.0.12.md) | v1.0.12 same-input contract refresh with explicit parser/process boundaries |
| [v1.0.12 validation matrix](../validation/results/cross-engine-matrix-v1.0.12.json) | Historical machine-readable ten-format contract matrix |
| [`2026-09-11-same-process-sdf-contract-v1.0.12.json`](2026-09-11-same-process-sdf-contract-v1.0.12.json) | v1.0.12 same-process SDF semantic contract |
| [`2026-09-11-same-process-v2000-mol-contract-v1.0.12.json`](2026-09-11-same-process-v2000-mol-contract-v1.0.12.json) | v1.0.12 same-process V2000 MOL semantic contract |
| [`2026-09-11-same-process-v3000-mol-contract-v1.0.12.json`](2026-09-11-same-process-v3000-mol-contract-v1.0.12.json) | v1.0.12 same-process V3000 MOL structural contract |
| [`2026-09-11-same-process-mol2-contract-v1.0.12.json`](2026-09-11-same-process-mol2-contract-v1.0.12.json) | v1.0.12 same-process MOL2 structural contract |
| [`2026-09-11-same-process-xyz-contract-v1.0.12.json`](2026-09-11-same-process-xyz-contract-v1.0.12.json) | v1.0.12 same-process XYZ coordinate contract |
| [`2026-09-11-same-process-extxyz-contract-v1.0.12.json`](2026-09-11-same-process-extxyz-contract-v1.0.12.json) | v1.0.12 same-process Extended XYZ coordinate contract |
| [`2026-09-11-same-process-pdb-contract-v1.0.12.json`](2026-09-11-same-process-pdb-contract-v1.0.12.json) | v1.0.12 same-process PDB coordinate/signature contract |
| [`2026-09-11-same-process-cdxml-contract-v1.0.12.json`](2026-09-11-same-process-cdxml-contract-v1.0.12.json) | v1.0.12 same-process CDXML signature contract |
| [v1.0.10 validation matrix](../validation/results/cross-engine-matrix-v1.0.10.json) | Historical 2026-09-10 ten-format matrix used by the fail-closed validator |
| [`2026-09-09-same-process-sdf-contract-v1.0.10.md`](2026-09-09-same-process-sdf-contract-v1.0.10.md) | Same-process schematic/RDKit SDF semantic contract |
| [`2026-09-09-same-process-sdf-contract-v1.0.10.json`](2026-09-09-same-process-sdf-contract-v1.0.10.json) | Machine-readable same-process SDF evidence |
| [`2026-09-09-same-process-v2000-mol-contract-v1.0.10.md`](2026-09-09-same-process-v2000-mol-contract-v1.0.10.md) | Same-process schematic/RDKit V2000 MOL semantic contract |
| [`2026-09-09-same-process-v2000-mol-contract-v1.0.10.json`](2026-09-09-same-process-v2000-mol-contract-v1.0.10.json) | Machine-readable same-process V2000 MOL evidence |
| [`2026-09-09-same-process-v3000-mol-contract-v1.0.10.md`](2026-09-09-same-process-v3000-mol-contract-v1.0.10.md) | Same-process schematic/RDKit V3000 MOL semantic contract |
| [`2026-09-09-same-process-v3000-mol-contract-v1.0.10.json`](2026-09-09-same-process-v3000-mol-contract-v1.0.10.json) | Machine-readable same-process V3000 MOL evidence |
| [`2026-09-09-same-process-mol2-contract-v1.0.10.md`](2026-09-09-same-process-mol2-contract-v1.0.10.md) | Same-process schematic/RDKit MOL2 semantic contract |
| [`2026-09-09-same-process-mol2-contract-v1.0.10.json`](2026-09-09-same-process-mol2-contract-v1.0.10.json) | Machine-readable same-process MOL2 evidence |
| [`2026-09-09-same-process-xyz-contract-v1.0.10.md`](2026-09-09-same-process-xyz-contract-v1.0.10.md) | Same-process schematic/RDKit XYZ frame contract |
| [`2026-09-09-same-process-xyz-contract-v1.0.10.json`](2026-09-09-same-process-xyz-contract-v1.0.10.json) | Machine-readable same-process XYZ evidence |
| [`2026-09-09-same-process-extxyz-contract-v1.0.10.md`](2026-09-09-same-process-extxyz-contract-v1.0.10.md) | Same-process schematic/RDKit Extended XYZ frame contract |
| [`2026-09-09-same-process-extxyz-contract-v1.0.10.json`](2026-09-09-same-process-extxyz-contract-v1.0.10.json) | Machine-readable same-process Extended XYZ evidence |
| [`2026-09-09-same-process-pdb-contract-v1.0.10.md`](2026-09-09-same-process-pdb-contract-v1.0.10.md) | Same-process schematic/RDKit PDB semantic contract |
| [`2026-09-09-same-process-pdb-contract-v1.0.10.json`](2026-09-09-same-process-pdb-contract-v1.0.10.json) | Machine-readable same-process PDB evidence |
| [`2026-09-10-same-process-sdf-contract-v1.0.10.json`](2026-09-10-same-process-sdf-contract-v1.0.10.json) | Re-run same-process SDF contract against the current v1.0.10 source extension |
| [`2026-09-10-same-process-v2000-mol-contract-v1.0.10.json`](2026-09-10-same-process-v2000-mol-contract-v1.0.10.json) | Re-run same-process V2000 MOL contract |
| [`2026-09-10-same-process-v3000-mol-contract-v1.0.10.json`](2026-09-10-same-process-v3000-mol-contract-v1.0.10.json) | Re-run same-process V3000 MOL contract |
| [`2026-09-10-same-process-mol2-contract-v1.0.10.json`](2026-09-10-same-process-mol2-contract-v1.0.10.json) | Re-run same-process MOL2 contract |
| [`2026-09-10-same-process-xyz-contract-v1.0.10.json`](2026-09-10-same-process-xyz-contract-v1.0.10.json) | Re-run same-process XYZ frame contract |
| [`2026-09-10-same-process-extxyz-contract-v1.0.10.json`](2026-09-10-same-process-extxyz-contract-v1.0.10.json) | Re-run same-process Extended XYZ contract with property boundary |
| [`2026-09-10-same-process-pdb-contract-v1.0.10.json`](2026-09-10-same-process-pdb-contract-v1.0.10.json) | Re-run same-process PDB contract with explicit lenient-parser boundary |
| [`2026-09-10-same-process-cdxml-contract-v1.0.10.json`](2026-09-10-same-process-cdxml-contract-v1.0.10.json) | Re-run same-process CDXML contract with explicit lenient-parser boundary |
| [`2026-09-10-same-process-contract-bundle-v1.0.10.md`](2026-09-10-same-process-contract-bundle-v1.0.10.md) | Reproduction notes and boundaries for the eight-format same-process bundle |
| [`2026-09-10-same-process-sdf-timing-v1.0.10.md`](2026-09-10-same-process-sdf-timing-v1.0.10.md) | Same-process equivalent SDF timing context against RDKit 2025.09.3 |
| [`2026-09-10-same-process-sdf-timing-v1.0.10.json`](2026-09-10-same-process-sdf-timing-v1.0.10.json) | Machine-readable paired SDF timing and semantic evidence |
| [`2026-09-10-mmff94-electrostatic-neighbor-list-v1.0.10.md`](2026-09-10-mmff94-electrostatic-neighbor-list-v1.0.10.md) | Opt-in MMFF94 electrostatic cutoff neighbor-list slice |
| [`../validation/results/mmff94-electrostatic-neighbor-list-v1.0.10.json`](../validation/results/mmff94-electrostatic-neighbor-list-v1.0.10.json) | Machine-readable opt-in MMFF94 electrostatic cutoff evidence |
| [`../validation/results/streaming-failure-taxonomy-v1.0.10.json`](../validation/results/streaming-failure-taxonomy-v1.0.10.json) | Bounded ten-format parser failure variant taxonomy |
| [`2026-09-07-streaming-equivalent.md`](2026-09-07-streaming-equivalent.md) | Same-input SDF ingestion with separately scoped Open Babel evidence |
| [`2026-09-04-streaming-formats.md`](2026-09-04-streaming-formats.md) | Original file-backed SDF/MOL/XYZ runner and non-equivalent RDKit reference |

### Same-input format contracts

| Format | Record |
|---|---|
| SDF | [`2026-09-08-streaming-cross-engine.md`](2026-09-08-streaming-cross-engine.md) |
| XYZ / Extended XYZ | [`2026-09-08-streaming-cross-engine-xyz.md`](2026-09-08-streaming-cross-engine-xyz.md) · [`2026-09-08-streaming-cross-engine-extxyz.md`](2026-09-08-streaming-cross-engine-extxyz.md) |
| V2000 / V3000 MOL | [`2026-09-08-streaming-cross-engine-mol.md`](2026-09-08-streaming-cross-engine-mol.md) · [`2026-09-08-streaming-cross-engine-v3000.md`](2026-09-08-streaming-cross-engine-v3000.md) |
| MOL2 | [`2026-09-08-streaming-cross-engine-mol2.md`](2026-09-08-streaming-cross-engine-mol2.md) · [`2026-09-08-streaming-cross-engine-openbabel-mol2.md`](2026-09-08-streaming-cross-engine-openbabel-mol2.md) |
| Open Babel supplemental | [`2026-09-08-streaming-cross-engine-openbabel.md`](2026-09-08-streaming-cross-engine-openbabel.md) · [`2026-09-08-streaming-cross-engine-openbabel-v3000.md`](2026-09-08-streaming-cross-engine-openbabel-v3000.md) |
| CML / CDXML | [`2026-09-08-streaming-cross-engine-cml.md`](2026-09-08-streaming-cross-engine-cml.md) · [`2026-09-08-streaming-cross-engine-cdxml.md`](2026-09-08-streaming-cross-engine-cdxml.md) |
| mmCIF / PDB | [`2026-09-08-streaming-cross-engine-mmcif.md`](2026-09-08-streaming-cross-engine-mmcif.md) · [`2026-09-08-streaming-cross-engine-pdb.md`](2026-09-08-streaming-cross-engine-pdb.md) |

### Gzip contracts

| Input | Record |
|---|---|
| SDF | [`2026-09-08-streaming-gzip-contract.md`](2026-09-08-streaming-gzip-contract.md) · [`2026-09-08-streaming-gzip-openbabel-sdf.md`](2026-09-08-streaming-gzip-openbabel-sdf.md) |
| XYZ / Extended XYZ | [`2026-09-08-streaming-gzip-xyz-contract.md`](2026-09-08-streaming-gzip-xyz-contract.md) · [`2026-09-08-streaming-gzip-extxyz-contract.md`](2026-09-08-streaming-gzip-extxyz-contract.md) |
| V3000 MOL / MOL2 | [`2026-09-08-streaming-gzip-openbabel-v3000.md`](2026-09-08-streaming-gzip-openbabel-v3000.md) · [`2026-09-08-streaming-gzip-openbabel-mol2.md`](2026-09-08-streaming-gzip-openbabel-mol2.md) |
| CML / CDXML | [`2026-09-08-streaming-gzip-openbabel-cml.md`](2026-09-08-streaming-gzip-openbabel-cml.md) · [`2026-09-08-streaming-gzip-openbabel-cdxml.md`](2026-09-08-streaming-gzip-openbabel-cdxml.md) |
| mmCIF / PDB | [`2026-09-08-streaming-gzip-openbabel-mmcif.md`](2026-09-08-streaming-gzip-openbabel-mmcif.md) · [`2026-09-08-streaming-gzip-openbabel-pdb.md`](2026-09-08-streaming-gzip-openbabel-pdb.md) |

## Artifacts and environment-sensitive records

| Record | Scope |
|---|---|
| [`2026-09-07-wasm-size-v1.0.9.md`](2026-09-07-wasm-size-v1.0.9.md) | v1.0.9 candidate WASM raw/gzip size, digest, toolchain, and commands |
| [`2026-09-09-wasm-size-v1.0.10.md`](2026-09-09-wasm-size-v1.0.10.md) | v1.0.10 current-candidate WASM raw/gzip size, digest, toolchain, and commands |
| [`2026-09-09-wasm-size-v1.0.10.json`](2026-09-09-wasm-size-v1.0.10.json) | Machine-readable v1.0.10 WASM artifact evidence |
| [`2026-09-06-wasm-size-v1.0.8.md`](2026-09-06-wasm-size-v1.0.8.md) | v1.0.8 candidate artifact snapshot |
| [`2026-09-06-wasm-size-v1.0.7.md`](2026-09-06-wasm-size-v1.0.7.md) | v1.0.7 tagged artifact snapshot |
| [`2026-09-09-clean-install-cold-start-v1.0.9.md`](2026-09-09-clean-install-cold-start-v1.0.9.md) | Isolated CPython 3.13 arm64 wheel build/install/import and cold-start evidence |
| [`2026-09-09-clean-install-cold-start-v1.0.10.md`](2026-09-09-clean-install-cold-start-v1.0.10.md) | v1.0.10 clean install, cold start, SMILES throughput, and peak RSS evidence |
| [`2026-09-09-clean-install-cold-start-v1.0.10.json`](2026-09-09-clean-install-cold-start-v1.0.10.json) | Machine-readable v1.0.10 Python evidence |
| [`2026-09-09-ensemble-diversity-v1.0.10.md`](2026-09-09-ensemble-diversity-v1.0.10.md) | Deterministic multi-seed ensemble reproduction and flexible-molecule diversity evidence |
| [`2026-09-09-ensemble-diversity-v1.0.10.json`](2026-09-09-ensemble-diversity-v1.0.10.json) | Machine-readable deterministic ensemble diversity evidence |
| [`2026-09-09-streaming-safety-v1.0.10.md`](2026-09-09-streaming-safety-v1.0.10.md) | Historical ten-format malformed, oversized, and gzip safety gate |
| [`2026-09-09-streaming-safety-v1.0.10.json`](2026-09-09-streaming-safety-v1.0.10.json) | Machine-readable streaming safety gate evidence |
| [`2026-09-09-3d-class-failure-rates-v1.0.10.md`](2026-09-09-3d-class-failure-rates-v1.0.10.md) | Historical 3D status-class failure rates for the 58-molecule gate |
| [`2026-09-09-3d-class-failure-rates-v1.0.10.json`](2026-09-09-3d-class-failure-rates-v1.0.10.json) | Machine-readable 3D class-level failure evidence |
| [`2026-09-09-3d-energy-sanity-v1.0.10.md`](2026-09-09-3d-energy-sanity-v1.0.10.md) | Finite and non-increasing force-field energy checks on the bounded 63-molecule pipeline gate |
| [`2026-09-09-3d-energy-sanity-v1.0.10.json`](2026-09-09-3d-energy-sanity-v1.0.10.json) | Machine-readable force-field energy sanity evidence |
| [`2026-09-09-symmetric-torsion-distance-v1.0.10.md`](2026-09-09-symmetric-torsion-distance-v1.0.10.md) | Local automorphism-aware torsion-distance invariants |
| [`2026-09-09-symmetric-torsion-distance-v1.0.10.json`](2026-09-09-symmetric-torsion-distance-v1.0.10.json) | Machine-readable symmetry-aware torsion-distance evidence |
| [`2026-09-09-symmetric-rmsd-oracle-v1.0.10.md`](2026-09-09-symmetric-rmsd-oracle-v1.0.10.md) | Automorphism-aware RMSD against an independent RDKit oracle |
| [`2026-09-09-symmetric-rmsd-oracle-v1.0.10.json`](2026-09-09-symmetric-rmsd-oracle-v1.0.10.json) | Machine-readable symmetry-aware RMSD oracle evidence |
| [`2026-09-09-rxn-atomic-number-h1-v1.0.10.md`](2026-09-09-rxn-atomic-number-h1-v1.0.10.md) | Bounded `[#N;H1]` reaction compatibility bridge |
| [`2026-09-09-rxn-atomic-number-h1-v1.0.10.json`](2026-09-09-rxn-atomic-number-h1-v1.0.10.json) | Machine-readable atomic-number H1 bridge evidence |
| [`2026-09-09-rxn-atomic-number-h1-v1.0.10.md`](2026-09-09-rxn-atomic-number-h1-v1.0.10.md) | Bounded `[#N;H1]` reaction compatibility bridge |
| [`2026-09-09-workspace-test-v1.0.10.md`](2026-09-09-workspace-test-v1.0.10.md) | Offline workspace-wide Rust unit, integration, and doctest gate |
| [`2026-09-09-workspace-test-v1.0.10.json`](2026-09-09-workspace-test-v1.0.10.json) | Machine-readable workspace test evidence |
| [`2026-09-09-node-wasm-contract-v1.0.10.md`](2026-09-09-node-wasm-contract-v1.0.10.md) | Node/WASM contract smoke with explicit stale-artifact version boundary |
| [`2026-09-09-node-wasm-contract-v1.0.10.json`](2026-09-09-node-wasm-contract-v1.0.10.json) | Machine-readable Node/WASM contract and artifact-version evidence |
| [`2026-09-04-mmff94-3d.md`](2026-09-04-mmff94-3d.md) | Experimental MMFF94, ETKDG, and 3D local microbenchmarks |

## Historical snapshots

These records are retained for provenance and trend context; their versions
and hardware must be read from the record before comparing them with current
results.

| Record | Scope |
|---|---|
| [`2026-09-03-competitive.md`](2026-09-03-competitive.md) | v1.0.1 six-operation competitive run |
| [`2026-09-10-official-rdkit-js-gate-v1.0.11.md`](2026-09-10-official-rdkit-js-gate-v1.0.11.md) | Historical v1.0.11 same-process official RDKit.js comparison |
| [`2026-08-23.md`](2026-08-23.md) | v0.18.0 accuracy, corpus, WASM, and CIP remeasurement |
| [`2026-07-17.md`](2026-07-17.md) | v0.4.29 throughput non-reproduction and descriptor accuracy |
| [`2026-06-25.md`](2026-06-25.md) | v0.4.20 baseline |

## Reproduction and reporting rules

Every new record must include:

- source revision and package versions;
- corpus identity and hash;
- hardware, OS, language/runtime, and build profile;
- exact operation boundary and configuration;
- warm-up, repetitions, aggregation, failure policy, and raw output location;
- correctness, ranking, or byte-equivalence checks relevant to the operation.

Do not relabel source-level A/B data as a published artifact result, compare a
streaming API with a materializing API without saying so, or generalize one
corpus to all chemistry workloads. For the canonical reproduction commands,
see [`docs/benchmark.md`](../docs/benchmark.md).
- [2026-09-09 canonical identity focused gate (JSON)](2026-09-09-canonical-identity-focused-v1.0.10.json) / [report](2026-09-09-canonical-identity-focused-v1.0.10.md)
- [2026-09-09 identity budget gate (JSON)](2026-09-09-identity-budget-gate-v1.0.10.json) / [report](2026-09-09-identity-budget-gate-v1.0.10.md)
- [2026-09-09 reaction and 3D focused gate (JSON)](2026-09-09-reaction-3d-focus-v1.0.10.json) / [report](2026-09-09-reaction-3d-focus-v1.0.10.md)
- [2026-09-09 RDKit TFD evidence boundary](2026-09-09-tfd-oracle-evidence-v1.0.10.md) / [machine result](../validation/results/tfd-oracle-evidence-v1.0.10.json)
- [2026-09-09 Ewald real-space cell-list parity](2026-09-09-ewald-cell-list-v1.0.10.md) / [machine result](../validation/results/ewald-cell-list-v1.0.10.json)
- [2026-09-09 orthorhombic periodic-neighbor cell-list parity](2026-09-09-periodic-neighbor-cell-list-v1.0.10.md) / [machine result](../validation/results/periodic-neighbor-cell-list-v1.0.10.json)
- [2026-09-09 UFF prepared-energy topology](2026-09-09-uff-prepared-topology-v1.0.10.md) / [machine result](../validation/results/uff-prepared-topology-v1.0.10.json)
- [2026-09-09 workspace unit/integration gate (JSON)](2026-09-09-workspace-unit-integration-v1.0.10.json) / [report](2026-09-09-workspace-unit-integration-v1.0.10.md)
