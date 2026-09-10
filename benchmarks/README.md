# chematic benchmark records

This directory contains dated, reproducible measurement records. Numbers are
scoped to the source revision, package versions, corpus, hardware, runtime,
and operation boundary written in each record. They are not universal speed,
accuracy, or compatibility claims.

## Start here

| Need | Start with |
|---|---|
| Understand the rules and how to report a result | [`docs/benchmark.md`](../docs/benchmark.md) |
| Compare current similarity search with RDKit | [`2026-09-09-similarity-search-v1.0.9.md`](2026-09-09-similarity-search-v1.0.9.md) |
| Reproduce the 1.10x hot-path gate | [`2026-09-05-hotpath-110.md`](2026-09-05-hotpath-110.md) |
| Check file-streaming contracts | [v1.0.10 validation matrix](../validation/results/cross-engine-matrix-v1.0.10.json) |
| Check current streaming safety gate | [`2026-09-09-streaming-safety-v1.0.10.md`](2026-09-09-streaming-safety-v1.0.10.md) |
| Check official RDKit.js comparison gate | [`2026-09-10-official-rdkit-js-gate-v1.0.11.md`](2026-09-10-official-rdkit-js-gate-v1.0.11.md) |
| Check official RDKit.js browser gate | [`2026-09-10-official-rdkit-js-browser-gate-v1.0.11.md`](2026-09-10-official-rdkit-js-browser-gate-v1.0.11.md) |
| Check WASM artifact size | [`2026-09-09-wasm-size-v1.0.10.md`](2026-09-09-wasm-size-v1.0.10.md) |
| Find older measurements | [Historical snapshots](#historical-snapshots) |

The current release line is v1.0.11. The newest checked-in performance
records are still versioned historical records where their headers say so;
the v1.0.10 release does not imply that an older measurement was rerun.

## Performance and scaling

### Search and fingerprint paths

| Record | Scope |
|---|---|
| [`2026-09-09-similarity-search-v1.0.9.md`](2026-09-09-similarity-search-v1.0.9.md) | 4,500-entry library / 500-query exact top-k comparison with RDKit; latency and ranking overlap are separate axes |
| [`2026-09-09-similarity-search-v1.0.9.json`](2026-09-09-similarity-search-v1.0.9.json) | Machine-readable similarity-search measurements and ranking checks |
| [`2026-09-09-wasm-rdkit-gate.md`](2026-09-09-wasm-rdkit-gate.md) | Same-corpus Node/WASM comparison with the installed official RDKit.js package |
| [`2026-09-09-wasm-rdkit-gate.json`](2026-09-09-wasm-rdkit-gate.json) | Machine-readable WASM comparison output and exact fingerprint parity count |
| [`2026-09-09-wasm-rdkit-paired.md`](2026-09-09-wasm-rdkit-paired.md) | Same-process paired Node/WASM timing follow-up |
| [`2026-09-09-wasm-rdkit-paired.json`](2026-09-09-wasm-rdkit-paired.json) | Machine-readable paired timing and fingerprint parity output |
| [`2026-09-10-official-rdkit-js-gate-v1.0.11.md`](2026-09-10-official-rdkit-js-gate-v1.0.11.md) | v1.0.11 same-process Node comparison against `@rdkit/rdkit@2025.3.4-1.0.0`, including artifact digests |
| [`../validation/results/competitive-benchmark-rdkitjs-2026-09-10-v1.0.11.json`](../validation/results/competitive-benchmark-rdkitjs-2026-09-10-v1.0.11.json) | Machine-readable official RDKit.js comparison result |
| [`2026-09-10-official-rdkit-js-browser-gate-v1.0.11.md`](2026-09-10-official-rdkit-js-browser-gate-v1.0.11.md) | Real Playwright Chromium comparison against the official package |
| [`../validation/results/competitive-browser-rdkitjs-2026-09-10-v1.0.11.json`](../validation/results/competitive-browser-rdkitjs-2026-09-10-v1.0.11.json) | Machine-readable browser comparison result |
| [`2026-09-10-generated-wasm-v3000-gate-v1.0.11.md`](2026-09-10-generated-wasm-v3000-gate-v1.0.11.md) | Generated Web-WASM V3000 metadata preservation gate |
| [`2026-09-05-prepared-index.md`](2026-09-05-prepared-index.md) | Reusable prepared fingerprint index on the pinned ten-molecule fixture |
| [`2026-09-05-tanimoto-parallel.md`](2026-09-05-tanimoto-parallel.md) | Serial/parallel dense Tanimoto matrix parity and scaling |
| [`2026-09-05-hotpath-110.md`](2026-09-05-hotpath-110.md) | Alternating source A/B gate for canonical SMILES, SDF, and parsing |
| [`2026-09-05-hot-path-follow-up.md`](2026-09-05-hot-path-follow-up.md) | Historical follow-up gate with rejected experiments and exact-output checks |

### Parsing, descriptors, and chemistry workloads

| Record | Scope |
|---|---|
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
| [`2026-09-09-streaming-safety-v1.0.10.md`](2026-09-09-streaming-safety-v1.0.10.md) | Bounded ten-format malformed, oversized, gzip, and generated parser-entry safety gate |
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
| [`2026-09-09-streaming-cross-engine-matrix-v1.0.10.json`](2026-09-09-streaming-cross-engine-matrix-v1.0.10.json) | Machine-readable current ten-format contract matrix |
| [v1.0.10 validation matrix](../validation/results/cross-engine-matrix-v1.0.10.json) | Current 2026-09-10 ten-format matrix used by the fail-closed validator |
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
| [`2026-09-09-streaming-safety-v1.0.10.md`](2026-09-09-streaming-safety-v1.0.10.md) | Current ten-format malformed, oversized, and gzip safety gate |
| [`2026-09-09-streaming-safety-v1.0.10.json`](2026-09-09-streaming-safety-v1.0.10.json) | Machine-readable streaming safety gate evidence |
| [`2026-09-09-3d-class-failure-rates-v1.0.10.md`](2026-09-09-3d-class-failure-rates-v1.0.10.md) | 3D status-class failure rates for the current 58-molecule gate |
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
