# chematic roadmap

> Status: strategic roadmap, revised 2026-09-04. v1.0.3 is the current patch release;
> the next work remains version-locked until a new release decision.

## North star

Make chematic the dependable, embeddable cheminformatics runtime for Rust,
Python, JavaScript/WASM, and AI agents: RDKit-compatible where compatibility
matters, smaller and easier to deploy where it does not, and more explicit
about uncertainty than either a silent fallback or a plausible-looking guess.

This is not a promise to reproduce every RDKit behavior. RDKit has a mature
ecosystem, broad language bindings, database integrations, and a long record
of production use. COSMolKit is an emerging Rust-native project with a similar
parity direction. chematic wins by making a narrower set of guarantees
reproducible, portable, inspectable, and useful in browser/agent/materials
workflows.

## Competitive strategy

| Dimension | Target position | Evidence required |
| --- | --- | --- |
| RDKit compatibility | Opt-in, named compatibility profiles for fingerprints, aromaticity, stereochemistry, standardization, and I/O | Version-pinned differential corpora; exact/within-tolerance rates; explicit unsupported outcomes |
| Native deployment | First-class Rust core, Python wheels, Node/WASM, no C/C++ toolchain for the common path | Clean-install matrix, startup/throughput/size benchmarks, API smoke tests |
| Correctness | Deterministic results, bounded work, typed failures, no silent fallback | Reproducible seeds, limits in results, failure and provenance fixtures |
| Explainability | Every important transformation can expose why an atom, bond, bit, parent, or score changed | Trace/report APIs and golden JSON contracts |
| Breadth | One molecule model spanning cheminformatics, reactions, 3D, crystals, and volumetric data | Cross-crate round trips and end-to-end examples |
| Agent/browser use | Safe, small, serializable operations usable through WASM and MCP | Size budget, resource limits, JSON schema, MCP conformance tests |

The comparison harness must never convert “unsupported” into “pass”. A
competitor adapter is useful only when its installation, version, corpus hash,
and operation status are recorded alongside the result.

## Use-case competition tracks

The practical competitors are use-case specific, not one universal toolkit:

| Track | Practical competitor | chematic objective |
| --- | --- | --- |
| Dataset processing and ML input | [sdfrust](https://github.com/HFooladi/sdfrust) | The safest high-throughput molecular dataset pipeline: streaming, metadata-preserving, columnar, resumable, and directly usable from Rust/Python/WASM |
| Fingerprints and molecular ML | [scikit-fingerprints](https://scikit-fingerprints.readthedocs.io/) | A Rust-native fingerprint engine with scikit-learn-compatible Python ergonomics, sparse/bit/count outputs, explainable features, and reproducible 2D/3D transforms |
| SMILES parsing and writing | [OpenSMILES specification](https://opensmiles.org/opensmiles.html) and conforming implementations | The most trustworthy OpenSMILES implementation: strict by default, explicit extensions, excellent diagnostics, bounded parsing, and stable round trips |

The goal is not to win by copying package names or making unsupported claims.
Each track needs a public conformance corpus, a clean-install/throughput
benchmark, and a documented “where to use which tool” migration guide.

## Competitive execution phases (priority order)

This program targets named workflows rather than a single undifferentiated
feature-count claim. RDKit remains the broadest and most mature reference;
Open Babel remains the format-conversion breadth reference; CDK remains the
reaction/SMARTS/QSAR reference; sdfrust remains the Rust dataset/ML reference;
and kekule remains a reference for polymer and molecular-modeling workflows.
COSMolKit is intentionally out of scope for this comparison program.

| Priority | Phase | Competitive wedge | Exit evidence | Weight |
| --- | --- | --- | --- | --- |
| P0 | Trust and measurement | Reproducible, typed, bounded behavior | Version-pinned matrix, corpus hashes, unsupported/failure separation, clean local validator | Critical / light |
| P1 | Interchange throughput | SDF/MOL/XYZ streaming with loss-preserving records | Same-input records/s, bytes/s, RSS, first-error latency, round-trip and malformed-file gates | Critical / heavy |
| P2 | Identity and ML primitives | Stable canonical keys, descriptors, fingerprints, explanations | Held-out differential corpus, cross-binding fixtures, stable-key fail-closed contract, bit explanations | Critical / heavy |
| P3 | Portable production surface | Rust/Python/Node/WASM parity with low deployment friction | Same fixture/config result across bindings, clean-install matrix, size/startup and batch manifests | High / heavy |
| P4 | Chemistry workflow depth | Reactions, SMARTS/MCS, standardization, scaffold and medicinal-chemistry reports | Curated reaction/query corpus, typed ambiguity outcomes, deterministic provenance reports | High / heavy |
| P5 | 3D and materials | Safe conformers, MMFF/UFF, crystals and volumetric interchange | Soundness gates, per-class quality/failure rates, unit/frame round trips, independent references | Medium / heavy |
| P6 | Ecosystem durability | Extensions, migration tooling, browser/agent workflows | Public compatibility dashboard, contributor corpus policy, downstream reproduction reports | Medium / external |

Execution rule: complete the highest-priority local slice that has a
reproducible gate before starting a lower-priority breadth feature. Heavy or
external work is prepared with contracts and fixtures first; it is not called
complete until the required measurement or external reproduction exists.

### P0 — Trust and measurement (current)

- [x] Keep competitor capability states separate from measured results and
  exclude COSMolKit from the active comparison scope.
- [x] Pin the current release, benchmark protocol, corpus hashes, and local
  validators; preserve historical measurements as historical records.
- [ ] Add a checked-in capability matrix for RDKit, Open Babel, CDK, sdfrust,
  kekule, and chematic with `supported`, `partial`, `unsupported`, and
  `not_measured` states per operation.
- [ ] Add a scorecard validator that rejects stale release versions, missing
  corpus/configuration metadata, and claims derived from unsupported rows.

### P1 — Interchange throughput and safety

- [ ] Extend the resumable benchmark to SDF/MOL V2000/V3000, XYZ, MOL2, CML,
  CDXML, PDB/mmCIF, and gzip, including malformed and oversized inputs.
- [ ] Implement bounded streaming batch APIs with cancellation, backpressure,
  deterministic ordering, and an explicit partial-result manifest.
- [ ] Measure chematic against RDKit and Open Babel only where both are
  installed under identical input/configuration; report sdfrust separately
  for dataset-oriented operations.

### P2 — Identity and ML primitives

- [ ] Finish canonical atom-order/E-Z invariance for the supported domain and
  keep `canonical_smiles_stable_key()` as the only dedup/cache recommendation.
- [ ] Inventory and stabilize descriptor/fingerprint shapes, sparse/count
  semantics, configuration, provenance, and explanation APIs.
- [ ] Add held-out parity reports for Morgan/ECFP, MACCS, topological, torsion,
  descriptor, and standardization operations across Rust/Python/WASM.

### P3 — Portable production surface

- [ ] Make the Rust/Python/Node/WASM contract suite consume one fixture schema
  and one versioned expected-result manifest.
- [ ] Publish clean-install, cold-start, throughput, peak-memory, and WASM
  size evidence with explicit platform/configuration metadata.
- [ ] Add browser/agent adversarial cases for limits, cancellation, malformed
  records, and JSON error stability.

### P4–P6 — Breadth after the trust lanes

- [ ] Add reaction/SMARTS/medicinal-chemistry breadth only after P0–P3 gates
  have current evidence and shared graph primitives are stable.
- [ ] Expand 3D/materials and polymer workflows behind soundness, unit, and
  topology contracts; never trade bounded failure for plausible output.
- [ ] Add extension points, migration guides, and public dashboards only after
  measurements are reproducible by a clean checkout.

### D-track — Dataset processing and ML pipeline

#### D0 — Dataset contract and benchmark

- [x] Define the first machine-readable competitive benchmark protocol with
  pinned operation/corpus metadata, fairness rules, explicit unsupported and
  failure statuses, and an offline validator. This records preparation only;
  current comparative measurements remain a separate heavyweight task.
- [x] Add resumable per-operation benchmark execution with atomic state,
  persistent logs, fail-fast status, and `--resume` support; interrupted
  measurements remain explicitly incomplete rather than being treated as
  results.
- [x] Add benchmark environment preflight that rejects a stale installed
  chematic package whose version differs from the target workspace release.
- [x] Execute the current six-operation local benchmark for chematic 1.0.1
  and record raw state, logs, corpus hashes, and operation-specific results in
  `benchmarks/2026-09-03-competitive.md`; Open Babel remains `not_installed`.
- [x] Add an SDF graph/property fast reader and serialization-only writer mode;
  measure them against RDKit with layout work reported separately in
  `benchmarks/2026-09-04-sdf-fast-path.md`.
- [x] Improve that SDF read/write fast path by at least 1.1× on the same
  benchmark: seven-process medians improve 1.26× for read and 1.33× for
  serialization-only write, while strict Z validation remains intact.
- [x] Optimize the canonical SMILES hot path and record repeated comparisons
  on two independent 5,000-molecule corpora; chematic leads RDKit by 2.5% and
  1.47× at the respective medians on the recorded macOS arm64 environment.
  These are scoped results, not a cross-platform or all-corpus claim.
- [x] Add local MMFF94 and 3D pipeline benchmarks covering prepared energy,
  one-shot energy, ETKDG generation, and L-BFGS minimization; record the
  measurements in the 1.0.3 changelog entry.

- Define a lossless record model for molecule, coordinates, SD properties,
  source location, parse status, diagnostics, units, and provenance.
- Benchmark SDF/MOL V2000/V3000, MOL2, XYZ, SMILES tables, and gzip input on
  both valid and malformed corpora. Measure records/s, bytes/s, peak RSS,
  allocations, first-error latency, and partial-result behavior.
- Compare against sdfrust and RDKit with identical files, hardware, and
  configuration; publish success, skip, error, and unsupported counts separately.

#### D1 — Streaming and hostile-file safety

- Add file-backed streaming readers/writers that do not require the complete
  dataset in memory; support bounded batches, backpressure, cancellation, and
  resumable offsets/checkpoints.
- Preserve property order or provide an explicit deterministic ordering policy;
  never silently discard duplicate fields, coordinates, stereo diagnostics, or
  malformed records.
- Apply S1 security limits to records, fields, lines, decompression, coordinates,
  and batch output. Return a manifest that identifies every rejected record.

#### D2 — Columnar and ML-ready representation

- Provide zero-copy or amortized-copy columnar export to NumPy and, where the
  dependency/license/size budget permits, Arrow/Parquet interoperability.
- Add typed labels, missing-value policy, multi-task targets, masks, molecule
  IDs, split IDs, and a stable schema version so dataset preparation is not an
  ad-hoc pandas script.
- Expose OGB/GNN-ready atom, bond, angle, dihedral, neighbor-list, coordinate,
  and mask tensors without changing the source molecule semantics.

#### D3 — Dataset transforms and reproducibility

- Ship streaming transforms for parse/standardize, deduplicate, filter,
  scaffold split, random/stratified split, salt/fragment policy, descriptor
  calculation, fingerprint calculation, and train/validation/test manifests.
- Make seeds, ordering, worker count, normalization, invalid-record policy,
  and source hashes explicit in a machine-readable dataset manifest.
- Support parallel execution without nondeterministic row order or hidden
  cross-process state; add restart-after-failure tests.

#### D4 — Dataset production gate

**Exit gate:** a 100M-record-class workflow can be processed in bounded memory,
resume after interruption, preserve or explicitly report every record decision,
and reproduce the same ML tensors and manifest from the same source hash.

### F-track — Fingerprint and molecular-ML platform

#### F0 — Fingerprint inventory and compatibility matrix

- Inventory every current fingerprint and descriptor representation by input
  requirements (2D/3D), output type (bit/count/float/sparse), dimensionality,
  seed, chirality/stereo policy, and explainability support.
- Prioritize parity gates for ECFP/Morgan, RDKit/RDK, atom pair, topological
  torsion, layered, pattern, MACCS, pharmacophore, Avalon, MHFP/SECFP, MAP4,
  ERG, and reaction fingerprints.
- Compare the full supported surface with scikit-fingerprints' feature catalog,
  while keeping “implemented” separate from “measured against the reference”.

#### F1 — One transformer API across Rust and Python

- Add a common configuration/result model with `fit/transform`-style Python
  ergonomics, direct Rust iterators, batch APIs, and WASM JSON equivalents.
- Support dense NumPy, SciPy sparse, packed bits, sparse counts, raw IDs,
  per-feature counts, and stable feature names/metadata without forcing one
  representation on every workload.
- Make parallelism, chunk size, ordering, invalid-input behavior, and memory
  budget explicit. Avoid one Python allocation per molecule in large batches.

#### F2 — Explainability and ML interoperability

- Return bit/count provenance: atom environments, paths, pharmacophore pairs,
  reaction-center changes, collisions, and the configuration that produced them.
- Add scikit-learn transformer adapters, feature unions, similarity/distance
  functions, nearest-neighbor indexes, and split/evaluation helpers without
  making scikit-learn a Rust-core dependency.
- Add feature-selection, count simulation, folding/collision diagnostics,
  normalization, and serialization contracts suitable for model persistence.

#### F3 — 2D/3D quality and performance

- Close parity and representation gaps before adding more algorithms; every
  compatibility mode gets version-pinned reference fixtures and a holdout.
- Add 3D fingerprints only after conformer validity, conformer provenance, and
  failure policy are explicit. Never silently generate or substitute a conformer.
- Benchmark cold start, throughput, peak memory, dense-vs-sparse crossover,
  parallel scaling, and WASM size against scikit-fingerprints/RDKit on the same
  dataset.

#### F4 — Fingerprint production gate

**Exit gate:** a user can replace a scikit-fingerprints pipeline with a
chematic transformer while preserving shape, ordering, labels, sparse semantics,
and reproducibility; parity and performance are reported per fingerprint rather
than as one aggregate score.

### M-track — OpenSMILES parser and writer

#### M0 — Specification conformance baseline

- Turn the OpenSMILES grammar and semantic sections into a versioned test matrix:
  organic subset, bracket atoms, isotopes, charges, hydrogens, chirality,
  aromaticity, ring closures, branches, dot-disconnected components, and bond
  symbols.
- Separate syntax acceptance, semantic validation, chemistry perception, and
  canonicalization. A parser must not silently reinterpret a valid token merely
  to make downstream chemistry convenient.
- Classify nonstandard behavior explicitly as opt-in extensions, with no hidden
  fallback from strict to permissive mode.

#### M1 — Diagnostics and bounded parsing

- Replace coarse parse errors with stable error codes, byte/character span,
  line/column, expected token class, and a short safe message.
- Add configurable limits for input bytes, atoms, bonds, components, ring labels,
  branch depth, bracket depth, explicit hydrogens, and canonical-search work.
- Guarantee no panic, stack overflow, runaway allocation, or nontermination on
  malformed/adversarial SMILES; integrate the corpus with S1/S2 fuzzing.

#### M2 — Stereo, aromaticity, and extension correctness

- Complete tetrahedral, allene, square-planar, trigonal-bipyramidal,
  octahedral, atropisomeric, and E/Z direction handling, including ring and
  aromatic bond carriers.
- Make aromaticity model, valence policy, implicit-H policy, isotope/charge
  normalization, CXSMILES fields, and unknown extensions configurable and
  visible in the parse result.
- Preserve information through parse → write → parse and canonicalization;
  when preservation is impossible, return a typed diagnostic rather than guess.

#### M3 — Differential and round-trip conformance

- Differential-test against OpenSMILES reference cases, RDKit, Open Babel, and
  other conforming implementations while distinguishing specification behavior
  from implementation convention.
- Maintain valid, invalid, ambiguous, and extension corpora with minimized
  failures. Test atom-order invariance, randomized writing, Unicode/line-ending
  handling, and large disconnected inputs.
- Provide a strict CLI/API mode for validators and a documented compatibility
  mode for migration from permissive toolchains.

#### M4 — OpenSMILES production gate

**Exit gate:** strict mode has complete corpus coverage for the declared
OpenSMILES scope, stable diagnostics, bounded resource use, and zero unexplained
parse/write round-trip corruption; every extension is named and opt-in.

## Release lanes

Each phase has its own acceptance gate. A phase is complete only when the code,
tests, documentation, and measured evidence land together. Roadmap candidates
are not release scope until their gate is closed.

- **Core lane:** graph model, parsing, perception, canonicalization, SMARTS,
  descriptors, fingerprints, reactions, and identity.
- **Compute lane:** 3D conformers, force fields, shape, simulation, crystals,
  and volumetric data.
- **Product lane:** Python/JS/WASM APIs, CLI, reports, MCP, deployment, and
  migration tooling.
- **Trust lane:** differential corpora, fuzzing, reproducibility, provenance,
  performance, and security gates. Every feature uses this lane.

## Security program

“穴を無くす”は証明できないため、ここでは攻撃面を減らし、再現可能な
検査で未修正リスクを見える状態にする。セキュリティ修正は通常の精度改善
と混ぜず、各 S-phase のゲートを満たしてからリリースに含める。

### S0 — Threat model and security baseline

**Goal:** 現在守るもの、信頼境界、未検査箇所を正確に把握する。

- Rust core、各ファイル parser/writer、SMILES/SMARTS/InChI、reaction、3D、
  crystal/volumetric、Python、Node/WASM、MCP、CLI の攻撃面 inventory を作る。
- attacker を「不正な入力を渡せる利用者」「大量入力を送れるサービス利用者」
  「悪意ある依存関係/CI変更者」「権限のあるローカル利用者」に分ける。
- confidentiality、integrity、availability、memory safety、supply chain
  の観点で脅威を記録し、accuracy bug と security bug の判定基準を定義する。
- `SECURITY.md` の現行バージョン、対応期間、連絡先、開示手順、実際に有効な
  GitHub 設定を同期する。未検証の “enabled” 宣言は残さない。
- dependency、GitHub Actions、release artifact、秘密情報、生成コードの
  provenance を一覧化し、各リリースに SBOM と検証記録を紐付ける。

**Exit gate:** 全公開入口に owner、入力形式、信頼境界、制限、検査方法、
既知リスクが記録され、`SECURITY.md` が現行リリースと一致する。

### S1 — Untrusted-input and resource safety

**Goal:** 壊れた/巨大な入力で panic、無限ループ、過剰な割り当て、長時間
停止を起こさない。

**Completed locally (v0.89.0 maintenance):**

- [x] MCP request/response JSON-RPC envelopes have 1 MiB bounds, including
  bounded stdio frame allocation and serialized response size.
- [x] MCP SMILES/SMARTS/MolJSON/PubChem paths have byte, depth, atom, bond,
  match, and upstream-response limits with explicit errors.
- [x] CLI input/output paths have bounded reads/writes, including stdin and
  file destinations; oversized results fail before output is written.
- [x] WASM format JSON, Extended XYZ, tautomer enumeration, and reactant
  execution paths have binding-level input, result-count, atom-count, and
  serialized-output limits.
- [x] Each local milestone updates the workspace version and documents its
  user-facing limit/error contract in `CHANGELOG.md` and
  `docs/error-and-limits.md`.
- [x] WASM MCS configuration and R-group core SMARTS inputs are byte-bounded
  before deserialization or query compilation.
- [x] WASM MMP pair results are capped before JSON materialization to avoid
  quadratic result expansion.
- [x] WASM R-group SMARTS queries are capped at 10,000 query atoms after
  compilation, in addition to their byte boundary.
- [x] The default workspace test gate is reproducibly green with long-running
  Experimental 3D and corpus-scale canonical measurements separated into
  documented explicit lanes; formatter, diff, manifest, and workspace clippy
  checks also pass.
- [x] WASM reaction and library products are checked against the 10,000-atom
  boundary before canonicalization and JSON materialization.
- [x] WASM library templates and fragment inputs are byte-bounded, and every
  fragment is atom-bounded before enumeration.
- [x] WASM UFF minimization and PDBQT conversion inputs are byte- and
  molecule-bounded before coordinate/charge JSON parsing.
- [x] WASM coordinate-driven depiction inputs are byte- and molecule-bounded
  before their lightweight coordinate parser runs.
- [x] WASM reaction, library, and MMP JSON payloads are capped at 16 MiB
  before crossing the JS boundary.
- [x] WASM PDBQT conversion rejects malformed coordinate and charge JSON
- [x] WASM validity-only SMILES and InChI/InChIKey helpers enforce byte and
  atom limits before expensive chemistry processing
- [x] WASM single-SMILES ECFP4/MHFP similarity helpers enforce the shared byte
  and atom limits before fingerprint generation
- [x] WASM pKa, ADMET, and BOILED-Egg SMILES helpers enforce the shared byte and
  atom limits before descriptor calculation
- [x] WASM reaction depiction, normalization, reaction-center, and balance
  helpers use bounded reaction parsing (bytes, components, atoms, and bonds)
  instead of silently replacing it with empty defaults.
- [x] WASM SMILES standardization and standardization-report helpers reject
  oversized input before canonicalization or report generation.
- [x] CLI single-molecule SMILES commands share explicit byte and atom limits
  before descriptor, fingerprint, similarity, substructure, standardization,
  or report work.
- [x] CLI format conversion applies a hard 64 MiB input ceiling and checks the
  parsed molecule atom ceiling independently of parser-specific defaults.
- [x] CLI single-SMILES boundary regressions cover oversized bytes and atom
  counts before chemistry work begins.

These checks cover the named paths above; parser-wide fuzzing and the
cross-binding adversarial contract suite are tracked separately below. The
remaining timeout/cancellation and host-policy work remains open until
measured locally.

**S2 audit progress (local):**

- [x] Unsafe inventory completed: only the optional native-InChI FFI module
  contains executable `unsafe`; default and WASM paths remain safe-only.
- [x] Native-InChI pointer/count/output ownership guards are documented in the
  security-surface inventory.
- [x] `scripts/check_unsafe_surface.py` enforces the reviewed native-InChI
  FFI allowlist as a repeatable local gate.
- [x] A standalone `cargo-fuzz` manifest now covers SMILES canonical round
  trips plus bounded XYZ and MolJSON parser inputs.
- [x] Fuzz smoke execution completed for `smiles`, `xyz`, and `moljson` with
  100 libFuzzer runs per target and no surviving crash.
- [x] A format-dispatch target now exercises the public text parser surface
  across chemistry, SMARTS, reactions, 3D/PDB, and the major Mol formats.
- [x] Long local campaigns completed for all four targets with no surviving
  crash: SMILES ~2.39M, XYZ ~2.35M, MolJSON ~3.63M, and dispatch ~1.78M runs.
- [x] Coverage-preserving corpus minimization completed: SMILES 1,985→1,644,
  XYZ 1,183→510, MolJSON 1,336→1,178, and dispatch 11,125→5,255 cases.
- [x] A dedicated sanitizer workflow and local runner now define AddressSanitizer,
  LeakSanitizer, and ThreadSanitizer test jobs on Linux nightly.
- [x] A dedicated scheduled/manual Miri workflow now exercises the core and
  SMILES parser library tests on Linux nightly, isolated from normal CI.
- [x] Fuzz execution, long campaigns, and corpus minimization are complete on
  the local macOS host, with parser panic regressions recorded.
- [x] Linux sanitizer execution and Miri completed on GitHub Linux after the
  runner was narrowed to the independently buildable core parser scope.
- [x] GitHub Linux sanitizer run 33568431894 completed AddressSanitizer and
  LeakSanitizer successfully on the core/parser scope.
- [x] The same GitHub run completed ThreadSanitizer successfully on the
  core/parser scope.
- [x] Miri run 33588355817 passed four focused high-risk remap/canonical tests;
  the earlier full-library run 33568432012 was cancelled after excessive
  runtime.

**Exit-gate audit (2026-09-03):** Fuzz harnesses, long campaigns, corpus
minimization, all three Linux sanitizer jobs, focused Miri tests, and the
shared adversarial contract suite for the complete common topology parser
surface are evidenced. The format-dispatch fuzz target remains the parser-wide
adversarial backstop for formats that do not have the same API in every
binding.

**Sanitizer execution audit (2026-09-02):** Local macOS runners cannot execute
the Linux nightly configuration. GitHub Linux run
`33568431894` completed AddressSanitizer, LeakSanitizer, and ThreadSanitizer
successfully on the core/parser scope. Focused Miri run `33588355817` also
passed on GitHub Linux.

**v1.0 candidate gate recheck (2026-09-03):** Python binding tests passed
locally with the pinned test environment (`653 passed, 72 skipped`), and the
Node/WASM Node and Web targets were rebuilt and passed all checked-in runtime
contracts. Focused Miri run `33711754307` passed and sanitizer run
`33711755623` passed AddressSanitizer, LeakSanitizer, and ThreadSanitizer on
GitHub Linux. The `atorvastatin_fragment` UFF rescue regression is fixed; its
focused test and the policy bridge group pass. The full default workspace
command `cargo test --workspace --all-targets --locked` now exits 0. The
Experimental 3D long-run tests and corpus-scale canonical measurements are
explicitly quarantined and reproducible with the commands in
`docs/v1.0-local-release-gate.md`; they are not represented as completed
long-run evidence by the default pass. The explicit long-run lanes subsequently
completed successfully: 3D 9/9 in 1467.02s, NCI 5k and descriptor census
passed in the corpus lane, and ChEMBL 4999 passed in a serialized rerun in
1504.50s. The machine-readable record is
`validation/manifests/v1.0.0-local-release-gate.json`.

**Binding execution update (2026-09-02):** The native-incompatible WASM test
is now explicitly gated to `wasm32`; the remaining 306 native WASM tests pass.
The release wheel builds at `0.89.0`, and the full Python 3.13 contract suite
passes (`806 passed`).

**Binding runtime update (2026-09-02):** The Node/WASM package was built with
`wasm-pack --target nodejs --release`, and all 10 checked-in JavaScript runtime
contract files pass, including format parity, pipeline limits/timeouts, RDKit
compatibility fixtures, stereo diagnostics, and typed-array accessors. The
`wasm-pack test --node` Rust harness itself has no registered
`#[wasm_bindgen_test]` cases; the Node package suite is the runtime evidence.
The browser workflow remains a smoke/adversarial supplement. The complete
common topology cross-binding contract suite is now separately executed from
the same `validation/cross_binding_contract.json` fixture by Rust, Python,
and Node-hosted WASM.

**Fuzz execution audit (2026-09-02):** The three bounded targets completed
60-second campaigns with no surviving crash. The format-dispatch campaign
then exposed and enabled fixes for SMARTS isotope/number overflow, PDBQT
Unicode byte slicing, CML Unicode attribute slicing, and chemical-formula
count overflow; each has a regression test. A post-fix 60-second dispatch
campaign completed without a crash, followed by corpus minimization. This is
local macOS/ASan evidence; Linux sanitizer and Miri evidence remain separate.

- 全 parser と writer に atoms/bonds/records、ネスト深度、文字列長、座標、
  voxel 数、反応生成数、match 数、再帰深度の上限を設ける。
- `unwrap`/`expect`/assert による公開 API の入力起因 panic を監査し、typed
  error または明示的な拒否に置き換える。既存の Hosoya/MCS/KET/CIF 系の
  DoS 修正を横断的な共通ポリシーへ昇格する。
- CPU budget、wall-clock timeout、allocation budget、cancellation を
  可能な API に伝播し、上限到達を成功結果に見せない。
- zip/圧縮入力、巨大な一行、重複レコード、循環/極端なグラフ、NaN/Inf、
  異常な isotope/charge/coordinate/unit を拒否または安全に報告する。
- malformed-input corpus、property test、coverage-guided fuzzing を全 parser、
  canonicalization、SMARTS/VF2、reaction、3D、format conversion に接続する。

**Exit gate:** fuzz/property/test harness が公開入力入口を一巡し、panic、
  unbounded allocation、timeout、プロセス停止が再現しない。全制限に negative
  test と利用者向け error contract がある。

### S2 — Memory safety and unsafe minimization

**Goal:** Rust のメモリ安全性を継続的に機械検査する。

- `unsafe` の全箇所を inventory 化し、各ブロックに safety comment、最小化、
  safe alternative の検討、専用 regression test を要求する。不要な `unsafe`
  は削除し、workspace lint で新規追加をレビュー対象にする。
- Miri、AddressSanitizer、LeakSanitizer、ThreadSanitizer、debug assertions、
  release binary の fuzz を CI または定期ジョブで回す。
- Python/PyO3、WASM、Node、C ABI/native-InChI feature の境界について、
  lifetime、buffer length、integer conversion、panic-unwind、threading を検査する。
- index、length、size multiplication、recursion、serialization の integer
  overflow を監査し、checked arithmetic と上限を標準化する。
- fuzz crash は最小 corpus と advisory/修正 commit に結び付け、再発防止 test
  を必須にする。

**Exit gate:** production feature set に未説明の `unsafe`、sanitizer failure、
  fuzz crash、未分類の panic がなく、境界 API の安全性テストが通る。

### S3 — Binding, MCP, CLI, and filesystem boundary hardening

**Goal:** library 内部だけでなく、外部から操作できる境界を安全にする。

- Python/WASM/Node/MCP/CLI で入力 schema、最大 request/response、timeout、
  concurrent jobs、出力サイズを共通化する。
- MCP の tool 引数を strict schema validation し、URL/path/format/engine 名を
  allowlist または安全な encoding で扱う。shell injection、path traversal、
  SSRF、意図しない network/file/process side effect がないことを test する。
- 明示的な file operation には path policy、symlink policy、overwrite policy、
  atomic write、permission/error handling を定義する。default は offline と
  least privilege にする。
- Python exceptions、WASM status JSON、MCP protocol errors、CLI exit code を
  監査し、stack trace、入力分子、秘密情報をログへ漏らさない。
- denial-of-service の rate/size/concurrency policy はホスト側設定と混同せず、
  chematic が保証する範囲を明記する。

**Exit gate:** 各 binding の adversarial contract test が通り、境界越しに
  panic、任意ファイル書込み、任意コマンド実行、意図しない通信、秘密情報の
  漏えいが起きないことを検証できる。

**Binding contract evidence (2026-09-02):**

- [x] A shared `validation/cross_binding_contract.json` fixture is consumed by
  the Rust, Python, and Node-hosted WASM contract tests, so canonical output
  and atom-count expectations cannot drift independently across bindings.
- [x] The same fixture now carries eight malformed cases covering every
  topology parser exposed by all four supported surfaces (SMILES, MOL V2000,
  MOL V3000, MOL2, CML, CJSON, MolJSON, and CDXML). Rust, Python, and
  Node-hosted WASM all assert rejection; the Rust suite is in the integration
  CI matrix.
- [x] Parser-wide adversarial coverage remains explicit in the 22-way
  `fuzz/fuzz_targets/formats.rs` dispatch target; binding parity is asserted
  only where all four surfaces expose the same parser contract.

- [x] Python 3.13 installed wheel contract suite: 806 passed.
- [x] WASM native unit suite: 306 passed after isolating wasm32-only paths.
- [x] Node/WASM runtime contract suite: all 10 checked-in `.test.mjs` files
  passed, including bounded pipelines, format parity, and typed arrays.
- [x] Re-run real browser execution on the rebuilt v1.0.0 WASM assets in the
  isolated Chromium/Firefox/WebKit CI workflow: caffeine computes HBA=6 and
  the page reports v1.0.0 with no console errors. Successful run:
  `33740014401`. Playwright is not installed locally, so this hosted CI run is
  the recorded browser evidence.
- [x] Added isolated Chromium/Firefox/WebKit browser smoke CI for the static
  demo: version badge, caffeine HBA=6, and page/console error absence.
- [x] GitHub Actions run `33611825746` executed the browser smoke suite on
  Chromium, Firefox, and WebKit successfully from the `ci/linux-gates-089`
  evidence branch.
- [x] Browser smoke also verifies the oversized-SMILES fail-closed path across
  all three engines.
- [x] Real-browser adversarial checks cover malformed/oversized SMILES,
  malformed SMARTS recovery, valid SMARTS recovery, and language-toggle state
  restoration.
- [x] Browser smoke now includes a 1 MiB-plus SDF input rejection assertion and
  verifies that no stale SVG is rendered after the fail-closed response.
- [x] Browser smoke covers malformed and empty SDF/MOL input error states;
  GitHub Actions run `33613146015` also passed the oversized-SDF assertion on
  Chromium, Firefox, and WebKit.
- [x] Browser smoke covers valid V2000 MOL/SDF loading, repeated records, and
  recovery from a malformed record after a successful load.
- [x] The valid SDF path is wired into the isolated Chromium/Firefox/WebKit
  matrix; local in-app browser execution verifies rendered output and
  descriptor reset, while matrix execution remains CI evidence.
- [x] Browser smoke verifies recovery from an oversized SMILES rejection to a
  valid calculation before continuing with SMARTS and SDF cases.
- [x] Browser smoke verifies whitespace-trimmed SMILES input after the
  oversized-input recovery path.
- [x] Browser smoke verifies valid SDF/MOL recovery after a malformed-record
  error, including restored rendering and descriptors.
- [x] Browser smoke verifies whitespace-only SMILES rejection and subsequent
  recovery to a valid calculation.
- [x] Browser smoke covers repeated valid calculations across aliphatic,
  aromatic, carboxylic-acid, and heterocyclic molecules without stale errors
  or descriptor state.
- [x] Browser smoke covers malformed ring-closure SMILES rejection and
  recovery after the oversized-input boundary.
- [x] Browser smoke covers empty and malformed-input similarity comparison,
  valid comparison recovery, and SVG output restoration.
- [x] Local browser adversarial pass exercised malformed/recovery paths for
  SMILES, SMARTS, SDF, Similarity, and Reaction without console errors.
- [x] Local browser format pass exercised example parsing and malformed-input
  recovery for Cube, OpenDX, mmCIF, PQR, QCSchema, ORCA Input/Output, and
  LAMMPS Data/Dump.
- [x] Expanded format/adversarial matrix passed in three-browser CI run
  `33616264568`: Chromium, Firefox, and WebKit all passed the 9-format
  example, malformed-input, recovery, and ORCA incomplete-result checks.
- [x] Expanded Reaction scheme/equation and exact oversized-SDF-boundary
  adversarial additions passed in three-browser CI run `33619623893` for
  Chromium, Firefox, and WebKit.
- [x] ORCA Output contract fixed: truncated or unrecognized fragments return
  `Ok` with `termination.kind = incomplete`; explicit ORCA error markers,
  resource-limit violations, and invalid numeric/geometry data are typed
  rejects. Covered by Rust and browser regression tests.

### S4 — Dependency, CI, and release supply chain

**Goal:** 正しいソースから正しい artifact を再現可能に届ける。

- Cargo/PyPI/npm の direct/transitive dependency と GitHub Actions を定期監査し、
  advisories、license、yanked release、unmaintained package を release gate にする。
- Actions を commit SHA または検証可能な固定版へ pin し、最小権限 token、環境分離、
  secret exposure 防止、PR からの publish 防止を確認する。
- clean-room build、reproducible metadata、artifact checksum、SBOM、署名/attestation、
  registry 上の version を検証し、tag や GitHub Release だけで成功扱いしない。
- publish workflow は dependency order、manual approval、rollback/documented
  yanking 手順、package ownership を持つ。公開鍵/token はログへ出さない。
- malicious fixture、build script、proc macro、generated file、downloaded oracle
  を untrusted code として隔離し、検証用依存を production path に混ぜない。

**Exit gate:** dependency/CI audit に未評価の high/critical finding がなく、同じ
  source revision から artifact、SBOM、checksum、provenance を再生成・照合できる。

**Supply-chain audit evidence (2026-09-02):**

- [x] `cargo audit` scanned 308 locked dependencies with zero vulnerability
  findings.
- [x] `cargo deny check` passed advisories, bans, licenses, and sources.
- [x] Copyright attribution is normalized to `Kentaro Tanabe (kent-tokyo)` across the
  license texts, Cargo/Python metadata, citation, README translations, and a
  root `NOTICE`; the SPDX license remains `MIT OR Apache-2.0`.
- [x] `cargo update` refreshed 78 compatible lockfile entries and the full
  workspace then passed `cargo check --workspace --locked`.
- [x] Dependabot #444 `jsonschema` update was integrated locally (`0.50.1` →
  `0.52.1`); workspace check and all 85 MCP unit tests pass.
- [x] PyPI/npm publish workflow action references were pinned to immutable
  commits; other pre-existing workflow references remain tracked below.
- [x] Security workflow action references were also pinned to immutable
  commits.
- [x] The audit's unmaintained `rustybuzz`/`ttf-parser` finding is documented
  as an explicit residual-risk acceptance with a re-evaluation trigger; the
  transitive duplicate versions remain tracked for future dependency cleanup.
- [x] `scripts/generate_sbom.py` produces a 310-package SPDX 2.3 SBOM from
  `cargo metadata --locked`.
- [x] `scripts/generate_provenance.py` records source revision and SHA-256
  hashes for the SBOM, WASM artifact, and Python wheel; OpenSSL detached
  signing and verification were exercised successfully with a disposable
  local key.
- [x] `.github/workflows/supply-chain-evidence.yml` generates SBOM/checksum
  artifacts and requests GitHub artifact attestation with SHA-pinned actions.
- [x] `SOURCE_DATE_EPOCH` fixes the SBOM/provenance timestamp; two independent
  local regenerations produced byte-identical SPDX and provenance files.
- [x] The newly added evidence workflow and the two pre-existing `master`
  toolchain references are pinned to immutable commit SHAs.
- [x] Remaining floating workflow action references for the checked-in CI,
  benchmark, validation, Pages, and publish workflows are pinned to immutable
  SHAs; `actionlint` passes across the complete workflow set.
- [x] Added a detached-signature verification helper that checks exact
  provenance bytes against an externally supplied public key; key custody and
  publication remain intentionally separate gates.
- [x] Added `scripts/check_workflow_pins.py` to prevent future mutable Action
  references from entering the repository; it verifies all checked-in pins.
- [x] Added a protected `Release key evidence` workflow and a portable
  fingerprint-plus-detached-signature verifier; disposable local-key testing
  passed without retaining private key material.
- [x] Production release-key custody was activated in the protected `release`
  environment; the main-branch evidence run `33699428867` verified the
  published public key and retained the signed provenance artifact for
  `v0.89.0`. Fingerprint:
  `f1147c10688e412d183cc6cc0f22017c67874327741815a971c40b362f06ac4e`.
- [x] `v0.89.0` was published from main commit `344e9dc2` with GitHub Release,
  crates.io, npm, and PyPI workflow success.
- [x] The repository-local S4/S5 review gate is now executable via
  `scripts/security_review_gate.py` and runs in the scheduled security
  workflow; it verifies the shared adversarial fixture, 22-way parser fuzz
  dispatch, immutable Action pins, and the reviewed unsafe-code allowlist.

### S5 — Adversarial verification and independent security review

**Goal:** 自分たちの正常系テストでは見つからない穴を継続的に見つける。

- corpus を parser、DoS、memory、binding、MCP、supply-chain、serialization、
  differential のカテゴリに分け、最小化済み regression として保存する。
- libFuzzer/AFL または同等の coverage-guided fuzzing、mutation test、長時間 soak、
  parallel race test、巨大入力 benchmark を定期実行する。
- pre-release ごとに threat-model review、dependency review、public API diff、
  error/log review、red-team checklist を実施する。
- 少なくとも major boundary と parser/serialization 層は、maintainer 以外の
  reviewer または外部監査で確認する。結果は finding、severity、再現手順、修正、
  residual risk に分けて公開する。
- security test が精度/互換性 test を通ったことだけで安全と判断しない。両方の
  gate と malformed/hostile input gate を要求する。

**Exit gate:** リリース候補に未分類の high/critical finding、再現可能な crash、
  未説明の resource exhaustion、未レビューの boundary change がない。

**Local S5 review evidence (2026-09-03):**

- [x] `docs/security-review/v0.89.0-local-review.md` records the local
  threat-model, dependency, public-boundary, error/log, and red-team review
  checklist plus explicit residual risks.
- [x] `scripts/security_review_gate.py` passed offline, including the shared
  adversarial binding fixture, parser-wide 22-way dispatch target, workflow
  pin audit, and unsafe-surface allowlist.
- [x] v1.0 independent-review packet prepared with scope, adversarial actions,
  acceptance criteria, and a sign-off record; preparation is not third-party
  review evidence.
- [ ] Independent maintainer or external review of major parser/serialization
  and binding boundaries; this remains an intentional S5 exit-gate blocker.

**v1.0 gate policy (2026-09-03):** The independent-review item, hosted CI
execution, and external oracle/benchmark campaigns are explicitly excluded
from the repository-local v1.0 release gate. They must not be described as
completed audit evidence; they remain post-release follow-up work. The v1.0
gate still requires the reproducible local S1-S4 checks, binding contracts,
focused Miri, sanitizer/fuzz regressions, and dependency/license review.

### S6 — Response, backport, and continuous maintenance

**Goal:** 発見後の被害と修正遅延を小さくする。

- severity、affected crates/bindings、修正期限、CVE/GHSA 判断、credit、embargo、
  disclosure、ユーザー通知の手順を `SECURITY.md` に固定する。
- supported versions を現行 release と実際に backport できる範囲に合わせ、各修正に
  regression test、CHANGELOG、advisory、影響範囲、migration note を付ける。
- security contact を定期確認し、private report の受領・再現・修正・公開を演習する。
- monthly dependency review、release ごとの attack-surface diff、quarterly fuzz/
  sanitizer refresh を運用化する。
- accuracy limitation、resource limit、security vulnerability を別の公開分類で
  管理し、利用者が危険度を誤認しないようにする。

**Exit gate:** security report の受付から修正版 artifact と advisory の公開までを
  記録付きで実行でき、supported-version 表・テスト・公開文書が同期している。

## Phase 0 — Measurement and compatibility contract

**Goal:** turn “better than the competitors” into falsifiable, versioned
claims.

- Maintain a feature matrix for RDKit, COSMolKit, and chematic with
  `supported`, `partial`, `unsupported`, and `not measured` as separate states.
- Extend the existing common JSONL comparison protocol to descriptors,
  stereochemistry, canonicalization, standardization, SMARTS/MCS, reactions,
  3D, and file formats.
- Publish small smoke corpora plus larger holdouts; include salts, charges,
  isotopes, metals, aromatic edge cases, stereochemistry, tautomer pairs,
  pathological graphs, and malformed input.
- Record engine version, feature/configuration, corpus hash, seed, runtime,
  status, error class, and provenance for every row.
- Define scorecards per operation. Never publish one composite “chemistry
  accuracy” number that hides parse failures or unsupported features.

**Exit gate:** a reproducible CI report can answer, per operation, where
chematic matches, differs, fails safely, or is not yet comparable.

## Phase 1 — Trustworthy molecular core

**Goal:** make the graph and identity layer safe enough to be a drop-in
building block.

- Close the remaining aromaticity, fused-ring, charge-aware kekulization,
  E/Z, CIP, allene, atropisomer, and stereo round-trip gaps.
- [x] Add fail-closed `canonical_smiles_stable_key()` and a regression fixture
  for cosmetic E/Z carrier spelling; the historical 275/5000 E/Z-only metric
  remains explicitly open, with the recovered Issue #11 corpus and a separate
  current diagnostic projection pinned as evidence.
- [x] Add a public RDKit-like aromaticity regression gate for purine and
  azulene; keep the compatibility-preserving Hückel default model-distinct.
- [x] Freeze the v0.89 support boundary for Python `RWMol`, loss-preserving
  CDXML editing, bounded polymer expansion, and partial RDKit/Morgan
  compatibility in `docs/compatibility-scope.md`.
- [x] Freeze aromaticity and CIP model boundaries: explicit `RdkitLike`
  aromaticity selection, opt-in `CipMode::Accurate`, and documented fused,
  non-alternant, symmetric-cage, and unresolved-tie residuals.
- [x] Freeze the proposed v1.0 compatibility contract: bounded CDXML/polymer
  scope, partial Python `RWMol`, fail-closed canonical identity, explicit
  aromaticity/CIP modes, and Experimental 3D/MMFF94 boundaries.
- Finish canonical SMILES invariance across atom order, equivalent spellings,
  isotope/charge/stereo state, and directional bond systems.
- Make standardization and Parent identity configurable, bounded, and
  explainable: fragments, salts, charges, zwitterions, tautomers, isotopes,
  stereo, and `super_parent` must return status and transformation reports.
- Add a strict SMARTS/MCS contract: query-preserving output, timeout/limit
  reporting, recursive-query coverage, and substructure self-verification.
- Expand exact compatibility profiles for RDKit Morgan/ECFP, RDK, layered,
  atom-pair, torsion, pattern, MACCS, and reaction fingerprints. Keep native
  fingerprints separate from compatibility fingerprints in names and docs.
- Complete exact/optional InChI behavior and document the boundary between
  the pure-Rust implementation and the native feature.
- Add property-based and fuzz testing for parsers, writers, canonicalization,
  SMARTS, and reaction input; every discovered crash becomes a minimized
  regression fixture.

**Exit gate:** no known silent corruption in the core paths; compatibility
claims have operation-level corpus evidence; failures are typed and bounded.

## Phase 2 — Production Python and batch platform

**Goal:** beat the installation and workflow friction of native-heavy stacks.

- Provide a stable Python facade with predictable naming, NumPy/pandas batch
  operations, streaming SDF/SMILES/CSV, multiprocessing-safe behavior, and
  clear migration examples from common RDKit calls.
- Add a small, composable CLI for parse, standardize, descriptors, fingerprints,
  substructure, similarity, reactions, report, and format conversion.
- Define memory, molecule-size, timeout, and parallelism policies for every
  batch operation; return partial-result manifests instead of losing the
  whole job.
- Benchmark cold start, throughput, peak memory, wheel size, and clean
  installation against RDKit and COSMolKit on the same corpus and hardware.
- Ship self-contained HTML/JSON/CSV reports with bit explanations, changed
  atoms, warnings, and provenance.

**Exit gate:** a new user can install and run a documented batch workflow on
Linux/macOS/Windows without a compiler, and the result is reproducible from a
single manifest.

## Phase 3 — 3D quality and conformer ensembles

**Goal:** make the experimental 3D lane useful without overstating parity.

- [x] Freeze the v0.89 3D/MMFF94 declaration boundary: legacy coordinate
  generation remains experimental, `embed_pipeline_v2` is opt-in and
  evidence-bearing, and incomplete MMFF94 coverage is an observable failure
  rather than an implicit fallback.
- [x] UFF rescue now rechecks and repairs declared tetrahedral/E-Z stereo after
  minimization; the `atorvastatin_fragment` catastrophic-blowup regression is
  covered by a passing policy-bridge test, with finite-coordinate and bond-
  length soundness checks retained. When the first stereo-preserving start does
  not yield a sound result, the rescue now tries two additional fixed seeds;
  the search is deterministic, bounded, and remains fail-closed.
- [x] Complete the default workspace test gate without a timeout by separating
  Experimental 3D long-run and corpus-scale canonical tests into explicit
  `--ignored` lanes; document both commands and record the distinction in the
  v1.0 local release manifest.
- [x] Execute the separated long-run lanes: all 9 Experimental 3D tests pass;
  NCI 5k, descriptor census, and serialized ChEMBL accuracy corpus checks pass
  with no failure or panic.
- [x] Freeze the long-run execution manifest with commands, candidate revision,
  durations, timeout policy, and the serialized ChEMBL rerun note.
- [x] Add an offline validator for the long-run manifest and include it in the
  local v1.0 release-gate command set.
- [x] Add an offline release-document consistency validator covering the fixed
  version, `chematic` product name, release-key secret name, and gate links.
- [x] Audit current-facing documentation examples for v1.0.2 and preserve old
  versions only where they are historical benchmark or release-note records.
- Stabilize the modern bounded embedding pipeline and clearly separate it from
  legacy APIs; remove or quarantine paths that can produce unsound geometry.
- Complete ETKDG-style knowledge and constrained embedding, including ring,
  macrocycle, stereochemical, metal, and difficult charge cases.
- Close MMFF94/UFF typing, charge, bonded-term, stretch-bend, torsion, and
  minimization gaps with independent oracle and soundness gates.
- Make ensemble generation deterministic, seedable, force-field aware, and
  provenance-rich; measure diversity/coverage, not only best-single-conformer
  RMSD.
- Add 3D descriptor, shape, pharmacophore, alignment, and similarity parity
  reports with symmetry-aware metrics.

**Exit gate:** every successful conformer passes bond/stereo/energy sanity
checks; quality and failure rates are reported per chemical class and force
field, with no unqualified “RDKit parity” claim.

## Phase 4 — Reactions and medicinal-chemistry workflows

**Goal:** provide useful design workflows around the same molecular model.

- Expand reaction SMARTS/SMIRKS parsing, matching, atom mapping, stereo/E/Z
  preservation, reagent/agent separation, and deterministic product
  enumeration.
- Add reaction fingerprints, reaction-center explanations, retrosynthesis
  templates, MMP/MMS, R-group decomposition, scaffold networks, BRICS, and
  diversity selection as composable APIs.
- Add explicit valence, mapping, aromaticity, and product-validity reports;
  never invent a product when a transform is ambiguous or unsupported.
- Validate against curated reaction corpora and hand-audited edge cases, with
  precision/recall and invalid-product rates rather than only output counts.

**Exit gate:** a reaction workflow is deterministic, explainable, and safe to
  use in batch design without silently dropping stereo, mapping, or products.

## Phase 5 — Materials and simulation bridge

**Goal:** extend the advantage beyond traditional small-molecule toolkits.

- Harden crystal/periodic structures: CIF/POSCAR symmetry, cell operations,
  disorder/occupancy policy, periodic bonds, and round-trip provenance.
- Add validated QCSchema, ORCA, Cube/OpenDX, PQR, mmCIF, LAMMPS, and related
  interchange workflows, including unit and coordinate-frame contracts.
- Provide Ewald/PME, force-field, charge, and molecular-dynamics interfaces
  behind explicit stability and energy-conservation gates.
- Add materials fingerprints, periodic-neighbor queries, volumetric analysis,
  and structure-to-molecule conversion where chemically justified.

**Exit gate:** cross-format round trips preserve declared fields and units;
  simulation features publish independent reference comparisons and stability
  limits before being labeled production-ready.

## Phase 6 — Browser, edge, and AI-native product surface

**Goal:** make serious chemistry practical where RDKit deployment is awkward.

- Keep the common WASM path within a published raw/gzip size budget and expose
  resource limits, cancellation, and typed JSON errors.
- Provide browser-safe parsing, standardization, fingerprints, substructure,
  descriptors, 2D depiction, reports, and selected 3D operations.
- Version the WASM/Node/Python/Rust schemas and generate conformance fixtures
  so the same input/configuration has the same result across bindings.
- Expand MCP tools with safe defaults, provenance, compact reports, and
  machine-readable uncertainty; include local/offline operation as a first
  class mode.
- Add a no-backend demo and notebook/agent examples that show practical
  workflows, not only API listings.

**Exit gate:** the same conformance suite passes across Rust, Python, Node,
  and WASM; browser execution remains responsive and bounded on adversarial
  inputs.

## Phase 7 — Ecosystem and durable advantage

**Goal:** make the project difficult to displace even when feature lists
  converge.

- Maintain a public compatibility dashboard and release notes that distinguish
  measured, inferred, and planned behavior.
- Offer stable extension points for descriptors, fingerprints, force fields,
  formats, reaction templates, and external engines.
- Build a contributor corpus policy: provenance, licensing, fixture minimization,
  oracle version pinning, and review requirements for chemistry changes.
- Publish migration guides and reference implementations for RDKit and
  COSMolKit users, including when chematic is not the right choice.
- Track real adoption metrics: successful clean installs, workflow completion,
  issue reproduction time, benchmark regressions, and downstream breakage.

**Exit gate:** releases are boring to upgrade, claims are independently
  reproducible, and users can choose a compatibility profile without coupling
  their application to internal implementation details.

## Near-term order

The next cohesive releases should follow this order:

1. **S0 security baseline**: synchronize `SECURITY.md` with v1.0.2, inventory
   every public input/binding boundary, and verify the repository's claimed
   security controls.
2. **S1 resource safety**: make parser, SMARTS/MCS, reaction, format, 3D, and
   MCP limits uniform, then add malformed-input and panic/DoS gates.
3. Phase 0 scorecard and corpus contract, including a maintained COSMolKit
   adapter and current RDKit version pin.
4. **D0/D1 + M0/M1 baselines**: establish the sdfrust dataset benchmark and
   OpenSMILES conformance/resource-safety corpus before adding convenience APIs.
5. **F0/F1 fingerprint inventory and transformer contract**: make output shape,
   sparse/count semantics, ordering, configuration, and provenance stable.
6. Phase 1 remaining identity/stereo/standardization gaps, prioritizing issues
   that affect fingerprints, substructure, and canonical output downstream.
7. **D2/D3 + F2**: deliver columnar/ML-ready dataset transforms and
   explainable fingerprint features with Python/Rust parity.
8. **S2/S3 boundary hardening**: audit unsafe and FFI/binding paths, then close
   Python, WASM, Node, MCP, CLI, and filesystem boundary tests.
9. **S4 supply-chain gate**: lock down dependencies, Actions, artifact
   provenance, SBOM, checksums, and registry verification before the next
   multi-registry release.
10. **M2/M3 + F3**: close stereo/extension behavior, differential conformance,
    2D/3D quality, and performance evidence.
11. Phase 2 Python/batch/CLI ergonomics plus clean-install and memory/throughput
   evidence.
12. Phase 3 3D soundness and ensemble-level measurement before adding more 3D
   features.
13. **D4/F4/M4 production gates**, then S5/S6 continuous security review and
    response drills, followed by Phase
   6 cross-binding conformance and agent/browser hardening.

Do not start a broad feature expansion when a trust-lane regression is open in
the same primitive. Fix the shared primitive, rerun the affected corpus and
holdout, then decide whether the result forms a release theme.

## Definition of “ahead”

chematic is ahead of a competitor only for a named use case and a measured
dimension. The claim must include the versions, corpus, configuration, hardware,
failure policy, and reproduction command. “Faster”, “more accurate”, “smaller”,
and “more compatible” are otherwise roadmap hypotheses, not product claims.
