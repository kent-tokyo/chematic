# MOL2 indexed reads and fixed4 writes — 2026-10-18

**Local source adoption gate passed.** Large-record parsing improves 1.688–3.829x; prepared writing improves 1.316–3.167x across nine profiles. All seven primary lanes win 24/24 blocks with paired median 95% lower bounds above one. All 27 source medians exceed one and no interval lies wholly below one. Four parse lanes (ethanol, chain-32, chain-64 and stereo fixtures) miss every-block superiority; this is not a strict speed win in every operation.

The [index-only packet](2026-10-18-mol2-index-only-numeric.json.gz) retains the rejected preceding candidate. The [combined packet](2026-10-18-mol2-index-fixed4-numeric.json.gz) is a distinct fresh source-only follow-up. No elapsed observations are pooled. Both complete decisions are independently recomputed by the checked-in stdlib auditor.

## Source and implementation

Runtime baseline main: 229881de8d1cc2639874eba3c5d416d2ceccc469. Measured baseline with only the common probe added: feb6c6f38cb9d84da9c1ac1abae60bbbffdc095e. Rejected index-only source: [2b5a11ab](https://github.com/kent-tokyo/chematic/commit/2b5a11ab). Measured combined source: [a23602ded7029d2a126c8e225b55a47b93d12685](https://github.com/kent-tokyo/chematic/commit/a23602ded7029d2a126c8e225b55a47b93d12685).

Later main 5993e574 adds tests/evidence only. After integration, all 924 archived runtime/harness files remain byte-identical to the combined measured source. Additional terminal N–H imine cases are covered by fresh integration tests. Timings remain assigned to the measured source, not a public artifact or arbitrary later main.

The editable MOL2 reader scanned all atom IDs twice for each bond. Above 32 atoms, it now sorts only the private ID-to-index vector in place and uses binary search. The input atom index is the secondary key, retaining the old first-match rule for duplicate IDs. Smaller records keep linear search. Atom/bond order, coordinates and metadata stay untouched; no array is allocated according to the maximum ID.

The writer still evaluates all original fields. Its x/y/z and partial-charge columns use the existing guarded push_fixed4 formatter, with general formatting for ties, large and nonfinite values. Integer/string padding, Unicode width, optional fields and status bits keep their original formatting contract. There is no new public API, dependency, cache, unsafe code or format policy.

## Combined source measurements

Median baseline/candidate ratios; all raw elapsed nanoseconds, counts, work units, order, warmup and host load observations are in the packet.

| Profile | Parse | Prepared write | Parse + write |
|---|---:|---:|---:|
| ethanol | 1.005255x | 1.659771x | 1.136013x |
| chain-32 | 1.005456x | 1.338977x | 1.088898x |
| chain-64 | 1.050514x | 1.316429x | 1.130140x |
| chain-256 | 1.200137x | 1.328780x | 1.245339x |
| chain-1024 | 1.762425x | 1.348142x | 1.623837x |
| chain-4096 | 3.828540x | 1.351399x | 3.049610x |
| sparse-chain-1024 | 1.688012x | 1.368750x | 1.590588x |
| reversed-chain-1024 | 1.752366x | 1.351574x | 1.619518x |
| mol2-stereo-fixtures | 1.000957x | 3.167467x | 1.044366x |

| Primary lane | Median | Paired median 95% CI | Faster blocks |
|---|---:|---|---:|
| chain-1024/parse | 1.762425x | [1.747677, 1.772645] | 24/24 |
| chain-4096/parse | 3.828540x | [3.797437, 3.864668] | 24/24 |
| sparse-chain-1024/parse | 1.688012x | [1.673588, 1.701800] | 24/24 |
| reversed-chain-1024/parse | 1.752366x | [1.744509, 1.756752] | 24/24 |
| ethanol/write | 1.659771x | [1.649815, 1.677392] | 24/24 |
| chain-1024/write | 1.348142x | [1.334227, 1.363480] | 24/24 |
| chain-4096/write | 1.351399x | [1.328889, 1.369439] | 24/24 |

## Rejected index-only epoch

The four read primaries improved 1.758x, 3.862x, 1.692x and 1.750x, all 24/24. Three prepared-write controls failed the unchanged 0.99 median floor:

| Control | Median | Paired median 95% CI |
|---|---:|---|
| chain-64/write | 0.989413214 | [0.980644096, 1.011749050] |
| chain-4096/write | 0.984838628 | [0.968789752, 0.995356963] |
| reversed-chain-1024/write | 0.988312645 | [0.985857660, 0.996485837] |

All four wholly negative intervals remain: chain-32/roundtrip, chain-4096/write, reversed-chain-1024/write and mol2-stereo-fixtures/write. The writer source was unchanged; no causal attribution follows. The read wins did not make this candidate acceptable. The follow-up changes writer code and retains the floor.

## Open Babel boundary

Open Babel 3.2.1 successfully processes all nine profiles, but each fails the fixed retained-field prerequisite. Residue labels differ: ethanol LIG1 becomes LIG11 on rewrite. Some profiles also differ in charge-method labels, IDs, bond representation or the indexed RDKit observation. These are differences under the strict observer, not a diagnosis that each representation difference is a chemistry error. It checks indexed graph, xyz within 1e-4, names/types/IDs, residue data, charge labels/values, status bits and extra sections. The independent graph observer is local RDKit 2025.09.3.

All 27 index-only external elapsed lanes stay **unranked**, including favorable ratios. Before any new candidate preflight/timing, the follow-up plan was amended to repeat comparator semantics but time only source arms. Its initial plan and amendment remain archived. No earlier external time is assigned to the combined candidate. There is no Open Babel speedup or complete-toolkit claim here. Other formats, SMARTS, descriptors, force fields, conformers, large-file RSS and public packages are outside this experiment.

## Protocol and accounting

Nine exposed development profiles contain 32 records: ethanol; zero-coordinate neutral chains of 32, 64, 256, 1024 and 4096 atoms; 1024-atom sparse maximum-range/reversed IDs; and the exact 24 existing MOL2 stereo fixtures named in the packet. Later-added fixtures are not inserted into this frozen set. Synthetic inputs retain explicit names, types, residue data and charge values. The benchmark separator between records is outside each MOL2 text.

Each profile has 24 balanced blocks, shuffled with seed 20261018; at least 250ms warmup in the measured mode, 100ms timed work and three complete batches per invocation. Parse/roundtrip create and destroy fresh vectors. Prepared writing reuses parsed records; both engines allocate returned strings. Startup/input loading are excluded; framing is inside parse/roundtrip. Full outputs are checked outside inner timing; work units and record counts remain checked.

Index-only uses all six three-arm permutations four times. Combined uses 12 AB/12 BA blocks and reuses the baseline executable by exact hash, with all fresh timing observations. Original source arms and the combined candidate are built from independent empty targets with Rust 1.97.0, opt3/full LTO/one CGU/one build job, incremental=false, strip=none. Host: Apple M4, Mac16,12, arm64, 10 CPUs, 16GiB, macOS 27.0.1. External scheduling is uncontrolled. Bootstrap: 10,000 resamples, median, seeds 0 source and 1 external, percentile endpoints at indices 249/9749.

- Index-only: 216 blocks, 1944 timed invocations, 12,566,442 complete batches, 13,454,656 processed records; gate false.
- Combined: 216 blocks, 1296 timed invocations, 12,775,699 complete batches, 13,891,360 processed records; gate true.

Both formal epochs complete without interruption, reruns, drops or replacement samples. A pre-build full Git archive ran out of disk before any measurement; the incomplete file was removed and source-only archives saved. Completed private compiler intermediates were removed after preserving source, executables and logs; original dirty docs and other worktrees stayed untouched.

Source/probe, binaries, input/preflight data and Open Babel/RDKit/NumPy package contents are sealed before/after profiles. System and package-external transitive libraries are outside the seal. Public packets contain allowlisted numbers, hashes, labels and relative public fixture references. Absolute paths, complete local manifests, SDK files, source archives, binaries and compound streams stay local.

## Validation and reproduction

- Both source contrasts match all 32 full indexed snapshots and rewritten streams exactly.
- 54 ID sets / 27,318 lookups match the frozen original first-match implementation, including zero, duplicate, descending, sparse, missing and maximum IDs.
- 40,981 complete atom rows match the original formatter on IEEE edges, ties/adjacent floats, random bit patterns/ordinary values, Unicode/oversize names and optional fields.
- Measured combined source: 1448 Rust tests, 0 failed, 0 ignored, 76 suites. Test-only main integration: 1449 passed, 0 failed, 0 ignored. Clippy with tests, formatting and refactor ratchet pass.
- Full-data audits recompute both decisions and reject corrupt record counts. The independent numeric auditor checks every observation, order, unit and interval; 12 corruption subcases cover counts, budgets, order, gate, intervals and semantic eligibility. Hosted CI is separate from these local results.

```bash
python3 scripts/check_mol2_index_numeric_packet.py benchmarks/2026-10-18-mol2-index-only-numeric.json.gz
python3 scripts/check_mol2_index_numeric_packet.py benchmarks/2026-10-18-mol2-index-fixed4-numeric.json.gz
python3 -m unittest scripts.tests.test_check_mol2_index_numeric_packet -v
```

For a new experiment, build the native example openbabel_mol2_throughput in separate clean checkouts of the exact baseline/candidate commits, with the recorded flags and empty CARGO_TARGET_DIR values. Copy the probes into a new journal directory as baseline-mol2 and candidate-mol2. Build tools/openbabel_mol2_throughput.cpp with C++17, -O3 and pkg-config --cflags --libs openbabel-3; copy it as openbabel-mol2. The portable driver needs the recorded RDKit/NumPy observer packages. From the final checkout:

```bash
export CHEMATIC_MOL2_EVIDENCE=/absolute/path/to/new-empty-journal-directory
export PYTHONDONTWRITEBYTECODE=1
python3 scripts/bench_openbabel_mol2_index.py --prepare-inputs benchmarks/2026-10-18-mol2-index-fixed4-numeric.json.gz
python3 scripts/bench_openbabel_mol2_index.py --preflight > "$CHEMATIC_MOL2_EVIDENCE/preflight.log" 2>&1
cargo test --locked --offline -p chematic-mol > "$CHEMATIC_MOL2_EVIDENCE/tests.log" 2>&1
# Continue only after successful completion of that test command.
python3 scripts/bench_openbabel_mol2_index.py --seal-after-tests "$CHEMATIC_MOL2_EVIDENCE/tests.log"
python3 scripts/bench_openbabel_mol2_index.py > "$CHEMATIC_MOL2_EVIDENCE/formal.log" 2>&1
```

The preparation mode recreates/verifies all nine input hashes; it was checked locally without collecting another performance epoch. The portable driver copies the source-only measured driver with explicit paths and preparation/seal modes; it does not replace either archived driver. Keep reproduction results separate. Original source archives, build scripts, snapshots, three-arm driver and runtime manifests remain locally retained.
