# Maintenance scripts

Scripts are retained when they are used by CI, reproduce a checked-in evidence
record, generate a committed fixture, or provide a documented release operation.
One-off diagnostics should be removed after their result is captured unless a
current issue or benchmark links to them.

## Find a script

| Task | Start here |
|---|---|
| CI and release gates | Search `.github/workflows/` for the script path; most gate entry points start with `check_`. |
| Evidence reproduction | Search `docs/`, `benchmarks/`, and `validation/` for the script path or basename. The linked record defines inputs and expected outputs. |
| Benchmarks | `bench_*` entry points and their dated records in `benchmarks/`. |
| Fixture generation | `generate_*`, `gen_*`, or `build_*` entry points; committed fixture data belongs under `validation/` or a test fixture directory. |
| Boundary fixture regeneration | `generators/generate_*_boundary_fixtures.py`; use the pinned RDKit version in the corresponding `validation/` record. |
| Historical diagnostics | `diagnostics/` contains retained reproductions with checked-in results or active references. |
| Same-process format contracts | `validation/check_same_process_*_contract.py`; the bundle gate is `validation/check_same_process_contracts.py`. |
| Earlier benchmark runners | `benchmarks/bench_*.py`; use the dated record for corpus and engine versions. |
| Shared code and tests | Importing scripts and `scripts/tests/`. A filename without a path reference may still be an imported helper or a discovered test. |

The [dated inventory](../docs/script-inventory-2026-10-20.md) records the
reference audit. `tools/script_root_baseline.list` lists the existing direct
children of this directory; `tools/check_script_layout.py` checks it in CI.

## Main entry points

| Purpose | Entry point |
|---|---|
| MOL2 Kekule bond output | `validation/validate_mol2_kekule.py` checks exact token changes, Rust/Python/WASM outputs, RDKit identity and explicit valence, and Open Babel cross-reads for [#2668](https://github.com/openbabel/openbabel/issues/2668). See the [reproduction instructions](../docs/mol2-kekule-bond-output.md). |
| Explicit PDB CONECT graph | `validation/validate_pdb_conect.py --wasm-module <node-module> --output <directory> [--issue-pdb <original-9581.pdb>]` compares indexed graphs/charges/XYZ/serials against RDKit 2026.09.1 with proximity disabled, using a newly built Python module and Node WASM. `validation/validate_pdb_conect_node.mjs` is its batch driver; see [the accepted profile](../docs/pdb-conect-graph.md). |
| Prepared MMFF dielectric updates | `validation/validate_mmff_dielectric.py --openbabel-probe <executable> --output <directory>` rechecks 56 pinned RDKit goldens, compares new conformers across screening changes and minimization, and records the original Open Babel #2811 stale-cache probe. See [the numerical contract](../docs/mmff-dielectric-updates.md). |
| ChemDraw byte-block dispatch | `generators/generate_chemdraw_fixtures.py` writes authored CDX/CDXML pairs; `validation/validate_chemdraw_block.py --wasm-module <module> --output <directory>` compares Python and Node/WASM graphs to RDKit 2026.09.1 XML reads and reproduces #9604. `validation/validate_chemdraw_block_node.mjs` is its Node driver; see [the profile](../docs/chemdraw-block.md). |
| Open Babel SDF batch throughput | `build_openbabel_sdf_speed_corpus.py --offset` generates fixed unfiltered V2000 slices and refuses format changes; `benchmarks/bench_openbabel_sdf_throughput.py` validates parse and rewrite semantics separately and records 21 native timing/RSS blocks. See [the original packet](../benchmarks/2026-10-14-openbabel-sdf-speed.md) and [the reader follow-up](../benchmarks/2026-10-14-openbabel-sdf-read-followup.md). |
| Open Babel V3000 block throughput | `generators/build_openbabel_v3000_speed_corpus.py` retains unfiltered force-V3000 slices; `benchmarks/bench_openbabel_v3000_throughput.py` compares pre-split native MOL API batches in 21 rotating blocks. Rejected outputs keep their row positions and fail the semantic gate. The [packet](../benchmarks/2026-10-14-openbabel-v3000-speed.md) separates diverse chemistry from controlled chain-size scaling. |
| Python speed frontier | `bench_python_speed_frontier.py` compares exact wheels, RDKit, and COSMolKit across first-use, repeated-use, and parse-inclusive loops; `--offset` selects a fixed corpus slice. `verify_python_descriptor_frontier.py` and `verify_python_morgan_frontier.py` retain float bits, errors, and fingerprints. See the [base packet](../benchmarks/2026-10-14-python-speed-frontier.md), [scalar/stereo follow-up](../benchmarks/2026-10-14-python-scalar-stereo.md), [core-model preparation](../benchmarks/2026-10-14-python-model-view-speed.md), and [sanitized-view reuse](../benchmarks/2026-10-14-python-sanitized-view-reuse.md). |
| Long scalar reuse batches | `bench_python_scalar_reuse.py` retains and validates every result of one million prepared TPSA/ASA calls per arm, independently of the short-batch frontier protocol. Both original and additional protocols are in [the follow-up packet](../benchmarks/2026-10-14-python-scalar-stereo.md). |
| Release/version checks | `check_release_*.py`, `check_publish_graph.py`, `bump_version.py` |
| Compatibility dashboard | `check_compatibility_profiles.py`, `generators/generate_compatibility_dashboard.py` |
| Benchmark index | `check_benchmark_index.py` |
| Fingerprint cross-binding parity | `validation/check_fingerprint_cross_binding.py --profile <operation> --output <new-report.json>` compares Rust, Python, and Node/WASM; select one of nine named profiles. Historical results remain under `validation/results/`. |
| Refactoring debt ratchet | `check_refactor_budget.py` checks `validation/refactor_budget_v1.json`; `--print-baseline` prints a reviewed replacement without writing files |
| Open Babel file-I/O comparison | `check_openbabel_file_io_contract.py`, `check_openbabel_file_io_gate.py`, `check_openbabel_file_io_semantics.py`, `benchmarks/bench_openbabel_file_io_paired.py`, `benchmarks/bench_openbabel_file_io_same_process.py`, and `check_openbabel_file_io_evidence.py` |
| Open Babel 3D MOL2 identity | `compare_openbabel_mol2_3d_stress.py` gates RDKit V2000 identity, feeds identical Open Babel-written MOL2 to both readers, and records per-row stereo identity. See the [measured report](../validation/openbabel_comparison/mol2-3d-stereo-2026-10-17.md). |
| Open Babel unmarked 3D SDF identity | `validation/openbabel_comparison/check_unmarked_3d_stereo.py` scores the same V2000/V3000 conformers in CheMatic and Open Babel; the [disjoint ChEMBL follow-up](../validation/openbabel_comparison/unmarked-3d-chembl-disjoint-2026-10-17.md) retains exact inputs, row outcomes, and a remaining Open Babel-only cage. |
| Parser security | `run_isolated_parser_security.py` |
| Browser comparison | `bench_browser_wasm_vs_rdkit_isolated.py` |
| Published browser replication | `check_published_browser_paired.py`, `check_v1030_published_browser_morgan_rows.py`, and `check_published_linux_browser_replication.py` validate the pinned Chromium and Ubuntu three-engine records |
| Sealed cohort | `prepare_sealed_accuracy_cohort.py` |
| Frozen trust holdouts | `validation/trust_holdouts/freeze.py` selects outcome-blind rows from pinned sources; `check_trust_holdouts.py` verifies the committed manifest. |
| V3000 external readers | `v3000_*_gate.py` |
| Stereo gates | `stereo_torture_suite_gate.py`, `stereo_spelling_invariance_gate.py` |
| Open Babel 0D MOL stereo source check | `compare_openbabel_0d_stereo.py` compares before and after CLI binaries on identical Open Babel-generated records; accepts JSONL or `.smi` corpora and records generation failures and source-gate exclusions |
| Paired Open Babel stereo preservation | `compare_openbabel_stereo_formats.py` runs CheMatic and Open Babel on the same RDKit-written V2000/V3000 inputs and retains every semantic outcome. |
| A6 source-wheel termination diagnostic | `public_package_3d_chematic.py` runs the fixed 265 rows; `summarize_a6_mmff94_termination.py` accounts for every outcome; `compare_a6_source_published_rows.py` checks same-host success retention without claiming the full quality gate |
| A6 platform gate (#739) | `pipeline_v2_vs_rdkit_common_scorer` scores one platform's 265 rows; `a6_platform_gate.py` requires 265 successes, all independently sound, stereo-clean and clash-free, and reports whether the coordinates equal another platform's row by row |
| Cross-binding Node dump | `binding_dump.mjs <mode>` |
| Published v1.0.30 artifact packet | `check_v1030_artifact_packet.py`, `check_v1030_published_chemistry_residuals.py`, `emit_rdkit_smarts_oracle.py`, `run_published_npm_chemistry_lane.mjs`, `run_published_npm_63op_outputs.mjs`, `run_reaction_compatibility_v2.py`, `reaction_template_map_gate.py`, `check_published_rust_63op_outputs.py`, and the version-pinned `tools/published_rust_gate*` lockfiles |
| Published-wheel output and speed comparison | `bench_python_op_matrix_vs_rdkit.py --outputs-only/--paired`, `check_published_python_version_outputs.py`, `bench_published_python_versions.py` |
| Isolated Python speed and process RSS | `bench_published_python_isolated_paired.py` (20 fresh-process AB/BA blocks), `check_published_python_isolated_paired.py` (archived six-lane integrity gate) |
| Isolated published-Rust speed and process RSS | `bench_published_rust_isolated_paired.py` builds two `--locked` crates.io graphs and runs fresh-process AB/BA blocks; `check_published_rust_isolated_paired.py` verifies archived outputs, paired intervals and RSS |
| Isolated published npm/WASM speed and process RSS | `bench_published_wasm_paired.mjs` pins v1.0.29/v1.0.30/RDKit.js tarballs and runs 20 fresh-process ABBA/BAAB blocks; `check_published_wasm_paired.py` verifies the archived output gate, raw block accounting, intervals and whole-process RSS |

## Organization rule

- Put new scripts in `scripts/checks/`, `scripts/validation/`,
  `scripts/benchmarks/`, `scripts/diagnostics/`, `scripts/generators/`, `scripts/release/`, or
  `scripts/lib/` according to purpose; put their tests in
  `scripts/tests/`. The `scripts/` root is a legacy path list, not the default
  destination for new scripts.
- When moving a legacy script, update every CI and documentation path, check
  imports and `__file__`-based repository paths, and remove its entry from
  `tools/script_root_baseline.list`. Do not rewrite archived measurements.
- Shared mechanics belong in one helper; callers select a named profile/mode.
- Evidence scripts write versioned JSON and keep raw failures visible.
- Historical measurements keep their original script contract in Git history.
- Historical generated-package reports are immutable evidence, not a check
  against the current source export surface.
- A script with no CI, documentation, issue, benchmark, or fixture reference is
  a deletion candidate.
- Generated packages, build output, and temporary corpora do not belong in this
  directory.

### Checked 3D distance-matrix validation

`validation/validate_checked_3d_matrix.py` compares the source-built checked and legacy
APIs with installed RDKit on 128 ordinary synthetic coordinate cases and eight
numeric boundary recipes. It records full matrices, finite/nonfinite categories,
source/runtime hashes and zero-valued failure counters. It does not execute
the oversized upstream heap-crash reproducer.

### HIN aromatic-bond validation

`validation/validate_openbabel_hin.py` compares the original synthetic HIN fixture through
a source-built chematic Python extension and a pinned installed Open Babel
binary. It records both MDL bond-type lists, source/artifact hashes and raw
outputs. This checks Open Babel #3000's bounded bond semantics, not general
chemical equivalence or performance.

### CCD molecular templates

`validation/validate_ccd.py --download` caches eight official RCSB component files, compares
CCD atom/bond/coordinate fields against pinned Gemmi 0.7.5, and records the local
Open Babel CIF/mmCIF diagnostic. It needs a source-built Python extension and
Gemmi in the oracle environment. Result hashes identify the downloaded bytes;
the hermetic CI fixture is independently authored. This validates single CCD
components, not a bulk dictionary or package release.

## Molecular CIF bond geometry

`validation/validate_cif_bonds.py` checks the source Python extension's explicit molecular
CIF bond output against Gemmi 0.7.5 and records bounded Open Babel CLI probes.
It uses original fixtures in `validation/openbabel_issues/cif-bonds/`.
Provide the compiled extension and the oracle environment on `PYTHONPATH`, then
run `python3 scripts/validation/validate_cif_bonds.py`. Native hermetic tests run with
`cargo test -p chematic-mol --test cif_bond_writer`; binding regressions are in
`crates/chematic-py/tests/test_cif_bond_writer.py`. See
[`docs/cif-molecular-bonds.md`](../docs/cif-molecular-bonds.md) for the projection
boundary and the distinction between Open Babel's `-ag` and `-xg` options.

## Rooted SMILES stereo

`validation/validate_rooted_stereo.py --cases /tmp/rooted-cases.json` compares original
constructed bridgehead configurations with RDKit 2026.09.1's accurate CIP
labeler. Generate cases with
`cargo run -p chematic-smiles --example rooted_stereo_probe`. Provide the
source-built Python extension and pinned RDKit on `PYTHONPATH`; optional
`--wasm-module` adds exact real Node/WASM versus Python checks. The output
contains a compressed complete packet and a source/runtime-hashed summary.
See [the checked API scope](../docs/rooted-smiles-stereo.md).
## XYZ element identity

`validation/validate_xyz_symbol_case.py` checks canonical/uppercase XYZ species against
RDKit 2026.09.1's periodic table and runs bounded Open Babel 3.2.1 controls.
Provide the source-built Python extension and pinned RDKit on `PYTHONPATH`;
optional `--wasm-module` adds real Node/WASM checks, and
`--baseline-pythonpath` probes a retained pre-fix source runtime separately.
The complete compressed packet and source/runtime-hashed summary preserve
the element/coordinate comparison boundary. See
[the format-local scope](../docs/xyz-element-symbols.md).
## CDXML stereo inference

`validation/validate_cdxml_stereo_options.py` compares 216 transformed non-directional
CDXML drawings with explicit directional RDKit 2026.09.1 controls. Supply a
fresh source extension and pinned RDKit on `PYTHONPATH`; `--wasm-module`
adds real Node/WASM, and `--baseline-pythonpath` diagnoses a retained pre-fix
source runtime separately. Full isomeric graph, CIP and default controls
have distinct counts. See [the inference boundary](../docs/cdxml-stereo-options.md).
## Hydrogen-count accounting

`validation/validate_hydrogen_count_contract.py` uses pinned RDKit 2026.09.1 graph
references for total nuclei and attached-H atom fields. Provide fresh source
chematic and RDKit on `PYTHONPATH`; `--wasm-module` adds real Node/WASM, and
`--baseline-pythonpath` probes an earlier source runtime separately.
Open Babel formula and serialized-output diagnostics have separate counts.
See [the accounting contract](../docs/hydrogen-count-accounting.md).
### Checked PDBQT output

`validation/validate_pdbqt_checked.py` checks fixed-column output, Python/WASM byte parity,
refusals, an Open Babel #2101 reproduction and independent Open Babel cross-read.
See `docs/pdbqt-checked-output.md` for runtime/probe setup and limits.
### Python stub packaging gate

`check_python_stubs.py` validates every source `.pyi` and the `py.typed` marker.
`--wheel "dist/chematic-*.whl"` additionally checks all matched archives without
importing native modules, requires the typed marker and compares all package
stub bytes with source. An unmatched wheel pattern is a failure. JSON reports
are optional via `--report`. See `docs/python-stub-packaging-gate.md` for scope.

`validation/validate_python_stub_gate.py --wheel /path/to/chematic.whl` checks the real CLI
against an existing wheel and eight temporary negative mutations. Optional
`--rdkit-stubs` and `--rdkit-version` record a read-only upstream syntax audit.
### Proof-backed tetrahedral cleanup

`validation/validate_tetrahedral_cleanup.py` compares Rust/Python/WASM decisions and
independently verifies anchored atom mappings, frames and RDKit-supplied E/Z
facts. Use the native `tetrahedral_cleanup_dump` example and actual Node/WASM
module; setup and precise scope are in `docs/tetrahedral-stereo-cleanup.md`.
