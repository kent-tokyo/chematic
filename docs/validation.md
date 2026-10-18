# Validation report

Updated 2026-10-18 for **v1.1.7**.
Every result is limited to its recorded artifact or source revision,
comparator, corpus, options, platform, and failure policy.

## Evidence boundary

- v1.1.7 is the current registry release.
- The v1.1.7 published-wheel chemistry workflow completed all 11 hosted jobs.
  The latest durable checked-in cross-platform packet remains v1.1.3.
- The broadest published-package RDKit performance packet is v1.0.42.
- Source and earlier package records are useful development evidence, but they
  are not relabelled as v1.1.7 package results.
- Typed refusals and abstentions stay in the denominator and are not matches.

## Current evidence

| Area | Result and limit | Record |
|---|---|---|
| Hosted Python artifacts | Exact v1.1.7 PyPI wheels completed 11 Linux, macOS, Windows, RDKit-comparison, A6, and autoconf jobs. The hosted run establishes pass/fail evidence; its durable hashes and counts still need a checked-in packet. npm/WASM and crates.io remain separate evidence. | [v1.1.7 hosted run](https://github.com/kent-tokyo/chematic/actions/runs/39078219154) · [v1.1.3 checked-in packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-v113-published-artifact-packet.md) |
| Unusual-input descriptors | The checked-in v1.1.3 packet confirms HBD on 4,072/4,072 rows and QED on all 4,071 comparable rows; RDKit rejects the remaining QED input. v1.1.7 completed its hosted cross-platform gate, but exact durable per-row counts should be quoted from a checked-in packet rather than inferred from workflow success. | [v1.1.7 hosted run](https://github.com/kent-tokyo/chematic/actions/runs/39078219154) · [v1.1.3 packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-v113-published-artifact-packet.md) |
| RDKit performance | Exact v1.0.42 CPython 3.13 wheels on hosted Linux x64 and macOS arm64 pass 10k TPSA, Labute ASA, and chiral-Morgan output gates. Parse-inclusive intervals favor chematic on both hosts; prepared macOS Labute loses and prepared Morgan crosses parity. | [published packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-11-v1042-published-rdkit-packet.md) |
| Reactions and workflows | The 83-row gate has 80 exact product graph/origin/map rows and three inputs invalid in both engines. The bounded RDKit workflow packet also has 39/39 largest-fragment, 254/254 Normalize after common-input filtering, and 16/16 labelled R-group rows. These are bounded fixtures, not universal parity. | [reaction follow-up](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-10-v1042-regression-followups.md) · [workflow packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-12-rdkit-workflows-expanded-source.md) |
| Current standardization source | Against RDKit 2026.09.1 on 4,069 comparable stress rows, current source has Normalize 4,065 exact, Cleanup 3,923 exact and Reionize 3,150 exact. Cleanup fixes 87 pre-integration residuals with zero regressions. Every retained residual is classified rather than promoted to a match. | [machine record](https://github.com/kent-tokyo/chematic/blob/main/validation/results/rdkit-standardization-current-source-2026-10-18.json) · [summary](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-18-trust-gates-and-residuals.md) |
| Browser Morgan | Published v1.0.30 npm artifacts have 250/250 direct and prepared bit parity in Chromium, Firefox, and WebKit, with favorable 20-block intervals. A separate 10k Chromium lane has 9,999 exact rows and one typed Fe refusal. | [three-browser packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-linux-three-browser-paired20.md) |
| COSMolKit comparison | Source comparisons cover 161 operations plus a 4,072-row unusual-SMILES corpus. Seeded ETKDG matches RDKit on 4,977/4,977 ChEMBL and 9,947/9,947 RDKit.js rows where RDKit succeeds. These are source results, not general toolkit superiority. | [0.5 comparison](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-10-cosmolkit-050.md) · [ETKDG packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-09-cosmolkit-parity-2.md) |
| 3D/MMFF94 | The current source reproduces RDKit 2026.09.1 topology bounds within `1e-12` on 768,664 frozen A6 cells. Equal-count 20-conformer quality is within 1 kcal/mol on 258/265 rows. Remaining tails, unused holdout, two-host timing, memory, and exact-package evidence stay open. | [bounds follow-up](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-12-a6-rdkit-2026-09-bounds-followup.md) |
| File interchange | A current-source CTfile stress rerun has 4,045 exact readbacks, ten stereo-only differences, one chematic readback failure and 16 RDKit oracle failures. Other records cover bounded V3000/SDF/MOL2/CML/CIF/CDXML semantics and typed loss reporting. | [current CTfile classification](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-18-trust-gates-and-residuals.md) · [V3000 packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-openbabel-v3000-speed.md) |
| Independent holdouts | After candidate `f5d6cd6f` was frozen, NCI NCISMA99 and PubChem CID 4,097–20,000 were acquired and hashed. MCS comparison (64 rows), SDF attributes (64), and stereo/3D (96) were selected by fixed SHA-256 rank without running an engine or inspecting outcomes. CI checks candidate/source ordering and file hashes. The rows are not yet candidate results. | [manifest](https://github.com/kent-tokyo/chematic/blob/main/validation/trust_holdouts/v1/manifest.json) |
| Rust execution coverage | The processed Linux Codecov lane reports 96.93% line coverage. Python and WASM coverage are not included. Coverage is not chemistry validation. | [coverage workflow](https://github.com/kent-tokyo/chematic/actions/workflows/coverage.yml) |

## Important bounded contracts

- **A0 core-eight:** the frozen candidate passed 2,000/2,000 development and
  one 8,000/8,000 sealed evaluation on eight declared descriptors. The sealed
  rows cannot be reused for tuning.
- **CIP:** versioned comparisons include exact labels and typed abstentions.
  Wildcard, lone-pair, and oracle-sensitive cases are never guessed.
- **SMARTS:** the current gate contains 310,000 cells. Historical package and
  later source counts differ; consult the exact dated record before quoting a
  number.
- **Canonical identity:** `canonical_smiles()` is not a universal identity
  key. `canonical_smiles_stable_key()` is fail-closed and bounded.
- **3D:** typing, energy, gradient, convergence, stereo, clash, conformer
  quality, time, and memory are independent gates.

## v1.1.7 safeguards

v1.1.7 adds checked scientific-file conversion and export APIs, custom
fingerprint and force-field inputs, source-preserving SDF/PDB/MOL2 workflows,
and typed compatibility helpers. Required source CI and Rust coverage passed;
the exact registry wheels completed the hosted cross-platform chemistry gate,
while a durable checked-in per-row artifact packet remains follow-up work.

## How to interpret results

A speed win requires equivalent declared work, acceptable equivalent output,
complete failure accounting, a paired 95% lower bound above 1.0, and recorded
source/package/host provenance. Parse-inclusive and prepared timings answer
different questions.

Formal reruns must use each record's exact versions, hashes, corpus, options,
and runtime. Earlier measurements remain available in the
[benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md), [complete inventory](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/archive-index.md),
and [versioned results](https://github.com/kent-tokyo/chematic/tree/main/validation/results/).

## Local consistency checks

```bash
python3 scripts/check_release_docs_consistency.py
python3 scripts/check_benchmark_index.py
python3 scripts/check_roadmap_disposition.py
python3 scripts/check_product_names.py
```

See [compatibility scope](compatibility-scope.md), [benchmark method](benchmark.md),
and the [accuracy plan](rdkit-accuracy-plan.md) for the support boundary.
