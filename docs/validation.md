# Validation report

Updated 2026-10-15 for **v1.1.4** plus current **Unreleased** changes.
Every result is limited to its recorded artifact or source revision,
comparator, corpus, options, platform, and failure policy.

## Evidence boundary

- v1.1.4 is the current registry release.
- The latest complete hosted Python chemistry packet is v1.1.3.
- The broadest published-package RDKit performance packet is v1.0.42.
- Later source records are useful development evidence, but they are not
  relabelled as v1.1.4 package results.
- Typed refusals and abstentions stay in the denominator and are not matches.

## Current evidence

| Area | Result and limit | Record |
|---|---|---|
| Hosted Python artifacts | Exact v1.1.3 PyPI wheels on Linux x86-64, macOS arm64, macOS x86-64, and Windows x86-64 reproduce the pinned reaction, 310k SMARTS, Morgan, descriptor, CIP, hybridization, MMFF94, MOL stereo, BioTransformer, and A6 counts. npm/WASM and crates.io are separate evidence. | [v1.1.3 packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-v113-published-artifact-packet.md) |
| Unusual-input descriptors | The v1.1.3 package packet confirms HBD on 4,072/4,072 rows and QED on all 4,071 comparable rows; RDKit rejects the remaining QED input. The v1.1.4 release adds further MW, Morgan, ASA, and format fixes, but its full cross-platform artifact packet is still open. | [v1.1.3 packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-v113-published-artifact-packet.md) · [QED follow-up](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-13-qed-stress-parity.md) |
| RDKit performance | Exact v1.0.42 CPython 3.13 wheels on hosted Linux x64 and macOS arm64 pass 10k TPSA, Labute ASA, and chiral-Morgan output gates. Parse-inclusive intervals favor chematic on both hosts; prepared macOS Labute loses and prepared Morgan crosses parity. | [published packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-11-v1042-published-rdkit-packet.md) |
| Reactions and workflows | The 83-row gate has 80 exact product graph/origin/map rows and three inputs invalid in both engines. The bounded RDKit workflow packet also has 39/39 largest-fragment, 254/254 Normalize after common-input filtering, and 16/16 labelled R-group rows. These are bounded fixtures, not universal parity. | [reaction follow-up](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-10-v1042-regression-followups.md) · [workflow packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-12-rdkit-workflows-expanded-source.md) |
| Browser Morgan | Published v1.0.30 npm artifacts have 250/250 direct and prepared bit parity in Chromium, Firefox, and WebKit, with favorable 20-block intervals. A separate 10k Chromium lane has 9,999 exact rows and one typed Fe refusal. | [three-browser packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-linux-three-browser-paired20.md) |
| COSMolKit comparison | Source comparisons cover 161 operations plus a 4,072-row unusual-SMILES corpus. Seeded ETKDG matches RDKit on 4,977/4,977 ChEMBL and 9,947/9,947 RDKit.js rows where RDKit succeeds. These are source results, not general toolkit superiority. | [0.5 comparison](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-10-cosmolkit-050.md) · [ETKDG packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-09-cosmolkit-parity-2.md) |
| 3D/MMFF94 | The current source reproduces RDKit 2026.09.1 topology bounds within `1e-12` on 768,664 frozen A6 cells. Equal-count 20-conformer quality is within 1 kcal/mol on 258/265 rows. Remaining tails, unused holdout, two-host timing, memory, and exact-package evidence stay open. | [bounds follow-up](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-12-a6-rdkit-2026-09-bounds-followup.md) |
| File interchange | Current source records cover bounded V3000/SDF/MOL2/CML/CIF/CDXML semantics and typed loss reporting. Diverse Open Babel parser rows with semantic differences are not ranked for speed. | [V3000 packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-openbabel-v3000-speed.md) · [SDF follow-up](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-openbabel-sdf-read-followup.md) |
| Rust execution coverage | The release workflow gate reports 96.5056% line coverage (153,415/158,970). Python and WASM coverage are not included. Coverage is not chemistry validation. | [coverage workflow](https://github.com/kent-tokyo/chematic/actions/workflows/coverage.yml) |

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

## Unreleased safeguards

Changes merged after v1.1.4 restore the v1.1.3 exhaustive-match surface of
`CmlError` and `CifError`, add `cargo-semver-checks`, pin known RDKit
2026.03.6/2026.09.1 ring-count disagreements, and update `rustls`. They are
not evidence about the published v1.1.4 artifacts.

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
