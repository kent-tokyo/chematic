# Validation report

Updated 2026-10-21 for the **v1.2.0 release**.
Every result is limited to its recorded artifact or source revision,
comparator, corpus, options, platform, and failure policy.

## Evidence boundary

- v1.2.0 is the current release. Its exact PyPI wheels passed the hosted Linux,
  macOS and Windows runtime smoke and 11-job chemistry workflows. This is
  package pass/fail evidence; source-only timing does not become package timing.
- One public v1.1.8 macOS arm64 wheel has completed the 10,000-row Labute ASA
  correctness and performance lane. Other operations and platforms do not
  inherit that result.
- The v1.1.8 release has not completed an equivalent broad, cross-platform
  public-wheel chemistry packet.
- The v1.1.9 published-wheel chemistry workflow completed all 11 hosted jobs.
  The latest durable checked-in cross-platform packet remains v1.1.3.
- The broadest published-package RDKit performance packet is v1.0.42.
- Source and earlier package records are useful development evidence, but they
  are not relabelled as v1.2.0 package results.
- Typed refusals and abstentions stay in the denominator and are not matches.

## Current evidence

| Area | Result and limit | Record |
|---|---|---|
| Hosted Python artifacts | Exact v1.2.0 PyPI wheels passed the three-platform runtime smoke and all 11 Linux, macOS, Windows, RDKit-comparison, A6, BioTransformer and autoconf chemistry jobs. All 25 wheels plus the sdist match their Actions-built SHA-256 values. The hosted runs establish package pass/fail evidence; the latest durable per-row packet remains v1.1.3. | [runtime smoke](https://github.com/kent-tokyo/chematic/actions/runs/39511354658) · [chemistry run](https://github.com/kent-tokyo/chematic/actions/runs/39511354905) · [v1.1.3 checked-in packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-v113-published-artifact-packet.md) |
| Current public performance | One exact v1.1.8 macOS arm64 wheel matches all 10,000 RDKit and COSMolKit Labute ASA values. It wins all 20 paired blocks over RDKit for first use, repeated use and parse-inclusive use; COSMolKit first use is a statistical tie. This is one host, wheel and operation. | [v1.1.8 public-wheel record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-19-v118-published-labute-asa-cosmolkit.md) |
| Unusual-input descriptors | The checked-in v1.1.3 packet confirms HBD on 4,072/4,072 rows and QED on all 4,071 comparable rows; RDKit rejects the remaining QED input. v1.1.9 completed its hosted cross-platform gate, but exact durable per-row counts should be quoted from a checked-in packet rather than inferred from workflow success. | [v1.1.9 hosted run](https://github.com/kent-tokyo/chematic/actions/runs/39249060732) · [v1.1.3 packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-v113-published-artifact-packet.md) |
| RDKit performance | Exact v1.0.42 CPython 3.13 wheels on hosted Linux x64 and macOS arm64 pass 10k TPSA, Labute ASA, and chiral-Morgan output gates. Parse-inclusive intervals favor chematic on both hosts; prepared macOS Labute loses and prepared Morgan crosses parity. | [published packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-11-v1042-published-rdkit-packet.md) |
| Reactions and workflows | The 83-row gate has 80 exact product graph/origin/map rows and three inputs invalid in both engines. The bounded RDKit workflow packet also has 39/39 largest-fragment, 254/254 Normalize after common-input filtering, and 16/16 labelled R-group rows. These are bounded fixtures, not universal parity. | [reaction follow-up](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-10-v1042-regression-followups.md) · [workflow packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-12-rdkit-workflows-expanded-source.md) |
| Current standardization source | Against RDKit 2026.09.1 on 4,069 comparable stress rows, source has Normalize 4,069 exact, Cleanup 4,053 exact and Reionize 4,064 exact. The ordered acid/base catalog fixes 25 Reionize rows with no regressions on this corpus. Five remaining Reionize strings differ because a tetrahedral tag retained in RDKit's in-memory result disappears when its SMILES is reparsed; they remain mismatches in the exact-output score. | [machine record](https://github.com/kent-tokyo/chematic/blob/main/validation/results/rdkit-standardization-source-2026-10-21-reionize-catalog.json) · [summary](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-21-rdkit-reionize-catalog-source.md) |
| Browser Morgan | Published v1.0.30 npm artifacts have 250/250 direct and prepared bit parity in Chromium, Firefox, and WebKit, with favorable 20-block intervals. A separate 10k Chromium lane has 9,999 exact rows and one typed Fe refusal. | [three-browser packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-linux-three-browser-paired20.md) |
| TPSA source follow-up | Frozen source `f8e9d7de` → `add2b6b4` improves first-use TPSA 1.037x/1.270x on two fixed macOS windows. All 42 source median ratios meet 0.99; the largest observed decrease is 0.845%. Scalar values/refusals remain source-exact over 42,216 cells. COSMolKit first-use deficits remain, and later main integration was not re-timed. | [full report and numeric packets](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-18-terminal-tpsa.md) |
| COSMolKit comparison | Source comparisons cover 161 operations plus a 4,072-row unusual-SMILES corpus. Seeded ETKDG matches RDKit on 4,977/4,977 ChEMBL and 9,947/9,947 RDKit.js rows where RDKit succeeds. These are source results, not general toolkit superiority. | [0.5 comparison](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-10-cosmolkit-050.md) · [ETKDG packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-09-cosmolkit-parity-2.md) |
| 3D/MMFF94 | The current source reproduces RDKit 2026.09.1 topology bounds within `1e-12` on 768,664 frozen A6 cells. Equal-count 20-conformer quality is within 1 kcal/mol on 258/265 rows. Remaining tails, unused holdout, two-host timing, memory, and exact-package evidence stay open. | [bounds follow-up](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-12-a6-rdkit-2026-09-bounds-followup.md) |
| File interchange | A current-source CTfile stress rerun has 4,045 exact readbacks, ten stereo-only differences, one chematic readback failure and 16 RDKit oracle failures. Other records cover bounded V3000/SDF/MOL2/CML/CIF/CDXML semantics and typed loss reporting. | [current CTfile classification](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-18-trust-gates-and-residuals.md) · [V3000 packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-openbabel-v3000-speed.md) |
| Independent holdouts | After the candidate was frozen, NCI NCISMA99 and PubChem CID 4,097–20,000 were acquired by fixed SHA-256 rank. Immutable v1.1.8 was scored once: all SDF columns are 64/64, MCS atom/bond counts are 31/64, and 3D embed/stereo status is 96/96 while conformer count is 95/96 and strict coordinate parity is 4/96. Success, mismatch, refusal, invalid input and internal error remain separate. | [score](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-19-trust-holdout-v1-score.md) · [row-level result](https://github.com/kent-tokyo/chematic/blob/main/validation/results/2026-10-19-trust-holdout-v1-v118.json) |
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

## v1.2.0 public artifacts and retained package evidence

The v1.2.0 public wheels passed their runtime and chemistry workflows. The
exact v1.1.8 macOS arm64 wheel remains the latest recorded public-package
performance lane for Labute ASA because v1.2.0 was not timed as a package.
Neither result supersedes the latest durable cross-platform per-row packet at
v1.1.3. Availability, chemistry pass/fail evidence, retained row evidence and
performance evidence remain separate claims.

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
