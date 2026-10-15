# chematic benchmark records

This directory stores dated benchmark evidence. Every result is scoped to the
named source or package, comparator, corpus, host, runtime and timed boundary.
The methodology and interpretation rules are in
[`docs/benchmark.md`](../docs/benchmark.md); the complete record inventory
is in [`archive-index.md`](archive-index.md).

## Current evidence

| Question | Record | Boundary |
|---|---|---|
| Does the V3000 SDF reader preserve molecular identity and stereo on the source-gated corpora? | [2026-10-15 V3000 recognition](2026-10-15-sdf-v3000-recognition.md) | RDKit 2025.09.3 generated 2D V3000 blocks that it read back unchanged. CheMatic and COSMolKit 0.5.0 each retain SMILES identity and atom count on 4,004 stress and 9,992 RDKit.js records. Source-changed and write-error rows are excluded; this is a source-built, format-specific parity result. |
| Does PDBQT parsing retain an over-width charge and reject a missing charge? | [2026-10-15 bounded Open Babel comparison](2026-10-15-openbabel-pdbqt-fields.md) · [raw record](2026-10-15-openbabel-pdbqt-fields.json) | Two local inputs against Open Babel 3.2.1. CheMatic retains the wider charge and `Cl` atom type, and refuses a missing charge. The wider charge exceeds the nominal fixed-column width; no general format claim. |
| Can scalar descriptors reuse an already sanitized model view? | [2026-10-14 sanitized-view reuse](2026-10-14-python-sanitized-view-reuse.md) · [slice 0 timings](2026-10-14-python-sanitized-view-reuse-slice0-paired20.json) · [slice 5000 timings](2026-10-14-python-sanitized-view-reuse-slice5000-paired20.json) · [source manifest](2026-10-14-python-sanitized-view-reuse-source.json) · [retention and validation](2026-10-14-python-sanitized-view-reuse-parity.json) | First-use TPSA improves 1.032x/1.068x and additional-slice ASA 1.065x over the prior source wheel in 20/20 blocks. All 27 eligible read-inclusive comparator lanes pass; initial deficits and small incremental regressions remain. Exact 14,072-row descriptor/Morgan retention, clone/edit invalidation and generic callback boundaries. One-host source wheels against RDKit 2026.09.1 and COSMolKit rc22, not general superiority. |
| Does dropping unused stereo preparation improve Python core-model callbacks? | [2026-10-14 core-model follow-up](2026-10-14-python-model-view-speed.md) · [slice 0 timings](2026-10-14-python-model-view-speed-slice0-paired20.json) · [slice 5000 timings](2026-10-14-python-model-view-speed-slice5000-paired20.json) · [source manifest](2026-10-14-python-model-view-speed-source.json) · [retention and validation](2026-10-14-python-model-view-speed-parity.json) | Exact previous source wheel, RDKit 2026.09.1 and COSMolKit rc22; two fixed unfiltered slices, seven operations and separate first-use/repeated-use/parse-inclusive boundaries. Matching outputs and all-block paired intervals gate each result. Source wheels on one host, not general superiority. |
| Does the V3000 reader improve while preserving its results? | [2026-10-14 V3000 packet](2026-10-14-openbabel-v3000-speed.md) · [slice 0 timings](2026-10-14-openbabel-v3000-speed-0-paired21.json) · [slice 5000 timings](2026-10-14-openbabel-v3000-speed-5000-paired21.json) · [controlled chains](2026-10-14-openbabel-v3000-speed-chains-paired21.json) · [source manifest](2026-10-14-openbabel-v3000-speed-source.json) · [verification](2026-10-14-openbabel-v3000-speed-verification.json) | Native MOL block APIs: diverse-slice reading improves 1.055x and 1.087x over the corrected baseline in 21/21 blocks. All six diverse Open Babel lanes remain unranked after semantic failures. Separate controlled linear chains pass all three comparator gates; prepared writing has mixed incremental results, including a 0.9% regression. One-host source evidence, not general superiority. |
| How do longer prepared scalar-call batches behave? | [scalar/stereo follow-up](2026-10-14-python-scalar-stereo.md) · [slice 0 long batches](2026-10-14-python-scalar-reuse-slice0-paired20.json) · [slice 5000 long batches](2026-10-14-python-scalar-reuse-slice5000-paired20.json) · [additional protocol manifest](2026-10-14-python-scalar-reuse-source.json) | One million retained and validated native float calls per operation, 20 rotating blocks per slice. Separate repeated-property throughput protocol added after short-batch outliers; original three-mode results are retained unchanged. |
| Does avoiding redundant first-use perception close the remaining Python gaps? | [2026-10-14 first-use follow-up](2026-10-14-python-first-use-followup.md) · [slice 0 timings](2026-10-14-python-first-use-followup-slice0-paired20.json) · [slice 5000 timings](2026-10-14-python-first-use-followup-slice5000-paired20.json) · [source manifest](2026-10-14-python-first-use-followup-source.json) · [parity and development](2026-10-14-python-first-use-followup-parity.json) | Additional-slice read-plus-ring-count improves 1.076x in 20/20 blocks; all 27 eligible parse-inclusive comparator lanes pass. Incremental results are otherwise mixed and first-use scalar/ring gaps remain. Exact descriptor/Morgan retention; broad stack-buffer experiment rejected. Source wheels, not toolkit-wide superiority. |
| Does scalar and stereo reuse close the prepared Python speed gaps? | [2026-10-14 scalar/stereo follow-up](2026-10-14-python-scalar-stereo.md) · [slice 0 timings](2026-10-14-python-scalar-stereo-slice0-paired20.json) · [slice 5000 timings](2026-10-14-python-scalar-stereo-slice5000-paired20.json) · [parity](2026-10-14-python-scalar-stereo-parity.json) · [source manifest](2026-10-14-python-scalar-stereo-source.json) | Isolated source candidate following the Morgan patch, two fixed 1,000-row slices, RDKit 2026.09.1 and COSMolKit rc22. Initial, repeated and parse-inclusive calls are separate; mismatching outputs suppress speed ranking. Not a published package, sealed holdout or toolkit-wide result. |
| Does the SDF reader follow-up pass both semantics and speed? | [2026-10-14 reader follow-up](2026-10-14-openbabel-sdf-read-followup.md) · [slice 0 timings](2026-10-14-openbabel-sdf-read-followup-slice0-paired21.json) · [slice 5000 timings](2026-10-14-openbabel-sdf-read-followup-slice5000-paired21.json) · [source manifest](2026-10-14-openbabel-sdf-read-followup-source.json) · [verification](2026-10-14-openbabel-sdf-read-followup-verification.json) · [full-corpus generation boundary](2026-10-14-openbabel-sdf-read-followup-full10k.json) | Candidate internal graphs agree on two unfiltered 1,000-row slices. All three additional-slice lanes beat Open Babel in 21/21 blocks; original Open Babel internal graph differences remain unranked. Allocation optimization has mixed incremental results. Native source evidence, one host, no general superiority claim. |
| How does diverse SDF batch throughput compare with Open Babel? | [2026-10-14 SDF source candidate](2026-10-14-openbabel-sdf-speed.md) · [21-block timings](2026-10-14-openbabel-sdf-speed-paired21.json) · [source manifest](2026-10-14-openbabel-sdf-speed-source.json) | Unfiltered 1,000-row V2000 2D SDF, macOS arm64, native loops and whole-process RSS. Rewrite semantics pass all rows; internal parse-stereo differences prevent a parser superiority claim. Source candidate, not a published package or toolkit-wide result. |
| Does compact stereo expansion and word packing close the Python speed gaps? | [2026-10-14 speed frontier](2026-10-14-python-speed-frontier.md) · [20-block timings](2026-10-14-python-speed-frontier-paired20.json) · [Morgan parity](2026-10-14-python-speed-frontier-parity.json) · [source manifest](2026-10-14-python-speed-frontier-source.json) | Source candidate on one macOS arm64 host: first-use Morgan improves 1.127–1.144x with 56,288 exact reference checks. Plain Morgan and MACCS beat both comparators; prepared chiral Morgan, TPSA/ASA and initial ring perception remain open. COSMolKit ring counts differ on 6/1,000 timing rows and are not ranked. |
| Which RDKit 2026.09.1 feature gap was closed next? | [Batch conformer best-RMSD follow-up](2026-10-14-rdkit-conformer-batch-rmsd.md) | Source candidate `1691ee5d`; Rust, Python and WASM expose same-molecule and cross-molecule batch APIs. A local macOS arm64 wheel agrees exactly with RDKit's existing three-pair batch result and within `3.4e-16` with repeated cross-molecule `GetBestRMS`. This is a focused API parity result, not a toolkit-wide claim. |
| Does public v1.1.3 reproduce its chemistry contract on hosted Linux, macOS and Windows? | [v1.1.3 published Python artifact packet](2026-10-14-v113-published-artifact-packet.md) · [machine-readable packet](../validation/results/published-wheel-chemistry-gates-v1.1.3.json) | Yes for the four exact PyPI wheels and the pinned reaction, SMARTS, Morgan, descriptor, CIP, hybridization, MMFF94, MOL stereo, BioTransformer and A6 gates. Platform outputs use the pinned Linux RDKit comparator where no matching native RDKit wheel exists. npm/WASM, crates.io, speed, and universal parity remain separate. |
| Did the v1.1.3 follow-up restore source compatibility, unusual-input parity, and ordinary-input speed? | [v1.1.3 evaluation follow-up](2026-10-13-v113-evaluation-fixes.md) · [parse](2026-10-13-v113-fix-v112-parse-paired20.json) · [TPSA](2026-10-13-v113-fix-v112-tpsa-paired20.json) · [Morgan](2026-10-13-v113-fix-v112-morgan-paired20.json) · [QED](2026-10-13-v113-fix-v113-qed-paired20.json) | Source candidate `90542885`: v1.1.2 exhaustive enum matches compile; MW/Morgan/rings are 4,072/4,072 and every finite ASA value agrees within `1e-9`; QED is 16.1% faster than public v1.1.3. Parse-inclusive v1.1.2 comparisons cross parity and therefore are ties, not speed wins. |
| Are the 50 public-v1.1.2 QED residuals on unusual SMILES closed? | [QED stress-corpus parity follow-up](2026-10-13-qed-stress-parity.md) · [result JSON](../validation/results/qed-stress-rdkit-source-2026-10-13.json) | Yes for the source candidate against RDKit 2026.03.6: 4,071/4,071 comparable rows agree, up from 4,021/4,071. RDKit itself raises on the remaining row. This is source-candidate reference parity, not a published-package or experimental-accuracy claim. |
| Does the source candidate reproduce RDKit 2026.09.1 topology bounds, and what quality tail remains? | [Bounds and conformer follow-up](2026-10-12-a6-rdkit-2026-09-bounds-followup.md) · [summary JSON](2026-10-12-a6-rdkit-2026-09-bounds-followup-summary.json) | Raw and triangle-smoothed bounds are within `1e-12` for all 768,664 cells. Equal-count quality is within 1 kcal/mol on 258/265 rows. Seven quality-tail rows and formal multi-host/public-artifact timing remain open; one loaded-host timing is not a speed claim. |
| Does shared, parallel RDKit-compatible ETKDG/MMFF94 close the A6 quality and speed gap? | [20-conformer source-candidate packet](2026-10-12-a6-multiconformer-etkdg-mmff.md) · [summary JSON](2026-10-12-a6-multiconformer-etkdg-mmff-summary.json) | Baseline 20-conformer quality is within 1 kcal/mol on 243/265 rows; bounded retry reaches 256/265. Two one-repeat baseline timings disagree on the winner and retry is slower, so this is a quality improvement and API result, not a speed-win claim. |
| Does full-hydrogen MMFF ranking improve A6 conformer quality without losing throughput? | [source candidate vs RDKit 2026.09.1](2026-10-12-a6-full-h-mmff-vs-rdkit.md) · [summary JSON](2026-10-12-a6-full-h-mmff-vs-rdkit.json) | On one macOS arm64 host and 265 molecules, equal-count 10-candidate medians favor CheMatic, 237/265 are within 1 kcal/mol, and a two-repeat same-work timing is 1.17x faster. Maximum residuals, 21-block timing, holdout, second-host, and public-artifact gates remain open. |
| What remains after the RDKit 2026.09.1 A6 bounds work? | [conformer-tail record](2026-10-13-a6-conformer-tail-development.md) · [v3 summary](2026-10-13-a6-conformer-tail-v3-summary.json) | Fixes a gradient early-exit bug and ring traversal order, and adds an optional rescue profile. On untouched v3, normal mode is effectively tied but 0.27% slower; rescue improves quality but is 2.1% slower. No simultaneous superiority claim. |
| Can the RDKit-compatible force fields avoid a redundant final energy pass? | [MMFF final-energy cache](2026-10-13-mmff-final-energy-cache.md) | Yes. All compared outputs are unchanged. The 100-row development median moves only slightly, from a 1.063x to 1.068x RDKit/CheMatic ratio on one host. |
| Does deterministic parallel conformer generation improve best-of-N throughput? | [source candidate vs RDKit 2026.09.1](2026-10-12-parallel-conformer-ensemble-vs-rdkit.md) · [CheMatic raw JSON](2026-10-12-parallel-conformer-ensemble-chematic.json) · [RDKit raw JSON](2026-10-12-parallel-conformer-ensemble-rdkit.json) | On one macOS arm64 host and three molecules, the parallel CheMatic path is exactly equal to its sequential output and 3.54x faster. Its median is 3.87x faster than RDKit's all-thread lane, but cross-engine ensemble quality is not equivalent-output-gated yet, so this remains candidate evidence. |
| Does lazy preparation remove the v1.1.1 parse-inclusive regression? | [published v1.1.1 vs source candidate](2026-10-12-v111-public-vs-lazy-fix.md) · [parse](2026-10-12-v111-public-vs-lazy-fix-parse-paired20.json) · [TPSA pipeline](2026-10-12-v111-public-vs-lazy-fix-tpsa-pipeline-paired20.json) · [Morgan pipeline](2026-10-12-v111-public-vs-lazy-fix-morgan-pipeline-paired20.json) · [TPSA precomputed](2026-10-12-v111-public-vs-lazy-fix-tpsa-precomputed-paired20.json) · [Morgan precomputed](2026-10-12-v111-public-vs-lazy-fix-morgan-precomputed-paired20.json) | Yes for the recorded macOS arm64 source candidate: parsing and both parse-inclusive lanes improve with equal outputs. Precomputed TPSA is unchanged; precomputed Morgan is about 2% slower, so it is not claimed as a win. Not a published-package or cross-platform result. |
| Are the former Labute ASA and first-use ring-count deficits closed? | [v1.1.0 source follow-up](2026-10-12-v1.1.0-source-labute-ring-speed.md) · [raw JSON](2026-10-12-v1.1.0-source-rdkit-labute-ring-speed.json) | Yes on the recorded macOS arm64 source-built wheel: both are 10,000/10,000 equivalent, with prepared medians of 1.71x and 1.22x versus RDKit 2026.09.1 and favorable paired intervals. Not a published-package or cross-platform claim. |
| How did the recorded v1.1.0 source snapshot compare with RDKit? | [v1.1.0-source / RDKit 2026.09.1 Python packet](2026-10-11-v1.1.0-main-vs-rdkit-2026.09.1.md) · [focused raw JSON](2026-10-11-v1.1.0-main-rdkit-python-macos-arm64.json) · [63-operation raw JSON](2026-10-11-v1.1.0-main-vs-rdkit-63op-paired20.json) | Source-built wheel at `44338e16`; four 10,000-row output gates plus 21 fresh-process blocks, and a 63-operation/20-block supplementary matrix on macOS arm64. Historical source evidence, not current-main or published-package evidence. |
| What is the latest published-package comparison? | [v1.0.42 / RDKit 2026.09.1](2026-10-11-v1042-published-rdkit-packet.md) | Exact PyPI/npm artifacts on hosted Linux x64 and macOS arm64. |
| What is the broad Python operation matrix? | [v1.0.30 63-operation matrix](2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json) | Historical published-wheel result against RDKit 2026.03.6; only output-equivalent lanes with favorable confidence intervals count as wins. |
| How does browser/WASM compare? | [v1.0.30 three-browser packet](2026-10-03-v1030-published-linux-three-browser-paired20.md) | Published npm artifacts, 250 exact Morgan rows and 20 blocks in Chromium, Firefox and WebKit. |
| How does similarity search compare? | [v1.0.13 similarity-search packet](2026-09-12-similarity-search-a3-v1.0.13.md) | Separates native/native, RDKit-compatible/RDKit and cross-profile overlap. |
| How does file interchange compare? | [Open Babel file-I/O program](../docs/openbabel-file-io.md) · [same-process candidate](2026-10-10-openbabel-file-io-same-process-v1.0.42.md) | Small pinned fixtures and source candidate; not broad-corpus or published-package evidence. |
| How do CDK and Indigo compare? | [bounded CDK/Indigo baseline](2026-10-10-cdk-indigo-bounded-baseline.md) | 24-row common matrix and 21 timing blocks; source-candidate smoke evidence. |
| How do reaction, parent selection, and R-group workflows compare with Indigo? | [Indigo workflow packet](2026-10-12-indigo-workflows-source.md) · [raw JSON](2026-10-12-indigo-workflows-source.json) | Source-built, one-host bounded lanes. CheMatic has higher exact counts for reaction and parent selection; R-group is tied. Common-correct medians favor CheMatic in all three lanes. This is not a toolkit-wide or published-package claim. |
| How do reaction, standardization, and labelled R-group workflows compare with RDKit? | [expanded RDKit workflow packet](2026-10-12-rdkit-workflows-expanded-source.md) · [raw JSON](2026-10-12-rdkit-workflows-expanded-source.json) | Source-built macOS arm64 lanes with 21 alternating blocks. Reaction has 80 semantic matches plus 3 jointly invalid rows; largest-fragment selection is 39/39; Normalize is 254/254 after common-input filtering; labelled R-group is 16/16. Every equivalent-output timing lane favors CheMatic, including prepared reaction and R-group lanes. Broader MolStandardize residuals remain explicitly bounded. |
| How does COSMolKit compare? | [RDKit/COSMolKit 0.5 comparison](2026-10-10-cosmolkit-050.md) | Source-built 161-operation and unusual-SMILES evidence; not a published-package claim. |
| What feature gap was closed after COSMolKit 0.5? | [Persistent force-field follow-up](2026-10-13-cosmolkit-persistent-force-field.md) | Python source-candidate API and one-host exploratory reuse timing. Exact energy tests pass; this is not a toolkit-wide or published-package comparison. |
| Which COSMolKit 0.5 Python interoperability gap was closed next? | [Fingerprint NumPy follow-up](2026-10-14-cosmolkit-fingerprint-numpy.md) | Source-candidate API review and packed-bit equivalence tests for dense NumPy output; no chemistry or toolkit-wide parity claim. |
| What is the current A6 MMFF distance-cache result? | [development record](2026-10-13-a6-mmff-distance-cache.md) · [once-only v4 holdout](2026-10-13-a6-mmff-distance-cache-v4.md) · [raw JSON](../validation/results/a6-mmff-distance-cache-holdout-v4-r3-2026-10-13.json) | RDKit 2026.09.1, macOS arm64. v4 is 95/100 within 1 kcal/mol and effectively tied in speed: CheMatic wins two of three blocks and has a 1.0038x median ratio. This is not a robust or published-artifact superiority claim. |
| Can a final-gradient gate fund MMFF rescue without losing the speed lead? | [gradient-rescue development record](2026-10-13-a6-mmff-gradient-rescue-development.md) · [raw JSON](../validation/results/a6-mmff-gradient-rescue-development-r3-2026-10-13.json) | On the already-inspected 265-row development set, the fixed profile reaches 261/265 within 1 kcal/mol and is faster in all three alternating blocks, with a 1.0179x median ratio. Fresh holdout, clean-artifact, 21-block, second-host, and memory gates remain open. |

## What counts as a speed win

A lane is called an equivalent-output speed win only when:

1. both engines perform the same declared work;
2. every compared output satisfies the lane's exact or numeric contract;
3. failures and unsupported rows are fully accounted for;
4. the lower bound of the paired 95% confidence interval is above `1.0`;
5. source, package, corpus, runtime and host provenance are recorded.

Parse-inclusive and prepared-input timings answer different questions. A
cached RDKit property lookup is not compared as though it included perception,
and a source-built result is not presented as a registry-package result.

## Reproduce the current Python packet

Build a wheel from the recorded revision and install it with the pinned RDKit
wheel in an isolated Python 3.13 environment. Then run:

```bash
python scripts/run_published_v1042_rdkit_python_packet.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --artifact-dir /path/to/exact-wheels \
  --output /tmp/chematic-rdkit.json \
  --expected-chematic 1.1.0 \
  --expected-rdkit 2026.9.1 \
  --schema rdkit-labute-ring-optimization/v1 \
  --blocks 41 \
  --chematic-revision f2d030fc6a2760fc9ad3f168f98713584c20828c
```

The runner checks all 10,000 outputs before emitting speed intervals and uses
41 alternating fresh-process blocks over the first 1,000 rows.

Validate the generated record with:

```bash
python scripts/check_published_v1042_rdkit_packet.py /tmp/chematic-rdkit.json \
  --expected-schema rdkit-labute-ring-optimization/v1 \
  --expected-chematic 1.1.0 \
  --expected-rdkit 2026.9.1 \
  --expected-revision f2d030fc6a2760fc9ad3f168f98713584c20828c \
  --expected-blocks 41 \
  --minimum-artifacts 1
```

## Record maintenance

- Add new claims as a Markdown summary plus machine-readable raw data.
- Keep superseded records; label them historical instead of rewriting them.
- Link every checked-in root-level `.md` and `.json` record from this file or
  [`archive-index.md`](archive-index.md).
- Validate the indexes with `python3 scripts/check_benchmark_index.py`.
- Do not infer universal accuracy or speed from one corpus, host or API lane.

The [complete inventory](archive-index.md) contains the full dated record
list, streaming/format matrices, 3D evidence, artifact audits and earlier
candidate measurements.
