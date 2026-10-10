# chematic benchmark records

This directory contains dated, reproducible measurement records. Numbers are
scoped to each record's source or package, corpus, host, runtime and operation.
The current release is **v1.0.41**. Its new RDKit/COSMolKit 0.5 results and the
v1.0.40 seeded ETKDG results are source-build evidence; published-package
reruns remain separate. Older records stay scoped to their stated artifacts,
comparators and hosts.

## Start here

| Need | Start with |
|---|---|
| Review the RDKit/COSMolKit 0.5 comparison | [`2026-10-10-cosmolkit-050.md`](2026-10-10-cosmolkit-050.md): 161 operations on ChEMBL 5k and RDKit.js 10k, a 4,072-row unusual-SMILES stress corpus, and equal-work timings; source-built evidence against RDKit 2026.03.1 and COSMolKit 0.5.0rc15, not a published-package claim |
| Review the output-preserving performance candidate | [`2026-10-10-rdkit-cosmolkit-performance-candidate.md`](2026-10-10-rdkit-cosmolkit-performance-candidate.md) and [raw JSON](2026-10-10-rdkit-cosmolkit-pipeline-candidate.json): clean-commit arm64 wheel, 21 rotating-order blocks, paired confidence intervals, direct ChEMBL/NCI ring and Morgan agreement; all seven parse-inclusive pipelines win, while five prepared-operation lanes remain open |
| Review the second RDKit/COSMolKit parity round and seeded ETKDG work | [`2026-10-09-cosmolkit-parity-2.md`](2026-10-09-cosmolkit-parity-2.md): source-built comparison against RDKit 2026.03.1 and COSMolKit 0.3.0; seeded ETKDG coordinates match RDKit on 4,977/4,977 ChEMBL and 9,947/9,947 RDKit.js rows where RDKit succeeds on the recorded Linux x86-64 host; broad writer, reader, alignment, MolHash and non-tetrahedral checks; not yet published-package evidence |
| Review the v1.0.38 RDKit/COSMolKit parity work | [`2026-10-08-cosmolkit-parity.md`](2026-10-08-cosmolkit-parity.md): named RDKit-compatible SMILES, InChI, fingerprint, stereoisomer, 2D-coordinate, MOL-block and force-field APIs; exact corpus denominators and remaining ETKDG/SMARTS-writer/PDB-writer boundaries; source-candidate evidence, not a published-package speed claim |
| Review the #754/#734 batch-24 follow-ups | [`2026-10-08-754-734-followups-batch24.md`](2026-10-08-754-734-followups-batch24.md): published v1.0.38 baseline (SMARTS 310,000 exact); property-based SMIRKS fuzz lane (implicit-H differences 125 → 28, refusals 52 → 18 on 2,000 templates); SMARTS `[ca]`/`[na]`/`[cr6]` read as elements fixed; element-change charge and hypervalent ring S follow RDKit; comparator platform policy; upstream probes; the 8 autoconf differences decided |
| Review #769 on the published packages and #779 | [`2026-10-09-issues-769-779.md`](2026-10-09-issues-769-779.md): the #769 valence check confirmed on published PyPI, npm and crates.io v1.0.38; opt-in `chematic-mcp --transport streamable-http` (2026-07-28 stateless dialect at `/mcp`, loopback by default, same protocol/server/tool layers as stdio); official conformance suite 0.2.0-alpha.12 run on all 40 2026-07-28 server scenarios: 82 checks pass, the 41 failures all need suite fixtures or unimplemented capabilities, and 3 real defects it found were fixed |
| Review #769 and #739 | [`2026-10-08-issues-769-739.md`](2026-10-08-issues-769-739.md): over-valent neutral O/F in SMILES is a typed `InvalidValence` error (SMIRKS templates exempt); A6 rows 53/246 reproduced on published v1.0.31 and on the source with glibc math; the `libm` change had only moved the penam inversion to other seeds (v1.0.37: 32 of 2,650 seeded runs fail); a stereo-safe run that ends in a final stereo violation now re-embeds (0 / 2,650); v1.0.37 A6 coordinates bit-identical on Linux aarch64 (qemu) and Windows (Wine); A6 gated on every CI platform |
| Review the #734/#754 batch-23 follow-ups | [`2026-10-08-734-754-followups-batch23.md`](2026-10-08-734-754-followups-batch23.md): octahedral metal complexes drawn with their pyridines 1.5 bonds out; rows with a clash as drawn 23 → 5 of 15,000 (RDKit 276), crossings 230 → 213 |
| Review the #734/#754 batch-22 follow-ups | [`2026-10-08-734-754-followups-batch22.md`](2026-10-08-734-754-followups-batch22.md): cheaper clash relief with byte-identical MOL blocks (652M → 589M instructions); long peptides' forks opened (narrow branch atoms 1,757 → 916 of 15,000 rows); two crossings fixed; remaining decided items restated |
| Review the #734/#754 batch-21 follow-ups | [`2026-10-08-734-754-followups-batch21.md`](2026-10-08-734-754-followups-batch21.md): Python parsers raise `ChematicInputError` (category / code / format); abstention matrix across Python and WASM; nucleic-acid edit sequences replayed identically in Rust, Python and WASM; published v1.0.37 depiction 4–13× faster than RDKit |
| Review the #734/#754 batch-20 follow-ups | [`2026-10-08-734-754-followups-batch20.md`](2026-10-08-734-754-followups-batch20.md): published-wheel chemistry gate workflow repaired (no RDKit on macOS/Windows, Linux compare with a pinned RDKit wheel, rule tables fetched by path); published v1.0.37 Linux wheel / npm / crate / sdist gates and A6; porphyrin template, straight CF2 chains, exterior gaps at ring-fusion atoms (clean 14,767 of 15,000); 3D quality tier: lowest of 10 stereo-safe MMFF94 conformers p90 3.69 kcal/mol above RDKit's best of 10 (single conformer 8.65) |
| Review the #734/#754 batch-19 follow-ups | [`2026-10-08-734-754-followups-batch19.md`](2026-10-08-734-754-followups-batch19.md): crossings as drawn 277 → 242 of 15,000 (RDKit 402), clashes 28 → 25 (RDKit 276), clean 14,757; adamantane-type cage template; wider turns for small substituents; 90° forks; no stacked atoms from the relief; SVG 23% and canonical SMILES 7% fewer instructions with identical output; fullerenes and long-peptide forks decided |
| Review the #734/#754 batch-18 follow-ups | [`2026-10-07-734-754-followups-batch18.md`](2026-10-07-734-754-followups-batch18.md): CI jobs run locally (binding inventory and MMFF94 source-wheel gate would have failed; refreshed), PR perf gate equal to main in instructions, `origin/main` merged with every gate unchanged, gem-disubstituted ring atoms (crossings 300 → 277 of 15,000; RDKit 402), canonical worst-case audit |
| Review the #734/#754 batch-17 follow-ups | [`2026-10-07-734-754-followups-batch17.md`](2026-10-07-734-754-followups-batch17.md): rows with a clash as drawn 53 → 26 of 15,000 (RDKit 276), clean 14,699; turns about ring double-bond atoms; stress majorization for overlapping fused/spiro systems; depiction 5–12× faster than RDKit; RDKit 2026.09.1 Python still not on PyPI (stopped early for the release) |
| Review the #734/#754 batch-16 follow-ups | [`2026-10-07-734-754-followups-batch16.md`](2026-10-07-734-754-followups-batch16.md): branch points drawn 120° apart (v1.0.36 drew every acyclic three-bond atom 150°/150°/60°); rows with a clash as drawn 155 → 53 of 15,000 (RDKit 276), clean 14,669 (RDKit 14,504); terminal wedge neighbours no longer stacked on ring atoms; Python `to_mol_block` writes 2D coordinates for every molecule; reaction SVGs and similarity maps draw stereo; canonical E/Z pruning: explicit-H rows 4111/4384 60/40 ms → 2.3/0.4 ms, 30,305 audit strings unchanged; A6 row 0036 analysed (decided) |
| Review the #734/#754 batch-15 follow-ups | [`2026-10-07-734-754-followups-batch15.md`](2026-10-07-734-754-followups-batch15.md): SVG/depict data draw declared stereo (RDKit reads the drawings as the input on 1,686/1,687 and 1,670/1,670 stereo rows; PyPI v1.0.36 286 and 317); SDWriter Å coordinates; batch-14 `write_mol` slowdown 6,960M → 625M instructions (batch 13 353M); bridged layout clashes 212 → 5 of 531; A6 264/265 converged (stalled minimization re-embeds); BTMR1063 pair 162.6 s → 1.83 s; equal-output reactions 0.72x v1.0.34 |
| Review the #734/#754 batch-14 follow-ups | [`2026-10-07-734-754-followups-batch14.md`](2026-10-07-734-754-followups-batch14.md): 2D layout rows without clash/crossing 10,619 → 14,592 of 15,000 (RDKit 14,504; bridged 48 → 245 of 531); 3D default 1,000 iterations (A6 263/265 converged); A6 conformers vs RDKit's best of 10 ETKDGv3 (155/265 within 1 kcal/mol, RDKit's first 125); canonical E/Z worst row 63.5 → 6.3 ms (corpus −40%); quinuclidinone C6 CIP as RDKit; SMARTS `^n` −17%, Kekulé `run_reactants` −10% |
| Review the #734/#754 batch-13 follow-ups | [`2026-10-07-734-754-followups-batch13.md`](2026-10-07-734-754-followups-batch13.md): MMFF94 per-term energies = RDKit on all 14,684 single-fragment rows of both corpora; opt-in RDKit interfragment rule (102/102 salts); RDKit BFGS minimizer ported (energies after 200 iterations within 2e-13 median; A6 converged 100 → 135, 265/265 sound); pseudo-asymmetric lone-pair CIP; canonical SMILES on rigid symmetric graphs seconds → ms; MMFF typing −22%/−0.6% vs v1.0.36; RDKit 2026.09.1 `[R<n>]` profile 310,000/310,000 |
| Review the #734/#754 batch-12 follow-ups | [`2026-10-06-734-754-followups-batch12.md`](2026-10-06-734-754-followups-batch12.md): published v1.0.36 Linux aarch64 wheel gives every x86-64 gate count; MMFF94 per-term energies within 1e-12 kcal/mol of RDKit on 262 rows (max 0.32 before); A6 265/265 sound; bridgehead amine CIP = RDKit's reading (520/520 spellings, 0 abstentions); Python `StereoLossWarning`; RDKit 2026.09.1 `[R<n>]` = relevant cycles; E/Z canonical ranks 572/572, 518/518, 1,090/1,090; MMFF typing −11% (aromatic input); canonical SMILES automorphism search 34 ms → 24 µs on a product order |
| Review the #734/#754 batch-11 follow-ups | [`2026-10-06-734-754-followups-batch11.md`](2026-10-06-734-754-followups-batch11.md): RDKit's canonical Kekulé structure ported (identical on all 12,594 aromatic rows); MMFF94 heavy atoms typed differently from RDKit 81/5/39 → 0/0/0 (exposed 10k, ChEMBL 5k, Kekulé-written 15k); RDKit ring order 8,830/8,830; `^n` 135,000/135,000; bryostatin ring E/Z drawn (MOL 1,670/1,670); accurate-mode imine E/Z = rdCIPLabeler; macOS/Windows published-wheel gates packaged for CI |
| Review the v1.0.36 published reruns and follow-ups | [`2026-10-06-v1036-published-reruns-and-followups.md`](2026-10-06-v1036-published-reruns-and-followups.md): published PyPI wheel, sdist, npm and crate give the source's BioTransformer corpus rows, 83-fixture 80+3, autoconf 75/83, 310k SMARTS, CIP/hybridization/MMFF/MOL numbers and WASM Node tests; LRU template cache; `write_mol` 12% fewer instructions; RDKit.js 2026.03.6 vs 2026.09.1 reactions: 8,627/8,706 rows the same, every change a MOL read-back stereo difference |
| Review the #734/#754 batch-10 follow-ups | [`2026-10-06-734-754-followups-batch10.md`](2026-10-06-734-754-followups-batch10.md): BTMR1031 typed (RDKit's single-sanitize products equal chematic's), cyclophosphazene P labelled as RDKit (ChEMBL 5k 4,186/4,186), RDKit ring-order port (8,825/8,830 and 4,929/4,929 identical) with MMFF typing differences 88 → 81 and 53 → 5, hybridization 220,015/220,015, MOL writer 1,686/1,687 and 1,669/1,670, reactions 0.94x v1.0.34 instructions |
| Review the #734/#754 batch-9 follow-ups | [`2026-10-06-734-754-followups-batch9.md`](2026-10-06-734-754-followups-batch9.md): reaction products independent of hash-map order (BTMR1032), ring sulfoxide CIP labels (bridgehead amines stay unlabelled: RDKit's reading changes with the spelling), hybridization 220,014/220,015 and 138,655/138,655 atoms, `^n` 134,999/135,000 cells, MMFF typing differences 89 → 88 and 103 → 53, `write_mol` 5% fewer instructions |
| Review the second BioTransformer follow-up batch | [`2026-10-05-biotransformer-corpus-followups-2.md`](2026-10-05-biotransformer-corpus-followups-2.md): source gives RDKit 2026.03.6's product sets for 1,051/1,051 implicit-H and 6,604/6,652 explicit-H pairs (macrolide E/Z rows closed); RDKit-written MOL stereo blocks read back 1,678/1,687 exposed and 1,663/1,670 ChEMBL-5k (squaraine writer fix, single-wedge parity); imine E/Z 1,949/1,949 vs RDKit CIPLabeler (ring-fusion N kekulization); CIP abstentions classified; corpus reaction loop 15.7 s → 1.15 s |
| Review the BioTransformer corpus follow-ups | [`2026-10-05-biotransformer-corpus-followups.md`](2026-10-05-biotransformer-corpus-followups.md): source gives RDKit 2026.03.6's product sets for 1,051/1,051 implicit-H and 6,601/6,652 explicit-H pairs, no refusals; remaining explicit-H classes are implicit-H equivalence, RDKit truncation, 3 macrocycle E/Z rows read through RDKit `RemoveHs` and 1 naphthalene kekulization; canonical writer fix for atom2-anchored direction stashes (92/272 changed, now 0/572) |
| Review the BioTransformer rule corpus and the published v1.0.35 reruns | [`2026-10-05-biotransformer-rule-corpus.md`](2026-10-05-biotransformer-rule-corpus.md): 974 public BioTransformer rules × 400 molecules vs RDKit 2026.03.6; where RDKit has a product, same product sets for 1,045/1,051 implicit-H pairs (v1.0.35: 775) and 6,453/6,652 explicit-H pairs (v1.0.35: 3,807), rules refused everywhere 241 → 5, remaining differences classified; published v1.0.35 PyPI/npm/crate 80 + 3 on the 83 reaction fixtures, sdist autoconf 75/83 |
| Review MOL stereo loss reporting, write_mol speed and the clean-commit A6 rerun | [`2026-10-05-mol-stereo-loss-and-clean-a6.md`](2026-10-05-mol-stereo-loss-and-clean-a6.md): RDKit reads 1,681/1,687 exposed and 1,663/1,670 ChEMBL-5k stereo rows back (macrocycle E/Z, bridge layouts, fused-ring layout fix), `MolStereoLoss` report and strict mode in Rust/Python/WASM, `write_mol` 77.9 ms vs 8.2 ms for v1.0.34 (9.5x slower; draft 190.1 ms), A6 265/265 from a clean commit; replaces those three rows of the record below |
| Review the MOL stereo writer, CIP, MMFF typing and hybridization follow-ups | [`2026-10-05-mol-stereo-cip-mmff-followups.md`](2026-10-05-mol-stereo-cip-mmff-followups.md): coordinate-and-wedge MOL writer (RDKit reads 1,581/1,687 exposed stereo rows back, none inverted), conjugated E/Z in the MOL reader, RDKit-model `hybridization_per_atom` (220,003/220,015 atoms), MMFF typing 451→89 differing heavy atoms, acyclic phosphorus and lone-pair CIP labels, published PyPI v1.0.34 310k opt-in gate, A6 external scorer, branch-vs-v1.0.34 timings |
| Review the stereo-integrity, SMARTS, A6 and CIP follow-ups | [`2026-10-04-stereo-integrity-smarts-a6-followups.md`](2026-10-04-stereo-integrity-smarts-a6-followups.md): order-invariance audit fixes (H round trip, `reionize`, `uncharge`, tautomer iminols, MOL aromaticity), RDKit hybridization for `^n`, independent ChEMBL 4,625-row SMARTS rerun (opt-in 143,375/143,375 exact on PyPI v1.0.34 and source), the 18 refusals shown to be RDKit atom-order dependent, A6 Linux 265/265 with `libm`, CIP label stability |
| Review the #734/#754 SMARTS/SMIRKS dialect audit | [`2026-10-04-xsmarts-autoconf-v1034.md`](2026-10-04-xsmarts-autoconf-v1034.md): xsmarts-autoconf 83 flags vs RDKit 2026.03.6 (v1.0.34 sdist build 65, unreleased source 75, the other 8 decided), sweep and tweak findings, the published v1.0.34 PyPI Linux wheel / npm / crates.io 83-row reaction gate (80 exact, three jointly invalid), npm formula and E/Z JSON, and the `random_smiles`/`write` stereo fix |
| Review the v1.0.30 published-artifact rerun | [`2026-10-02-v1.0.30-published-artifact-gates.md`](2026-10-02-v1.0.30-published-artifact-gates.md): Python/npm/Rust 10k chemistry and 310k SMARTS, 63-operation output differential and 83 reaction strata; P0.1 accounting complete, strict parity and P0.2/P1 still open |
| Review the current-source reaction 83 graph/origin/map audit | [`2026-10-03-source-reaction-83-checked-graph-origin-map.md`](2026-10-03-source-reaction-83-checked-graph-origin-map.md): 76 exact, three typed unsupported, one refusal, three jointly invalid on an opt-in source profile; not published-package evidence |
| Review the opt-in Python SMARTS 310k source-wheel gate | [`2026-10-03-python-source-smarts-optin-310k.md`](2026-10-03-python-source-smarts-optin-310k.md): 309,982 exact match sets, 18 typed unsupported, zero wrong confident; locally built wheel, not a published-package result |
| Review post-merge reaction binding gate | [`2026-10-03-reaction-83-postmerge-binding-gate.md`](2026-10-03-reaction-83-postmerge-binding-gate.md): same 83-row Rust outcome, full WASM graph/status test, Linux/macOS source-wheel CI results, and 3/3 new-product-map supplement matches on current Rust source; not published-package evidence |
| Review checked Python reaction atom origins and template maps | [`2026-10-03-reaction-83-python-provenance-source.md`](2026-10-03-reaction-83-python-provenance-source.md): local source extension reproduces 76 exact graph/origin/map rows, three typed unsupported, one diagnosed refusal and three jointly invalid; published and CI source-wheel reruns remain open |
| Review published v1.0.31 MMFF94 typing | [`2026-10-03-a6-published-v1031-typing-census.md`](2026-10-03-a6-published-v1031-typing-census.md): 9,774 comparable molecules, 2,908 heavy-atom type differences; full type/status breakdown unchanged from v1.0.25 |
| Review the local source-wheel MMFF94 type gate | [`2026-10-03-a6-source-wheel-mmff94-typing.md`](2026-10-03-a6-source-wheel-mmff94-typing.md): 10k inputs, 451 differing heavy atoms and eight H atoms vs RDKit 2026.03.6; typing only, not published or 3D-quality parity |
| Review v1.0.31 MMFF94 400-iteration diagnostic | [`2026-10-03-a6-v1031-iteration400-diagnostic.md`](2026-10-03-a6-v1031-iteration400-diagnostic.md): same published wheel and 265 molecules; convergence 100→164 with independent sound/stereo/clash counts unchanged, but 101 remain non-converged; no default change |
| Classify the 101 A6 non-converged rows | [`2026-10-03-a6-v1031-nonconvergence-classification.md`](2026-10-03-a6-v1031-nonconvergence-classification.md) and [row-level JSON](../validation/results/a6-mmff94-v1031-nonconvergence-paired-200-400.json): 97 cap-hit, four early stops with unchanged iteration and high residual force; missing torsion parameters on six overlapping rows |
| Review the completed P0.1 artifact-audit decision | [`2026-10-03-v1030-published-p0-acceptance-policy.md`](2026-10-03-v1030-published-p0-acceptance-policy.md): all published inputs accounted for, 200 SMARTS failures and four npm API gaps retained; P0.2/P1 remain open |
| Review isolated v1.0.30 Python speed and process RSS | [`2026-10-03-v1.0.30-isolated-python-time-memory.md`](2026-10-03-v1.0.30-isolated-python-time-memory.md): 20 fresh-process paired blocks on HBA and compatible Morgan, with parse/prepared/precomputed lanes and output-gated claims |
| Review published Rust v1.0.29/v1.0.30 speed and process RSS | [`2026-10-03-v1029-v1030-published-rust-isolated-time-memory.md`](2026-10-03-v1029-v1030-published-rust-isolated-time-memory.md): six 20-block fresh-process lanes; Morgan is near parity and HBA output changed, so no equivalent-output HBA speed claim |
| Review published npm/WASM vs RDKit.js on Node | [`2026-10-03-v1029-v1030-rdkitjs-node-isolated-paired20.md`](2026-10-03-v1029-v1030-rdkitjs-node-isolated-paired20.md) and [raw JSON](2026-10-03-v1029-v1030-rdkitjs-node-isolated-paired20.json): 250 bit-identical Morgan rows, four lanes, 20 ABBA/BAAB fresh-process blocks each; Node-only, not browser speed |
| Review published npm/WASM vs RDKit.js in Chromium | [`2026-10-03-v1030-rdkitjs-published-chromium-paired20.md`](2026-10-03-v1030-rdkitjs-published-chromium-paired20.md), [raw JSON](2026-10-03-v1030-rdkitjs-published-chromium-paired20.json), and [checked summary](2026-10-03-v1030-rdkitjs-published-chromium-paired20-summary.json): 20 alternating fresh-process batch-timed blocks; fingerprint output-gated, browser/host scope only |
| Review first-use versus reused-object browser Morgan | [`2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20.md`](2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20.md), [raw JSON](2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20.json), and [checked summary](2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20-summary.json): 20 alternating fresh-process pairs, separate warm-up objects, packed Morgan output gate |
| Review published Chromium row-by-row Morgan output | [`2026-10-03-v1030-published-chromium-morgan-row-parity.md`](2026-10-03-v1030-published-chromium-morgan-row-parity.md) with [250 direct](2026-10-03-v1030-published-chromium-morgan-row-direct.json), [250 prepared](2026-10-03-v1030-published-chromium-morgan-row-prepared.json), [10k direct](2026-10-03-v1030-published-chromium-morgan-10k-direct.json), and [10k prepared](2026-10-03-v1030-published-chromium-morgan-10k-prepared.json): 250/250 exact on speed corpus, 9,999 exact plus one typed Fe refusal on independent 10k corpus |
| Review published Linux three-browser replication | [`2026-10-03-v1030-published-linux-three-browser-paired20.md`](2026-10-03-v1030-published-linux-three-browser-paired20.md): Ubuntu 24.04, 20 paired blocks per engine, 250/250 direct and prepared Morgan bit parity in each browser; speed and memory claims scoped to each engine |
| Linux Chromium raw evidence | [timing](published-browser-chromium.json), [checked summary](published-browser-chromium-summary.json), [direct bits](published-browser-chromium-row-direct.json), [prepared bits](published-browser-chromium-row-prepared.json) |
| Linux Firefox raw evidence | [timing](published-browser-firefox.json), [checked summary](published-browser-firefox-summary.json), [direct bits](published-browser-firefox-row-direct.json), [prepared bits](published-browser-firefox-row-prepared.json) |
| Linux WebKit raw evidence | [timing](published-browser-webkit.json), [checked summary](published-browser-webkit-summary.json), [direct bits](published-browser-webkit-row-direct.json), [prepared bits](published-browser-webkit-row-prepared.json) |
| Raw published-Rust paired lanes | [HBA parse-inclusive](2026-10-03-rust-v1029-v1030-hba-parse-inclusive-paired20.json), [HBA first-use](2026-10-03-rust-v1029-v1030-hba-prepared-first-use-paired20.json), [HBA precomputed](2026-10-03-rust-v1029-v1030-hba-precomputed-paired20.json), [Morgan parse-inclusive](2026-10-03-rust-v1029-v1030-morgan-parse-inclusive-paired20.json), [Morgan first-use](2026-10-03-rust-v1029-v1030-morgan-prepared-first-use-paired20.json), [Morgan precomputed](2026-10-03-rust-v1029-v1030-morgan-precomputed-paired20.json) |
| Raw isolated Python paired lanes | [v29/v30 HBA](2026-10-03-v1.0.29-vs-v1.0.30-isolated-hba-parse-inclusive.json), [v29/v30 Morgan](2026-10-03-v1.0.29-vs-v1.0.30-isolated-morgan-parse-inclusive.json), [RDKit/v30 HBA](2026-10-03-v1.0.30-vs-rdkit-isolated-hba-parse-inclusive.json), [RDKit/v30 Morgan parse-inclusive](2026-10-03-v1.0.30-vs-rdkit-isolated-morgan-parse-inclusive.json), [prepared first-use](2026-10-03-v1.0.30-vs-rdkit-isolated-morgan-prepared-first-use.json), [precomputed hot-call](2026-10-03-v1.0.30-vs-rdkit-isolated-morgan-precomputed.json) |
| Review the published-crate Rust 63-operation output differential | [`../validation/results/v1.0.29-to-v1.0.30-published-rust-63op-diff.json`](../validation/results/v1.0.29-to-v1.0.30-published-rust-63op-diff.json): 63 operations and 210,410 rows per release, exact to corresponding published Python wheel; only HBA/bundle change across pinned versions; output only |
| Review the v1.0.30 published-wheel 63-operation paired timing | [`2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json`](2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json): 20 alternating blocks with output checks and paired intervals; 20 bounded full-agreement wins, perception/memory scope open |
| Review the opt-in checked reaction source profile | [`2026-10-02-reaction-checked-source-profile.md`](2026-10-02-reaction-checked-source-profile.md): 83 stratified cases with typed unsupported/refusal accounting; unpublished candidate, no product-map provenance parity |
| Review source SMARTS residual classification | [`2026-10-03-source-smarts-310k-profiles.md`](2026-10-03-source-smarts-310k-profiles.md): published 200 match-set vs 43 Boolean residuals reproduced; opt-in ring profile fixes 154 match sets and 37 Booleans, with 29/6 wrong-confident and 18 typed refusals; source-only |
| Review the opt-in `[kN]` follow-up | [`2026-10-03-source-smarts-ring-size-follow-up.md`](2026-10-03-source-smarts-ring-size-follow-up.md): bounded symmetrized ring-size evaluation fixes ten more match sets; 19 wrong-confident match sets, six Boolean residuals and 18 typed refusals remain; source-only |
| Review mixed-size bridge-ring fallback | [`2026-10-03-source-smarts-mixed-ring-fallback.md`](2026-10-03-source-smarts-mixed-ring-fallback.md): six more opt-in match sets and one Boolean corrected; 13/5 wrong-confident and 18 typed refusals remain; source-only |
| Review `[k6]` mixed-size ring-set follow-up | [`2026-10-03-source-smarts-k6-mixed-size.md`](2026-10-03-source-smarts-k6-mixed-size.md): one more opt-in match set corrected; 12/5 wrong-confident and 18 typed refusals remain; source-only |
| Review dense-cage ring-count follow-up | [`2026-10-03-source-smarts-dense-cage-rings.md`](2026-10-03-source-smarts-dense-cage-rings.md): three fullerene-like rows corrected; six match-set and two Boolean wrong-confident cells remain on one Fe-containing row, plus 18 typed refusals; source-only |
| Review Fe organometallic ring-view follow-up | [`2026-10-03-source-smarts-organometallic-ring-view.md`](2026-10-03-source-smarts-organometallic-ring-view.md): all original 43 Boolean residuals and 183/200 match-set residuals corrected; zero wrong-confident cells, 18 typed refusals; source-only |
| Classify the 18 opt-in ring-count refusals | [`2026-10-03-smarts-ring-refusal-census.md`](2026-10-03-smarts-ring-refusal-census.md): six doubly charged polycycles, 17 formerly wrong match sets and one formerly correct set; no parity claim |
| Diagnose refused ring families | [`2026-10-03-smarts-ring-family-diagnostic.md`](2026-10-03-smarts-ring-family-diagnostic.md): six actual ring-set comparisons plus 1,000 fixed-seed RDKit atom permutations each; three source ring families differ from original RDKit and five RDKit cases are order-dependent |
| Review published-crate reaction atom origins and template maps | [`2026-10-02-v1.0.30-published-artifact-gates.md#rust-product-atom-provenance`](2026-10-02-v1.0.30-published-artifact-gates.md#rust-product-atom-provenance): 83 exposed rows, 74 graph-and-origin matches; of these 73 also match template-map identity and one styrene row is map-only different |
| Review newly created product-atom map evidence | [`2026-10-03-v1.0.30-reaction-new-product-maps.md`](2026-10-03-v1.0.30-reaction-new-product-maps.md): three supplemental cases on the pinned published Rust crate and RDKit, original 83 rows unchanged |
| Review distinct products vs raw reaction matches | [`2026-10-03-v1.0.30-reaction-multiplicity.md`](2026-10-03-v1.0.30-reaction-multiplicity.md): all 83 published Python/Rust/npm rows classified against RDKit; fewer raw embeddings alone are not a missing distinct product |
| Review the published v1.0.29 Python accuracy packet and unreleased HBA correction | [`2026-10-02-v1.0.29-python-accuracy.md`](2026-10-02-v1.0.29-python-accuracy.md): 10k chemistry rows, 310k SMARTS cells, 57 reactions and 5k operation matrix; published and source evidence separated |
| Review the post-v1.0.25 output-identical speed branch (`perf/speed-4`) | [`2026-09-25-perf-speed4-output-identical.md`](2026-09-25-perf-speed4-output-identical.md): `912c4b7d` vs `5269d69f` and follow-up `98f57746`; 49 operations x 46,736 molecules with eight base-only nondeterministic error-text differences. Shared 2-vCPU VM source evidence, not a package claim. |
| Review the previous output-identical speed branch | [`2026-09-24-perf-speed3-output-identical.md`](2026-09-24-perf-speed3-output-identical.md) |
| Review the M4 source-candidate SSSR/cache/VF2 browser rerun | [`2026-09-24-rdkitjs-m4-sssr-cache.md`](2026-09-24-rdkitjs-m4-sssr-cache.md) |
| Understand the rules and how to report a result | [`docs/benchmark.md`](../docs/benchmark.md) |
| Review current public-package fingerprint and 3D evidence | [`2026-09-23-public-package-fingerprint-3d-v1.0.20.md`](2026-09-23-public-package-fingerprint-3d-v1.0.20.md) |
| Review published v1.0.26 MMFF94 stereo-safe quality (#637) | [`2026-09-26-mmff94-public-v1.0.26.md`](2026-09-26-mmff94-public-v1.0.26.md): two fixed 265-molecule wheel runs with complete raw rows and common scoring; single-host boundary |
| Review separate A6 quality gates | [`2026-10-03-a6-separated-quality-matrix.md`](2026-10-03-a6-separated-quality-matrix.md): typing, convergence, stereo, clash, same-coordinate energy and independent conformer quality kept distinct; historical artifacts do not clear v1.0.31 |
| Review published v1.0.31 MMFF94 quality rerun | [`2026-10-03-a6-published-v1031-mmff94-quality.md`](2026-10-03-a6-published-v1031-mmff94-quality.md): current wheel 265/265 independently sound, stereo-clean and clash-free; only 100/265 converged, full A6 remains open |
| Review MMFF94 source-wheel termination diagnostic | [`2026-10-03-a6-source-mmff94-termination.md`](2026-10-03-a6-source-mmff94-termination.md): typed optimizer outcomes on a local dev-profile 265-row candidate; three typed failures require optimized-wheel confirmation before adoption |
| Review source MMFF94 aromatic-state propagation | [`2026-10-03-a6-mmff94-aromatic-propagation-source.md`](2026-10-03-a6-mmff94-aromatic-propagation-source.md): paired 10k source typing differences 2,964→459; complex rings and published-package A6 gates remain open |
| Review the post-freeze MMFF94 stereo-safe source speed gate | [`2026-09-23-mmff94-stereo-safe-performance.md`](2026-09-23-mmff94-stereo-safe-performance.md) |
| Review MMFF94 per-term same-coordinate energy (#637) | [`2026-09-25-mmff94-per-term-energy.md`](2026-09-25-mmff94-per-term-energy.md) |
| Historical v1.0.19-source MMFF94 same-coordinate energy | [`2026-09-23-mmff94-current-source-energy.md`](2026-09-23-mmff94-current-source-energy.md) |
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

The current release line is v1.0.39. Older records remain versioned historical
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
| [`2026-10-02-v1.0.30-published-artifact-gates.md`](2026-10-02-v1.0.30-published-artifact-gates.md) | Published v1.0.30 Python/npm/Rust accuracy packet and bounded paired-speed diagnostics; P0/P1 still open |
| [`2026-10-04-rdkit-2026-09-1-npm-rebaseline.md`](2026-10-04-rdkit-2026-09-1-npm-rebaseline.md) | Official RDKit.js 2026.03.6 → 2026.09.1 on fixed exposed 10k: 9,999/9,999 CheMatic Morgan bits; graph-checked CIP has 9,988 exact, five abstentions, six pre-existing E/Z mismatches, one unproven row in both versions; published CheMatic WASM SMARTS has 194 → 195 mismatches and 31 index-unproven cells per 310k lane; 12 old/new SMARTS cells changed; 20-run browser records; Python 2026.09.1 remains pending |
| [`2026-10-04-rdkit-2026-09-1-python-baseline.md`](2026-10-04-rdkit-2026-09-1-python-baseline.md) | Published Python 2026.03.6 baseline on the same exposed 10k: 9,989 CIP exact + five abstentions and six E/Z mismatches after correcting the cis/trans comparator, 9,999 Morgan exact + one typed refusal, 200/310k SMARTS mismatches; seven Python boundary timing samples; 2026.09.1 Python wheel still unavailable |
| [`2026-10-04-rdkit-2026-09-1-native-source-rebaseline.md`](2026-10-04-rdkit-2026-09-1-native-source-rebaseline.md) | Pinned independent C++ source builds on the same exposed 10k: old/new canonical, CIP and Morgan unchanged; 12 `[R2]`/`[R3]` SMARTS cell changes across six rows exactly match npm/WASM; old C++ matches the old distributed Python oracle on canonical/CIP/Morgan; not a distributed binary or new Python binding |
| [`2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json`](2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json) | Published v1.0.30 Python versus RDKit 2026.03.6 on 63 operations, 20 counterbalanced repeats with raw blocks and 95% paired intervals; not an end-to-end or cross-host claim |
| [`2026-10-02-reaction-checked-source-profile.md`](2026-10-02-reaction-checked-source-profile.md) | Opt-in unpublished reaction profile: 76/83 semantic matches, 3 typed unsupported, 1 typed refusal, 3 jointly invalid; zero wrong-confident rows in bounded corpus |
| [`2026-10-02-v1.0.30-python-op-matrix.json`](2026-10-02-v1.0.30-python-op-matrix.json) | Published v1.0.30 Python 5k/63-operation three-repeat diagnostic (not counterbalanced) |
| [`2026-10-02-v1.0.29-vs-v1.0.30-paired-parse.json`](2026-10-02-v1.0.29-vs-v1.0.30-paired-parse.json), [`paired-morgan`](2026-10-02-v1.0.29-vs-v1.0.30-paired-morgan.json), [`paired-hba`](2026-10-02-v1.0.29-vs-v1.0.30-paired-hba.json) | Published-wheel version comparison, 20 alternating blocks each; HBA outputs differ |
| [`Morgan parse-inclusive`](2026-10-02-v1.0.29-vs-v1.0.30-paired-morgan-parse-inclusive.json), [`Morgan precomputed`](2026-10-02-v1.0.29-vs-v1.0.30-paired-morgan-precomputed.json), [`ring first-use`](2026-10-02-v1.0.29-vs-v1.0.30-paired-ring-prepared.json), [`ring parse-inclusive`](2026-10-02-v1.0.29-vs-v1.0.30-paired-ring-parse-inclusive.json) | Additional 20-block published-wheel mode-separated runs; all intervals cross or approach parity |
| [`2026-10-02-v1.0.30-vs-rdkit-paired-parse.json`](2026-10-02-v1.0.30-vs-rdkit-paired-parse.json), [`paired-morgan`](2026-10-02-v1.0.30-vs-rdkit-paired-morgan.json), [`paired-hba`](2026-10-02-v1.0.30-vs-rdkit-paired-hba.json) | Published v1.0.30 versus RDKit 2026.03.6, 20 alternating blocks each; parse work differs |
| [`2026-10-02-v1.0.30-vs-rdkit-paired-ring-lanes.json`](2026-10-02-v1.0.30-vs-rdkit-paired-ring-lanes.json) | Four 20-block ring/perception lanes; output mismatch means no equivalent-output win |
| [`2026-10-02-v1.0.29-python-accuracy.md`](2026-10-02-v1.0.29-python-accuracy.md) | Published Python wheel vs RDKit 2026.03.6: exposed 10k/310k chemistry packet, 57 reaction fixtures, HBA residual and source-only correction; not a cross-binding or broad speed claim |
| [`2026-10-02-v1.0.29-python-op-matrix.json`](2026-10-02-v1.0.29-python-op-matrix.json) | Raw Apple M4 published-wheel 5k/63-operation diagnostic with per-repeat values, artifact hashes and output agreement; timing is not counterbalanced |
| [`2026-09-24-perf-speed3-output-identical.md`](2026-09-24-perf-speed3-output-identical.md) | v1.0.23 `6f8645ef` vs perf branch `6bd3fae0`: byte-identical outputs (29 operations x 46,736 molecules, 0 differing rows), paired Rust cold timing, and the Python matrix vs RDKit 2026.03.6 re-run on both revisions; RDKit agreement unchanged. 2-vCPU VM, not a package claim |
| [`2026-09-24-perf-digest-diff-v1.0.23-20b2cbc3.json`](2026-09-24-perf-digest-diff-v1.0.23-20b2cbc3.json) | Output-identity digest summary (all operations) for the record above |
| [`2026-09-24-perf-digest-diff-v1.0.23-6bd3fae0-smarts.json`](2026-09-24-perf-digest-diff-v1.0.23-6bd3fae0-smarts.json) | Output-identity digest summary (SMARTS-affected operations) after the warm-up fix |
| [`2026-09-24-perf-digest-time-pair-v1.0.23-6bd3fae0.json`](2026-09-24-perf-digest-time-pair-v1.0.23-6bd3fae0.json) | Paired Rust cold timings for the record above |
| [`2026-09-24-python-op-matrix-vs-rdkit-v1.0.23-6f8645ef.json`](2026-09-24-python-op-matrix-vs-rdkit-v1.0.23-6f8645ef.json) | Python matrix, v1.0.23 wheel |
| [`2026-09-24-python-op-matrix-vs-rdkit-branch-20b2cbc3.json`](2026-09-24-python-op-matrix-vs-rdkit-branch-20b2cbc3.json) | Python matrix, superseded branch wheel `20b2cbc3` (single-query `has_substructure` regression, fixed in `6bd3fae0`) |
| [`2026-09-24-python-op-matrix-vs-rdkit-branch-6bd3fae0.json`](2026-09-24-python-op-matrix-vs-rdkit-branch-6bd3fae0.json) | Python matrix, branch wheel `6bd3fae0` |
| [`2026-09-24-rdkit-agreement-branch-6bd3fae0.json`](2026-09-24-rdkit-agreement-branch-6bd3fae0.json) | RDKit agreement counts for `6bd3fae0` (identical to the accuracy-branch record) |
| [`2026-09-25-perf-digest-diff-v1.0.24-issue-650.json`](2026-09-25-perf-digest-diff-v1.0.24-issue-650.json) | Issue #650 candidate vs v1.0.24: canonical SMILES, largest fragment, standardize, InChI, CIP, and stereocenter outputs are byte-identical on seven benchmark corpora (280,416 rows, 0 differing). The new atom-order and fragment-source contracts are covered by Rust/Python tests; this digest does not measure the intentionally corrected plain-writer ring-closure spelling. |
| [`2026-09-25-perf-digest-diff-912c4b7d-5269d69f.json`](2026-09-25-perf-digest-diff-912c4b7d-5269d69f.json) | Base-to-first speed-4 output digest |
| [`2026-09-25-perf-digest-diff-912c4b7d-98f57746.json`](2026-09-25-perf-digest-diff-912c4b7d-98f57746.json) | Base-to-final speed-4 output digest |
| [`2026-09-25-perf-digest-time-pair-912c4b7d-5269d69f.json`](2026-09-25-perf-digest-time-pair-912c4b7d-5269d69f.json) | Paired Rust timing for the first speed-4 candidate |
| [`2026-09-25-perf-digest-time-pair-912c4b7d-98f57746.json`](2026-09-25-perf-digest-time-pair-912c4b7d-98f57746.json) | Paired Rust timing for the final speed-4 candidate |
| [`2026-09-25-python-op-matrix-vs-rdkit-912c4b7d.json`](2026-09-25-python-op-matrix-vs-rdkit-912c4b7d.json) | RDKit 2026.03.6 Python operation matrix for the base |
| [`2026-09-25-python-op-matrix-vs-rdkit-branch-5269d69f.json`](2026-09-25-python-op-matrix-vs-rdkit-branch-5269d69f.json) | RDKit 2026.03.6 Python operation matrix for the first speed-4 candidate |
| [`2026-09-25-python-op-matrix-vs-rdkit-branch-98f57746.json`](2026-09-25-python-op-matrix-vs-rdkit-branch-98f57746.json) | RDKit 2026.03.6 Python operation matrix for the final speed-4 candidate |
| [`2026-09-25-rdkit-agreement-912c4b7d.json`](2026-09-25-rdkit-agreement-912c4b7d.json) | RDKit agreement counts for the base |
| [`2026-09-25-rdkit-agreement-branch-5269d69f.json`](2026-09-25-rdkit-agreement-branch-5269d69f.json) | RDKit agreement counts for the first speed-4 candidate |
| [`2026-09-25-rdkit-agreement-branch-98f57746.json`](2026-09-25-rdkit-agreement-branch-98f57746.json) | RDKit agreement counts for the final speed-4 candidate |
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
