# Configured Morgan borrowing and single-ring scratch — 2026-10-18

**The combined source candidate passes the declared adoption gate. The borrowing-only candidate fails and remains recorded separately.** First-use chiral Morgan improves by **1.061750600×** and **1.039879497×** in the two 1,000-row windows. Both paired 95% intervals are wholly above one; both win 24/24 blocks. All 42 source medians meet the unchanged 0.99 floor. This does not mean every operation improves: offset-5,000 repeated Labute ASA is 0.995588003×, with its entire interval below one.

The combined candidate meets the strict external criterion in **18/18 RDKit** and **17/18 COSMolKit** Morgan lanes. Offset-0 first-use chiral Morgan remains **0.875397757× COSMolKit**, with 0/24 wins. Complete toolkit superiority, all-input parity, other hosts, and published-package performance are not established.

Baseline: `c9b63c17bebe914208b97e544cf8bece3e120e94`. Borrowing-only candidate: `2c448750b5ff8bcf045d50f89117139baa453d54`. Combined candidate: `f6f83d440558641c599ce51f23db221eed8a6feb`. The integration also retains the rejected stereo-projection evidence from PR #1027; that projection's runtime changes are excluded. The source archives, native wheels, complete row outputs, build and validation logs remain local. The public packets contain only timing/accounting values and hashes, without molecular outputs, SDK payloads or local filesystem paths.

## Changes and preserved behavior

Configured Morgan uses the existing `with_rdkit_parity_view` identity proof to borrow a molecule whose RDKit aromaticity view would otherwise be an identical clone. Nonidentity molecules still use the existing shared view. Hydrogen suppression, coordination rejection, aromaticity error propagation, Morgan invariants, and the legacy-stereo/accurate-CIP fallback boundary are unchanged. The fixed ECFP4 entry point is unchanged.

For exactly one parser-recorded ring closure, the existing RDKit ring-graph reachability check stores visited bits in `SmallVec<[u64; 2]>` and pending atoms in `SmallVec<[AtomIdx; 16]>`. Both retain heap fallback. Query-bond rejection, excluded/dative bonds, neighbor order, and counting semantics are unchanged. No dependency or public API changes.

## Separate experiments and adoption

Both plans were fixed locally before their candidate timings. The second plan was fixed before the ring-scratch code change. It explicitly follows the first candidate's failure; the first observations are retained rather than replaced. A real ring-count implementation change distinguishes the second candidate. No elapsed observations are pooled between epochs and no performance-selected reruns occur.

Each epoch requires both chiral first-use source CI lower bounds > 1, all 42 unrounded source medians >= 0.99, complete journals, and exact source output agreement. The external criterion independently requires exact outputs, CI lower > 1 and 24/24 faster blocks for each lane.

The borrowing-only epoch fails offset-0 `num_rings/parse_inclusive`: **0.988262691×**, CI **[0.985367911, 1.011172424]**, 10/24 wins. Its interval crosses one; the failed median floor is not proof that borrowing caused a ring-count slowdown. Both chiral primaries pass, but that does not waive the floor.

| Offset | Operation | Mode | Comparator | Paired median | Paired 95% CI | Faster blocks | Strict external win |
| ---: | --- | --- | --- | ---: | --- | ---: | --- |
| 0 | morgan2_chiral | first_use | baseline | 1.064428492 | [1.061426111, 1.068160725] | 24/24 | — |
| 5000 | morgan2_chiral | first_use | baseline | 1.042560599 | [1.036233568, 1.046960652] | 24/24 | — |
| 0 | num_rings | parse_inclusive | baseline | 0.988262691 | [0.985367911, 1.011172424] | 10/24 | — |

| Epoch | Complete arm journals | Accounting cells | Timed full batches | Timed rows | Decision |
| --- | ---: | ---: | ---: | ---: | --- |
| Borrowing only | 192 | 2880 | 189594 | 189594000 | reject |
| Borrowing + ring scratch | 192 | 2880 | 189868 | 189868000 | pass |

The combined minimum source median is **0.991886993×** (0/labute_asa/first_use). All 42 floors pass. The one source interval wholly below one is retained here and in the full table; a small permitted regression is not described as an improvement.

| Offset | Operation | Mode | Comparator | Paired median | Paired 95% CI | Faster blocks | Strict external win |
| ---: | --- | --- | --- | ---: | --- | ---: | --- |
| 5000 | labute_asa | repeated_use | baseline | 0.995588003 | [0.991114111, 0.999436254] | 7/24 | — |
## Correctness and existing limits

The borrowing-only source check set passes 1,230 Rust tests (45 suites, zero failures, five existing ignored); the combined source check set passes **1,670 Rust tests** (43 suites, zero failures, five existing ignored). These are different package selections, not 440 new tests. Both pass workspace all-target Clippy, formatting and seven runner tests. Both isolated source wheels pass **329 selected Python tests**.

Morgan's normal 10,000 and stress 4,072 rows match baseline/candidate/RDKit exactly for R2/R3 with and without chirality: **56,288 fingerprint cells per epoch**. All 40,000 normal COSMolKit cells match; the stress corpus retains six refusals per operation and seven additional chiral mismatches at each radius. They are not filtered or counted as agreements.

The new borrowing tests compare all consumed atom/bond/ring/invariant/stereo inputs against the former owned-view boundary, including repeated calls, clones, hydrogen suppression, nonidentity/error cases, and graph/stereo edits. Ring tests keep the old full counting implementation as a reference on 10,000 rows, exhaustive five-atom graphs and edited graphs. Added single-closure rings and comb-shaped side chains cross 64/128-atom visited-word boundaries and spill the pending stack.

The combined wheel's additional **14,072 ring-count rows match baseline exactly**, including all large-graph fallbacks. Against RDKit 2026.09.1, all 10,000 normal ring counts agree; four existing stress discrepancies remain unchanged:

| Zero-based stress row | Baseline and candidate rings | RDKit rings |
| ---: | ---: | ---: |
| 605 | 22 | 70 |
| 608 | 10 | 11 |
| 1014 | 16 | 70 |
| 3774 | 13 | 14 |

The two formal windows have exact Morgan outputs across all four engines. All source TPSA/ASA/ring/MACCS values match in every timed batch. These checks establish preservation of the selected source behavior; the four ring discrepancies and COSMolKit refusals/mismatches prevent an all-input parity claim.

## Measurement boundary and provenance

One Apple M4 arm64 host (Mac16,12, 10 CPUs, 16 GiB, macOS 27.0.1), CPython 3.13.6, Rust 1.97.0; installed RDKit 2026.09.1, COSMolKit 0.5.0 and NumPy 2.5.3. These are source wheels, not measurements of a public chematic release artifact.

Both first-epoch source wheels are built from separate empty targets. The second epoch reuses the **exact hash-verified c9b63c17 baseline wheel and source archive**, and builds the changed combined candidate from a new empty target. It reuses no elapsed values. Native-inchi, opt3/full LTO/one CGU/one build job/no incremental/no stripping and the CPython interpreter are identical; only target/output paths differ. Installed native hashes match their wheels. All 1,631 archived source files per arm match build inputs. All 1,953 sealed runtime files and measured source files are unchanged after each formal epoch. No own builds, tests, profilers or source edits run during formal timing.

The public runner is [bench_python_morgan_linear.py](../scripts/bench_python_morgan_linear.py). Each epoch uses 1,000 rows at offsets 0 and 5,000 of `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`. Each profile runs every permutation of baseline/candidate/RDKit/COSMolKit once: 24 paired blocks, each arm appearing six times in each position, shuffled by seed 20261017+offset.

Operations are R2/2048 nonchiral, R2/2048 chiral and R3/2048 nonchiral Morgan; source controls are TPSA, Labute ASA, RDKit number of rings and MACCS. Each has first-use, repeated-use and parse-inclusive modes. First-use reparses fresh objects outside each timed batch and calls the requested operation once per object. Repeated-use primes once outside timing and reuses prepared objects. Parse-inclusive constructs fresh molecules inside every timed batch. Warmup uses separate objects.

Each operation accumulates complete independently timed batches for at least 50 ms (Morgan) or 20 ms (controls). It normalizes and checks every batch outside timing, preserving all raw durations, batch counts, rows and output digests. Native public return values are timed; conversion to canonical on-bits is outside timing. Whole-process RSS includes imports, preparation and validation and is not ranked. This accumulated-operation metric is not process cold start or a contiguous steady-state run.

Ratios are comparator batch time / candidate batch time paired within blocks. The unrounded median of ratios and 10,000 paired bootstrap resamples with fixed seeds decide the gates. Printed values below are rounded only for display.

## Public numeric packets

- [2026-10-18-morgan-borrowed-view-only-numeric.json.gz](2026-10-18-morgan-borrowed-view-only-numeric.json.gz): 617,489 compressed bytes; SHA-256 `e7ae50a2a302f8cf0e8a30351e980f9a013e84015c5025de6f06449387d5ccb0`.
- [2026-10-18-morgan-borrowed-view-ring-scratch-numeric.json.gz](2026-10-18-morgan-borrowed-view-ring-scratch-numeric.json.gz): 616,893 compressed bytes; SHA-256 `6ce0dac27f6603598b1216a815af8a56440a8a545c12cc63fc10442401bd0bcd`.

Recompute timing counts, every interval and each declared decision with:

```sh
python3 scripts/check_morgan_numeric_packet.py benchmarks/2026-10-18-morgan-borrowed-view-only-numeric.json.gz
python3 scripts/check_morgan_numeric_packet.py benchmarks/2026-10-18-morgan-borrowed-view-ring-scratch-numeric.json.gz
python3 -m unittest scripts.tests.test_check_morgan_numeric_packet -v
```

The checker verifies numeric consistency and hashes, not chemistry from hashes alone; full row outputs remain in the local audited journals. CI runs both new packet checks and corruption controls. The complete combined 78-lane table follows. The borrowing-only packet also retains all 78 comparisons, including its failed floor.

## Complete combined results

| Offset | Operation | Mode | Comparator | Paired median | Paired 95% CI | Faster blocks | Strict external win |
| ---: | --- | --- | --- | ---: | --- | ---: | --- |
| 0 | morgan2_2048 | first_use | baseline | 1.170419990 | [1.164507487, 1.174065594] | 24/24 | — |
| 0 | morgan2_2048 | first_use | rdkit | 3.909490167 | [3.887637357, 3.934115252] | 24/24 | yes |
| 0 | morgan2_2048 | first_use | cosmolkit | 2.455133161 | [2.443103100, 2.464622874] | 24/24 | yes |
| 0 | morgan2_2048 | repeated_use | baseline | 1.005404535 | [0.999672199, 1.011799111] | 16/24 | — |
| 0 | morgan2_2048 | repeated_use | rdkit | 5.431870108 | [5.382335646, 5.459598003] | 24/24 | yes |
| 0 | morgan2_2048 | repeated_use | cosmolkit | 3.278043969 | [3.269849176, 3.296593203] | 24/24 | yes |
| 0 | morgan2_2048 | parse_inclusive | baseline | 1.131719782 | [1.127981683, 1.139091606] | 24/24 | — |
| 0 | morgan2_2048 | parse_inclusive | rdkit | 9.765490692 | [9.717380168, 9.789949911] | 24/24 | yes |
| 0 | morgan2_2048 | parse_inclusive | cosmolkit | 8.610269476 | [8.558524476, 8.635920683] | 24/24 | yes |
| 0 | morgan2_chiral | first_use | baseline | 1.061750600 | [1.057355485, 1.065482007] | 24/24 | — |
| 0 | morgan2_chiral | first_use | rdkit | 1.375882802 | [1.368727849, 1.384726871] | 24/24 | yes |
| 0 | morgan2_chiral | first_use | cosmolkit | 0.875397757 | [0.873502150, 0.880731862] | 0/24 | no |
| 0 | morgan2_chiral | repeated_use | baseline | 0.996198942 | [0.990697619, 1.003156814] | 10/24 | — |
| 0 | morgan2_chiral | repeated_use | rdkit | 5.319931202 | [5.298827677, 5.356855176] | 24/24 | yes |
| 0 | morgan2_chiral | repeated_use | cosmolkit | 3.207670991 | [3.189718273, 3.232030408] | 24/24 | yes |
| 0 | morgan2_chiral | parse_inclusive | baseline | 1.064844216 | [1.060316466, 1.070118305] | 24/24 | — |
| 0 | morgan2_chiral | parse_inclusive | rdkit | 4.762485871 | [4.749461660, 4.780093022] | 24/24 | yes |
| 0 | morgan2_chiral | parse_inclusive | cosmolkit | 4.176595902 | [4.159281306, 4.201583490] | 24/24 | yes |
| 0 | morgan3_2048 | first_use | baseline | 1.126509231 | [1.124032306, 1.128949577] | 24/24 | — |
| 0 | morgan3_2048 | first_use | rdkit | 3.899860643 | [3.876678621, 3.929842061] | 24/24 | yes |
| 0 | morgan3_2048 | first_use | cosmolkit | 2.529846313 | [2.517501006, 2.537190068] | 24/24 | yes |
| 0 | morgan3_2048 | repeated_use | baseline | 1.004855480 | [0.996196420, 1.016101744] | 15/24 | — |
| 0 | morgan3_2048 | repeated_use | rdkit | 5.040643951 | [5.009469394, 5.089000844] | 24/24 | yes |
| 0 | morgan3_2048 | repeated_use | cosmolkit | 3.175315193 | [3.156639860, 3.185196648] | 24/24 | yes |
| 0 | morgan3_2048 | parse_inclusive | baseline | 1.108415640 | [1.106960221, 1.115212178] | 24/24 | — |
| 0 | morgan3_2048 | parse_inclusive | rdkit | 9.010255709 | [8.992324733, 9.053423805] | 24/24 | yes |
| 0 | morgan3_2048 | parse_inclusive | cosmolkit | 7.844019069 | [7.810394893, 7.876345604] | 24/24 | yes |
| 0 | tpsa | first_use | baseline | 0.996370037 | [0.991608467, 1.005656660] | 10/24 | — |
| 0 | tpsa | repeated_use | baseline | 0.999011322 | [0.996310275, 1.003154219] | 11/24 | — |
| 0 | tpsa | parse_inclusive | baseline | 1.001561050 | [0.985371337, 1.011708005] | 12/24 | — |
| 0 | labute_asa | first_use | baseline | 0.991886993 | [0.981714118, 1.001035027] | 8/24 | — |
| 0 | labute_asa | repeated_use | baseline | 0.997301029 | [0.992276084, 1.005768741] | 10/24 | — |
| 0 | labute_asa | parse_inclusive | baseline | 1.007486631 | [1.001642434, 1.013341722] | 18/24 | — |
| 0 | num_rings | first_use | baseline | 1.001589810 | [0.991352859, 1.011367616] | 12/24 | — |
| 0 | num_rings | repeated_use | baseline | 0.998244044 | [0.989760667, 1.001835091] | 9/24 | — |
| 0 | num_rings | parse_inclusive | baseline | 1.010835525 | [1.004886975, 1.016100287] | 19/24 | — |
| 0 | maccs | first_use | baseline | 1.002979513 | [1.001681377, 1.005038893] | 19/24 | — |
| 0 | maccs | repeated_use | baseline | 1.004903243 | [1.003220362, 1.006604064] | 23/24 | — |
| 0 | maccs | parse_inclusive | baseline | 1.006626723 | [1.000655332, 1.008303943] | 20/24 | — |
| 5000 | morgan2_2048 | first_use | baseline | 1.046030410 | [1.039119155, 1.048039008] | 24/24 | — |
| 5000 | morgan2_2048 | first_use | rdkit | 1.666704440 | [1.652747177, 1.685644648] | 24/24 | yes |
| 5000 | morgan2_2048 | first_use | cosmolkit | 1.039678578 | [1.035672985, 1.058355495] | 24/24 | yes |
| 5000 | morgan2_2048 | repeated_use | baseline | 1.000397904 | [0.996804165, 1.009206369] | 12/24 | — |
| 5000 | morgan2_2048 | repeated_use | rdkit | 5.092697763 | [5.053986088, 5.122166880] | 24/24 | yes |
| 5000 | morgan2_2048 | repeated_use | cosmolkit | 3.017306884 | [2.998497477, 3.034861740] | 24/24 | yes |
| 5000 | morgan2_2048 | parse_inclusive | baseline | 1.033053202 | [1.024238569, 1.039914244] | 24/24 | — |
| 5000 | morgan2_2048 | parse_inclusive | rdkit | 4.780047400 | [4.740657919, 4.828188432] | 24/24 | yes |
| 5000 | morgan2_2048 | parse_inclusive | cosmolkit | 4.618966784 | [4.588684920, 4.650012987] | 24/24 | yes |
| 5000 | morgan2_chiral | first_use | baseline | 1.039879497 | [1.033708276, 1.046253917] | 24/24 | — |
| 5000 | morgan2_chiral | first_use | rdkit | 1.690210058 | [1.676061707, 1.696050410] | 24/24 | yes |
| 5000 | morgan2_chiral | first_use | cosmolkit | 1.058131603 | [1.050327694, 1.065977590] | 24/24 | yes |
| 5000 | morgan2_chiral | repeated_use | baseline | 0.996490872 | [0.986986945, 1.007509909] | 11/24 | — |
| 5000 | morgan2_chiral | repeated_use | rdkit | 5.073988637 | [4.995360404, 5.104268754] | 24/24 | yes |
| 5000 | morgan2_chiral | repeated_use | cosmolkit | 2.995597172 | [2.989201958, 3.008614694] | 24/24 | yes |
| 5000 | morgan2_chiral | parse_inclusive | baseline | 1.031433351 | [1.029478486, 1.048631026] | 24/24 | — |
| 5000 | morgan2_chiral | parse_inclusive | rdkit | 4.795942846 | [4.768851265, 4.837577686] | 24/24 | yes |
| 5000 | morgan2_chiral | parse_inclusive | cosmolkit | 4.638650006 | [4.603065610, 4.653753367] | 24/24 | yes |
| 5000 | morgan3_2048 | first_use | baseline | 1.037427979 | [1.032771265, 1.044841823] | 24/24 | — |
| 5000 | morgan3_2048 | first_use | rdkit | 1.883156270 | [1.877515069, 1.891739382] | 24/24 | yes |
| 5000 | morgan3_2048 | first_use | cosmolkit | 1.220381422 | [1.213502492, 1.226258382] | 24/24 | yes |
| 5000 | morgan3_2048 | repeated_use | baseline | 0.998969217 | [0.994549627, 1.004450610] | 9/24 | — |
| 5000 | morgan3_2048 | repeated_use | rdkit | 4.795227412 | [4.751867978, 4.850804993] | 24/24 | yes |
| 5000 | morgan3_2048 | repeated_use | cosmolkit | 3.015805354 | [2.975022758, 3.034177604] | 24/24 | yes |
| 5000 | morgan3_2048 | parse_inclusive | baseline | 1.034619549 | [1.029126737, 1.040271602] | 24/24 | — |
| 5000 | morgan3_2048 | parse_inclusive | rdkit | 4.712128537 | [4.688964357, 4.744537547] | 24/24 | yes |
| 5000 | morgan3_2048 | parse_inclusive | cosmolkit | 4.491568686 | [4.477888463, 4.507644681] | 24/24 | yes |
| 5000 | tpsa | first_use | baseline | 1.003469500 | [0.998265666, 1.010625984] | 15/24 | — |
| 5000 | tpsa | repeated_use | baseline | 1.003085554 | [0.995089511, 1.015561647] | 14/24 | — |
| 5000 | tpsa | parse_inclusive | baseline | 0.999978268 | [0.992520620, 1.021217210] | 12/24 | — |
| 5000 | labute_asa | first_use | baseline | 0.999778496 | [0.990127911, 1.005276711] | 11/24 | — |
| 5000 | labute_asa | repeated_use | baseline | 0.995588003 | [0.991114111, 0.999436254] | 7/24 | — |
| 5000 | labute_asa | parse_inclusive | baseline | 1.005076749 | [0.993522311, 1.011412514] | 14/24 | — |
| 5000 | num_rings | first_use | baseline | 1.003828973 | [0.993580115, 1.013593447] | 14/24 | — |
| 5000 | num_rings | repeated_use | baseline | 1.002099616 | [0.991770932, 1.006551371] | 13/24 | — |
| 5000 | num_rings | parse_inclusive | baseline | 1.004599117 | [0.996603351, 1.010848026] | 14/24 | — |
| 5000 | maccs | first_use | baseline | 1.003140787 | [0.998422797, 1.006737479] | 15/24 | — |
| 5000 | maccs | repeated_use | baseline | 1.002244730 | [0.997627710, 1.007433228] | 14/24 | — |
| 5000 | maccs | parse_inclusive | baseline | 1.000543902 | [0.994901429, 1.004334485] | 12/24 | — |
