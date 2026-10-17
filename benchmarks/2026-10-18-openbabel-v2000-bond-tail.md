# Open Babel-oriented V2000 bond-tail experiment — 2026-10-18

**Adoption gate failed; runtime candidate excluded from main.** Combining the existing V2000 bond-type, stereo-code and newline appends improves SDF prepared writing by paired ratios **1.053468087 and 1.073880333** against the pinned baseline. A 3D control's parse median nevertheless falls below the predeclared 0.99 floor: **0.988703242**. Both primary intervals pass, but the regression floor fails; neither the documented storage interruption nor strong external ratios changes that decision.

Measured baseline: `832450783b3faff6ad455a4d514cadf6837b597d`. Measured rejected candidate: [`f0a6e1abe99ce28eadb74bbf386b3d75521eb772`](https://github.com/kent-tokyo/chematic/commit/f0a6e1abe99ce28eadb74bbf386b3d75521eb772) ([source diff](https://github.com/kent-tokyo/chematic/compare/832450783b3faff6ad455a4d514cadf6837b597d...f0a6e1abe99ce28eadb74bbf386b3d75521eb772)); retained on `codex/openbabel-bond-tail-candidate-20261018`. This evidence-only change contains no writer implementation or new runtime tests. Later main advances and published packages were not measured. Earlier accepted zero-atom-field and substituent-stack changes are already part of this baseline; do not multiply their older ratios into these measurements.

Ratios are comparator time / candidate time, computed within each paired block and then summarized. Ratios above one favor the candidate. Unrounded medians and paired bootstrap intervals determine the gate; dividing aggregate medians can give a different answer, especially on this shared host.

| Input | Mode | Comparator | Paired median ratio | Paired 95% CI | Faster blocks |
| --- | --- | --- | ---: | --- | ---: |
| sdf-0 | write | baseline | 1.053468087 | [1.041505863, 1.060281655] | 19/24 |
| sdf-5000 | write | baseline | 1.073880333 | [1.059771969, 1.079352391] | 22/24 |

All three source intervals wholly below one are retained below. The final row fails the 0.99 median floor; the others remain regressions even though they fit that tolerance. These observations do not establish that the changed writer causes a parse regression: no parse algorithm changes in this candidate, and this shared-host experiment does not isolate code layout, allocator, scheduling or competing-work effects.

| Input | Mode | Comparator | Paired median ratio | Paired 95% CI | Faster blocks |
| --- | --- | --- | ---: | --- | ---: |
| v3000-0 | roundtrip | baseline | 0.990867320 | [0.985748730, 0.997280109] | 6/24 |
| v3000_3d-branched-spiral | parse | baseline | 0.992266831 | [0.984093134, 0.996390702] | 6/24 |
| v3000_3d-branched-offset | parse | baseline | 0.988703242 | [0.974267672, 0.997585515] | 6/24 |

All five semantically eligible fresh Open Babel 3.2.1 SDF lanes meet interval lower bound above one and 24/24 faster blocks:

| Input | Mode | Comparator | Paired median ratio | Paired 95% CI | Faster blocks |
| --- | --- | --- | ---: | --- | ---: |
| sdf-0 | write | openbabel | 7.221626713 | [7.145771157, 7.310005872] | 24/24 |
| sdf-0 | roundtrip | openbabel | 7.417029605 | [7.354760398, 7.467666215] | 24/24 |
| sdf-5000 | parse | openbabel | 9.197055305 | [9.134850851, 9.287288529] | 24/24 |
| sdf-5000 | write | openbabel | 10.524843500 | [10.480555205, 10.624402468] | 24/24 |
| sdf-5000 | roundtrip | openbabel | 9.264953754 | [9.207913442, 9.326426383] | 24/24 |

`sdf-0/parse` stays **unranked**. Fresh independent RDKit 2025.09.3 canonical comparison differs for Open Babel at rows 90 and 745; both source parses match all 1,000 reference rows, and all three tools' rewritten graph/atom/bond/coordinate/name/property/unknown-stereo checks pass all 1,000 rows. This does not diagnose an Open Babel reader defect. Numeric observations for the unranked lane remain in the packet. There is no new external ranking for V3000, other formats, SMARTS, descriptors, force fields or conformer generation, and no universal-superiority claim. The external ratios describe the rejected source candidate, not adopted main or a published artifact.

The rejected implementation adds a compile-time ASCII table for the existing seven-byte V2000 bond tail. Bond endpoints, aromatic/Kekule order selection, wedge direction, either markers, order mapping (including quadruple, dative and query orders), all atom fields, chemistry rules, public APIs and the separate 3D-conformer writer remain unchanged. Removing the helper and restoring its three appends reconstructs the complete old writer file byte for byte. No unsafe code, cache or dependency is added.

Local correctness before timing: **1,442 Rust tests passed, 0 failed, 0 ignored, 74 suites**, including every 13 supported bond-order variants × 7 stereo codes against both the old appends and Rust formatting. All-target Clippy, fmt, refactor budget, normal commit hook and 10 runner tests pass. Across all 10 profiles and **4,384 records**, both source parse snapshots and complete warmed rewrite streams are exactly equal. Full output streams and their hashes remain local; a numeric auditor alone cannot independently reconstruct the unexported chemistry bytes. This evidence PR's hosted checks validate integration separately from those measured-source tests.

The plan predates candidate timing. Two unfiltered 1,000-record SDF windows use all six baseline/candidate/Open-Babel permutations four times each. Eight source-only controls use 12 AB and 12 BA blocks each: two normal V3000 windows, long chains, continued lines and four 3D coordinate profiles. Each measured invocation warms in its own mode for at least 500ms, then measures at least 100ms and three full batches. Fresh parse/roundtrip objects are created; prepared write reuses parsed inputs. Startup/input loading are excluded. Rust reuses its output String and Open Babel returns its native allocated output. Whole-process RSS is unranked. Work-unit counts are checked in inner repetitions; complete streams are checked after the same warm/timed budget in preflight, rather than hashed in every inner repeat.

There was one **operational disk pause** after 220 complete blocks, at 4/24 blocks of `v3000_3d-branched-offset`. The 1GiB capacity guard stopped before the next block; no arm or block was partial. The original failure log and partial report were retained. After checking use and deleting only stopped, unopened incremental compilation caches, precisely the 20 still-unmeasured preregistered blocks continued. All **1,464 original invocation journals** were rehashed unchanged, and 120 new invocations were added. No observations were filtered, repeated or replaced; candidate, inputs, order, timing budgets, per-block 1GiB guard and adoption criteria stayed fixed. An extra conservative continuation-start capacity check also stopped before any timed call, with zero added observations; that log is retained. The original source/runtime/input seal was rechecked unchanged before and after continuation. An external workspace Cargo build was observed during disk recovery and later finished; this was a shared host, with no guarantee of global quiescence. There were no own builds, source edits or profiling during measured blocks. The capacity recovery and continuation were explicit operational amendments after partial observations; this is not an uninterrupted epoch. The gate remains **failed**.

All **240 blocks / 1,584 invocations / 40,562 complete timed batches / 20,589,632 processed records** are retained. No earlier elapsed observations were pooled. Both source arms were built from independently empty targets with Rust 1.97.0, opt3/full LTO/one CGU/one job, incremental compilation disabled and no stripping, on macOS 27.0.1 arm64 (Apple M4, 10 CPUs, 16GiB). Comparator executable and Open Babel 3.2.1 runtime/data seals were rechecked; only its executable was reused, with all fresh timings. Input bytes are unchanged from the already-public `9cfba0ad` artifact references in the packet. Sources, six native executables, raw streams, journals, build logs and provenance archives are retained locally. Completed private build targets were removed after saving binaries; original dirty docs and other worktrees were preserved.

The [2026-10-18-openbabel-v2000-bond-tail.json](2026-10-18-openbabel-v2000-bond-tail.json) numeric packet contains all invocations, count/normalization data, all 36 source/external lanes including the unranked observations, input/output hashes, known labels, versions, public input references and an explicit operational-interruption record. It excludes local paths, process/environment inventories, compound streams, SDK files, binaries and source archives. A separate independent full-data audit and the stdlib numeric auditor both recompute a failed gate; a corrupted record-count negative control is rejected.

```sh
python3 scripts/check_openbabel_numeric_packet.py benchmarks/2026-10-18-openbabel-v2000-bond-tail.json
```

All 30 source lanes follow, including the negatives:

| Input | Mode | Comparator | Paired median ratio | Paired 95% CI | Faster blocks |
| --- | --- | --- | ---: | --- | ---: |
| sdf-0 | parse | baseline | 0.995058728 | [0.984039626, 1.001966502] | 10/24 |
| sdf-0 | write | baseline | 1.053468087 | [1.041505863, 1.060281655] | 19/24 |
| sdf-0 | roundtrip | baseline | 1.009022526 | [1.000604105, 1.016115880] | 17/24 |
| sdf-5000 | parse | baseline | 0.994606116 | [0.985926033, 1.001627759] | 9/24 |
| sdf-5000 | write | baseline | 1.073880333 | [1.059771969, 1.079352391] | 22/24 |
| sdf-5000 | roundtrip | baseline | 1.012915074 | [1.006730795, 1.021585444] | 18/24 |
| v3000-0 | parse | baseline | 0.998234822 | [0.996597673, 1.002152913] | 11/24 |
| v3000-0 | write | baseline | 0.996178441 | [0.989760915, 1.001720498] | 9/24 |
| v3000-0 | roundtrip | baseline | 0.990867320 | [0.985748730, 0.997280109] | 6/24 |
| v3000-5000 | parse | baseline | 1.002468123 | [0.994587764, 1.007477198] | 14/24 |
| v3000-5000 | write | baseline | 0.997594165 | [0.996602444, 1.002624463] | 9/24 |
| v3000-5000 | roundtrip | baseline | 0.999867396 | [0.994731001, 1.005407793] | 11/24 |
| v3000-chains | parse | baseline | 1.004941143 | [1.000378755, 1.013782868] | 18/24 |
| v3000-chains | write | baseline | 0.998799988 | [0.986298081, 1.013278428] | 11/24 |
| v3000-chains | roundtrip | baseline | 1.001513652 | [0.993761347, 1.005015590] | 13/24 |
| v3000-continued | parse | baseline | 1.000512204 | [0.997043554, 1.012878413] | 14/24 |
| v3000-continued | write | baseline | 0.999045188 | [0.991939190, 1.011262736] | 11/24 |
| v3000-continued | roundtrip | baseline | 1.000979336 | [0.992935383, 1.005379245] | 13/24 |
| v3000_3d-spiral | parse | baseline | 0.998158877 | [0.994405792, 1.008819487] | 10/24 |
| v3000_3d-spiral | write | baseline | 1.001348132 | [0.990012386, 1.005617885] | 13/24 |
| v3000_3d-spiral | roundtrip | baseline | 1.001693754 | [0.997423977, 1.005164887] | 13/24 |
| v3000_3d-branched-spiral | parse | baseline | 0.992266831 | [0.984093134, 0.996390702] | 6/24 |
| v3000_3d-branched-spiral | write | baseline | 1.001498224 | [0.993412319, 1.018149070] | 13/24 |
| v3000_3d-branched-spiral | roundtrip | baseline | 0.994094480 | [0.990169227, 1.003460196] | 9/24 |
| v3000_3d-offset | parse | baseline | 1.000714247 | [0.995747353, 1.002997998] | 13/24 |
| v3000_3d-offset | write | baseline | 0.994488022 | [0.983677748, 1.005356753] | 8/24 |
| v3000_3d-offset | roundtrip | baseline | 1.000070941 | [0.996866896, 1.007067278] | 12/24 |
| v3000_3d-branched-offset | parse | baseline | 0.988703242 | [0.974267672, 0.997585515] | 6/24 |
| v3000_3d-branched-offset | write | baseline | 1.000319968 | [0.990205811, 1.010088081] | 12/24 |
| v3000_3d-branched-offset | roundtrip | baseline | 0.992629876 | [0.987994082, 1.000014391] | 8/24 |
