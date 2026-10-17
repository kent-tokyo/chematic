# Open Babel-oriented E/Z substituent allocation experiment — 2026-10-17

**Local adoption gate passed.** Replacing each E/Z end's temporary substituent Vec with a two-atom stack buffer improves SDF writing by **4.22% and 4.52%** against the current-main baseline. Every source median remains at least 0.99. One SDF parse lane nevertheless regresses about 0.8%, with its paired interval wholly below one; the gate deliberately permits up to a 1% median loss, and this is not an all-operations speedup.

Baseline: `c84c68f297f95562d1d29c75234f2190cd5f9d75`. Measured candidate: `5a536791b6c7e7565516b4abfeadc9a4b309d888`. Later commits add only the numeric packet, auditor, report and index. Neither the earlier zero-field V2000 writer change nor rejected Morgan/degree-only candidates is included. Do not multiply ratios from those experiments.

Ratios mean **comparator time / candidate time**; above one favors the candidate. Unrounded paired medians and intervals determine the decision, rather than a ratio of aggregate medians.

| Input | Mode | Comparator | Paired median ratio | Paired 95% CI | Faster blocks |
| --- | --- | --- | ---: | --- | ---: |
| sdf-0 | write | baseline | 1.042244207 | [1.029112228, 1.047185683] | 23/24 |
| sdf-5000 | write | baseline | 1.045162746 | [1.034975337, 1.054549305] | 22/24 |

The full source regression with interval upper bound below one follows. It remains visible despite being within the predeclared 0.99 floor.

| Input | Mode | Comparator | Paired median ratio | Paired 95% CI | Faster blocks |
| --- | --- | --- | ---: | --- | ---: |
| sdf-5000 | parse | baseline | 0.992036659 | [0.984703823, 0.996312095] | 5/24 |

The five fresh Open Babel 3.2.1 SDF lanes that pass semantic prerequisites also meet paired interval lower bound above one and 24/24 faster blocks:

| Input | Mode | Comparator | Paired median ratio | Paired 95% CI | Faster blocks |
| --- | --- | --- | ---: | --- | ---: |
| sdf-0 | write | openbabel | 6.497845529 | [6.430567898, 6.564847497] | 24/24 |
| sdf-0 | roundtrip | openbabel | 7.228635352 | [7.207976586, 7.248574350] | 24/24 |
| sdf-5000 | parse | openbabel | 9.168115709 | [9.136092294, 9.223939038] | 24/24 |
| sdf-5000 | write | openbabel | 9.223213177 | [9.147793321, 9.291402038] | 24/24 |
| sdf-5000 | roundtrip | openbabel | 8.959451079 | [8.932122115, 8.985604065] | 24/24 |

`sdf-0/parse` is **unranked**: the independent RDKit 2025.09.3 canonical comparison differs for Open Babel at rows 90 and 745. Both source parses pass all 1,000 rows; both tools' rewritten field checks pass all 1,000 rows. This does not isolate the discrepancy to an Open Babel reader bug. Its numeric observations remain in the packet without a speed claim. No fresh external ranking is made for V3000/3D, SMARTS, other formats, descriptors, force fields or conformer generation. These are exposed development inputs and local native source executables on one macOS 27.0.1 arm64 host, not published-wheel or universal-superiority evidence. Complete superiority over Open Babel remains unproven.

The implementation still visits every neighbor, counts even more than two qualifying substituents, excludes the opposite endpoint and dative bonds identically, and then computes cumulation and implicit H exactly as the old code did. Both ends are evaluated eagerly. Ring-size rejection, lazy exact equivalence classes, distinctness and final bond ordering are unchanged. No degree-only shortcut, unsafe code, new dependency, numeric formatting, writer metadata, chemistry rule or public API change is introduced. The classifier is moved into a private submodule to keep the parent within the refactor ratchet. The actual production iterator changes from 446 to 411 instructions; this proves code realization, not speed by itself.

Correctness: full frozen old-classifier comparisons cover **14,072 corpus rows × three declaration maps**, plus **784** generated high-degree/dative/cumulated cases, empty input and explicit boundary molecules. Changed-crate Rust tests: **1,390 passed, 0 failed, 0 ignored, 68 suites**. All-target Clippy, fmt, refactor ratchet and normal commit Clippy hook pass; steady/SDF/V3000 runner tests: **10 passed**. Both source parse snapshots and complete warmed rewrite bytes match for all **4,384** records across ten profiles. Python binding/whole-workspace hosted validation is separate from these local checks.

The experiment was registered before candidate timing. Two unfiltered 1,000-record SDF windows use all six baseline/candidate/Open-Babel permutations four times each. Eight source-only controls use twelve AB and twelve BA blocks each: normal V3000 windows, long chains, continued lines, and four coordinate profiles. Orders are fixed in advance. Every native invocation warms in its measured mode for at least 500ms, then measures at least 100ms and three complete batches. Timed parse/roundtrip work creates fresh objects; prepared write reuses parsed objects. Input file loading/startup are excluded. Rust write reuses its output String while Open Babel returns its native allocated output. RSS is whole-process and is not ranked. Inner repeats validate work-unit counts; complete output bytes are checked after the same warm/timed budget in preflight, not hashed on every inner repetition.

All **240 blocks / 1,584 native invocations / 40,867 complete timed batches / 20,582,008 processed records** were retained. No elapsed epoch was repeated, filtered or pooled. The original baseline binaries were reused only after their c84c68f2 empty-target build/source/binary proof and seals were rechecked; no old elapsed observations were reused. Candidate starts from an independent empty target with the same Rust 1.97.0, opt3/full LTO/one CGU/one job, no incremental compilation or stripping. All seven probe files match exactly between arms. No own builds, source edits or profiling ran during formal timing. Source, native binaries, SDK/observer runtimes, inputs and preflight streams are hash-checked before/after each profile. Only completed, unopened, private regenerable build targets were cleaned; original dirty docs and other worktrees remain intact.

The [small numeric packet](2026-10-17-openbabel-substituent-pair.json) contains **all 1,584 raw numeric invocation records**, actual repeat/record/work-unit counts, warm/timed durations, all 36 unrounded source/external comparisons and input/output hashes. It contains no local paths, full source/runtime inventories, process listings, compound/output streams or binaries. Input bytes are unchanged from the referenced already-public artifact files at commit `9cfba0ad`; the generation step verifies every decompressed SHA-256. Full source/build/runtime seals, raw output streams and original observer records are retained locally for review. No previous rejected CIP packet is included in this branch or numeric dataset.

Recompute every paired interval, count and gate from the numeric packet with standard-library Python:

```sh
python3 scripts/check_openbabel_numeric_packet.py benchmarks/2026-10-17-openbabel-substituent-pair.json
```

The original full-data independent audit also passes. The numeric-only audit cannot itself reconstruct molecular output bytes; output equality is backed separately by the retained full streams, hashes and checked-in frozen-reference tests. The packet is an audit of this local experiment, not a release attestation. All source lanes follow for inspection:

| Input | Mode | Comparator | Paired median ratio | Paired 95% CI | Faster blocks |
| --- | --- | --- | ---: | --- | ---: |
| sdf-0 | parse | baseline | 0.995367127 | [0.991094737, 1.000875306] | 8/24 |
| sdf-0 | write | baseline | 1.042244207 | [1.029112228, 1.047185683] | 23/24 |
| sdf-0 | roundtrip | baseline | 1.005974048 | [1.001620803, 1.009188370] | 19/24 |
| sdf-5000 | parse | baseline | 0.992036659 | [0.984703823, 0.996312095] | 5/24 |
| sdf-5000 | write | baseline | 1.045162746 | [1.034975337, 1.054549305] | 22/24 |
| sdf-5000 | roundtrip | baseline | 1.005047499 | [1.002091355, 1.009258348] | 20/24 |
| v3000-0 | parse | baseline | 1.001556106 | [0.998012658, 1.010213590] | 13/24 |
| v3000-0 | write | baseline | 1.018449755 | [1.013154660, 1.022905876] | 21/24 |
| v3000-0 | roundtrip | baseline | 1.003740464 | [1.000939044, 1.008524407] | 18/24 |
| v3000-5000 | parse | baseline | 1.007380314 | [1.004097478, 1.011059242] | 19/24 |
| v3000-5000 | write | baseline | 0.998793188 | [0.996243569, 1.000615148] | 10/24 |
| v3000-5000 | roundtrip | baseline | 1.010634678 | [1.004249219, 1.017001845] | 20/24 |
| v3000-chains | parse | baseline | 1.009200654 | [1.004333461, 1.014340161] | 20/24 |
| v3000-chains | write | baseline | 0.997941696 | [0.992797399, 1.004937056] | 11/24 |
| v3000-chains | roundtrip | baseline | 1.005257746 | [1.002334940, 1.011985172] | 17/24 |
| v3000-continued | parse | baseline | 1.001422622 | [0.993624474, 1.010540161] | 12/24 |
| v3000-continued | write | baseline | 0.995482794 | [0.990438436, 1.000585040] | 9/24 |
| v3000-continued | roundtrip | baseline | 1.005354507 | [1.001434464, 1.010984078] | 17/24 |
| v3000_3d-spiral | parse | baseline | 1.001446125 | [1.000001219, 1.007518543] | 17/24 |
| v3000_3d-spiral | write | baseline | 0.998097286 | [0.990216731, 1.004667977] | 11/24 |
| v3000_3d-spiral | roundtrip | baseline | 1.000331934 | [0.998480277, 1.002270296] | 13/24 |
| v3000_3d-branched-spiral | parse | baseline | 1.007932829 | [1.000963430, 1.016067552] | 17/24 |
| v3000_3d-branched-spiral | write | baseline | 1.002171114 | [0.995444467, 1.009185508] | 14/24 |
| v3000_3d-branched-spiral | roundtrip | baseline | 1.000681060 | [0.998181284, 1.007457924] | 14/24 |
| v3000_3d-offset | parse | baseline | 0.998876478 | [0.997218040, 1.008871244] | 11/24 |
| v3000_3d-offset | write | baseline | 0.993697660 | [0.983840797, 1.011764291] | 10/24 |
| v3000_3d-offset | roundtrip | baseline | 1.000416977 | [0.997429866, 1.002588582] | 13/24 |
| v3000_3d-branched-offset | parse | baseline | 1.001007736 | [0.996499783, 1.008392392] | 13/24 |
| v3000_3d-branched-offset | write | baseline | 1.003200500 | [0.994032958, 1.011696262] | 14/24 |
| v3000_3d-branched-offset | roundtrip | baseline | 1.001883953 | [0.997158283, 1.008219082] | 15/24 |

