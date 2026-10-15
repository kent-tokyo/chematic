# Inline E/Z ring-search queue — 2026-10-15

The production optimization is adopted in the attached source candidate. This does not establish toolkit-wide superiority, current main or a release.

Decision: Adopt the inline E/Z ring-search queue as the fourteenth unmerged source candidate. ChEMBL parse improves 1.145x [1.136,1.147] and roundtrip 1.108x [1.105,1.118]; NCI parse improves 1.180x [1.174,1.183] and roundtrip 1.149x [1.143,1.156]. All four intended lanes are faster in all 21 blocks. Prepared write remains approximately 1.000x in both slices; it is not a new writer algorithm. All 4963770 independent old/new BFS queries match, including 90005 above the 32-atom inline boundary; all existing 506592 public E/Z comparisons pass. Native snapshots/rewrite bytes match on all 2064 original inputs, 64 continued chains and 64 additional branched records; 1139 Rust tests (one existing ignored), ten runner gates, Clippy, formatting and unchanged refactor budget pass. No initial source-relative CI is wholly below 1, so the saved adaptive policy selects no follow-up. Retain weak negatives: original-chain parse 0.991x [0.984,1.029] and continued parse 0.994x [0.963,1.008], neither a clear regression nor a gain. A separate post-hoc exploratory large-neighborhood plan was saved before generating/timing its 64 valid branched graphs: parse 1.012x [1.010,1.025], 16/21 faster; write 0.992x [0.947,1.011]; roundtrip 1.007x [0.994,1.014]. None is a strict all-block gain and no clear regression is demonstrated; all observations remain without pooling or recursive retesting. Controlled-chain Open Babel ratios are parse 127.521x, write 14.821x, roundtrip 109.917x, all 21 faster with semantic agreement. These unchanged no-double-bond controls are not attributed to the new ring queue. Diverse Open Babel profiles remain unranked for graph/stereo/rewrite differences. All five matrices completed without checkpoint interruption, lost or recollected observations. All thirteen accepted predecessors and rejected experiment evidence are preserved. This is source-version 1.1.4 evidence without complete toolkit-wide superiority, current-main verification, hosted CI, merge or publication.

## Change and correctness

The small-ring breadth-first search uses one inline queue of up to 32 atom indices as both the visited sequence and the level frontiers. Captured begin/end positions preserve the exact depth boundary and adjacency traversal order. Before adding a 33rd visited atom, the candidate seeds a HashSet from the queue and uses hashed membership thereafter. Large neighborhoods therefore retain hashed lookup complexity. An exhausted frontier returns false without unused empty iterations. No chemical predicate, coordinate arithmetic, tolerance, input filtering, public API, dependency or unsafe code is changed. The preceding inactive E/Z guard and every fallback classification/diagnostic remain unchanged.

The immediate baseline contains all thirteen accepted source candidates, including the inactive E/Z topology guard. Rejected prototypes are excluded. Its exact native executable is reused only after all 10,060 current baseline source/input hashes and its binary hash match the accepted source13 build. Baseline build/test logs and origin are retained. The candidate is built freshly in a dedicated target; no timing samples are reused or pooled. Every build and CPU sample completes before formal timing.

The independent BFS oracle freezes the original traversal; only its return values add a visited-count diagnostic. Every bond of all 14,072 existing pinned inputs is checked at 15 thresholds (0 through 12, 33 and 65), totaling 4,963,770 queries. Results match exactly: 4,873,765 queries stay within 32 visited atoms and 90,005 exceed that boundary, exercising both storage paths. True and false outcomes both occur. Existing 506,592 frozen public E/Z comparisons also pass again with coordinate failures, unspecified-bond masks, annotations, parsed/Kekule representations and repeated calls. Independent native snapshots and rewritten bytes match on all 2,064 fixed MOL inputs and all 64 continued chains. These oracle comparisons reuse existing licensed input records. A separate project-authored synthetic diagnostic is described below.

Scoped release Rust tests: 1,139 passed, 0 failed, 1 existing ignored. Full MOL and perception libraries plus six integration suites cover E/Z/wedges, 3D, square planar, semantic preservation and source-ID lookup. Ten runner gates, release all-target Clippy for MOL/perception, formatting and unchanged refactor budget (487 production files, 8,063 functions, 18 file allowances, 88 function allowances) pass. Complete logs remain in the packet.

Fresh source13 CPU profiles over the fixed 1,000-input batches are retained. The harness also performs untimed preparation; samples diagnose the executed read path, not pure timed latency. Inclusive counts overlap and are not additive; LTO can hide symbols. Sampling confirms ring-search work as a target, not a predicted speedup.

## Protocol

Fixed unfiltered profiles: 1,000 records at offsets 0 and 5,000, plus 64 controlled linear chains with lengths 16 through 1,024. One exposed development corpus on one macOS arm64 host, without a sealed holdout or toolkit-wide coverage. Open Babel 3.2.1 libraries/plugins and RDKit 2026.09.1 observer are pinned by hashes. Each profile uses 21 rotating three-arm fresh-process blocks, five batch repeats and all parse/prepared-write/roundtrip modes.

Ratios are medians of per-block paired reference/candidate elapsed ratios; above 1 is faster. Milliseconds are separate marginal medians, so their quotient need not equal the paired ratio. Brackets are paired median-ratio 95% intervals using 10,000 deterministic bootstrap resamples. A strict win requires semantic agreement, CI lower > 1 and all 21 blocks faster. Original and follow-up observations remain separate, without pooling or substitution. Movement in unchanged paths may include host noise and binary layout; it is not proof of a new algorithm.

Timers exclude startup, loading/framing and plugin discovery. Preparation retains a full parsed batch and performs one warm write. Parse evaluates new batches with default diagnostics, prepared-write creates strings from prepared records, roundtrip parses and writes fresh records. RSS includes preparation/plugins and is whole-process high-water usage. Output, source, input, executable and native library hashes are checked.

An external checkpoint/artifact wrapper reserves 16 MiB of owned storage, releases it on ENOSPC and retries writing the same computed report for at most 60 seconds. The frozen native timers, observer/bootstrap functions, selection policy and measured observations are unchanged. Logs retain any such event; collection-completion.json records interruptions if any.

## Original full matrix

| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |
|---|---|---:|---:|---|---|
| 0 | parse | 78.522 | 68.795 | 1.145 [1.136, 1.147]; 21/21 **strict win** | unranked (semantic mismatch) |
| 0 | write | 13.733 | 13.821 | 1.000 [0.993, 1.007]; 10/21 | unranked (semantic mismatch) |
| 0 | roundtrip | 93.665 | 84.263 | 1.108 [1.105, 1.118]; 21/21 **strict win** | unranked (semantic mismatch) |
| 5000 | parse | 32.560 | 27.596 | 1.180 [1.174, 1.183]; 21/21 **strict win** | unranked (semantic mismatch) |
| 5000 | write | 3.707 | 3.691 | 1.000 [0.990, 1.011]; 11/21 | unranked (semantic mismatch) |
| 5000 | roundtrip | 36.243 | 31.509 | 1.149 [1.143, 1.156]; 21/21 **strict win** | unranked (semantic mismatch) |
| chains | parse | 46.288 | 39.756 | 0.991 [0.984, 1.029]; 10/21 | 127.521 [125.004, 130.277]; 21/21 **strict win** |
| chains | write | 9.136 | 8.214 | 1.014 [0.981, 1.036]; 12/21 | 14.821 [14.412, 15.901]; 21/21 **strict win** |
| chains | roundtrip | 56.062 | 47.841 | 0.998 [0.981, 1.027]; 10/21 | 109.917 [107.823, 112.002]; 21/21 **strict win** |

## Semantic exclusions

Diverse Open Babel lanes remain unranked whenever indexed graphs, stereo or rewritten fields differ. No records or stereo markers are removed to obtain a ranking. Controlled-chain agreement establishes only that checked workload. Source-only snapshots include indexed atoms/bonds, stereo order, exact coordinate bits, names/metadata and geometry/stereo diagnostics; rewrite bytes are also required to match.

| Profile | Open Babel parse graph mismatches | Parse rejected | Rewrite rejected | Rewrite graph mismatches | Unknown-stereo mismatches |
|---|---:|---:|---:|---:|---:|
| 0 | 3 | 0 | 153 | 1 | 1 |
| 5000 | 137 | 0 | 0 | 137 | 137 |
| chains | 0 | 0 | 0 | 0 | 0 |

## Separate continuation diagnostic

All numeric atom/bond rows of the 64 chains are physically continued after the id. Original and continued source snapshots/rewrites must be byte-identical in both arms. This source-only diagnostic has 21 alternating fresh-process blocks, five repeats and all modes. It does not rank Open Babel or replace an original profile.

| Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks |
|---|---:|---:|---|
| parse | 80.921 | 82.104 | 0.994 [0.963, 1.008]; 9/21 |
| write | 10.581 | 10.717 | 0.991 [0.946, 1.049]; 10/21 |
| roundtrip | 92.551 | 93.157 | 0.993 [0.977, 1.014]; 9/21 |

## Adaptive regression checks

The rule saved before timing selects an entire profile for one separate 21-block source-only all-mode follow-up when any initial source-relative CI is wholly below 1. Initial outcomes remain; newly appearing follow-up negatives do not cause recursive retries.

```json
{
  "policy": "One separate 21-block two-source-arm follow-up per selected profile, all three modes; retain initial results; no pooling or recursive retries.",
  "profiles": {}
}
```

## Additional large-neighborhood diagnostic

This exploratory source-only diagnostic was declared after gains were observed in the two initial fixed slices, before generating or timing its inputs. It neither replaces nor pools with any original matrix and does not compare Open Babel. All outcomes remain, without recursive follow-up.

The project-authored recipe creates 64 valid acyclic carbon graphs: a central C=C with equal balanced binary subtrees at both ends, eight depths and eight fixed atom permutations, spanning 6 through 1,022 atoms. All graphs sanitize in pinned RDKit; no record is removed. Production ring search visits more than 32 atoms on 32 records and stays inline on 32. Full provenance, license, minimization, exact recipe and pinned oracle are in large-neighborhood-corpus-manifest.json.

Before timing, native snapshots and rewrite bytes match between both source arms for every record. RDKit verifies the generated versus native canonical/rewritten graph, name and exact coordinates. The diagnostic uses 21 alternating two-arm fresh-process blocks, five repeats and all modes through the unchanged source-only driver. It is exposed synthetic evidence, not a sealed holdout or a guarantee for every large graph.

| Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks |
|---|---:|---:|---|
| parse | 25.019 | 24.689 | 1.012 [1.010, 1.025]; 16/21 |
| write | 4.275 | 4.379 | 0.992 [0.947, 1.011]; 10/21 |
| roundtrip | 29.292 | 29.191 | 1.007 [0.994, 1.014]; 14/21 |

## Reproduction and status

[The evidence packet](2026-10-15-openbabel-ring-queue-artifacts/README.md) retains exact measured binaries, every matrix/snapshot/rewrite, failures if any, build/test/profile logs, full frozen source/input archive and patch. Archive readback verifies all 10061 candidate hashes and reconstructs all 10060 baseline hashes. An independent tiny Git repository checks patch reconstruction. Source reconstruction is not bit-reproducibility of rebuilt binaries.

Last verified shared Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`; refs were not refreshed. Source version is 1.1.4, without merge, hosted CI or publication. Immediate-baseline improvements do not establish recovery of every earlier regression; older writer-induced read regressions, Python Morgan/MACCS observations and other unresolved lanes remain separate evidence.
