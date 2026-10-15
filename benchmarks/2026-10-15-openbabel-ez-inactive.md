# Skip unrequested E/Z work — 2026-10-15

The production optimization is adopted in the attached source candidate. This does not establish toolkit-wide superiority, current main or a release.

Decision: Adopt the inactive E/Z topology guard as the thirteenth unmerged source candidate. NCI parse improves 1.423x [1.417,1.432] and roundtrip 1.373x [1.370,1.392]; ChEMBL parse improves 2.027x [2.019,2.035] and roundtrip 1.874x [1.862,1.890]. All four intended lanes are faster in all 21 blocks. Source snapshots and rewrites match exactly for all 2064 original inputs and 64 continued chains; 506592 frozen-API comparisons, 1138 Rust tests (one existing ignored), ten runner gates, Clippy, formatting and unchanged refactor budget pass. No initial source-relative interval is wholly below 1, so the saved adaptive policy selects no follow-up. Retain weaker observations: controlled-chain parse 0.998x [0.983,1.001], 7/21 faster; continued parse 0.993x [0.990,1.005], 8/21. Neither establishes a gain or clear regression. Prepared writer movement is not a new writer algorithm. Controlled-chain Open Babel ratios are parse 125.751x, write 14.467x, roundtrip 106.948x, all 21 faster with semantic agreement; these mostly unchanged lanes are not attributed to this guard. Diverse Open Babel profiles remain unranked for graph/stereo/rewrite differences. All four matrices completed without checkpoint interruption; no timing samples reused or pooled. Existing twelve candidates and all rejected prototype evidence are preserved. This is source-version 1.1.4 evidence, without complete toolkit-wide superiority, current-main verification, hosted CI, merge or publication.

## Change and correctness

E/Z classification first checks whether every literal double bond is already excluded by the existing topology rules. A ring of fewer than eight atoms is unconditionally NotRequested. A terminal end qualifies only if there is no other double bond at either endpoint and the bond is not explicitly unspecified, preserving the earlier cumulene/explicit-unspecified diagnostics. Empty double-bond sets retain the previous fast path. Only when every bond is proven NotRequested are aromaticity, full classification and carrier-conflict work skipped. If any possible request remains, the original complete classifier and carrier resolution run unchanged for all bonds, in their original order. No chemistry, coordinate arithmetic, tolerance, input filtering, public API, dependency or unsafe code is changed. Work is removed when unused; it is not moved to preparation.

The immediate baseline contains all twelve accepted source candidates, including direct V3000 integer output. Rejected prototypes are excluded. Its exact native executable is reused only after all 10,059 current baseline source/input hashes and the binary hash match the accepted source12 build. Baseline build/test logs and origin are retained. The candidate is built freshly in a dedicated target; no timing samples are reused or pooled. Every build and CPU sample completes before formal timing.

A frozen copy of the pre-change public E/Z routine and its private classifier checks exact diagnostics and every atom/bond field, adjacency, stereo annotations, bond directions/anchors and metadata. All 14,072 existing pinned inputs are checked in parsed and attempted Kekule representations (failed Kekule assignments retain the parsed form), three coordinate modes (finite deterministic geometry, missing coordinates, nonfinite coordinates), three unspecified-bond masks (none/all/even) and repeated calls. Alternate rows also carry existing direction annotations and an atom tag. This gives 506,592 exact API comparisons. Independent native snapshots and rewritten bytes must additionally match on all 2,064 fixed MOL inputs and 64 continued chains. No new molecular input records are added.

Scoped release Rust tests: 1,138 passed, 0 failed, 1 existing ignored. Full MOL and perception libraries plus six integration suites cover E/Z/wedges, 3D, square planar, semantic preservation and source-ID lookup. Ten runner gates, release all-target Clippy for MOL/perception, formatting and unchanged refactor budget (487 production files, 8,063 functions, 18 file allowances, 88 function allowances) pass. Complete logs remain in the packet.

Baseline CPU profiles over the fixed 1,000-input batches are retained. The read harness also performs untimed preparation; samples diagnose the executed read path, not pure timed latency. Inclusive symbol counts overlap and are not additive. They identify aromaticity within E/Z as a target, not a predicted speedup.

## Protocol

Fixed unfiltered profiles: 1,000 records at offsets 0 and 5,000, plus 64 controlled linear chains with lengths 16 through 1,024. One exposed development corpus on one macOS arm64 host, without a sealed holdout or toolkit-wide coverage. Open Babel 3.2.1 libraries/plugins and RDKit 2026.09.1 observer are pinned by hashes. Each profile uses 21 rotating three-arm fresh-process blocks, five batch repeats and all parse/prepared-write/roundtrip modes.

Ratios are reference elapsed / candidate elapsed; above 1 is faster. Brackets are paired median-ratio 95% intervals using 10,000 deterministic bootstrap resamples. A strict win requires semantic agreement, CI lower > 1 and all 21 blocks faster. Original and follow-up observations remain separate, without pooling or substitution. Movement in unchanged paths may include host noise and binary layout; it is not proof of a new algorithm.

Timers exclude startup, loading/framing and plugin discovery. Preparation retains a full parsed batch and performs one warm write. Parse evaluates new batches with default diagnostics, prepared-write creates strings from prepared records, roundtrip parses and writes fresh records. RSS includes preparation/plugins and is whole-process high-water usage. Output, source, input, executable and native library hashes are checked.

An external checkpoint/artifact wrapper reserves 16 MiB of owned storage, releases it on ENOSPC and retries writing the same computed report for at most 60 seconds. The frozen native timers, observer/bootstrap functions, selection policy and measured observations are unchanged. Logs retain any such event; collection-completion.json records interruptions if any.

## Original full matrix

| Profile | Lane | Baseline ms | Candidate ms | Baseline / candidate [95% CI]; faster blocks | Open Babel / candidate [95% CI]; faster blocks |
|---|---|---:|---:|---|---|
| 0 | parse | 159.214 | 78.590 | 2.027 [2.019, 2.035]; 21/21 **strict win** | unranked (semantic mismatch) |
| 0 | write | 13.922 | 13.820 | 1.008 [1.002, 1.018]; 15/21 | unranked (semantic mismatch) |
| 0 | roundtrip | 176.264 | 94.035 | 1.874 [1.862, 1.890]; 21/21 **strict win** | unranked (semantic mismatch) |
| 5000 | parse | 46.289 | 32.431 | 1.423 [1.417, 1.432]; 21/21 **strict win** | unranked (semantic mismatch) |
| 5000 | write | 3.781 | 3.679 | 1.033 [1.025, 1.037]; 21/21 **strict win** | unranked (semantic mismatch) |
| 5000 | roundtrip | 49.697 | 36.157 | 1.373 [1.370, 1.392]; 21/21 **strict win** | unranked (semantic mismatch) |
| chains | parse | 32.720 | 32.814 | 0.998 [0.983, 1.001]; 7/21 | 125.751 [124.556, 127.019]; 21/21 **strict win** |
| chains | write | 6.891 | 6.607 | 1.035 [1.021, 1.057]; 18/21 | 14.467 [14.385, 14.712]; 21/21 **strict win** |
| chains | roundtrip | 39.783 | 39.386 | 1.007 [0.995, 1.009]; 14/21 | 106.948 [106.437, 107.885]; 21/21 **strict win** |

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
| parse | 50.230 | 50.353 | 0.993 [0.990, 1.005]; 8/21 |
| write | 6.763 | 6.430 | 1.054 [1.041, 1.066]; 20/21 |
| roundtrip | 57.157 | 56.794 | 1.007 [1.005, 1.010]; 18/21 |

## Adaptive regression checks

The rule saved before timing selects an entire profile for one separate 21-block source-only all-mode follow-up when any initial source-relative CI is wholly below 1. Initial outcomes remain; newly appearing follow-up negatives do not cause recursive retries.

```json
{
  "policy": "One separate 21-block two-source-arm follow-up per selected profile, all three modes; retain initial results; no pooling or recursive retries.",
  "profiles": {}
}
```

## Reproduction and status

[The evidence packet](2026-10-15-openbabel-ez-inactive-artifacts/README.md) retains exact measured binaries, every matrix/snapshot/rewrite, failures if any, build/test/profile logs, full frozen source/input archive and patch. Archive readback verifies all 10060 candidate hashes and reconstructs all 10059 baseline hashes. An independent tiny Git repository checks patch reconstruction. Source reconstruction is not bit-reproducibility of rebuilt binaries.

Last verified shared Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`; refs were not refreshed. Source version is 1.1.4, without merge, hosted CI or publication. Immediate-baseline improvements do not establish recovery of every earlier regression; older writer-induced read regressions, Python Morgan/MACCS observations and other unresolved lanes remain separate evidence.
