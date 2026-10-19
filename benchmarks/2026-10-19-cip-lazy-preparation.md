# RDKit CIP speed audit and lazy preparation — 2026-10-19

CIP assignment was a confirmed current-source weakness. Removing unused
Kekulé/MANCUDE preparation, initializing E/Z ranking on demand, avoiding the
unused Python E/Z lookup pass, and skipping legacy tetrahedral rankings that
accurate labels or explicit abstentions supersede improves this Python API.
Every fallback and assignment order remains covered by an eager reference.

This is a source-wheel result on one macOS arm64 host. Both source wheels
identify as 1.1.8; neither measurement identifies a new published package.
Standalone CIP, canonical SMILES and SMARTS writing still trail RDKit here.

## Retained evidence

- [Numerical packet](2026-10-19-cip-lazy-preparation-numeric.json.gz): all 24 blocks, all 42 comparisons, source/native/wheel identities, screening and correctness counts.
- [All comparison lanes](2026-10-19-cip-lazy-preparation-lanes.csv): includes ties and every slower lane.
- [Runner](../scripts/bench_python_cip_frontier.py), [full CIP output audit](../scripts/verify_python_cip_frontier.py), [independent packet checker](../scripts/check_python_cip_frontier_packet.py).

Baseline wheel was built at `23f774e1b087e500d0a3c325517b1a4a5bc3a760`;
all 1,663 inventoried source files were verified against starting main
`3679321e5af4d8fffe7c3e54aff9639954e7f491` before editing. We reused that
binary and took fresh elapsed samples. Candidate runtime source is
`076ae03b915a890757acc001e8341428c1b4fa3b`; source hashes include its two new
oracle test modules. Compiler, release flags, wheels and imported native
module hashes are in the packet. Both builds enable native InChI and use
one build job, no incremental compilation, release debug 0, and no stripping.

The first hosted Test check rejected file/function growth under the existing
refactor-budget ratchet. Final source moves the existing inline tests to
`cip_tests.rs` and moves fallback explanation into API documentation, without
changing executable logic or loosening limits. A subsequent release build is
recorded under `post_measurement_build`: its native extension and every Python
runtime wheel member are byte-identical to the measured wheel. Only CycloneDX
SBOM metadata and its RECORD checksum changed. The checker validates both source
inventories and the runtime equivalence; the original measurement plan and all
elapsed observations remain unchanged. The final source revision is recorded
in that attestation. Refactor-budget and all 1,081 Rust tests pass locally.

## Protocol and output gates

Fresh sequential processes run baseline, candidate and RDKit in all six
orders four times each, shuffled with seed 20261019. The plan, including
script hashes, identities, scope and adoption policy, was saved before the
first formal observation. No blocks or input rows were removed. Three earlier
two-engine screening blocks cover thirteen operations; these are
diagnostics without a complete output gate and are retained separately.

Formal input is exactly rows 0–999 of the 10,000-row normal v1 corpus.
CPython 3.13.6 and RDKit 2026.09.1 run on macOS 27.0.1 arm64, 10 CPUs.
First-use means newly parsed molecules; repeated-use performs one untimed
call on those molecules; parse-inclusive executes parsing and the operation
inside the same timed loop. It is not a sum of separate timers. Warm-up uses
a disjoint object. Molecules and native return values are retained in every
engine; normalization and output hashing occur outside timing.

`Mol.cip_stereo(mode="accurate")` returns diagnostic dictionaries. RDKit's
`rdCIPLabeler.AssignCIPLabels` mutates molecule properties and returns `None`.
The timed calls retain those native API behaviors; label extraction happens
after timing. These are public Python API timings, not identical internal
kernels. Other operation mappings come from the linked runner and corpus table.
Cached QED and Morgan repeated-use lanes include their existing object caches.
Process peak RSS includes imports and untimed normalization; it is not a
per-operation memory result.

All seven formal operations match RDKit on all 1,000 rows (QED absolute
tolerance 1e-9; strings, labels and bit positions exact). Old/new values also
match. Independently, all 14,072 normal+stress rows preserve complete CIP
dictionaries, optional flags, bond metadata, unresolved results and legacy
results in both raw and RDKit-model representations. The complete serialized
old/new output files have the same SHA-256, recorded in the packet. No engine
refused any CIP audit row. Normal raw CIP labels match RDKit on all 10,000 rows.

Stress raw CIP labels match RDKit on 4,014/4,072 rows; RDKit-model conversion
matches 4,070/4,072. These 58 raw and two model differences predate this change
and remain unchanged. All stress CIP speed comparisons stay unranked. No
filter of matching stress rows is used to establish a speed win.

## CIP measurements

Speedup is comparator time / candidate time. Values above 1 favor candidate.
Brackets are deterministic paired-median bootstrap 95% intervals (10,000
resamples); parentheses count winning blocks. Times below are median ms per
1,000 inputs, while ratios are medians of paired block ratios.

| Boundary | Baseline ms | Candidate ms | RDKit ms | Baseline / candidate | RDKit / candidate |
|---|---:|---:|---:|---|---|
| first_use | 44.208 | 25.386 | 11.007 | 1.742× [1.733, 1.751] (24/24) | 0.433× [0.430, 0.437] (0/24) |
| repeated_use | 43.821 | 25.049 | 10.752 | 1.750× [1.743, 1.754] (24/24) | 0.428× [0.425, 0.431] (0/24) |
| parse_inclusive | 49.852 | 30.913 | 62.896 | 1.616× [1.609, 1.625] (24/24) | 2.040× [2.027, 2.061] (24/24) |

CIP passes the preregistered adoption policy: exact old/new output retention
and source first-use/parse-inclusive confidence bounds above 1. The standalone
CIP API still takes about 2.31 times RDKit's time; the full parse-and-CIP workflow
is 2.04 times faster. Toolkit-wide superiority is not established.

## Other measured operations

Every cell retains the paired ratio, interval and win count. Source control
ratios and all raw elapsed samples are in the packet/CSV. None of the eighteen
source control intervals has an upper bound below 1; several include 1 and
should be treated as ties. Small observed gains are binary-level measurements,
not evidence that unrelated algorithms changed.

| Operation | RDKit / candidate, first | repeated | parse-inclusive |
|---|---|---|---|
| canonical_smiles | 0.967× [0.965, 0.972] (0/24) | 0.962× [0.956, 0.964] (0/24) | 2.447× [2.432, 2.462] (24/24) |
| smarts_write | 0.912× [0.907, 0.920] (0/24) | 0.912× [0.904, 0.915] (0/24) | 2.747× [2.733, 2.750] (24/24) |
| qed | 3.338× [3.331, 3.348] (24/24) | 12284.962× [11915.691, 12521.360] (24/24) | 3.542× [3.539, 3.548] (24/24) |
| morgan2_chiral | 1.565× [1.558, 1.571] (24/24) | 5.355× [5.326, 5.370] (24/24) | 4.316× [4.297, 4.329] (24/24) |
| maccs | 1.817× [1.808, 1.821] (24/24) | 1.860× [1.852, 1.864] (24/24) | 2.004× [1.999, 2.007] (24/24) |
| fp_layered | 1.589× [1.579, 1.593] (24/24) | 1.610× [1.601, 1.618] (24/24) | 1.682× [1.675, 1.685] (24/24) |

Canonical SMILES and SMARTS writing remain slower when parsing is excluded
(approximately 3.4% and 9.6% more first-use time). Their parse-inclusive
workflows remain faster on this fixed corpus. Those standalone gaps need
further work; this change addresses CIP's redundant computation.

## Verification

- Rust: 1,081 passed, one existing ignored test across the three changed crates.
- Independent eager references cover 14,072 low-level assignments, E/Z on raw and Kekulé forms, and the complete accurate-mode fallback merge; tiny budgets, noncarbon centres and lazy initialization have dedicated controls.
- Python: 417 existing tests passed against the isolated candidate wheel.
- Clippy with warnings denied and formatting passed.
- The packet checker recalculates every interval and validates all orders and identities. Twelve corruption/control tests reject changed timing, partial runs, dropped correctness rows, stale identity, invented speedups ranking of unequal stress results, and a changed rebuilt runtime. Hosted CI runs these checks on the full checkout.

The unchanged ordinary/edge-case chemistry results and substantial CIP gain
justify this bounded source change. Published artifacts, other hardware and
other workloads require separate evidence.
