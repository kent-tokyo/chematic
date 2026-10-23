# v1.2.1 priority follow-ups (2026-10-24)

Candidate chemistry commit: `29b3fdc9ddf17c26cdb65b73be101e8e82046841`.
RDKit compatibility oracle: 2026.09.1. The machine summary is
[JSON](2026-10-24-v121-priority-followups.json).

## Results

| Lane | Result | Boundary |
|---|---:|---|
| Benzene/toluene MCS | 315.05 µs → 3.9332 µs median (80.1×) | Same-host Criterion; same maximum atom/bond counts; not a broad FMCS speed claim |
| Cleanup source path | 4,069/4,069 | Comparable rows in `stress-v1.smi` |
| Cleanup public Python path | 4,069/4,069 | Same pinned RDKit oracle; three incompatible rows remain typed errors |
| Compatibility baseline | v1.2.0 | Changed Rust crates, Python API/stubs, WASM declarations, shared JSON contracts |

The MCS speed recovery uses a safe absolute atom/bond upper bound. Search stops
only after a valid mapping reaches the theoretical maximum. The mixed chiral
three-input fixture is now input-order invariant in chematic. RDKit 2026.09.1
returns either 2 atoms/1 bond or 3 atoms/2 bonds for that fixture depending on
input order, so reproducing that instability is not a compatibility target.

Cleanup now carries explicit writer provenance through the public Python
object and pickle format. The public path therefore uses the same
RDKit-compatible serialization contract as the dedicated Rust path instead of
silently selecting the generic writer.

## Frozen holdout v2

The candidate was committed before the external sources were acquired. v2
excludes all 288 source IDs used by v1, selects rows by fixed SHA-256 ranks,
contains no expected outputs, and was committed before scoring. The exact
macOS arm64 wheel SHA-256 is
`fcba25a934c021ccc8ec5cdc5ded01191ef7484f369a99e624c39529333481f6`.

| Lane | Success | Mismatch/refusal boundary |
|---|---:|---|
| MCS, 64 rows | 54 completed; 29 atom-count and 31 bond-count exact | 9 typed refusals, 1 oracle-invalid row; broad FMCS parity is not achieved |
| SDF attributes/stereo, 64 rows | 64/64 graph, title, properties and stereo | No mismatch or refusal |
| 3D, 96 rows | 96/96 both embedded; 96/96 stereo graph; 95/96 conformer count | Only 4/96 coordinate sets met the strict 1e-4 Å pair-distance contract |

The 3D result distinguishes successful embedding from coordinate parity.
Matching success and stereo does not imply that the generated conformer
ensembles match RDKit. The row-level result is
[retained verbatim](../validation/results/trust_holdout_v2_macos_arm64.json);
inspected v2 rows must not be reused as sealed data or tuned against.

The Linux x86-64 workflow rebuilds the frozen candidate commit and uploads a
second-host row packet. Until that hosted job completes, the scored result is
one-host evidence.

## Verification performed before freezing

- MCS focused Rust tests: 52 passed.
- Cleanup Rust tests: 21 passed; RDKit writer tests: 154 passed.
- Focused Python Cleanup/pickle tests: 14 passed.
- Binding compatibility unit tests: 3 passed.
- Related clippy checks and `git diff --check` passed.

These focused checks and holdout scores do not imply universal superiority to
RDKit or complete compatibility across unsupported FMCS and 3D options.
