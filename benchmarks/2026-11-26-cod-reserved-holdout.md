# Reserved COD scaffold holdout

The holdout partition is reserved for a future frozen candidate. No conformer
generation or scoring was run on the holdout during preparation. It supplies ordinary
COD structures and does not replace a separate macrocycle holdout or prove
that structures were absent from every unrecorded source.

From the same COD archive used for the two previous 12,000-record samples,
select 12,000 SDF member paths after skipping the first 24,000 in SHA-256
path order. The archive SHA-256 is
`e1d6191fd665dbcc8ce5a24731fefe6ecf6460e0d25f5e290e264646543a16a2`.
The resulting SDF SHA-256 is
`e39322d1eba60bc041fa380f043146e95247d4209d000374591582e6e2acd117`.
The exact member paths and selection rule are in
`validation/manifests/cod-reserved-2026-11-26-source.json`. No paths overlap
either previous 12,000-record sample.

`platinum_scaffold_split.py` with RDKit 2026.09.1 excluded the recorded
training and evaluation manifests, the two previous COD splits, the older
COD split and the macrocycle split. The four external split inputs are stored
compressed in `validation/manifests/cod-reserved-2026-11-26-exclusions/`;
the split JSON records their original paths and uncompressed hashes. The
result contains 5,193 holdout records from 4,991 scaffold groups, 1,345
development records, 5,413 excluded records and 49 parse failures. No
holdout scaffold hash overlaps any of the four earlier splits. Eight invalid
exclusion SMILES were audited conservatively as atom-number subgraphs with
`platinum_invalid_exclusion_audit.py`; no eligible overlap was found.
The split and audit are recorded in
`validation/manifests/cod-reserved-2026-11-26-split.json` and
`validation/manifests/cod-reserved-2026-11-26-invalid-audit.json`.

The holdout graph groups were assigned before generating or scoring any of
their conformers. An RDKit atom-ring basis member of at least 12 atoms takes
precedence as `macrocycle`; otherwise a molecule with at least eight RDKit
strict rotatable bonds is `flexible`, and the remainder is `ordinary`. The
manifest `cod-reserved-2026-11-26-groups.json` records 142, 630 and 4,421
records in these groups, respectively. The ring-basis rule may miss a larger
cycle in a fused system. The user-specified acceptance margin is applied to
the lower bound of each group's paired scaffold-cluster 95% bootstrap
interval, using 10,000 draws and seed 42. Converged-MMFF recovery additionally
uses only molecules independently audited as parameterizable by both engines;
support-audit failures are reported separately. The margin is -2 percentage
points. This prepared partition is not yet an acceptance result.

The timed comparison uses one thread per engine, at most 20 candidates,
seed schedule `42 * (attempt + 1)`, a 0.2-second worker deadline after
initialization, geometry screening and strict declared-stereo scoring.
Recovery requires a candidate at or below 1 Å heavy-atom symmetry-aware RMSD
with all supported geometry checks passing. Timeouts and failures remain in
the relevant group denominator. Report the geometry-check coverage,
stereo-unevaluable counts, worker cleanup time and the exact wheel hashes
alongside the intervals. Both Linux and macOS runs are required; lock the
candidate source and wheel hashes before scoring the holdout.

Preparation used RDKit 2026.09.1 to match the existing split hashes. In the
local Python 3.13 macOS arm64 environment, the 2026.09.2 published wheel
raised `RuntimeError: bad any cast` for Murcko decomposition even for
ethylbenzene. The wheel's embedding and geometry operations were usable for
the earlier pilot; the scaffold issue needs independent investigation before
assuming equivalent decomposition between versions.

The lowest 20 SHA-256(name) records in the separate development partition
were evaluated on macOS 27.0.1 arm64 with the published chematic 1.2.5 and
RDKit 2026.09.2 wheels and, in a separate run, an optimized source candidate.
The fixed 0.2-second worker budget and full-geometry 1 Å screen matched the
pilot protocol. Published chematic recovered 9/20, RDKit 6/20, and the
candidate 9/20; requiring converged MMFF gave 5/20, 3/20 and 5/20,
respectively. The published chematic run had two MMFF setup errors; the
candidate returned unoptimized coordinates for both. Candidate and RDKit
workers were forcibly stopped at the deadline on all 20 inputs. These are
development results, from a small sample, and do not establish parity.
The candidate library SHA-256 is
`208bb5955389f810b40749ac82055537b8954544e27b2e9fa12e356a38e51592`;
the input, script snapshot, and per-record outputs are in
`validation/results/2026-11-26-3d-wheel-pilot/`.

The lowest 100 names by SHA-256 in the same development partition were then
run with the later candidate, `--strict-stereo`, and the same candidate limit,
deadline, geometry screen and seed schedule. Its optimized library SHA-256
was `af8b3a615fc163232d6d277f11bf5829d7334b94eef5adf986dfad5f37af3893`.
The candidate and RDKit returned at least one conformer for 94 and 93 inputs,
and recovered 64 and 53 inputs respectively under the full screen. Requiring
converged MMFF gave 38 and 33. One candidate input raised a generation error
for aromatic sanitization; no RDKit engine errors occurred. Workers were
forcibly stopped on 87 and 92 inputs, including those with earlier candidates.
All 100 inputs had complete geometry-parameter coverage. This deliberately
small development subset, not previously scored in this work, supports diagnosis only and
does not replace the reserved holdout comparison. Its raw output, input and
scorer snapshots are in `validation/results/2026-11-26-3d-wheel-pilot/`.

The same 100 development inputs were subsequently rerun using an installed
candidate wheel and, separately, the published 1.2.5 wheel. Candidate wheel
SHA-256 was
`93819473dc08f988072157852d4247a82ae893685be284f9e498b80e88d39cbb`;
the published wheel hash is recorded in the pilot report. Both used strict
stereo and the same per-engine deadline on the same macOS host. Published
chematic recovered 58/100 and had 11 generation errors; the candidate
recovered 64/100 and had one generation error. RDKit recovered 53/100 in
each run. The two chematic runs each had 38 converged-MMFF recoveries, and
RDKit had 33. Six records gained recovery in the candidate and none lost it;
all six gains involved an unoptimized MMFF-unavailable candidate. All 100
had complete geometry parameter coverage. The candidate wheel was built from
the development worktree with LTO disabled and 16 codegen units. These are
development results, not the holdout acceptance measurement.

An independent support audit found 90 of these 100 inputs in the common
MMFF-parameterizable, matching-graph domain; the other ten were unsupported
by both engines. On the common 90, the candidate and RDKit had 38 and 33
converged recoveries. A paired scaffold-cluster bootstrap of the candidate
wheel run, with 10,000 draws and seed 42, gave a recovery difference of
+11.0 percentage points (95% interval +5.05 to +17.65) over all 100 and a
converged-recovery difference of +5.56 points (interval +1.12 to +11.11)
over the common 90. These are small development-set intervals, not acceptance
evidence. The audited split and support records are saved alongside the raw
results.

## Reserved evaluation: source candidate on GitHub-hosted runners

The frozen candidate commit was `cbefae5d4eefa504185df29c744def7baac60be7`.
This is a source-built development candidate, not the published 1.2.5 wheel.
The run used RDKit 2026.09.2, at most 20 candidates, a hard 0.2-second
per-engine generation deadline, strict declared-stereo and geometry screens,
and a 1 Å reference RMSD threshold. The split and subgroup assignments above
were fixed before these outcomes were inspected. The 5,193 records and their
failures remain in the denominator. This is one reserved COD-derived sample,
not a representative measure of all molecules or all conformer workflows.

| Platform and subgroup | Records | chematic recovered | RDKit recovered | Paired difference (95% scaffold-cluster interval) |
| --- | ---: | ---: | ---: | ---: |
| Linux, all | 5,193 | 2,915 | 2,697 | +4.20 pp (+3.54 to +4.88) |
| Linux, ordinary | 4,421 | 2,840 | 2,639 | +4.55 pp (+3.80 to +5.32) |
| Linux, flexible | 630 | 53 | 40 | +2.06 pp (+0.64 to +3.51) |
| Linux, macrocycle | 142 | 22 | 18 | +2.82 pp (-1.41 to +7.14) |
| macOS, all | 5,193 | 2,764 | 2,295 | +9.03 pp (+8.17 to +9.89) |
| macOS, ordinary | 4,421 | 2,705 | 2,258 | +10.11 pp (+9.16 to +11.08) |
| macOS, flexible | 630 | 40 | 23 | +2.70 pp (+1.43 to +4.10) |
| macOS, macrocycle | 142 | 19 | 14 | +3.52 pp (-0.71 to +8.45) |

The independent MMFF support audit found 4,838 records with support in both
engines and matching graphs, 354 outside that domain and one graph mismatch;
there were no unknown-support records. On that common domain, converged-MMFF
recovery was 2,023 versus 1,878 on Linux and 1,939 versus 1,648 on macOS.
The predeclared lower-bound threshold of -2 percentage points was met for
both all-record recovery and common-domain converged recovery in every
subgroup on both platforms. Absolute recovery for flexible and macrocycle
records remains low under this short budget; satisfying that relative test
does not establish useful absolute performance or general RDKit parity.
The intervals are paired scaffold-cluster bootstrap intervals (10,000 draws,
seed 42), not intervals for an independently sampled population.

The Linux jobs ran on a mixture of Azure and AWS x86-64 runners: 3,927 and
1,266 records respectively. Each record's two engines ran on the same runner,
but hardware and kernel variation prevents treating the aggregate as a
single-host speed measurement. The macOS runner platform string was
`macOS-26.7.2-arm64-arm-64bit-Mach-O` for all shards; individual hosted
runners are still distinct machines. Operating systems use the same molecules,
so their results are not independent datasets.

The original [Actions run](https://github.com/kent-tokyo/chematic/actions/runs/45829103857)
completed all 16 measurement shards. Its Linux report rejected the differing
platform strings; the archived shards were subsequently merged with an
explicit same-OS, same-architecture mixed-platform option that records the
platform per record. The macOS support audit and summaries were reproduced
locally using the exact downloaded macOS candidate wheel. Raw merged records,
support audits, shard manifest and summaries are in
`validation/results/2026-11-26-reserved-cod-timed/`. The wheel SHA-256 values
were `c6c30540957cd97bf3c19a9524263bb467ca48f4611fca72580de481a4712333`
(Linux) and `c1cde89335fc67eb505e3db8843876e228fd9431925055bc82292250bba65b5c`
(macOS). The scorer SHA-256 was
`3ae6bdc4941926c552fbb6f488e4afba7c73bae396889bbfa4558cf6297bc109`.

## Subsequent flexible development probe

After examining the reserved outcomes, 20 flexible molecules were selected
from the separate development partition by the lowest SHA-256 of name, using
the same ring and rotatable-bond classification rule. This subset may overlap
earlier development probes. It was run once per budget on the local macOS
arm64 host with RDKit 2026.09.2 and the earlier candidate wheel SHA-256
`93819473dc08f988072157852d4247a82ae893685be284f9e498b80e88d39cbb`.
Other settings matched the strict 20-candidate protocol above.

| Deadline | chematic recovered | RDKit recovered | chematic median candidates | RDKit median candidates |
| --- | ---: | ---: | ---: | ---: |
| 0.2 s | 2/20 | 1/20 | 2 | 1 |
| 1.0 s | 5/20 | 6/20 | 10 | 7.5 |

This small, unreplicated development probe suggests that the short deadline
limits search coverage; it does not establish a general scaling curve or a
precision ranking. Input and per-record data are stored in
`validation/results/2026-11-26-flexible-development/`. Future changes informed
by these results need a new independent holdout; the reserved cohort above
has now been observed and must not be reused as an untouched acceptance set.

### Partial timeout results on the development subset

The same 20 development molecules were compared using the public Python
`generate_conformers_timed` wrapper, with 20 requested attempts and a
0.5-second deadline including worker startup. The opt-in `return_partial`
profile checkpoints completed, stereo-checked, corrected-gradient attempts;
the default mixed batch discards all results on timeout. Calls alternated in
order on a macOS 27.0.1 arm64 host (Python 3.13.6, one thread per configured
engine). The Rust extension SHA-256 was
`af8b3a615fc163232d6d277f11bf5829d7334b94eef5adf986dfad5f37af3893`;
this is a local source extension, not a published wheel.

All 20 calls in both modes timed out. The default returned zero conformers
for all 20 inputs; the partial profile returned at least one conformer for
all 20, with 110 conformers from 112 completed attempts in total. Median
observed elapsed time including cleanup was 0.5032 and 0.5035 seconds,
respectively. The small, non-isolated development run was not an accuracy or
same-algorithm performance comparison: the profiles use different gradient
schedules, and neither RMSD nor independent geometry was measured here. It
only demonstrates retained results under this protocol. Input, results, code
hashes and environment are in
`validation/results/2026-11-26-flexible-development/timed-partial-0.5.json`;
the runner is `scripts/benchmarks/timed_partial_development.py`.
