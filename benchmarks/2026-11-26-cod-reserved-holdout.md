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
