# Platinum Diverse holdout protocol

This protocol is fixed before running the holdout. Its split is
`validation/results/platinum-diverse-scaffold-split-audited-2026-11-21.json`;
the corresponding invalid-exclusion graph audit must have no overlaps. The
known-exposure exclusion inputs and their hashes are recorded in that split.
Other historical exposure cannot be ruled out. The 1,606 holdout molecules
are selected without reference to conformer outcomes.

Compare a pinned source-built chematic wheel against RDKit 2026.09.2. Record
the source commit, wheel hashes, operating system, and tool versions. For
each molecule generate at most 20 conformers with seed 42, corrected
embedding, MMFF94 at most 200 iterations, and chematic's optional bounded
400-iteration near-linear nitrogen retry. Apply the same predefined stereo
and heuristic geometry screens to each engine. A recovery requires at least
one conformer within 1 Å heavy-atom symmetry-aware RMSD of the crystal
structure, retaining specified stereochemistry and passing the screens.
Include failures in the full denominator and report them by engine. Record
convergence-required recovery separately.

The fixed-candidate comparison is assessed against a prespecified 2 percentage
point noninferiority margin using the lower end of a paired 95% bootstrap
interval. Report the raw difference and interval even when the margin is not
met. The interval resamples molecules and does not account for correlation
among molecules with a common scaffold; include a scaffold-cluster interval
before making a parity claim. This analysis does not establish chemical
accuracy for molecules outside the dataset or force-field parameter coverage.

Independently repeat on the same holdout at 0.2 seconds per molecule per
engine, at most 20 attempts, seed schedule `42 * (attempt + 1)`, with the
same screens and bounded retry. Timing includes embedding, minimization and
retry; RMSD and geometry scoring are outside the timed region. Engine order
alternates. Report completion counts, recovery, convergence-required recovery,
paired intervals and errors. Time-budget results are specific to the host and
load. The source-built wheel is not evidence about a published release; run
and report the same verification with published macOS and Linux wheels after
release. Do not conclude general RDKit parity from this benchmark alone.

The macOS arm64 source wheel was built before opening holdout results from
source commit `d30d186` (the subsequent benchmark-only commits do not change
`crates/`). Its SHA-256 is
`018eeb2fd87f6b0958610b4f3e2974de33a50bd31cf0a9a82324a7f01954ec35`.
The holdout probe and reporting scripts are pinned by their Git commits and
the probe records its own SHA-256 in its output.

## macOS arm64 source-wheel results

The fixed-candidate run completed all 1,606 holdout records using RDKit
2026.09.2 and the pinned source wheel above. Under the prescribed stereo and
heuristic geometry screens, chematic recovered 1,012 molecules within 1 Å
and RDKit recovered 1,007. The difference is +0.31 percentage points; the
paired molecule-bootstrap 95% interval is −0.68 to +1.31 points, and the
scaffold-cluster interval is −0.74 to +1.31 points. Both lower bounds exceed
the preregistered −2 point noninferiority margin. Requiring a converged
matching conformer yields 703 and 696 recoveries respectively.

The results record 32,100 minimized conformers for each engine. No
specified-stereo violations were detected. The heuristic screen flagged
angle violations in 75 chematic and 88 RDKit conformers, and no bond-length
or heavy-atom clash flags. All eligible molecular bonds and angles had
reference parameters for the screen. For one input, `840_3I8P_A`, both engines returned no conformers. An
independent rerun confirmed that RDKit embedded zero conformers before its
optimization call raised `ValueError('Bad Conformer Id')`. Both were counted
as failures. The screen is a
practical flagging rule rather than a proof of chemical validity.

At 0.2 seconds per engine and molecule, the first run recovered 803/1,606
for chematic and 707/1,606 for RDKit, with paired scaffold-cluster 95%
interval +4.61 to +7.41 points. Removing unused local build intermediates to
recover disk space overlapped this run, so that run may reflect the workload.
The independently repeated run recovered 796 and 683 respectively; its
scaffold-cluster interval is +5.58 to +8.52 points. Requiring converged
matching conformers in that repeat yielded 623 and 530, with a cluster
interval of +4.57 to +7.02 points. There were no reported per-row errors
or declared stereo violations in either timed run. These are measurements
on one macOS arm64 host; they do not establish universal performance.

The gzip-compressed per-record results are
`validation/results/platinum-holdout-macos-fixed-2026-11-22.json.gz`,
`validation/results/platinum-holdout-macos-time-2026-11-22.json.gz`, and
`validation/results/platinum-holdout-macos-time-repeat-2026-11-22.json.gz`.
Their hashes are `c506d6c197293e65faceed40dd74c89e7b3427dc3c46b571600e63ce04c4bdb1`,
`1a4b77521f2c37ba5a3998a2f0c7fec3b54dc609b55f942b41f3aa12d9f54874`, and
`1faa6558ae9851a59baf63817ddcfa672d3dd162d60c7b86c64084d3d11b9708`,
respectively. Machine-readable molecule and cluster summaries are alongside
them. The 400-iteration retry and corrected embedding are opt-in APIs; these
results do not represent the default compatibility API. This source-built
wheel is not the already-published 1.2.2 wheel. Linux and published-package
holdout evidence remain separate requirements.

## Common-coordinate force-field follow-up

Using the same source wheel, the holdout molecules were separately embedded
by RDKit with seed 42, and MMFF94, MMFF94s and UFF were evaluated on the
same explicit-H coordinates. This follow-up uses molecules already selected
for RDKit MMFF support, so it is not an independent estimate of general
parameter coverage. Of 1,606 molecules, 1,605 produced coordinates and
all three fields were constructed for both engines; the remaining molecule
was `840_3I8P_A`, for which RDKit embedding failed. The 95th percentile of
absolute energy error was at most `1.71e-13 kcal/mol`, and that of individual
gradient-component error was at most `1.55e-13 kcal/mol/Å`.

The maximum MMFF94 and MMFF94s energy error was **309.36 kcal/mol**, on
`PP9_1Z53_A`; the maximum component gradient error was 22.37 kcal/mol/Å.
The error came entirely from electrostatics. The numeric atom types and all
other force-field terms matched in the diagnostic. RDKit selects the first
ring containing a type-76 atom in its [MMFF charge implementation](https://github.com/rdkit/rdkit/blob/Release_2026_09/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp#L3220). RDKit assigned the two
type-76 nitrogen formal charges independently on their respective five-member
rings; the candidate shared the charge across a larger ring. Their MMFF
partial-charge totals were −6 and −5 respectively, although the molecule's
input formal charge is −4. Substituting RDKit's partial charges into the
candidate reproduced the RDKit energy and gradient at the measured geometry;
this is a diagnostic override, not a correction or a claim that RDKit's
charge total is physically correct. A general, validated ring-order fix is
still required. Because this input has now been inspected, it cannot be used
as untouched holdout evidence for changes made in response to it.

Maximum UFF energy and gradient-component errors were `6.71e-12 kcal/mol`
and `0.00629 kcal/mol/Å`. The raw per-record follow-up is
`validation/results/platinum-holdout-macos-force-fields-2026-11-22.json.gz`
(SHA-256 `2eb33a5e5c908e86e4e2fde743b33ee687aaf8c1ac611a0edee18ff9ca3340f0`);
the adjacent JSON summary includes sampled finite-difference errors. The
run resumed from a validated prefix after a disk-space interruption.

## Linux x86-64 source-wheel results for the preceding candidate

[Hosted run 44917477313](https://github.com/kent-tokyo/chematic/actions/runs/44917477313)
built merge commit `615187d54fd59c0d572e25732cfb03bd88414f4b`, combining
base `f97aa95253f2dc38b0d98a58fde1f646febe8aca` and head
`9dd693ef8874833fba5d1ceeb029111710ed08cd`. Its CPython 3.12 Linux wheel
SHA-256 is `c6252725d3e632c485ca6289597c4450aed1e74242598fcb1c0922a0deb16c0e`.
The fixed-candidate recovery was 1,012/1,606 for chematic and 1,005/1,606 for
RDKit; the paired scaffold-cluster 95% interval for the difference was −0.56
to +1.44 points. Converged matching recovery was 699 and 701 respectively.

At 0.2 seconds, recovery was 735 and 679 respectively, with a
scaffold-cluster interval of +2.42 to +4.59 points. Requiring a converged
matching conformer yielded 574 and 540, interval +1.23 to +3.03 points.
No timed rows reported errors or specified-stereo violations. These are
single-host source-wheel results, not published-wheel evidence. They predate
the type-76 charge change and do not independently validate that change.

Raw records are `validation/results/platinum-holdout-linux-fixed-2026-11-22.json.gz`
(SHA-256 `eb9ad8fdadfeec2776b87f9b99d0cd9b0efe7f3b48f4274635c1832143f2a4a3`)
and `validation/results/platinum-holdout-linux-time-2026-11-22.json.gz`
(SHA-256 `eede0b2905ad049781d6a937c0fafaf62c64374de641882cf7101096932dce04`).

## Separate macrocycle dataset for the revised candidate

The [CONFORGE supplementary archive](https://zenodo.org/records/8137603)
contains `prime_dataset.sdf`, a collection of 235 experimental reference
structures used in the study's macrocycle assessment. The archive SHA-256 is
`30a5218e80e442a7e5f9129433d307b7b11dd569930d33244fd9fbd26fbf173b`, and
the extracted SDF SHA-256 is
`28fc27d7bd33c7df5ab4ff63acb320cd5216ae91262d1288cd1fe18a045adbdb`.
All 235 records parsed with RDKit 2026.09.1. They have 208 distinct names
and 206 distinct isomeric SMILES; multiple records can represent the same
molecule. They include atoms in rings of size 10–40 and 11–119 heavy atoms;
this is a distinct, challenging population and is not interchangeable with
the Platinum result.

The recorded scaffold split first excludes all 2,859 Platinum structures,
including the now-exposed holdout, and the other known development manifests
listed in `validation/results/prime-scaffold-split-2026-11-22.json`. It leaves
31 development, 188 holdout, and 16 excluded records (28, 166, and 14
distinct names respectively). Eight invalid exclusion
inputs were also checked conservatively by atom-and-connectivity subgraph;
`validation/results/prime-invalid-exclusion-graph-audit-2026-11-22.json`
records zero overlaps. Unknown historical exposure remains possible. The
whole-Platinum graph manifest is generated by
`scripts/benchmarks/platinum_exposure_manifest.py` from the pinned Platinum
SDF, and is recorded as
`validation/manifests/platinum-diverse-2017-exposure.json`.

Before evaluating the 188 holdout structures, run the 31 development records
and record support, errors, geometry coverage, runtime and fixed-candidate
recovery. Freeze the revised source and wheel hashes, and fix any operational
limits and acceptance criteria before opening the holdout. The smaller
holdout may be insufficient to establish a 2-point noninferiority margin; a
wide interval is an inconclusive result. The existing Platinum holdout has
already been inspected and is development evidence for the revised code.
