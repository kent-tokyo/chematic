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
ring containing a type-76 atom in its [MMFF charge implementation](https://github.com/rdkit/rdkit/blob/Release_2026_09_2/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp#L3220). RDKit assigned the two
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

## Additional Platinum structures for the revised candidate

The [complete Platinum 2017 dataset](https://comp3d.univie.ac.at/fileadmin/user_upload/p_comp3d/datasets/platinum-dataset-2017-01-sdf.zip)
contains 4,548 structures, including the previously evaluated Diverse subset.
The ZIP SHA-256 is `4370821a41ac452145dc97ec20e442b6d3ea838d4f923607aec00c94e47030ee`;
the extracted SDF SHA-256 is
`715fa18108340f31f09c68f62a9484970b49c78653cbf950dbc411e3ab1b53d8`.
Its recorded split excludes scaffolds in the entire Diverse and Prime inputs
and the other known development manifests. It contains 129 development,
430 holdout and 3,989 excluded records. The conservative audit of the eight
invalid exclusion graphs found no candidate overlap. Both Prime and this
split remain bounded by the known exposure list and cannot rule out other
historical exposure. The two new candidate sets are scaffold-disjoint because
the complete-Platinum split excludes every Prime scaffold.

As with Prime, run the development records first. Freeze the revised source,
wheels, evaluation, error handling and acceptance criteria before scoring the
holdout. Confirm the number of scorable records and interpret uncertainty in
view of the 430-record denominator; a wide interval cannot establish the
2-point margin. The published-package checks remain a separate step.

The type-76 fix's ten exposed regression molecules contain 36 type-76 atoms.
Among them, no atom has equally short containing rings with different numbers
of type-76 nitrogens. They therefore do not test the equal-length cycle-order
boundary. Relevant-cycle enumeration is capped at 10,000 cycles; above that
cap the code retains its earlier SSSR fallback, which is not established as
RDKit 2026.09 parity for every graph.

The revised candidate's [Linux development run](https://github.com/kent-tokyo/chematic/actions/runs/44926444020)
built merge commit `fdb86b4598785a68ee0eb12c9ace49a1a98790fc` from base
`f97aa95253f2dc38b0d98a58fde1f646febe8aca` and head
`d3d320839993ce18437c496483f62ba032b2d83e`. Its source-wheel SHA-256 was
`5dc978885dab49252e20c5f955f9aa23a1bf74e3249129edd245900454fcaa2c`.
All 129 development records were processed by both engines, and all bond
and angle references were available. The screened 1 Å recovery was 68 for
chematic and 63 for RDKit; the paired molecule-bootstrap interval was +0.78
to +7.75 percentage points. The 2,580 minimized candidates per engine had
no detected specified-stereo violations; the heuristic screen flagged angle
violations in 23 and 24 conformers, respectively. This is development data,
not evidence on the remaining 430 records. The wheel was not retained as an
artifact, so it must not be treated as an executable artifact for holdout
scoring. Raw development results are in
`validation/results/platinum-complete-linux-development-2026-11-22.json.gz`
(SHA-256 `1681d406cfba385cb2a6466e5ab92010b36c4623a64885c4cdf853f0bc653f60`).

The revised candidate's [Prime development run](https://github.com/kent-tokyo/chematic/actions/runs/44926079891)
built Linux merge commit `8451329616ebac9d4089d77caa094097e653fb2d` from base
`f97aa95253f2dc38b0d98a58fde1f646febe8aca` and head
`4b24f16884d4d57d676a37803c578e041a1714ea`. The source-wheel SHA-256 was
`d9d011803b68df28309b353645b030f09c639ff87bf585fb1b9ce34349abbccc`.
All 31 records were processed by both engines with complete bond and angle
screen references. Twenty candidates per engine gave 11 screened recoveries
within 1 Å for chematic and 8 for RDKit; the paired molecule-bootstrap
interval was −3.23 to +22.58 points. The molecular and scaffold dependence
of repeated records is not modeled by this development summary. The results
are highly uncertain and neither engine recovered most references. Across
620 minimized candidates each, no specified-stereo violations were detected;
the angle screen flagged one chematic conformer and no RDKit conformers.
The fixed-candidate elapsed time per engine reached approximately 54 seconds
for chematic and 189 seconds for RDKit on individual records on this host.
Those batch times are exploratory and do not establish a speed comparison.
The executable wheel was not retained. Raw records are
`validation/results/prime-linux-development-2026-11-22.json.gz` (SHA-256
`15b45a33a7f8015f1579c103b55178197b0f5248ca9133bfff84298281fe9230`).

## Preregistered evaluation of the revised Linux candidate

After reviewing the two development splits, the candidate wheel is frozen
before scoring either new holdout. [Candidate build 44927371343](https://github.com/kent-tokyo/chematic/actions/runs/44927371343)
produced Linux CPython 3.12 wheel SHA-256
`046fd444e4131fbc5fabc4c1d14f5c9fe64b6525f42158abc961496f741a7898`
from merge commit `ac055be0c0262b32bfe52f7a9dd1cb3607cd1919`, parents
`f97aa95253f2dc38b0d98a58fde1f646febe8aca` and
`eea47ed329b1f67b7066005aadbcb2d8ebd1a4ab`. The executable wheel, rather
than a later rebuild, is the input to the evaluation.

Evaluate all 430 complete-Platinum holdout records and all 188 Prime holdout
records, separately, with RDKit 2026.09.2 and the same 20-candidate,
seed-42, maximum-200-iteration MMFF94, corrected-embedding and optional
400-iteration retry protocol used above. Score symmetry-aware heavy-atom
RMSD (up to 100,000 graph matches), specified stereochemistry and the same
heuristic geometry screen. Keep all records and count errors as failures;
report engine-specific errors and convergence-required recovery. The split
SHA-256 values are `00a3871e2157745dd1a90efeaa0c357912c5e04ec97a668b8c926c9815e63b0b`
and `52cb71f7b9d45b6085251fdd72c689d9f99f3f281ef65098e2446748e42cd330`.
Use paired molecule and scaffold-cluster bootstrap intervals. The fixed-count
noninferiority criterion for each population is a lower 95% interval bound
above −2 percentage points under both resampling methods; otherwise the test
is inconclusive or fails that criterion, not proof of equivalence. Repeated
Prime records may describe one molecule, making cluster intervals important.

Separately measure both engines using a 0.2-second per-engine deadline on
complete-Platinum and a 2-second deadline on Prime, both with at most 20
attempts and the original seed schedule. Prime's longer budget is selected
before scoring because its development molecules are substantially slower.
Report completion, recovery and convergence, molecule and cluster intervals,
errors and host details. The same −2-point lower-bound criterion applies to
each fixed-time comparison. Measurements are limited to the recorded host,
retry settings and populations. The geometry checks are flags, not a proof
of physical validity. Unknown historical exposure, unsupported chemistry,
equal-length cycle order and the relevant-cycle cap remain boundaries.

Before either new holdout result is opened, run the same two protocols on a
separately frozen macOS arm64 CPython 3.12 source wheel from the same candidate
build and merge commit. Its SHA-256 is
`5e1f34dab867a41a97fe66ea42d01f27199fd88325c00ae36ddf229bb9406920`.
The Linux and macOS runs use the same molecules and are platform reproductions,
not independent molecular samples. Report each platform separately. Neither
source wheel is evidence for a published package.

The revised Linux candidate's complete-Platinum fixed-count job processed all
430 records without reported engine errors or stereo violations. Chematic
recovered 189 and RDKit 191; the paired molecule 95% interval was −2.56 to
+1.63 points and the scaffold-cluster interval was −2.52 to +1.61 points.
The prespecified −2-point noninferiority criterion was not met; the interval
does not prove a difference in either direction. Converged matching recovery
was 97 and 103. The raw fixed result SHA-256 is
`01dda14cc25da574b9736363fb85f8c44e95a61aae136d78d7e951039ba756fc`, and
its committed gzip SHA-256 is
`e7392a2d6c55093457bfa40d684dcba28163991378e70e6aa2445de250781e13`.

The time-budget step did not begin: its shell variable name collided with
Bash's special `SECONDS` variable, causing the requested 0.2 seconds to be
parsed as zero. The recovery run uses a different variable name, the exact
same frozen wheel and source/SDF hashes, and the fixed-result artifact above.
It applies the previously registered 0.2-second protocol without selecting
or excluding any records based on the fixed-count outcome. Its host is a
separate hosted runner, so the timing is an independent run on the same
population, not a continuation of the original process.

## Additional reference candidate after the inconclusive result

Because the revised candidate's 430-record Linux fixed-count interval did not
meet the preregistered margin, an additional, post-hoc validation population
is being prepared from the [processed Crystallography Open Database archive](https://github.com/dkoes/conformer_analysis#raw-data).
This does not change or replace the result above. The preparation workflow
selects the 12,000 SDF members with lowest SHA-256 of their archive paths,
without examining conformer outcomes; it records the archive and sample
hashes. Before scoring, exclude the full 4,548-record Platinum dataset, the
full Prime dataset, and the other known exposure manifests by scaffold, audit
invalid exclusions, set the candidate limit, split, scoring and acceptance
rules, and retain a separate development partition. The archive is 510,146,560
bytes and is downloaded only on the hosted runner. Any later result must be
reported as an additional test initiated after the inconclusive 430-record
result, not as if it were part of the original preregistration.

The [Linux time recovery](https://github.com/kent-tokyo/chematic/actions/runs/44929831053)
completed all 430 records with no reported per-row errors or specified-stereo
violations. At 0.2 seconds chematic recovered 113 and RDKit 101; the paired
molecule interval was +0.93 to +4.88 points and the scaffold-cluster interval
was +0.91 to +4.82 points. Both exceed the preregistered −2-point lower-bound
criterion. Converged matching recovery was 75 and 70, with scaffold-cluster
interval −0.47 to +2.80 points. These do not change the fixed-candidate
noninferiority outcome. The compressed row data SHA-256 is
`31c751ad296878d0b3f5c35b31ec50760725d77d173d5af5b42303f53fff411d`.

The revised macOS candidate's complete-Platinum fixed-count run also handled
all 430 records with no reported engine errors or specified-stereo violations.
Chematic recovered 191 and RDKit 192; the molecule interval was −2.33 to
+1.86 points and the scaffold-cluster interval was −2.40 to +2.02 points.
The −2-point noninferiority criterion was not met. Requiring convergence
produced 98 and 105 recoveries, with scaffold interval −3.03 to −0.45 points.
The raw fixed-result SHA-256 is
`4a85751faf74d9599a4a96d459bce43d270302a049cb8442f4da8a21d83996d1`; the
committed gzip SHA-256 is
`84c429ebe6c31087214506538868dff4cc453d0c6559539e2e046adf920e28a8`.
Its time step hit the same shell-variable collision before scoring and is
recovered separately from the exact fixed-result artifact and frozen macOS
wheel, using the already-registered 0.2-second protocol.

[Preparation run 44930127544](https://github.com/kent-tokyo/chematic/actions/runs/44930127544)
identified the source archive SHA-256 as
`e1d6191fd665dbcc8ce5a24731fefe6ecf6460e0d25f5e290e264646543a16a2` and
the 12,000-record sample SDF SHA-256 as
`6ad24cf292a62561e55d3556395b841440f9d0d8588c7f2b85c8fb7e4675af74`.
The full list of selected archive members is preserved in
`validation/results/cod-reference-sample-manifest-2026-11-22.json.gz`.
The sample contained 55 records that RDKit 2026.09.1 could not parse; they
are reported separately and cannot be assigned scaffolds. Of the remaining
records, 3,008 match known exposed scaffolds. Within the other records,
200 development and 2,000 holdout scaffold groups were chosen by hashes of
partition and scaffold, retaining every selected group's records. One holdout
scaffold also matched the conservative invalid-exclusion graph audit and was
removed. The final split contains 223 development records and 2,263 holdout
records; the rest remain unused. The initial and final audit records and split
are stored under `validation/results/cod-*2026-11-22.json*`.

Run the 223 development records first. Measure MMFF and UFF support separately
for each engine, so unsupported molecules and graph mismatches cannot be
mistaken for equal successful coverage. Before evaluating the 2,263 holdout
records, fix the scorer, eligible population, handling of one-sided support,
comparison denominator, fixed and timed budgets, artifact hashes, and
bootstrap/noninferiority rule. The earlier 430-record criterion remains unmet
regardless of this additional assessment.

The [macOS time recovery](https://github.com/kent-tokyo/chematic/actions/runs/44931051084)
completed all 430 records without reported per-row errors or specified-stereo
violations. Chematic recovered 104 and RDKit 82 within 0.2 seconds; the
paired molecule interval was +2.79 to +7.67 points and the scaffold-cluster
interval was +2.75 to +7.58 points. Both exceed the registered −2-point
lower-bound criterion. Converged matching recovery was 72 and 59, with
scaffold interval +0.91 to +5.33 points. The compressed row data SHA-256 is
`e18f9025b7046fcb02f231e9786567c03342faf86ca1076eef048724422ad751`.
These timed results do not alter the fixed-candidate criterion's failure.
