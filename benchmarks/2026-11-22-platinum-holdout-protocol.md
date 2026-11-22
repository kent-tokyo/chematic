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

The revised candidate's Prime fixed-count runs each processed all 188 records
without reported engine errors or specified-stereo violations, with complete
bond and angle screen references. Linux recovered 59 for chematic and 61 for
RDKit; its molecule interval was −4.26 to +2.13 points and scaffold interval
−4.79 to +2.29 points. macOS recovered 60 and 61, with molecule interval
−3.72 to +3.19 points and scaffold interval −4.37 to +3.14 points. Neither
platform met the prespecified −2-point fixed-count noninferiority criterion.
Converged matching recoveries were 40/40 on Linux and 39/38 on macOS. Both
engines recovered roughly one third of the experimental reference records;
this population remains difficult for both implementations.

The Linux fixed raw SHA-256 is
`fc45178879a365f9781a601c16b8228836845f67563ec683baa3668715fdc8df`
(gzip SHA-256 `8ad6d4da273a415bca40aa4e1848741c157dc77e7af299f17a8066a06906c7d7`);
the macOS raw SHA-256 is
`3ce64efde873fcde6a265a94d1346a7e78adaf210dee36661ad79feab0ffc2d5`
(gzip SHA-256 `88611204e5d9f871bdeca1c7d62f809ff23a7a5ca61045c0d09011c9345a466c`).
Their time steps also stopped before scoring due to the `SECONDS` collision.
Separate recovery jobs use the corresponding frozen wheel, fixed-result hash
and the preregistered 2-second protocol on the same 188 records.

The [Linux COD development run](https://github.com/kent-tokyo/chematic/actions/runs/44931051201)
used the frozen revised Linux wheel. Of 223 records, both engines constructed
MMFF for 210, neither constructed it for 12, and only chematic constructed it
for one. UFF construction succeeded on both for 222 records and on chematic
alone for one. All 223 molecular graphs matched. These are constructor checks
with synthetic coordinates, not proof that every assigned parameter is
chemically correct. On the common MMFF domain, the screened fixed-candidate
recovery was 171/210 for chematic and 173/210 for RDKit; the paired
molecule-bootstrap 95% interval was −2.38 to 0 points. No specified-stereo
violations were detected among the 4,184 recorded minimized conformers per
engine. The 223-record full-denominator summary is also retained, but it must
not be interpreted as equal successful coverage for the 13 RDKit-unsupported
records. Raw development result gzip SHA-256 is
`4eb69617a5598e63126d374196bb98007ee8fd5632920bd6fff39fbe11499285`.

### Preregistered additional COD holdout protocol

Before scoring any of the 2,263 selected COD holdout records, freeze the
Linux and macOS candidate wheels and the sample hash above. First run the
same support probe on all 2,263 records, independently recording each engine's
MMFF and UFF construction result, graph mismatches and errors. Publish all
support counts and asymmetric cases. The conformer-recovery denominator is
the subset for which both engines construct MMFF and the molecular graphs
match; unsupported, one-sided or mismatched cases are excluded from that
performance denominator but remain in the separately reported full support
population. No reference geometry or generated conformer outcome influences
this eligibility filter. `scripts/validation/cod_eligible_split.py` verifies
the support record identity and derives a subset of the previously clean
exclusion audit without introducing any records.

On the common domain, use the same 20-candidate, seed-42, MMFF94 iteration
and stereo/geometry-screen protocol as the earlier holdouts. Also use the
same 0.2-second per-engine budget and maximum 20 attempts. For each platform
and each budget, report molecule and scaffold-cluster paired bootstrap 95%
intervals, candidate failures and convergence-required recovery. A lower
bound above −2 percentage points under both resampling methods is the
prespecified noninferiority criterion. Report both platforms separately; they
are repeated measurements of the same molecules. Compare MMFF94, MMFF94s and
UFF energies and gradients on common generated coordinates, reporting
construction failures and numerical outliers rather than only percentiles.
The maximum 10,000 relevant-cycle enumeration fallback and equal-size ring
ordering remain unproven boundaries. This additional COD assessment was
initiated after the inconclusive 430-record test and cannot retroactively
make that original criterion pass.
The support-probe SHA-256 is
`606d5d3f8b2a8cb17a3fe5534f9a1b9f09257eb4e8a9ae16c32e41fc7d474dcc`; the
common-domain selector SHA-256 is
`d2825bd9951a939333bf5f8ad926883f2801318f1b66ba68b28c68634d757d15`.
The decompressed split SHA-256 is
`138fff6f10426bf4041faae25fc07717e436eb7f965e2df65c5983cde65e3449`.

A supplementary identity audit was added after the COD jobs started because
their original graph check compares atom numbers and bond connectivity/order,
but does not itself compare per-atom formal charge or aromaticity. An audit of
all 2,263 selected holdout records on the local macOS wheel from source commit
`d30d186` found matching atomic number, formal charge and aromaticity at every
atom. The parsing and hydrogen-addition sources are unchanged between that
wheel and the revised candidate; the force-field code differs. One record
contains isotopic atoms, and the two implementations produced the same
isotope-aware molecular formula; atom-by-atom isotope identity was not
independently checked. This retrospective input audit did not change the
registered eligibility or scoring rules. Its records are in
`validation/results/cod-holdout-atom-identity-audit-2026-11-22.json.gz`
(SHA-256 `9b8c17f2afbd8c3746f102014c3a2517117ec132c1305bd3a14a1e3da5c95010`).

The [Prime Linux time recovery](https://github.com/kent-tokyo/chematic/actions/runs/44931816542)
completed all 188 records without reported errors or specified-stereo
violations. At the preregistered 2-second budget both engines recovered 59
records. The molecule interval for their difference was −2.66 to +2.66 points
and the scaffold-cluster interval was −3.61 to +3.26 points, so the
preregistered −2-point noninferiority criterion was not met. Converged
matching recovery was 40 for each engine. The gzip row data SHA-256 is
`69a206dfdd6b08df49efec4f948da1763334d17da623f208b1284f930050f11d`.

A separate regression fixture contains 30 deterministic atom-order variants
of the ten already-exposed type-76 molecules, with per-atom RDKit 2026.09.2
types and partial charges. It tests atom-order sensitivity on both frozen
candidate wheels. It is development evidence and does not cover the untested
case of equally short relevant rings with different type-76 membership.

The [Prime macOS time recovery](https://github.com/kent-tokyo/chematic/actions/runs/44931816522)
completed all 188 records without reported errors or specified-stereo
violations. At 2 seconds chematic recovered 60 and RDKit 58; the molecule
interval was −2.13 to +4.26 points and the scaffold-cluster interval was
−2.69 to +4.71 points. The preregistered −2-point noninferiority criterion
was not met. Converged matching recovery was 39 and 38. The gzip row data
SHA-256 is `d9f74cb2f6f6b70c8dc6bd9a7fef5730ebc7387a64ee086d46a7afd94909a283`.

A separate RDKit 2026.09.2 charge-sum audit over the 1,606 exposed Platinum
records found two cases where MMFF partial charges differ from the input
formal charge by more than `1e-6` elementary charge: `NR4_3L4Y_A` (formal 0,
MMFF sum −1) and `PP9_1Z53_A` (formal −4, MMFF sum −6). All 1,606 were
processed. This is a property of the measured RDKit reference behavior and is
not evidence that its charge assignment is physically valid. The audit is
reproducible with `scripts/validation/rdkit_mmff_charge_conservation.py`;
its compressed per-record data SHA-256 is
`b286358190270d1d541c1bd2d1b4c5d05c6fbe42c525cada7010f37c252d4f8e`.
The [atom-order regression run](https://github.com/kent-tokyo/chematic/actions/runs/44933256225)
passed all 30 variants on both frozen Linux and macOS candidate wheels.

To investigate the low Prime recovery without tuning on its exposed holdout,
a separate development-only run generates up to 100 candidates for the same
31 development records on the frozen Linux candidate and RDKit. Both-missed
records in the 20-candidate development run numbered 19; simple RMSD-only
and geometry-screened recovery counts were identical for both engines. The
100-candidate run is exploratory: check whether the first 20 candidates
actually match the prior run before interpreting a difference as a nested
budget effect, and do not treat it as evidence on the Prime holdout.

## Aromatic silicon follow-up (experimental)

A supplementary constructor audit on the selected COD records, using the
older local macOS wheel, found one molecule for which RDKit constructs MMFF
but chematic rejects aromaticity re-perception: an aromatic negatively charged
silicon connected to positively charged nitrogens. RDKit assigns a double
bond to the silicon; the current core candidate rule excludes all anionic
silicon from matching. This support observation has exposed that COD input,
so any revised implementation must use another untouched evaluation set.

The targeted valence rule was first recorded as an experimental patch in
`validation/experiments/silicon-anion-kekulization.patch`. In an isolated
copy of `chematic-core`, a synthetic silaanion regression and all 171 core
unit tests passed. [Experimental CI](https://github.com/kent-tokyo/chematic/actions/runs/44935610099)
then confirmed atom types, charges, and common-coordinate MMFF94/MMFF94s
energy and gradients against RDKit for two synthetic examples and the exposed
COD input. The rule is now included in the PR source with regression tests.
The earlier frozen wheels and their COD holdout measurements predate it; a
new candidate and a separate untouched validation set are required.


The [original COD Linux run](https://github.com/kent-tokyo/chematic/actions/runs/44932574044)
completed on the earlier frozen wheel, after its support audit had exposed a
holdout input. On 2,093 common-MMFF inputs, fixed 20-attempt recovery was
1,669/1,664 for chematic/RDKit; the paired molecule 95% interval for the
difference was −0.38 to +0.86 percentage points and the scaffold-cluster
interval was −0.38 to +0.87 points. At 0.2 seconds recovery was 1,333/1,269;
the respective intervals were +2.10 to +4.06 and +2.07 to +4.06 points.
These bounds exceed the prespecified −2-point threshold on this host, but do
not establish overall parity: the fixed run returned 41,440/41,436 candidates
with 75/67 specified-stereo violations, while the timed run had 22/15. The
newer stereo changes are not included. This result cannot be considered a
fully unexposed holdout after the earlier support investigation. The original
macOS job reached its six-hour limit after writing a prefix of 1,808 out of
2,093 fixed-score records. Its saved checkpoint SHA-256 is
`72d2ccabd7bba29105d0cc6fda9894a9826a075397a09665df8330ad9ffe4211`.
The remaining 285 records are resumed with the same wheel and scoring
configuration; the join rejects mismatched metadata, incomplete records and
non-macOS-arm64 hosts. The combined fixed result spans two identified hosts,
and the timed comparison runs separately on a third host. No macOS result is
claimed until all records are complete and validated.

On shared RDKit-embedded coordinates, the same Linux artifact reports 2,071
force-field comparisons and 22 row errors. Maximum absolute differences were
57.67 kcal/mol for MMFF energy, 121.10 for MMFF gradient components and 75.12
for UFF gradient components; maxima are not typical-case estimates. The
corresponding raw fixed, timed and force-field JSON SHA-256 hashes are
`8e15c8afc4e7a242184b1b82f10376836850e616a57493bbdaafa3bf02d57190`,
`028cdc503828489ba632c3ce8a09c83ad0402215f9366c1d817bffb93b83652c`, and
`886bde8880e77d3ed30455d04084a4ec65755d4326314f2452a1b4c100d4b075`.
The artifact retains the raw records and summaries. Results are for the older
wheel and do not predict the separately running updated-wheel assessment.

The first holdout's UFF gradient outliers were subsequently exposed as
development cases. On three phosphorus-containing structures, the earlier
wheel differed from RDKit 2026.09.2 by 47.75–74.22 kcal/mol/Å in a gradient
component, while RDKit agreed with a central difference of the energy. The
implementation retained a P/As/Sb/Bi inversion-gradient sign from RDKit
2026.03.1. RDKit corrected that sign in
[issue #9229](https://github.com/rdkit/rdkit/issues/9229), released in
2026.03.2. A subsequent source change updates the sign and tests against
2026.09.2 plus numerical derivatives. These three exposed structures are
development diagnostics, not fresh holdout evidence; the change is not in the
wheel currently being scored on the second COD holdout.

## Second COD scaffold split for the silicon-corrected candidate

After the first COD support audit exposed a holdout molecule, a second split
was selected from the original sample without scoring conformer performance.
Whole scaffold groups were ordered by SHA-256 of `round2:<partition>:` plus
the scaffold hash, taking 200 development and 2,000 holdout groups from the
previously unused partitions. The 2,486 previously selected development or
holdout records and the known Platinum and Prime scaffolds remain excluded.
The initial split has 225 development and 2,290 holdout records (uncompressed
JSON SHA-256 `52b37893b77224a4dd187b3d3c13ecae2898a8f0c09af0bf3ceb5be6bf80c330`).
The conservative invalid-input graph audit found two possible overlaps in
the holdout and excluded their scaffold groups before scoring. The final
split contains 225 development and 2,288 holdout records; its uncompressed
JSON SHA-256 is `b12af9cb05683735ab8fad725dcf38e0f8e8f51bdef8e4aa3db22e87590109c3`.
No scaffold occurs in more than one of the new development, new holdout,
and previously evaluated partitions. Of the 12,000 sampled archive members,
55 failed parsing and are outside these partitions. The selection and both
invalid-input audits are retained under `validation/results/`; the selector
is `scripts/validation/select_cod_scaffolds_round2.py`. These public COD
structures may have been used elsewhere; independence from unknown prior
training or tuning cannot be established. This additional dataset is a
post-hoc follow-up, not a replacement of the original preregistered tests.

The [100-candidate Prime development run](https://github.com/kent-tokyo/chematic/actions/runs/44933792400)
completed 31 records using the earlier frozen Linux wheel. Recovery increased
from 11 to 14 for chematic and 8 to 11 for RDKit. For both engines, the first
20 RMSDs, statuses, stereo checks and geometry checks matched the earlier
20-candidate run on every row, so this is a nested comparison. Among 3,100
minimized candidates per engine, each engine had one specified double-bond
stereo mismatch, on different structures. These results are exploratory and
do not establish holdout performance or stereo safety for larger ensembles.

The [silicon-corrected Linux candidate build](https://github.com/kent-tokyo/chematic/actions/runs/44936328481)
uses PR merge commit `b4852cefede06d4a864cec87684d5fea7d316304` with parents
`f97aa95253f2dc38b0d98a58fde1f646febe8aca` and
`567bf5e43eaa3d073c7be389dc2032940407e335`; its crates tree is
`b08ecd24da790430dad6eb0ff3a509981f78f8bf`. Its Python 3.12 Linux wheel
SHA-256 is `35e6bdf3a3c37e082c7d83d31414c37023b6c535687ceb28c997cc81c37a5b8d`.
A development-only comparison is configured for the 225 round-two development
records, plus a Linux reproduction of the exposed Prime stereo discrepancy.
No round-two holdout conformer outcome has been inspected.

The [Linux round-two development run](https://github.com/kent-tokyo/chematic/actions/runs/44936872973)
completed all 225 selected records with the silicon-corrected, pre-stereo-fix
wheel. Atomic-number, formal-charge and aromaticity identity matched on all
225; no isotopic atoms occurred. Both engines constructed MMFF for 212 and
neither did for 13. At 20 attempted candidates, with all 225 records in the
denominator, chematic and RDKit recovered 161 and 163 references; the paired
molecule bootstrap 95% interval for the difference was −4.00 to +2.22
percentage points. This development result is not a holdout noninferiority
claim. Lowest-energy recovery was 104 and 109. Of 4,180 minimized candidates
per engine, the stereo checker reported violations for 40 chematic and 38
RDKit candidates, all on the same two records and all with unfinished status.
The two input formal charges were −51 and −2; these records are retained in
the reported counts. The raw RDKit embedded candidates had no reported
stereo violations; this profile did not collect raw chematic candidates.
Source hashes of the uncompressed row, support and atom-identity files are
`5b58fa73ade8a92ce14b0d0253e1dbfd11ba6f4e52f73cdd0a98b186dc08770f`,
`53626b44dbcda1832f4744a42100fd92d59579f2b50d212d7cf43df1ea8e20c1`, and
`8860ef9ceaada101a93eae60fe95493193bfccf4cae6bfb7f01b3ea79b3cc202`.
Compressed evidence is retained as `validation/results/cod-linux-development-v2-*`.

A retrospective SMILES-direction audit on those 225 development records
identified 12 with an unspecified double bond adjacent to a directional bond
already used by a different specified alkene. The unmodified stereo verifier
reported an ambiguous direction for two inspected examples. A development
regression fixture records all 12 structures and the number of specified E/Z
bonds assigned by RDKit 2026.09.2; after the marker attribution change, the
macOS Rust test matches those counts without unsupported entries. The fixture
is development evidence, not an independent holdout or a complete stereo
validation of the dataset.

The macOS arm64 Python 3.12 wheel from the silicon-corrected merge commit has SHA-256
`b5c01642b5c545835071dedc808ad0d930865512444dd0f3b80e58c7d43b0201`.
The [independent Linux stereo diagnostic](https://github.com/kent-tokyo/chematic/actions/runs/44937409465)
reproduced the exposed Prime development mismatch at candidate index 43:
the specified E double bond had a −0.49-degree dihedral (Z-like) after the
MMFF stage, with nonconverged status 1. The same case was not reproduced in
the earlier local macOS run. A [second Linux diagnostic](https://github.com/kent-tokyo/chematic/actions/runs/44937817100)
measured the same seed's embedded dihedral at −177.49 degrees before MMFF,
confirming that the specified E bond flipped during the MMFF stage on that
runner. The checked profile now attempts a stereo-valid embedded fallback,
marks it unfinished, and omits a candidate when both geometries fail its
stereo check. The subsequent Linux regression and independent RDKit checks are
recorded below.

## Prospective second COD holdout protocol after stereo changes

The [stereo-corrected wheel build](https://github.com/kent-tokyo/chematic/actions/runs/44950293868)
used merge commit `9db581ac6d2cab75a9fb5809fae1a5ee60870054`, parents
`f97aa95253f2dc38b0d98a58fde1f646febe8aca` and
`91f53c1a60cd599bf92884c7db3a17b32c71776b`, crates tree
`1c8c06ddd74f06762e969948f2c08042aafd484f` and benchmark tree
`00a8c3615c849288c7d31f718a502bbd40da4fdc`. The CPython 3.12 Linux wheel
SHA-256 is `3db87234071e820b55c2acfce0b61d42f2edf870094f840d4d22ff9acd725957`;
the macOS arm64 wheel SHA-256 is
`438528ccacc94551233a35b3be3b7f914eb71f77f6e00ed0af9ad6c89d993037`.
These are source-built candidates, not published packages. The
[Linux regression](https://github.com/kent-tokyo/chematic/actions/runs/44950706502)
passed the targeted Rust tests and the independent RDKit check of all 44
generated Prime development candidates, with no specified-stereo mismatches.
The [COD development rerun](https://github.com/kent-tokyo/chematic/actions/runs/44952391922)
then completed on this wheel. Chematic and RDKit recovered 160 and 163 of
225 references, respectively; lowest-energy recovery was 103 and 109. Among
4,180 candidates per engine, 20 chematic and 38 RDKit candidates were flagged
by the RDKit stereo checker. The remaining chematic violations all came from
one −2 charged development molecule; every flagged candidate was unfinished.
These are development observations, not independent test results.

The follow-up investigation found that the affected double bond approached
linear geometry during optimization. A first five-degree guard also rejected
an existing Prime development conformer, where another substituent was far
from the axis. The revised guard keeps the marked substituent when another arm
is usable, and treats an end as unevaluable when every arm is within five
degrees of the axis. This angle is an implementation rule, not a demonstrated
replication of RDKit's stereo-perception threshold. The
[Linux regression](https://github.com/kent-tokyo/chematic/actions/runs/44957585587)
passed the Rust tests and an independent RDKit check of all 44 candidates for
the exposed Prime case and all 20 candidates for the charged COD case. These
two cases do not establish stereo correctness across the full development set.

The [updated wheel build](https://github.com/kent-tokyo/chematic/actions/runs/44958776980)
used merge commit `e4efaf29afd6bfa3689d1187dfc391343d7ceac0` and crates tree
`c953928e32e912a94b5bd4462c9b49e31f5461c8`, with benchmark tree
`00a8c3615c849288c7d31f718a502bbd40da4fdc`. The CPython 3.12 Linux wheel
SHA-256 is `e0f95eb569e14a5f69df04b566acbe0e86b86832afa18da079eedff6a3e82a78`,
and the macOS arm64 wheel SHA-256 is
`1a952db5bd458a4072df4b92f617b818d45373be2a000d6cfd9aafe7bf771001`.
The [225-record development rerun](https://github.com/kent-tokyo/chematic/actions/runs/44959301502)
completed. Atom identity matched all 225 records; both engines had MMFF
parameters for 212 and neither for 13. Chematic returned 4,080 candidates,
with no RDKit-detected specified-stereo violations; RDKit returned 4,180, with
38 violations. Both recovery counts remained 160/163 out of 225; lowest-energy
recovery was 103/109. The molecule-bootstrap interval for the former
difference was −4.44 to +1.78 percentage points, so this development result
does not establish noninferiority. Chematic returned no candidates for eight
of the 212 supported inputs, including five for which RDKit returned 20 each;
RDKit reported a `Bad Conformer Id` error on the other three. The same eight
inputs had zero chematic candidates in the previous wheel's development run.
The zero-violation count is therefore conditional on returned candidates and
must not hide those generation failures.

For those eight already-exposed records, a separate
[single-seed pipeline diagnostic](https://github.com/kent-tokyo/chematic/actions/runs/44962800236)
using the same frozen wheel reported unsupported stereochemical coordination
on all five records for which RDKit returned 20 candidates: two with
tetrahedral P bearing a double-bonded ligand, two sulfoxide S centers and one
square-planar S declaration. One sulfoxide record also had a separate
tetrahedral violation in this diagnostic. This pipeline and seed are distinct
from the corrected ensemble, so these observations explain an encountered
limitation without proving the cause of every failed ensemble attempt.
The result is retained in
`validation/results/cod-zero-candidates-diagnostic-2026-11-22.json`.
The subsequent phosphorus verifier change is a separate development
candidate and is not part of the frozen-wheel holdout evaluation below.

The uncompressed result SHA-256 is
`7419b7cbc20f5514699ba44ecb2c2eeb9a7f3624acaa2a3619a4d32aa72acaec`.
The retained compressed result, support, and identity evidence is
`validation/results/cod-linux-{development,support,identity}-v4-2026-11-22.json.gz`.
These wheels are source-built candidates, not published packages.

No conformer outcome from the 2,288 round-two holdout records has been
inspected. After the stereo regression passes, freeze and fingerprint new
Linux x86-64 and macOS arm64 wheels, their exact source and benchmark trees,
and the RDKit 2026.09.2 oracle. First repeat the 225 already-exposed
development records for the final candidate, including constructor coverage,
input identity, stereo violations, errors and convergence. Do not use their
results as holdout evidence.

On each platform, audit input identity and MMFF/UFF constructor support on all
2,288 holdout records before scoring, retain one-sided support and parsing
failures in a separate coverage table, and derive the common MMFF domain using
the existing eligibility and exclusion scripts. Stop before scoring if the
input identity audit reveals a mismatch requiring a protocol change. On that
common domain, use the same geometry and specified-stereo screening,
20-attempt seed-42 fixed protocol with 200 MMFF iterations and up to 400
retry iterations, and a separate 0.2-second per-engine, 20-attempt deadline.
Count runtime errors as failures, report convergence-required results and
stereo violations, and compare both paired-molecule and scaffold-cluster 95%
bootstrap intervals. The lower bound must exceed −2 percentage points for
both resampling methods on each platform and each budget. Also report
specified-stereo violations for every returned candidate; any violation by the
checked chematic profile is an open correctness failure. Score force fields on
common coordinates separately. This extra post-hoc test cannot make the
earlier Platinum or Prime criteria pass retrospectively and cannot establish
independence from unrecorded prior exposure.

To keep individual scoring jobs bounded, partition the eligible records into
four shards by the integer value of each scaffold SHA-256 modulo four. Every
record in a scaffold goes to the same shard. Run both engines for each record
on the same runner, preserve the host and result hashes for each shard, and
recombine all four before the molecule and scaffold-cluster bootstrap. A
missing, duplicate, misassigned, or differently configured shard makes the
aggregate invalid. Runtime results are therefore a paired comparison across
four separately identified hosts per platform, not one homogeneous host.
The timed scorer alternates engine order within each shard, restarting the
alternation at each shard boundary; per-record order remains in the results.
Report recovery separately for each shard as well. The combined intervals
condition on these particular runners and do not quantify uncertainty from
host-to-host variation.
This partitioning is fixed before any second-holdout conformer score.
Before scoring, a local integration check split and reassembled the 225
already-exposed COD development rows: every row was identical, and the
existing summarizer retained recovery counts of 161 and 163 across 200
scaffolds. A separate check reassembled all 430 previously-exposed Platinum
Linux time rows identically and retained the 385-scaffold recovery difference
of +2.79 percentage points. These are tests of the merge procedure using old
results, not new holdout measurements. Synthetic tests also verify rejection
of missing, duplicated, misassigned or differently configured shards.
