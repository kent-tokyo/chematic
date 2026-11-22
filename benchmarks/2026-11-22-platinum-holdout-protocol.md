# Platinum Diverse holdout protocol

**Stereo measurement limitation:** Results recorded below before the strict
verifier was introduced check tetrahedral atoms only where RDKit assigned a
CIP label, together with detected E/Z bonds. They do not establish that every
declared tetrahedral or non-tetrahedral center was verified. In particular, a
zero reported count is not a comprehensive stereo guarantee. The strict
verifier and its separate unknown/contradicted counts are being evaluated on
development inputs before any new holdout scoring.

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
and the timed comparison runs separately on a third host. At that checkpoint no macOS result was claimed; the completed result is recorded below.

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
wheel currently being scored on the second COD holdout. The combined
[Linux regression](https://github.com/kent-tokyo/chematic/actions/runs/44965076304)
passed the UFF contract and boundary tests, then compared those three
structures against RDKit 2026.09.2. The maximum energy and gradient differences
were zero on that runner; the largest phosphorus finite-difference discrepancy
was 3.62e−8 kcal/mol/Å. These are three exposed development cases, not a
general UFF guarantee. Results are in
`validation/results/uff-inversion-development-2026-11-22.json`.
For a broader development-only recheck, 51 already-exposed first-COD records
were selected if their earlier MMFF94 energy difference exceeded 1e−4
kcal/mol, their MMFF94 finite-difference discrepancy exceeded 1, or their UFF
gradient-component difference exceeded 1. The input list, rule and hashes of
the earlier source files are retained in
`validation/results/cod-first-field-outliers-development-2026-11-22.json`.
These cases were selected because of observed errors and are not a random
sample or independent test set.
The [v5 Linux recheck](https://github.com/kent-tokyo/chematic/actions/runs/44966536889)
completed all 51 without row errors. UFF energies and analytic gradients
matched RDKit on these inputs; the largest of the six sampled finite-difference
gradient discrepancies per molecule was 6.09 kcal/mol/Å. MMFF94 still had
maximum energy and gradient differences of 57.67 kcal/mol and 121.10
kcal/mol/Å, respectively, and a sampled finite-difference discrepancy as
large as 279.13 kcal/mol/Å. These maxima are not representative percentiles.
Per-record output is retained as
`validation/results/cod-first-fields-v5-development-2026-11-22.json.gz`
(uncompressed SHA-256
`df7958d0379175bccf9b9d228ee9cfa0f27916918d9e7db730d90470b4a4a6f6`).
Some of the MMFF energy outliers exposed a missing one-neighbor carbon
type-60 branch, which has a subsequent development fix and regression job;
that fix is not included in the v5 results. In the subsequent
[Linux regression](https://github.com/kent-tokyo/chematic/actions/runs/44967041251),
three exposed structures matched RDKit 2026.09.2 in atom types, partial
charges, MMFF energy and gradient on the tested coordinates. The previous
wheel's energy differences on those structures were approximately 27.84,
33.14 and 57.67 kcal/mol. The combined phosphorus and UFF development checks
also passed; this is targeted evidence, not a general MMFF claim. The result
is `validation/results/mmff-terminal-carbon-development-2026-11-22.json`.
In the [v6 51-case recheck](https://github.com/kent-tokyo/chematic/actions/runs/44969183955),
the number with MMFF94 energy error above 1e−4 kcal/mol fell from 19 to 11;
the maximum fell from 57.67 to 41.52 kcal/mol. UFF analytic gradients still
matched RDKit on all 51, while the largest sampled finite-difference discrepancy
remained 6.09 kcal/mol/Å. The uncompressed result SHA-256 is
`6c7a9e723c672b43391ea3dbf0dc4042359ef6f6a471e9c51428a308d4d04ea1`;
the per-record compressed result is
`validation/results/cod-first-fields-v6-development-2026-11-22.json.gz`.
Another exposed case had matching atom types and charges but a 41.52 kcal/mol
energy difference on the original Linux coordinates. Its type-75 phosphorus
was missing from the special single-bond classification, changing the angle and
stretch-bend parameters. A subsequent
[Linux regression](https://github.com/kent-tokyo/chematic/actions/runs/44968222238)
passed a registry-consistency test and matched RDKit energy and gradient on
that development structure; the result is
`validation/results/mmff-phosphaalkene-development-2026-11-22.json`. This
change is not included in the frozen v6 wheel or the running v4 holdout.
In the [v7 51-case recheck](https://github.com/kent-tokyo/chematic/actions/runs/44971010367),
seven inputs still had MMFF94 energy error above 1e−4 kcal/mol, with a maximum
of 32.94 kcal/mol. All 51 UFF energies and analytic gradients matched RDKit
on the tested coordinates. The compressed per-record output is
`validation/results/cod-first-fields-v7-development-2026-11-22.json.gz`
(uncompressed SHA-256
`a58c583e3e8204753deab5f373fe3261caec792fa176436ecef99177ddb742f7`).

Two of the seven MMFF outliers had directional single bonds adjacent to an
N-oxide; treating those bonds as single during MMFF typing corrected the
oxygen type. In the [subsequent regression](https://github.com/kent-tokyo/chematic/actions/runs/44970866778),
both exposed molecules matched RDKit 2026.09.2 atom types, charges, MMFF
energies and gradients at the tested coordinates. The combined earlier
regressions also passed. Evidence is
`validation/results/mmff-directional-nitrone-development-2026-11-22.json`.
This is another targeted development fix and is not included in v7 or the
running v4 holdout.
The [v8 51-case recheck](https://github.com/kent-tokyo/chematic/actions/runs/44972185553)
left five inputs with MMFF94 energy error above 1e−4 kcal/mol, down from seven
in v7; the maximum remained 32.94 kcal/mol and the maximum gradient-component
difference was 101.80 kcal/mol/Å. UFF energies and analytic gradients matched
RDKit on all 51, with the same 6.09 kcal/mol/Å maximum sampled finite-difference
discrepancy. This selected set remains development evidence. Per-record output
is `validation/results/cod-first-fields-v8-development-2026-11-22.json.gz`
(uncompressed SHA-256
`e45dbb23a0ccb87f00467ae9b6eaffa5734a115184c2b1fff2e4c49000e6593f`).

Among the five remaining v8 energy outliers, three exposed development
inputs have the same ring sets as RDKit 2026.09.2 but a different iteration
order in chematic's 2026.03.6-based ring port. Supplying the RDKit ring order
directly to the MMFF aromaticity view changed the problematic aromatic atom
and bond states in each of these three cases. This is a diagnostic, not yet
a tested force-field correction. RDKit 2026.09.2 uses RingDecomposerLib
relevant cycles for its default symmetrized ring information; chematic's
relevant-cycle enumeration also returns a different order for two of the
cases, so replacing the ring set or sorting by size alone is insufficient.

The largest sampled v5 finite-difference discrepancy occurred on an exposed
structure with a 179.998-degree bond angle. In the separate
[Linux diagnostic](https://github.com/kent-tokyo/chematic/actions/runs/44969777281),
the RDKit-compatible native and RDKit gradients on the worst sampled component
were both −166.015 kcal/mol/Å. Central differences changed from 113.118 at a
1e−5 Å step to 106.337 at 1e−7 Å; the existing opt-in consistent-plane
gradient was 106.338. Most of the analytic discrepancy was in the out-of-plane
term. These are one geometry and one component near a singular configuration,
not a broad claim that all default gradients are inaccurate. The compatible
optimizer still uses the RDKit-compatible gradient; the consistent-plane
gradient is currently diagnostic only. Full term and step-sweep results are
`validation/results/mmff-near-linear-gradient-development-2026-11-22.json`.

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
On the same combined Linux development wheel, the two phosphorus-containing
records that previously returned no candidate returned 20 each. RDKit found
no specified-stereo violations among those 40 candidates; 11 and 20,
respectively, had unfinished status. The results are in
`validation/results/phosphoryl-candidate-development-2026-11-22.json`.
The sulfur cases remain outside this change.

The phosphorus and UFF fixes were selected from already-exposed development
inputs and the first COD evaluation, before inspecting any conformer scores
from the running second COD holdout. A new candidate must be built from the
tested source, fingerprinted on both platforms, and checked again on all 225
development inputs before any prospective follow-up on that holdout. Any
such follow-up is related to the already-running v4 assessment and must report
the complete set of candidate results, rather than select a favorable one.
The [v5 build](https://github.com/kent-tokyo/chematic/actions/runs/44965786491)
uses merge commit `acc79a611d5c3bac571a7b1a6f274700b20f534e`, crates tree
`ebe322bf287fd58e9731bee836f682b19efd7a7e`, and benchmark tree
`00a8c3615c849288c7d31f718a502bbd40da4fdc`. Its Linux CPython 3.12 wheel
SHA-256 is `e113e7f74b6e5018a03913836824c203d2366e34189264c0d85291fbd663b7b3`;
its macOS arm64 wheel SHA-256 is
`402918d1553f4760d87eb5e3f1602f91bc3c37c164c8b01ac7cb048f4deeb464`.
The [v5 development comparison](https://github.com/kent-tokyo/chematic/actions/runs/44966536884)
matched atom identity on 225 records, with 212 common-MMFF inputs and 13
unsupported by both engines. Chematic/RDKit recovered 161/163 references out
of 225; the paired molecule-bootstrap interval was −4.00 to +2.22 percentage
points. Among 4,120/4,180 returned candidates, specified-stereo violations
were 0/38. Chematic returned no candidates for six inputs: three for which
RDKit returned 20 candidates and three for which RDKit reported a conformer
error. This is an exploratory development result, not noninferiority evidence.
Its uncompressed JSON SHA-256 is
`9c4b6fec1eeb8798fc58c5a8f262f0842fc2a78e30215bb2c2bb232e9f7af9b9`, and the
compressed result, support and identity records are retained as
`validation/results/cod-linux-{development,support,identity}-v5-2026-11-22.json.gz`.
The
[v6 build](https://github.com/kent-tokyo/chematic/actions/runs/44967752750)
includes the terminal-carbon fix and uses merge commit
`76a6d3c8ef9642155821d6810a1a0d28c6bfb146`, crates tree
`c23b1b1353602e84e5923d7b507f0aafacd05cfa`; its Linux and macOS CPython
3.12 wheel SHA-256 hashes are
`c5f0941fa318d22b280d172585c944052c6aafa74bdab995ca29d35ef04a5aea` and
`08179d1a33f4af6296c064831eecbd85de66c16529b63309b3d539a98587e23a`.
It was fixed before inspecting any second-COD holdout conformer scores and
requires its own development rerun. The
[v7 build](https://github.com/kent-tokyo/chematic/actions/runs/44969286762)
also includes the phosphorus-bond fix, and was fixed before inspecting those
scores. It uses merge commit `cd8a7f266a40531fa203dc7ec3434bd22b9f5834`
and crates tree `4871f3f53cb1ceb3940de32e744e223ae3d641b1`; its Linux and
macOS CPython 3.12 wheel SHA-256 values are
`cb9a9ca17f377693ebbbc8ba176e7da8defc8185d2e7cd3abdc77f4d5d7ec295` and
`accff6f48fa6567e94197dab8793e0e0e94325f12d22e628344b27fa992fba7d`.
The directional-bond fix was likewise selected from already-exposed inputs,
before inspecting any second-COD holdout conformer results, and a further
candidate including it was [frozen separately](https://github.com/kent-tokyo/chematic/actions/runs/44971533548).
Its merge commit is `98a270e4ba39c7f7db61a0e1980bfa29fd55f615`, crates tree
`46563b2448e2176184d826187ac4b570ad51d52d`, and Linux and macOS CPython 3.12
wheel SHA-256 values are
`8305fd27765c4d10670fa2a0900f3b86b286d35a711f3f2eecb5c412f96613ac` and
`b37bfac64e9ce9550eea1e36eb0e7554294b71f7e9efb22be470da44a5ef9608`.
The [225-record v8 Linux development evaluation](https://github.com/kent-tokyo/chematic/actions/runs/44972185356)
completed: identity matched for all 225, both MMFF implementations supported
212 and neither supported 13. Reference recovery was 161 / 163
(chematic / RDKit; denominator 225; paired molecule 95% interval −4.00 to
+2.22 pp). The engines returned 4,120 / 4,180 minimized candidates with
0 / 38 detected stereo violations. Six supported records had no chematic
candidates; three of those had no RDKit candidate count because RDKit
encountered an error. This remains development evidence. Compressed raw
files are `validation/results/cod-linux-{development,identity,support}-v8-2026-11-22.json.gz`;
uncompressed development SHA-256 is
`3f9503d09caf42ccbdfde3706f96cb343d37a03ab58f538181225b64d45c5522`.

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

### Second COD holdout: v4 fixed-count and timed results

The [v4 evaluation](https://github.com/kent-tokyo/chematic/actions/runs/44962477871)
completed both platforms' scoring and aggregation before its force-field jobs.
Of 2,288 scaffold-separated holdout records, 2,098 were MMFF-supported by both
engines, 181 by neither, and nine by chematic only. All 2,288 passed the identity
checks on both platforms. The following recovery denominators are the 2,098
common-supported records, not the entire holdout; none were dropped for a
missing conformer. These results are from the v4 wheels, before the subsequent
P, UFF, terminal-carbon, phosphorus-bond and directional-bond fixes.

| Platform and budget | chematic / RDKit recovered | Paired molecule 95% interval | Paired scaffold 95% interval |
| --- | ---: | ---: | ---: |
| Linux, 20 candidates | 1,582 / 1,633 | −3.38 to −1.48 pp | −3.53 to −1.40 pp |
| macOS, 20 candidates | 1,586 / 1,625 | −2.81 to −0.91 pp | −3.00 to −0.81 pp |
| Linux, 0.2 s per engine | 1,195 / 1,148 | +1.19 to +3.29 pp | +1.10 to +3.36 pp |
| macOS, 0.2 s per engine | 1,162 / 1,067 | +3.29 to +5.72 pp | +3.27 to +5.76 pp |

Thus the preregistered −2 pp lower-bound criterion failed on both platforms
for fixed count and passed on both for the timed condition. The fixed run
returned no chematic candidates for 63 inputs on each platform; 52 of those
had RDKit candidates, while RDKit returned an error for the other 11. These
11 errors have no recorded RDKit candidate count, rather than a count of zero.
Linux checked 40,532 / 41,551 chematic / RDKit minimized candidates and found
0 / 65 stereo violations; macOS checked 40,538 / 41,548 and found 0 / 70.
Only 18,829 / 19,535 of these Linux candidates and 18,781 / 19,563 of the
macOS candidates respectively reported convergence. Timed chematic candidates
had zero detected stereo violations on each OS; RDKit had seven on each. These
checks do not prove stereochemistry for cases the verifier cannot evaluate.
Common-coordinate force-field comparisons are still running. Raw fixed and
timed artifacts and summaries are preserved under `validation/results/` with
`cod-second-holdout-` prefixes; uncompressed fixed SHA-256 values are
`912be2665293e6ce76e7b349ff720cc3127158c6610cc7c60faf51d2b70e9f8c`
(Linux) and `19506c7988bfdbd901936f1c5cbd0c7400985725fdff0817da4b0d49dc7e51e5`
(macOS). The v8 candidate was frozen before these results were inspected, but
any comparison on this same holdout is a related follow-up, not a fresh
independent test. This was itself a post-hoc follow-up after first-COD exposure.

The v4 common-coordinate force-field jobs subsequently completed on both
platforms. Of 2,098 records, RDKit reference embedding failed on 22 Linux and
21 macOS records; all remaining records were comparable for MMFF94, MMFF94s,
and UFF. For MMFF94 the maximum absolute energy differences were 271.60 and
275.45 kcal/mol and maximum gradient-component differences were 513.78 and
558.40 kcal/mol/Å on Linux and macOS respectively. For UFF the corresponding
maximum energy differences were 35.77 and 36.07 kcal/mol and gradient
component differences 633.65 and 820.85 kcal/mol/Å. The p95 of energy and
gradient differences was zero on Linux and below 2e−13 on macOS; those
percentiles do not remove the significant outliers. Maximum sampled numerical
finite-difference discrepancies for MMFF94 were 141,918 and 302,089 kcal/mol/Å;
their individual causes have not yet been established. These results concern
the old v4 wheels, not the later fixes or the sulfoxide change. Raw results and
summaries have `cod-second-holdout-{linux,macos}-fields` prefixes under
`validation/results/`; uncompressed raw SHA-256 values are
`325606dfb7a7e9b43de49507ae58e99cced799f52d54bc1c7aace54d6f198939`
and `163543e8cba0ef5ae41d9d23f129e04369cce79143426110d1ca3d88de690d9e`.

The ring-order investigation is pinned to the actual RDKit 2026.09.2 build:
[its RingFamilies CMake configuration](https://github.com/rdkit/rdkit/blob/Release_2026_09_2/External/RingFamilies/CMakeLists.txt)
selects RingDecomposerLib `v1.1.3_rdkit` (archive MD5
`e9a0bcdda8b921a35e812b9888a9a874`) and replaces `qsort` with a stable
sort in `RDLcycleFams.c` and `RDLdimacs.c`. A compatible port must account for
that stable ordering as well as the biconnected-component, ring-family and
cycle iteration order. This source identification does not itself fix or
validate the MMFF discrepancies.
The pinned RingDecomposerLib implementation first performs biconnected
component decomposition, enumerates cycle families using an ordering involving
vertex degree, and stably sorts candidate families by cycle weight; its cycle
iterator advances through component, URF, RCF and path order. The existing
`relevant_cycles.rs` instead collects candidates globally using atom-index
ordering and expands relevant candidate families. Matching the set therefore
does not guarantee matching the order. Compare
[`RDLcycleFams.c`](https://github.com/rareylab/RingDecomposerLib/blob/v1.1.3_rdkit/src/RingDecomposerLib/RDLcycleFams.c)
and [`RDLhandler.c`](https://github.com/rareylab/RingDecomposerLib/blob/v1.1.3_rdkit/src/RingDecomposerLib/RDLhandler.c).

### Exposed second-COD follow-up probes

After inspecting v4 holdout results, 29 previously zero-candidate records with
specified neutral three-coordinate S(=O) stereo and RDKit candidates were
selected for a development-only sulfoxide probe. The updated source passed the
[Linux stereo regression](https://github.com/kent-tokyo/chematic/actions/runs/44973738608).
In the [candidate probe](https://github.com/kent-tokyo/chematic/actions/runs/44974156820)
the source merge commit was `6379b78989c7b56d37123635d86fd1f59cc3a92d`;
all 29 returned 20 candidates, compared with zero in v4. RDKit detected no
stereo violations in the 580 returned candidates, of which 176 were unfinished.
These cases informed the change and are not an independent holdout. The raw
result is `validation/results/cod-sulfoxide-candidates-development-2026-11-22.json.gz`,
with uncompressed SHA-256
`30316e0d0f4cf36d0c319854607b2405b529ea500923d5ccd882218efab39e76`.

A separate [frozen-v8 field probe](https://github.com/kent-tokyo/chematic/actions/runs/44974275953)
selected 57 exposed records with v4 MMFF94 or UFF energy error above 1e−4,
gradient-component error above 1, or sampled finite-difference discrepancy
above 1, on either OS. One RDKit reference embedding failed; on the other 56,
the count with MMFF94 energy error above 1e−4 fell from 18 in Linux v4 to 11
in v8, and UFF gradient-component error above 1 fell from 39 to 4. These are
Linux re-embeddings at the same seed, not a promise of bit-identical coordinates
across runners. The v8 maximum MMFF94 energy error remained 254.53 kcal/mol,
and the maximum MMFF94 and UFF sampled finite-difference discrepancies were
141,918 and 113,682 kcal/mol/Å. The largest MMFF energy outlier has multiple
negative formal charges; the finite-difference outlier is a different molecule.
The raw result is `validation/results/cod-second-fields-v8-development-2026-11-22.json.gz`,
with uncompressed SHA-256
`0f077e7eebaae5824b05cda0d8c0e73dd2f92832530e5ff90588060f08f0f4a5`.
Neither follow-up is independent evidence for v8 or a later candidate.

For the exposed finite-difference outlier YEBDEXSGHHFXMC-FNCQTZNRSA-N, a
[frozen-v8 diagnostic](https://github.com/kent-tokyo/chematic/actions/runs/44975074988)
recreated the geometry independently on each OS. Its angle (atoms 10–9–0) was
179.99832° on Linux and 179.99879° on macOS. For the sampled component, the
MMFF94 compatibility gradient was 39,499.68 and −11,138.73 kcal/mol/Å,
matching RDKit; the central differences at 1e−9 Å were 46,440.73 and
−13,637.32. The optional consistent-plane gradients were 46,440.79 and
−13,637.15. UFF analytic and 1e−9 Å difference were similarly close on these
coordinates, while 1e−5 Å gave dramatically different values. These are two
near-singular geometries, and small-step subtraction has its own rounding error;
they do not establish general gradient correctness. The checked optimizer
still uses its RDKit-compatible gradient, not the optional consistent-plane
one. Both raw records are in
`validation/results/near-linear-second-cod-{linux,macos}-v8-development-2026-11-22.json`.

For two exposed cage structures (XTULIPKREBMLNM-XHGUTPRISA-N and
OXBMFCGRQJNAOY-TWXUZQDTSA-N), v8 MMFF94 types and charges matched RDKit but
energy errors were 42.17 and 41.62 kcal/mol. Term-level diagnosis traced the
difference to angle and stretch-bend terms: one SSSR omitted a relevant
four-membered ring, affecting four angles on each molecule. An intermediate candidate used relevant cycles for angle parameterization in
both the prepared field and other MMFF paths, falling back to SSSR if
enumeration exceeded 10,000 cycles.
The [Linux regression](https://github.com/kent-tokyo/chematic/actions/runs/44976034714)
passed the RDKit-derived fixed-coordinate angle and stretch-bend energy and
gradient test, the existing `chematic-ff` library tests, and a source-built
wheel comparison on ETKDG coordinates. In the latter comparison both inputs
had zero energy and gradient-component difference on that runner; the tested
merge commit was `86fe06426ec53651c1f584e64f0e9d5d3bbbe563` and results are
`validation/results/mmff-cage-angle-development-2026-11-22.json`. This is a
targeted development result; the fallback, performance cost, broader corpus,
macOS and independent holdout have not been validated for this change.

The [v8/v9 constructor microbenchmark](https://github.com/kent-tokyo/chematic/actions/runs/44977893135)
ran on one Linux runner using verified source-built wheels and synthetic fixed
coordinates for 212 previously exposed COD development records. Four separate
processes ran in v8–v9–v9–v8 order; each performed a warm-up then seven timed
field constructions per molecule. No construction errors occurred. Pooling
repeats within each version, the median per-molecule v9/v8 ratio was 1.060,
and its p95 was 1.114. The sums of per-molecule median times were 86.25 ms
(v8) and 91.51 ms (v9). The run is short, includes runtime noise, and measures
warm repeated construction, not cold startup or end-to-end conformer speed.
Results are `validation/results/mmff-construction-v8-v9-development-2026-11-22.json.gz`.
The v9 Linux wheel SHA-256 was
`bffb1cd496eef80f842ad26e9c036ba7cd46075204aa99374d1909a3a1642002`;
it was built from the already tested merge commit
`86fe06426ec53651c1f584e64f0e9d5d3bbbe563` (crates tree
`fe6e1e1a7631d81e847e4be0c019d6f571a9834b`). The macOS wheel SHA-256 is
`605e1e0c0f256c9216b6f954a8a1cb6f8edc32b8375f71dd23308867182e828f`;
both platform wheels were downloaded and their hashes verified.
The 225-record v9 development conformer evaluation remains running.

A later intermediate change cached the bounded MMFF parameter ring set on the
input molecule and invalidated the cache on edits. The [Linux cache validation](https://github.com/kent-tokyo/chematic/actions/runs/44978178977)
passed the core and force-field library tests, a cache reuse/invalidation test,
the cage angle regression and the two exposed shared-coordinate cage checks
(zero observed energy and gradient-component differences). Its source merge
commit was `d038846ec894036371ae34cf961a9e76886642c9`; the Linux wheel SHA-256
was `0ecc715e6bf5915f24d3d7237eff5fc378c22d4c2fe9a13f1f340405ff84b69f`.
On 212 exposed development molecules, the same warm-construction A-B-B-A
microbenchmark reported a median cached/v9 time ratio of 0.940 (p95 0.969),
with sums of per-molecule medians 65.65 ms and 70.16 ms respectively. This
is a separate runner from the v8/v9 measurement and does not establish the
cached candidate's speed relative to v8 or end-to-end conformer performance.
Raw timings and cage diagnostics are in
`validation/results/mmff-ring-cache-development-2026-11-22.json.gz`.

The [v9 Linux development conformer evaluation](https://github.com/kent-tokyo/chematic/actions/runs/44977687329)
completed on the frozen pre-cache, pre-aromaticity-fix wheel. All 225 atom
identities matched; both implementations supported MMFF for 212 records and
neither for 13. Reference recovery was 161 / 163 (chematic / RDKit, denominator
225; paired molecule 95% interval −4.00 to +2.22 pp). The 4,160 / 4,180
returned candidates had 0 / 38 detected stereo violations. Four supported
inputs had zero chematic candidates, compared with six in the v8 development
run. These exposed development results are not a noninferiority test. Raw data
are `validation/results/cod-linux-{development,identity,support}-v9-2026-11-22.json.gz`;
the uncompressed development result SHA-256 is
`9612ca09c4f3608e1cc06861e72c17ccc912d4dfca5ef41144a9e46a48937653`.

Two exposed macrocycles (SAWYAYVQRLIVEW-UHFFFAOYSA-P and
FFDBWJMOMIIWCE-JNLZMNOBSA-N) had UFF discrepancies because the shared
RDKit-parity view lost aromatic furan bonds when a larger ring was processed
first. Relevant cycles are now used in both ranked perception and the SMILES
sanitizer's aromaticity pass, with the existing path retained for the unranked
shortcut and on the 10,000-cycle overflow. This does not reproduce the exact
RingDecomposerLib family order in all graphs. The
[Linux validation](https://github.com/kent-tokyo/chematic/actions/runs/44980080896)
passed 274 perception tests, 447 SMILES tests, the force-field library tests,
and fixed-coordinate RDKit 2026.09.2 aromatic-bond, UFF energy and gradient
regressions. On a source-built wheel (merge commit
`6fae31e8341a11913a20e46f740105044a6a5b59`, wheel SHA-256
`e5b9211c1d72d37bf1bfc583d1771219c56afa0572bca66a8d05195067868c08`), both
macrocycles matched RDKit energy and gradient exactly on that runner's shared
ETKDG coordinates. Among the 56 successfully embedded cases of the 57-case
exposed second-COD outlier follow-up, the UFF energy-error count above 1e−4
kcal/mol and gradient-component error count above 1 kcal/mol/Å each fell from
4 on v8 to 2. MMFF94 energy-error count fell from 11 to 9 and gradient-error
count from 10 to 8, including the later cage-angle fix; its maximum energy
error remained 254.53 kcal/mol. One RDKit embedding failed in both runs.
This is development-only evidence with re-embedded coordinates, not an
independent holdout or a general correctness claim. Raw result:
`validation/results/cod-second-fields-aromatic-development-2026-11-22.json.gz`,
uncompressed SHA-256
`7f4d4b3d62519ada8c5950f6d960ae006f9428e4d9cfab8264768b109fe1fde6`.

RDKit 2026.09.2 [`getMMFFAngleType`](https://github.com/rdkit/rdkit/blob/Release_2026_09_2/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp)
uses local bond adjacency for three- and four-membered angle classification.
A four-cycle may be absent even from the union of minimum cycle bases, for
example when its edge set is a sum of triangles. The current source replaces
the intermediate relevant-cycle and cache implementations with this local
topology check; the public `rings` parameter is retained but no longer
consulted. A direct local Rust build passed the two earlier cage tests and an
oracle-checked chorded-cycle classification probe. Full source-built CI,
shared-coordinate outlier comparison and construction timing are in progress.
The 10,000-cycle bound still applies to aromaticity and other separate uses.

### First COD macOS checkpoint resume

The [checkpointed macOS resume](https://github.com/kent-tokyo/chematic/actions/runs/44964491252)
completed the remaining 285 fixed-count records and merged them with the
original 1,808, after checking the frozen wheel, original checkpoint, split,
input identities and scored names. This is the earlier pre-stereo-fix wheel,
not v9 or the latest source. Among 2,093 common-supported records, fixed
20-candidate recovery was 1,662 / 1,662 (chematic / RDKit); the paired
molecule interval was −0.669 to +0.669 pp and scaffold-cluster interval was
−0.670 to +0.662 pp. Those numerical bounds meet −2 pp, but the fixed
candidates had 66 / 70 detected stereo violations among 41,441 / 41,437
returned candidates; therefore the stereo gate fails. The split had also been
partly exposed before scoring, so it is not independent noninferiority
evidence. The original and resumed records ran on separate macOS arm64 hosts;
this is not an uninterrupted timing measurement. Timed and force-field checks
subsequently completed on a third host, as recorded below. The merged and resumed portions are retained as
`validation/results/cod-first-holdout-macos-{fixed-resumed,resumed-285}-2026-11-22.json.gz`,
with uncompressed SHA-256
`621de9dcf4aa71ed90418e4d2cbe0072d65cdfa883e39f0b6ddea6f43f45123e` and
`f0b533f09ae2c65378641bcd10f437f553541e1413d691a224eec6f283751437`;
the fixed and stereo summaries are alongside them.

The [Linux local-angle regression](https://github.com/kent-tokyo/chematic/actions/runs/44981578055)
subsequently passed the core and force-field library tests, the chorded-cycle
classification and the prior cage and macrocycle regressions. On shared ETKDG
coordinates the two exposed cage structures again had zero observed MMFF94
energy and gradient-component difference. In the 56 successfully embedded
cases of the same 57-record development follow-up, the counts above the
previous thresholds remained MMFF 9 energy / 8 gradient and UFF 2 / 2; one
RDKit embedding failed. An A-B-B-A warm-construction microbenchmark on 212
exposed development inputs measured median local-angle/v9 ratio 0.906 and p95
0.940; sums of per-molecule medians were 63.33 ms / 70.52 ms. This short
single-runner microbenchmark does not establish end-to-end speed. The tested
merge commit was `55c5db999467dff3e30619b7a9e0467f55cc87bb`, with Linux wheel
SHA-256 `55c448bc30cea6ebffdd930125730f61c8da094d24fd60033b4bad7437139c52`.
Raw timing and field outputs are retained in
`validation/results/mmff-local-angle-development-2026-11-22.json.gz`.

An exposed MMFF outlier with formal charge −23 contains two sp-hybridized
cumulene carbons in a six-membered ring. RDKit 2026.09.2 requires each ring
carbon or nitrogen to be sp2 before accepting an MMFF aromatic ring; the old
port checked degree instead and incorrectly accepted that ring. The current
source uses sanitized hybridization, retained through Kekulization. The
[Linux regression](https://github.com/kent-tokyo/chematic/actions/runs/44981910098)
passed the RDKit-derived atom-type, charge, energy and full-gradient test,
the local-angle and UFF macrocycle tests, and the force-field library suite.
Its source-built wheel (merge commit `1d868c13991773c4fdf9f834389e5cf71eba8ea2`,
SHA-256 `21f482a615d2c4c2e7f5f224932a611d05cfc10bceb14c42cd3c51b03f099702`)
matched RDKit energy and all gradient components exactly on that runner's
re-embedded coordinates for the outlier. Across 56 successfully embedded
records in the same 57-case exposed development follow-up, MMFF94 energy
errors above 1e−4 fell from 9 to 8 and gradient-component errors above 1
from 8 to 7 compared with the preceding aromaticity candidate. The maximum
energy error fell from 254.53 to 41.64 kcal/mol and maximum gradient-component
error from 513.78 to 141.19 kcal/mol/Å; one RDKit embedding failed in both
runs. UFF counts remained 2 / 2. These are selected development inputs, not a
general or independent guarantee. Raw data:
`validation/results/cod-second-fields-hybrid-development-2026-11-22.json.gz`,
uncompressed SHA-256
`463ce55606ff29af018c86452d7a80fabe06b2bc5bb1d0f2ae4d12fc9f29c311`.

The next MMFF aromatic-typing candidate enumerates relevant cycles, with the
previous SSSR implementation retained as a fallback above the 10,000-cycle
limit. This is not a reproduction of every RDKit ring ordering. The
[Linux regression](https://github.com/kent-tokyo/chematic/actions/runs/44982865443)
passed the force-field library and previous targeted tests, and all eight
previously exposed ring-typing outliers matched RDKit MMFF94 energies and full
gradients on the runner's shared coordinates. Its tested merge commit was
`12ca0abae3040563d503c7e13a059e8c91a808cc`; the Linux wheel SHA-256 was
`2b5cd7a9c7766b6e9ccc97e52c4522c7aa6f85c5476ac3f4a283f048f5af1a6e`.
Among 51 first-split selected cases, one MMFF94 energy error exceeded 1e-4
and one gradient-component error exceeded 1. Among 56 successfully embedded
of 57 second-split selected cases, the counts were two and one respectively;
one embedded in RDKit for the remaining case. These exposed, selected sets
are development checks, not independent acceptance tests. Some sampled
finite-difference errors also remain large and require investigation. The
complete runner outputs are retained in
`validation/results/cod-first-fields-relevant-development-2026-11-22.json.gz`
and `validation/results/cod-second-fields-relevant-development-2026-11-22.json.gz`.

RDKit's UFF builder assigns SP2D, SP3D and SP3D2 hybridization from
non-tetrahedral stereo tags. The UFF port now accounts for these tags when
selecting labels and angular terms. The [Linux UFF regression](https://github.com/kent-tokyo/chematic/actions/runs/44983432313)
passed the library suite, 11 fixed-coordinate cases and the two previously
exposed phosphorus UFF comparisons. The tested merge commit was
`17958cde22c1ce6ec1a44d006b0f1a857401a9e4`; Linux wheel SHA-256 was
`1558a13111912671ebbc0b7a9a5ca408b3c9373059277eba5a801bf05b7a0c8c`.
The job checked both selected field sets, but its artifact glob omitted the
JSON outputs; the passing targeted assertions remain in the job log. A later
workflow corrects the artifact path. Local exploratory checks also found the
MMFF torsion builder incorrectly included SP2D central atoms, and that the
empirical torsion rule used MMFFProp aromatic flags instead of the separate
MMFFArom type set. The [MMFF non-tetrahedral Linux regression](https://github.com/kent-tokyo/chematic/actions/runs/44983825231) passed its fixed-coordinate cases, the force-field library suite and selected shared-coordinate outlier assertions. The subsequent aromatic-torsion comparison is described below. These exposed checks do not establish parity.

The remaining selected first-split MMFF energy/gradient outlier has six bonds
whose MMFF aromatic status differs from RDKit. The relevant-cycle sets match,
but their iteration order differs, and the fused-ring aromaticity pass is
order dependent. A local trial sorting fused-ring centers after peripheral
rings failed the existing `DEZSLFWYSCQUHD-QUMXLFQESA-N` atom-type regression;
that heuristic was not adopted. Reproducing RingDecomposerLib's ordering, not
just its cycle set, remains an open compatibility problem.

The finite-difference outliers were also checked locally against RDKit
2026.09.2 itself, using ETKDGv3 seed 42 and six sampled Cartesian components.
For `YEBDEXSGHHFXMC-FNCQTZNRSA-N`, RDKit UFF's largest sampled difference
between its analytic gradient and central energy difference at step 1e-5 was
253,964.55, with analytic component magnitude 11,053.00; RDKit MMFF94's
corresponding maximum was 302,088.59. These macOS-embedded coordinates differ
from the Linux CI coordinates. This establishes that the discrepancy can also
occur in RDKit on this exposed structure; it neither explains its cause nor
establishes stability or derivative accuracy for either implementation.

The [aromatic-torsion Linux regression](https://github.com/kent-tokyo/chematic/actions/runs/44984248482)
passed the force-field library and all targeted tests and shared-coordinate
assertions. The tested merge commit was
`22578e3829bfee7aaab03b4c1bd829e580c3a8fe`, with Linux wheel SHA-256
`9b43183d8fd7889520a1b166fcfb1cbb7e59c44b893910f6855b9ffeafaaa8ac`.
Among the selected 57 second-split development records, one RDKit embedding
failed; all 56 remaining records had reported zero energy and maximum
component gradient difference for MMFF94, MMFF94s and UFF. Among 51 selected
first-split records, one MMFF94/MMFF94s outlier remains, with maximum energy
difference 4.958 kcal/mol and gradient-component difference 49.682 kcal/mol/Å;
UFF differences were reported as zero. These selected and exposed inputs are
not independent acceptance evidence. Full outputs:
`validation/results/cod-first-fields-aromatic-torsion-development-2026-11-22.json.gz`
and `validation/results/cod-second-fields-aromatic-torsion-development-2026-11-22.json.gz`.

Linux and macOS candidate wheels were separately built from the tested merge
commit `22578e3829bfee7aaab03b4c1bd829e580c3a8fe` in the
[v11 wheel run](https://github.com/kent-tokyo/chematic/actions/runs/44984753739).
The Linux wheel SHA-256 is
`0052c7ad06ce6fe92ddf919fc0e2a15869b7d683c530f794c62a052edf0d1603`;
the macOS arm64 wheel SHA-256 is
`fa33b33382589fd46124131525e5efabc06150bca33eba7153a4f551c9d985d6`.
These binaries are not the Linux wheel used in the earlier one-off regression
job, despite sharing its tested source. New development comparisons use the
v11 wheel hashes; neither wheel has passed a fresh independent holdout.

The first COD macOS [timed and force-field job](https://github.com/kent-tokyo/chematic/actions/runs/44964491252)
completed on the earlier frozen wheel, separate from both fixed-count hosts.
At 0.2 seconds, recovery was 1,280 / 1,134 out of 2,093 for chematic / RDKit;
the paired scaffold-cluster 95% interval was +5.80 to +8.17 percentage points.
The returned ensembles contained 16 / 8 detected stereo violations among
13,592 / 7,438 candidates, so the chematic stereo gate fails despite that
recovery interval. There were no row errors. These results are from an earlier
partly exposed split and predate the later fixes; they do not establish a
current or independent parity result.

On shared RDKit coordinates for the same 2,093 cases, RDKit embedding failed
for 23. For the remaining 2,070, maximum MMFF94 energy and gradient-component
differences were 57.67 kcal/mol and 121.10 kcal/mol/Å; UFF's maximum gradient
difference was 74.36 kcal/mol/Å. The per-record timed and force-field results
are `validation/results/cod-first-holdout-macos-time-2026-11-22.json.gz` and
`validation/results/cod-first-holdout-macos-fields-2026-11-22.json.gz`, with
uncompressed SHA-256
`3379bdc09929e04002eb1a01b03317f1ef92a21e38728249b7e9304cf738045c` and
`883f8cbfa5dfb628556aefbf32dcd35c3326ad391391367be7070e47ef1cb1b7`.
Machine-readable molecule, scaffold and field summaries are alongside them.

The v11 [225-record COD development comparison](https://github.com/kent-tokyo/chematic/actions/runs/44985345006)
completed both platforms and both API profiles under the earlier CIP-based
stereo check. All 225 input identities matched; 212 had MMFF support in both
engines and 13 in neither. The following recoveries use all 225 in the
denominator and require the existing geometry and stereo screens:

| Host and profile | Chematic / RDKit recovered | Chematic / RDKit detected stereo violations | Chematic / RDKit candidates checked |
| --- | --- | --- | --- |
| Linux compatible | 162 / 163 | 42 / 38 | 4,180 / 4,180 |
| Linux corrected | 161 / 163 | 0 / 38 | 4,160 / 4,180 |
| macOS compatible | 161 / 163 | 38 / 39 | 4,180 / 4,180 |
| macOS corrected | 161 / 163 | 0 / 39 | 4,160 / 4,180 |

In each compatible run, three supported inputs yielded no candidates for
both engines; in each corrected run, chematic yielded none for one additional
input while RDKit produced 20. That input has a declared square-planar sulfur
center; the corrected path conservatively rejects stereo it cannot verify. A
local five-conformer diagnostic on an earlier macOS wheel and RDKit 2026.09.2
found that RDKit's 3D reassignment did not recover that square-planar tag for
any of those generated or RDKit-optimized geometries. This small diagnostic
is not a comparison of the final v11 wheels.

These are exposed development inputs, not an independent acceptance test.
The old checker did not inspect specified centers without CIP labels or
non-tetrahedral stereo, so the zero counts in the table are limited to what
it checked. All four raw outputs and summaries are in
`validation/results/cod-{linux,macos}-development-v11-{compatible,corrected}-2026-11-22.json.gz`
and corresponding `-summary-` files. Uncompressed raw SHA-256 values, in
table order, are
`1395f0534ee70e1fd782b19d3306805c808fcc0fe608b5ae9ccda89ac8c571b7`,
`e81fd9c1eba4c8fa218308ad1e312939b8c3e7ec306684711f2341449d9448bb`,
`4f924a25251c3f6d550474a55a2b9b3498e4ee035aae943e542dc6e0db0433fe`, and
`58671c5408896f26b3c39453f29d50e185d47fcf8e3b9ff81dbfc093b4eda32f`.
The v11 development support and atom-identity audit outputs were byte-identical
on Linux and macOS and are retained as
`validation/results/cod-support-v11-development-2026-11-22.json.gz` and
`validation/results/cod-identity-v11-development-2026-11-22.json.gz`, with
uncompressed SHA-256
`53626b44dbcda1832f4744a42100fd92d59579f2b50d212d7cf43df1ea8e20c1` and
`8860ef9ceaada101a93eae60fe95493193bfccf4cae6bfb7f01b3ea79b3cc202`.

A separate local RDKit 2026.09.2 structural audit of the 225 COD development
inputs found 259 specified tetrahedral centers, 58 specified double bonds and
one square-planar center. Of the tetrahedral centers, 28 in nine records had
no assigned CIP label. This audited input declarations without generating or
scoring conformers. The per-record output is
`validation/results/cod-development-declared-stereo-audit-local-2026-11-22.json.gz`,
uncompressed SHA-256
`2a9f2c6250f7cb9721fac77b86d2396405555f1fef340a5aa5e9bf82e77dfafc`.

The Linux corrected-profile arm of the [strict-stereo development repeat](https://github.com/kent-tokyo/chematic/actions/runs/44990053917)
completed on the frozen v11 wheel. Recovery remained 161 / 163 of 225 for
chematic / RDKit, and all 225 per-record RMSD arrays for both engines matched
the corresponding earlier v11 Linux run. Of 4,160 chematic candidates, zero
had contradicted or unevaluable declared stereo. Of 4,180 RDKit candidates,
19 had contradictions, 56 had unevaluable declarations, and 58 had either
(the categories overlap). Seventeen / sixteen records had no returned
candidate, including the 13 unsupported by both engines. These are exposed
development data, not a fresh noninferiority result. The raw output is
`validation/results/cod-linux-development-v11-corrected-strict-2026-11-22.json.gz`,
uncompressed SHA-256
`2f59294f7b22a6f257dc284b8c14fc2fede18c289c72b425d510f729c1c11fbd`;
recovery and strict-stereo summaries are alongside it. The Linux structural
audit of all 225 development inputs is
`validation/results/cod-development-declared-stereo-audit-linux-2026-11-22.json.gz`,
uncompressed SHA-256
`4f3390c46d4e89023f5d4f1959a38e9265c844e7ea451b48276446052ed6ce7b`.
The macOS corrected-profile arm also recovered 161 / 163; all 225 per-record
RMSD arrays and convergence-status arrays matched the earlier macOS v11 run.
Among 4,160
chematic candidates there were zero contradicted or unevaluable declarations;
among 4,180 RDKit candidates, 19 were contradicted, 57 unevaluable and 59
had either. Seventeen / sixteen records returned no candidate, including 13
unsupported by both. The raw macOS output is
`validation/results/cod-macos-development-v11-corrected-strict-2026-11-22.json.gz`,
uncompressed SHA-256
`a52c8071108b1221b21f0117fb5cb2e351094db0eab31ad575ac7de245b4aeb0`.
Its structural audit is
`validation/results/cod-development-declared-stereo-audit-macos-2026-11-22.json.gz`,
uncompressed SHA-256
`660a44b16ee8371af2f74d8e41e2e4f3d755baa64cba474308e83fac0f7ccacf`.
The compatible-profile arms also completed. Linux recovery was 162 / 163;
macOS was 161 / 163. In 4,180 returned chematic candidates, Linux had 21
contradicted, 59 unevaluable and 62 with either; macOS had 18, 55 and 58.
For 4,180 RDKit candidates Linux had 19, 56 and 58 respectively, and macOS
19, 57 and 59. The categories overlap. All 225 per-record RMSD arrays and
status arrays matched their earlier platform-specific compatible v11 runs.
Raw results are
`validation/results/cod-linux-development-v11-compatible-strict-2026-11-22.json.gz`
and `validation/results/cod-macos-development-v11-compatible-strict-2026-11-22.json.gz`,
with uncompressed SHA-256
`b9c4712e5972c00674a042db1ffbea4fc70b92e9d45523fea60fd2dda07ba6f1` and
`30d017e53528e33729f7a79f96a241c01435ff5cf40852a002e05071ca57dfa9`.
These are exposed development comparisons, not a fresh acceptance test.

The Linux arm of the [fallback strict-stereo COD development repeat](https://github.com/kent-tokyo/chematic/actions/runs/44994764474)
recovered 161 / 163 references among 225 records, exactly as the v11
corrected-profile run. All 225 chematic candidate counts and RMSD arrays
matched the earlier run, as did all RDKit RMSD arrays. Chematic had zero
contradicted or unevaluable declarations among 4,160 candidates. The same
three RDKit `Bad Conformer Id` errors occurred in both runs; 13 additional
records were unsupported by both. This is exposed development evidence and
shows no recovery gain from fallback in this sample. The raw result is
`validation/results/cod-linux-development-fallback-strict-2026-11-22.json.gz`,
uncompressed SHA-256
`71c65b728c0cf6a065f2d1c4dea8555eca87cfa7dffb67c8ab024f718fba663b`,
with recovery and strict-stereo summaries alongside it. The macOS arm also
recovered 161 / 163, with zero contradicted or unevaluable declarations in
4,160 chematic candidates. All 225 chematic candidate counts and chematic and
RDKit RMSD arrays matched the earlier macOS corrected-profile run; the same
three RDKit conformer errors occurred. The raw macOS result is
`validation/results/cod-macos-development-fallback-strict-2026-11-22.json.gz`,
uncompressed SHA-256
`568983c0bc413a23c00a56de188c774177b01f449edc76cbee75edf1ca4a9b82`.
These comparisons used the candidate before the subsequent phosphorus/sulfur
verifier changes.

For development diagnosis, 52 already-exposed second-COD holdout records were
selected solely because the older Linux v4 run returned zero chematic and at
least one RDKit candidate. The selection script verifies the prior result and
original split hashes, then labels these records as development; it does not
create a new independent holdout. In the [v11 corrected follow-up](https://github.com/kent-tokyo/chematic/actions/runs/44996384331),
all 52 records were processed, 45 had at least one chematic candidate, and
recovery was 36 / 40 (chematic / RDKit). Among 900 / 1,016 candidates, the
strict verifier found chematic 0 contradictions and 0 unevaluable, and RDKit
19 contradictions and 36 unevaluable. This selected sample cannot support a
general performance estimate. The raw result is
`validation/results/cod-exposed-empty-linux-corrected-strict-2026-11-22.json.gz`,
uncompressed SHA-256
`83c31eae906112bc82840e713dbf3ad90bb9255e18eb4ec28dcbde8e5372e146`;
the exact selection and summaries are alongside it. The fallback arm also
recovered 36 / 40, with seven zero-candidate records and no chematic stereo
contradictions or unknowns among 900 returned candidates. All 52 candidate
counts and chematic and RDKit RMSD arrays matched the corrected arm. The
fallback raw result is
`validation/results/cod-exposed-empty-linux-fallback-strict-2026-11-22.json.gz`,
uncompressed SHA-256
`94834e1cd261610c2f2fc98d95b7f75a1999f292c3d129e94fc115b3707d303c`.
The seven remaining records include non-tetrahedral and sulfur/phosphorus
stereo declarations; the reason for each zero result requires diagnosis.

A [same-coordinate stereo diagnostic](https://github.com/kent-tokyo/chematic/actions/runs/44996942332)
regenerated up to five RDKit 2026.09.2 optimized geometries for each of the
seven exposed zero-candidate records, checking atom and bond identity first.
All seven had at least one unevaluable declaration in the native verifier.
Four contained externally evaluable, satisfied examples: a four-coordinate
S(VI) center (5/5), a three-coordinate S=N center (3/5), a four-coordinate
phosphorus center with two double bonds (5/5), and a three-coordinate
sulfonium center (4/5). The other three contained square-planar,
trigonal-bipyramidal or degenerate double-bond declarations; these remain
unevaluable in this diagnostic. The RDKit geometry input and native reports
are `validation/results/cod-exposed-stereo-rdkit-input-2026-11-22.json.gz`
and `validation/results/cod-exposed-stereo-diagnostic-2026-11-22.json.gz`,
with uncompressed SHA-256
`d6029b3350c3ea6f33f14b4713bf6b6be1830a10299ac730c4bb5a31045a167d` and
`8227f938e47562df9a889438f8609fb7aacb0521eabb4fba4ca04b12c87cfe3a`.
These reports identify a verifier coverage gap on RDKit geometries; they do
not establish why each chematic embedding attempt failed.

The verifier was then extended only for four-neighbor phosphorus with one
double bond (existing behavior), four-neighbor P(-)(=N)(=S) or neutral S(=N)(=O), three-neighbor neutral S with one S=O or S=N bond, and
three-neighbor S(+) with single bonds. Unsupported non-tetrahedral geometry
remains unevaluable. The [cross-platform development tests](https://github.com/kent-tokyo/chematic/actions/runs/44997866129)
checked 17 externally satisfied RDKit geometries from the four exposed cases
and 20 additional RDKit-generated geometries from four smaller controls;
RDKit independently confirmed that each original matched its declaration and
that reflection contradicted the target center. Native checks passed for all
37 originals and mirrors on both hosts, along with existing stereo and
corrected-ensemble tests. A neutral four-coordinate sulfur radical outside the
specified topology remains unevaluable. These fixtures demonstrate bounded
verifier behavior, not general sulfur or phosphorus stereochemistry parity;
wheels and a fresh exposed-case repeat are needed to assess generated candidates.

A [follow-up with the first phosphorus/sulfur verifier candidate](https://github.com/kent-tokyo/chematic/actions/runs/44998708441)
used the same 52 exposed records, the fallback profile and a Linux wheel built
from merge `5020f56de0f7bd03e4fcf87f46bab09b4d5956f6` (SHA-256
`e21d849fbac5795bccf6d99d35242a6e698639fdb2a2cf3bac78863a5d023b63`).
The number of zero-candidate records fell from seven to three, with 980
chematic and 1,016 RDKit candidates. Chematic / RDKit recovery rose from
36 / 40 to 39 / 40, while no chematic candidate had contradicted or
unevaluable declared stereo. Exactly four chematic RMSD arrays changed: the
four sulfur/phosphorus records previously diagnosed as unevaluable; three of
those four were within 1 Å of their references. All RDKit RMSD arrays matched
the previous run. The three remaining zero-candidate cases declare
square-planar, trigonal-bipyramidal or degenerate double-bond stereo. The raw
result is `validation/results/cod-exposed-ps-linux-fallback-strict-2026-11-22.json.gz`,
uncompressed SHA-256
`62db0ce8ca010e418182005c7640cdf9bb7da97335ee70cd201fdfe3bf9ad4ff`;
summary files are alongside it. This selected, exposed sample is not a
noninferiority result. The two-double-bond support has since been narrowed to
the exact tested S(=N)(=O) and P(-)(=N)(=S) topologies and requires a further
repeat on the final candidate.
