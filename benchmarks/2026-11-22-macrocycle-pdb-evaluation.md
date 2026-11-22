# PDB macrocycle reference preparation

The [MacroCycle Benchmark](https://github.com/dptech-corp/MacroCycle-Benchmark)
reports 1,861 protein-bound macrocycle cases, with reference ligand SDF files
and an index. The repository's curation, index, and documentation are
attributed to the MacroCycle Benchmark contributors under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); the repository
states that the PDB-derived structure data follow the PDB CC0 policy. Cite
the original experimental structures when discussing individual cases. The
index and split below are derived from that compilation.

This preparation is separate from the earlier Prime benchmark. Protein-bound
reference geometries and the Prime crystal population need not have the same
distribution. A result on this dataset cannot retroactively make the earlier
Prime noninferiority test pass, and unknown historical exposure to these
public structures cannot be ruled out.

## Pinned index and graph-only split

The repository commit is `6afbb010770f8adaadb9db7f83f7a86538757a62`.
Its `benchmark_index.csv` Git blob is
`bb13d7f3c6cfff0f0e3fc06e85b1f38b5c4c8908`, and the exact CSV SHA-256 is
`cf408bbcbeed78fe8862f8f1ccf48a2714f18dcaad5b9ac1fa6c482f4d7155dc`.
The file contains 1,861 distinct case IDs and reported primary ring sizes of
8–40 atoms. Use the index's *unprotonated* isomeric SMILES and the matching
`ligand_unprotonated.sdf`; the pH-prepared derivative is a different input.
The reference structures and metadata are provided by the source; these
preparation checks do not independently establish their chemical correctness.

`scripts/validation/macrocycle_reference_split.py` uses RDKit 2026.09.1,
matching the version used to create prior scaffold hashes. It excludes every
scaffold present in the supplied prior Platinum Diverse, complete Platinum,
Prime and COD split manifests, including unused records in the COD sample.
This conservative exclusion is based only on index SMILES and previous
manifests. Remaining non-chiral canonical Murcko scaffold hashes are assigned
by the first 32 bits modulo five: zero is development and all other values are
holdout. All records sharing a scaffold remain together. The script records
each exclusion file hash, index row, case ID, scaffold hash and input-SMILES
hash; it does not use reference coordinates or conformer outcomes.

| Partition | Cases | Distinct isomeric SMILES | Scaffold groups |
| --- | ---: | ---: | ---: |
| Development | 363 | 324 | 230 |
| Holdout candidate | 1,243 | 1,067 | 802 |
| Excluded known scaffold | 255 | 179 | 88 |

The split was produced twice with identical uncompressed JSON SHA-256
`cd50ed3f5460e15b7e8a8dc0e33075884490d0e3c090d7a1ee7d39a81bd361c2`.
Its compressed copy is
`validation/results/macrocycle-pdb-scaffold-split-2026-11-22.json.gz`
(SHA-256 `e7ed8dff46d366cce15f90cc0504c1044b6e5544a26acd4fb9de2f6db2713795`).
The partitions have no shared scaffold hashes. Eight prior exclusions could
not be sanitized when the earlier manifests were created. The graph-only
`scripts/validation/macrocycle_invalid_exclusion_audit.py` compares each
eligible index scaffold as an atom-number subgraph against those eight
unsanitized structures, ignoring bond orders. It found no possible overlap
among 1,032 eligible scaffold groups; the audit is
`validation/results/macrocycle-pdb-invalid-exclusion-audit-2026-11-22.json`
(SHA-256 `1293b264be555e2802edb0a1067ec5e35fcfb5cc979935fdc2cf51ca9902c855`).
The same conservative exclusion audit was repeated after preparing the SDF
graphs and again found no possible overlaps among 1,032 eligible scaffolds.
For an input-preparation pilot only, the first development case,
`8ebk__Ligand_phospholipids_THR_B_4`, was fetched from the pinned commit.
Its source SDF Git blob is `a5d104cb80fac9bcd4220ef9b13c8274d26b0e55`
and SHA-256 is `aaac0c30b0df766e492011ba8bcfa1c1356d87982088afd29f13f7f7a4cc8ee4`.
RDKit 2026.09.1 parsed one 164-heavy-atom record with one coordinate set;
its canonical isomeric graph matched the index. A synthetic archive test of
`scripts/validation/prepare_macrocycle_reference.py` retained every source
SDF byte after the title line. This is one development input check, not a
conformer evaluation or validation of the other cases.

On this macOS arm64 CPython 3.13 environment, RDKit 2026.09.2 raised
`RuntimeError: bad any cast` for Murcko decomposition even on cyclohexane;
the isolated 2026.09.1 environment completed it. This observation is limited
to these local installations and is not a cross-platform RDKit claim.

## Required checks before performance scoring

The [preparation run](https://github.com/kent-tokyo/chematic/actions/runs/44955399937)
checked all 1,861 source SDF records against the index, including formal
charges, isotopes and declared stereochemistry. It found no input failures and
retained 363 development, 1,243 holdout-candidate and 255 excluded records.
The source archive SHA-256 was
`91e39f1a90f81982bf5781c54e7899952120f25a72a4d296a90bb709331fefcb`;
the prepared SDF SHA-256 was
`71bf05a54cd402146e8892acd88969ec52fd45d2a02c4194b7aa905fe900222d`,
and the prepared split SHA-256 was
`ffb256b0a9dbc779a5c6b3dd2ac96deec6174dc25408dc1d0691efcef8c4c4ff`.
The run artifact retains the inputs, source license, citation and audit.
Preparation does not establish MMFF support or conformer accuracy.

An exploratory [development pilot](https://github.com/kent-tokyo/chematic/actions/runs/44960741979)
completed with the frozen Linux wheel described in the Platinum holdout
protocol. It selects the first 32 development cases by SHA-256 of case name,
uses 20 attempted candidates per engine, seed 42, 200 MMFF iterations and up
to 400 retry iterations. Unsupported cases and errors remain in its output.
This pilot neither scores the candidate holdout nor establishes an acceptance
threshold. Larger structures may require a separate runtime protocol. A
pilot run is bounded to 100 minutes and uploads its checkpoint and exit status
even if it times out; this run completed all 32. An earlier
unbounded [pilot](https://github.com/kent-tokyo/chematic/actions/runs/44959489605)
was cancelled before it produced results so the bounded run could retain
partial evidence.

RDKit reported MMFF parameter support for all 32 selected records; chematic
reported no runtime error, though its zero-candidate case does not establish
successful force-field construction. At the
20-attempt budget, chematic and RDKit recovered 2 and 4 references under the
1 Å heavy-atom RMSD, specified-stereo and geometry screen; lowest-energy
recovery was 1 and 1. Across 603/607 returned candidates, neither engine had
a specified-stereo violation. Chematic returned no candidate for one case;
RDKit returned two for it. Another case produced 3/5 candidates. The paired
molecule-bootstrap interval for the recovery difference was −15.63 to 0
percentage points; this small development sample is not an independent
noninferiority test. Among cases with at least one returned candidate, the
median best heavy-atom RMSD, before stereo or geometry screening, was 2.85 Å
for chematic (31 cases) and 2.91 Å for RDKit (32 cases). Before MMFF, two RDKit
embedded ensembles contained a candidate within 1 Å; raw chematic candidates
were not recorded. Protein-bound
conformations are a distinct population,
and these numbers do not imply that either engine reproduces solution-state
conformer distributions.

The uncompressed result SHA-256 is
`4fedaddc479dd19ac55ac5008dff4d3b4d8119f8b55b03bfd8e21ca135686860`.
The compressed per-record results are retained as
`validation/results/macrocycle-pdb-linux-development-pilot-2026-11-22.json.gz`
(SHA-256 `dadd3f2a8d91dcd037c7540c9ea023353771bf7de55a31591a8f79fc39675bd3`).

Before scoring, report unsupported chemistry and preserve experimental coordinates.
If a mismatch affects the candidate holdout, resolve the inclusion rule and
whether a fresh holdout is needed before viewing any conformer outcomes.
Record source and prepared-data hashes, wheel identities, error handling,
candidate and time budgets, scorer, noninferiority threshold and confidence
interval method before scoring. Use development cases first. At the time this
split and index audit were fixed, no holdout SDF coordinate or conformer
outcome had been inspected. Reference preparation checks coordinate finiteness
without scoring, and no performance claim is made for these partitions.

## Proposed fixed candidate evaluation, before holdout scoring

The candidate population is the 1,243 prepared, scaffold-disjoint holdout
records in the pinned split described above. Protein-bound conformations are
assessed as a separate population; this test does not establish solution-state
sampling or Prime-set equivalence. The structural files, split and preparation
identity audits remain unchanged. Before starting a conformer evaluation,
record exact tested source and Linux/macOS wheel hashes and run the chemistry
support audit without generating conformers. Resolve any input-identity
mismatch before proceeding; if doing so exposes conformer outcomes, create a
fresh scaffold-disjoint holdout rather than silently reusing that record.

For each engine, generate up to 20 candidates with seed 42 using the same
stereo and heuristic geometry checks and the same 1 Å symmetry-aware
heavy-atom RMSD criterion used by the COD fixed-count protocol. Use the
corrected chematic profile with the bounded 400-iteration retry; retain all
failures in the 1,243-record primary denominator. Report the matched-MMFF
support subset separately and report parameterization, embedding,
minimization, no-candidate, nonconvergence, and stereo failures separately.
Apply a paired 95% percentile scaffold-cluster bootstrap interval to the
difference in recovery rates, using 10,000 resamples and seed 42. The
proposed noninferiority gate is a lower interval bound
above −2 percentage points, together with zero detected specified-stereo
violations or unevaluable declared stereo among returned chematic candidates.
Report the raw recovery and
converged-only recovery for both engines even if the noninferiority gate passes.
A low recovery for both engines remains a limitation of the reference and
budget, not proof of adequate absolute quality.

A separate timed comparison uses 0.2 seconds per molecule per engine, up to
20 attempts, the seed schedule `42 * (attempt + 1)`, and the same screens;
its measured time includes embedding, minimization and any retry but excludes
RMSD scoring. Apply the same paired, scaffold-cluster 95% interval and −2
percentage-point margin. Publish host/OS, elapsed-time distributions,
completed attempt counts and errors; report Linux and macOS separately. The
support audit and candidate/time budgets are proposed here before observing
any holdout conformer outcomes; this paragraph does not itself report a test
result or claim acceptance.

Before any holdout conformer scoring, the
[v11 Linux chemistry audit](https://github.com/kent-tokyo/chematic/actions/runs/44985669374)
checked all 1,243 holdout records against the prepared input and canonical
explicit-H graph. All 1,243 graphs matched and both implementations
constructed MMFF fields for all 1,243. Comparing the equivalent
`UFFHasAllMoleculeParams`-style checks, 1,242 records had all UFF atom
parameters in both engines and one had missing parameters in both. The missing
label is `S_5+4`; successful construction of a UFF field alone does not imply
all parameters are present. These are support checks, not conformer outcomes.
The per-record inputs and results are retained in
`validation/results/macrocycle-holdout-support-v11-2026-11-22.json.gz` and
`validation/results/macrocycle-uff-support-v11-2026-11-22.json.gz`. The Linux audit used the v11 wheel and pinned structural files. A
[macOS v11 audit](https://github.com/kent-tokyo/chematic/actions/runs/44987000227)
produced byte-identical support outputs and the same counts.

For a possible fixed-count holdout run, the 1,243 names are partitioned by
SHA-256(name) order into 16 round-robin shards before any conformer outcome
is inspected. Each derived split is a subset of the already audited holdout;
its derived audit retains the parent split and audit hashes. The manifests pin
the source, wheel, scorer, geometry checker, sharding program and fixed-count
summarizer. Their uncompressed SHA-256 values are
`ef294fbef650d7839e7ebf78878654e5c31c2b500d0047fa473ba52744c33059`
for Linux and
`1c5926f87645f9c092941ce04e30dbf30b1fb7a5ca3d8747f5c9da00ed542dec`
for macOS. They are stored in
`validation/results/macrocycle-shard-manifest-v11-linux-2026-11-22.json.gz`
and `validation/results/macrocycle-shard-manifest-v11-macos-2026-11-22.json.gz`.
The earlier manifests were superseded before scoring to pin the summary
helper and require the stricter declared-stereo verifier described below;
the values above identify the replacement manifests. The
sharding utility rejects missing, duplicated or unexpected records,
different protocol settings, mismatched input audits or scorer fingerprints,
and wheel hashes that differ from the manifest. Synthetic integration checks
exercised preparation, scoring, sealing, merging and summary without scoring
any of these holdout records. Timed-run sharding and development results are
still needed before deciding whether to run the candidate holdout.

The timed evaluation can derive the same shard membership only from a
verified, complete fixed-count baseline. Within each shard, cases run in
SHA-256(name) order and engine order alternates. Each timed shard verifies the
wheel, baseline, scorer, geometry code, 0.2-second budget and 20-attempt
limit; merging rejects missing or duplicate cases and mixed inputs or
settings. The timed sharding helper SHA-256 is
`c5f5158e493baf9857ba05fb6d0cb6993397e9187ca27eaf8dd6fdbaa0d7e7ff`,
the timed scorer SHA-256 is
`d777845f55fabe450b0919329b4aa95f5df3ef5bb34012edc0bf821ba09a7ef6`,
and the timed scaffold-summary SHA-256 is
`3ac65dead9bae503c078a8943235bc3dc75cac18647396543b40802872e60c5f`.
A synthetic two-molecule test exercised the fixed-to-timed pipeline and
summary; no macrocycle holdout conformer result was used. A timed manifest
is generated only after the complete fixed-count baseline exists.

The [v11 Linux development pilot](https://github.com/kent-tokyo/chematic/actions/runs/44985344904)
completed the same 32 previously exposed cases with the fixed v11 wheel.
Under the 20-candidate, stereo and geometry screened, 1 Å criterion, chematic
recovered 2 and RDKit 4; both engines recovered the same two cases and RDKit
recovered two additional cases. Chematic produced 603 candidates and RDKit
607, with zero detected stereo violations in either set. Chematic returned no
candidate for one case; RDKit returned at least one for every case. These
aggregate outcomes match the earlier pilot and remain too small and exposed
to establish noninferiority. Full output:
`validation/results/macrocycle-pdb-linux-development-v11-2026-11-22.json.gz`,
uncompressed SHA-256
`821e4cbb342bbb94db5f4a44a62560d17001ead9e1a7f1190f0cb58ab3d14708`.

A [follow-up development experiment](https://github.com/kent-tokyo/chematic/actions/runs/44987000248)
selected the two exposed cases recovered by RDKit but not chematic in the
20-candidate v11 pilot and generated 100 candidates per engine. Both engines
recovered both references, with zero stereo violations detected by the earlier CIP/E/Z check among 200
candidates per engine. The first candidates passing the same 1 Å, stereo and
geometry screen were numbers 38 / 12 for chematic / RDKit on
`8bi2__Ligand_peptide_QR7_A_701` and 35 / 5 on
`6exq__Ligand_organic_C2W_B_303`. The first 20 RMSDs in each 100-candidate
run matched the earlier pilot. This targeted observation suggests the extra
budget helps these cases but does not estimate performance on the remaining
molecules. The two-case derived development split SHA-256 was
`94dd509bbc9bfdc13c07f244eb5ea695c80599a7c38dc4694a8fcf843740e4ab`;
results are in `validation/results/macrocycle-pdb-linux-development-v11-100-2026-11-22.json.gz`,
uncompressed SHA-256
`5af57fccb40f3b1ec60810ae0fd7202b982429c7341cb4db4a6b8a7ae5b5ea16`.


A local macOS structural preflight with RDKit 2026.09.2, before scoring any holdout conformer, found 10,141
specified tetrahedral centers and 132 specified double bonds in the 1,243
records; 26 specified tetrahedral centers in 12 records lacked RDKit CIP
labels. Four additional raw tetrahedral declarations in three records were
not listed by RDKit as potential stereo; all four have four distinct ranked
neighbors. The earlier scorer checked atoms only when CIP labels existed and
therefore did not check those centers or non-tetrahedral declarations. Its
reported zero violations must not be interpreted as comprehensive stereo
validation. The new strict verifier checks declared raw tetrahedral parity
against RDKit's 3D assignment, checks eligible unrecognized four-neighbor
tetrahedra by signed volume, and checks recognized non-tetrahedral
permutations when unambiguous. Degenerate or otherwise unevaluable cases are
marked unknown. Recovery requires both zero contradictions and zero unknowns; the
strict summary reports them separately. Structural preflight output is
`validation/results/macrocycle-holdout-declared-stereo-audit-2026-11-22.json.gz`,
uncompressed SHA-256
`4a446fec8e08e8a753294927ae900b589f6783882acb2e8b194c8680409d1430`.
The strict verifier SHA-256 is
`3b09595c02ede46bccb1aa8b45e62b3e038fcbb30cdb5a974ef8ef0f1acd0132`.
The 225-record COD and 32-record macrocycle development runs must be repeated
with this verifier before interpreting strict stereo or recovery outcomes.

A local RDKit 2026.09.2 check compared the old and strict stereo detectors on
30 exposed COD development molecules, using an embedded geometry and its
mirror for each. No old finding was lost; two geometries gained a detected
contradiction and two gained an unknown. On 50 development molecules,
raw signed-volume parity agreed with RDKit's 3D tag for 236 comparable
center/geometry pairs; 16 were conservatively unknown. This is a small
calibration set, not an independent accuracy estimate. Reproducible per-case
selection and counts are retained in
`validation/results/stereo-verifier-crosscheck-local-2026-11-22.json`.
The final pre-score manifests supersede the previous strict-verifier hashes
to include finite-overflow handling; no holdout conformer scores were viewed.

The [cross-platform verifier check](https://github.com/kent-tokyo/chematic/actions/runs/44990054068)
passed on Linux and macOS with RDKit 2026.09.2 and the frozen v11 wheels.
Both independently reproduced the local development-set findings: zero
old-only findings in 60 geometries, two additional contradictions and two
additional unknowns; signed-volume parity agreed with RDKit's 3D tags for
236 comparable center/geometry pairs, with 16 raw-parity cases marked
unknown. The 11 synthetic unit and integration checks also passed on both
hosts. Per-platform outputs are retained as
`validation/results/stereo-verifier-crosscheck-{linux,macos}-2026-11-22.json`.
These checks do not validate every stereo class or establish holdout recovery.

The raw-parity calibration script was then tightened to count only centers
that RDKit itself recognizes, so a fallback cannot be compared with itself.
The local result remained 236 comparisons with zero disagreements and 16
unknowns. The [tightened cross-platform rerun](https://github.com/kent-tokyo/chematic/actions/runs/44990230618)
passed on Linux and macOS with the same counts; its outputs are
`validation/results/stereo-verifier-crosscheck-v2-{linux,macos}-2026-11-22.json`.
This refinement was made before any macrocycle holdout conformer scoring.

The same local structural audit on the 363 macrocycle development records
found 3,025 specified tetrahedral centers and 27 specified double bonds;
eight tetrahedral centers in four records lacked CIP labels. No raw declared
atom was omitted by RDKit's potential-stereo list in that development set.
This is input inspection only. The per-record audit is
`validation/results/macrocycle-development-declared-stereo-audit-local-2026-11-22.json.gz`,
uncompressed SHA-256
`fb7d8b05c33f2f385777d4e3483b1264496d5cccaef6170509352c6663fe546d`.

The [strict-stereo repeat](https://github.com/kent-tokyo/chematic/actions/runs/44990435587)
of the same two 100-candidate development cases also recovered 2 / 2 for
both engines. Among 200 returned candidates per engine it found zero
contradicted and zero unevaluable declared stereo elements; the first passing
candidate indices remained 38 / 12 and 35 / 5 for chematic / RDKit. This is
still a deliberately selected two-case development experiment. Its raw
results are
`validation/results/macrocycle-pdb-linux-development-v11-100-strict-2026-11-22.json.gz`,
uncompressed SHA-256
`9b53e9ba5e06f1fa8c06d538ae63efcfd06898c0a353b47368d11fc42f7abb1c`;
recovery and strict-stereo summaries are alongside them.

The [strict shard-pipeline smoke test](https://github.com/kent-tokyo/chematic/actions/runs/44990731715)
passed on Linux and macOS using the respective frozen v11 wheels. Each host
regenerated and byte-compared the 16-shard manifest from the audited real
inputs without scoring those inputs. Each then ran two synthetic molecules
through fixed-count and timed generation, sealing, merging, recovery summary
and strict-stereo summary. All four synthetic engine/profile combinations
returned 40 candidates on each host, with no contradictions or unknowns.
These checks establish that this small pipeline example works, not that the
1,243-record experiment or physical-device acceptance has passed. Outputs:
`validation/results/strict-shard-pipeline-smoke-{linux,macos}-2026-11-22.json`.

The [strict-stereo 32-case Linux repeat](https://github.com/kent-tokyo/chematic/actions/runs/44990053908)
completed using the frozen v11 wheel and the strict verifier. Recovery remained
2 / 4 for chematic / RDKit. Among 603 / 607 returned candidates, neither
engine had a contradicted or unevaluable declared stereo element; chematic
returned no candidate for one case. All 32 per-case RMSD arrays for each
engine were identical to the earlier v11 run. The 363-record development
input audit also reproduced the local structural counts. Results are in
`validation/results/macrocycle-pdb-linux-development-v11-32-strict-2026-11-22.json.gz`,
uncompressed SHA-256
`9d6ad9ca14d1539d24433e89c3a20c7d8c4d3ab5f8c060bdace1a14238e23ab1`,
with recovery and strict-stereo summaries alongside. The Linux structural
audit is `validation/results/macrocycle-development-declared-stereo-audit-linux-2026-11-22.json.gz`,
uncompressed SHA-256
`27aa5575e410eee2dfa38043d5f935a51f527e4705b60d974945104dfdf13b8f`.
These remain exposed development results.

The exposed 208-atom `5vk0__Ligand_phospholipids_ALA_L_4` case was
examined with the frozen v11 wheels in a separate [stage diagnostic](https://github.com/kent-tokyo/chematic/actions/runs/44991392140).
On both hosts, compatibility embedding produced two candidates and both
remained after MMFF; RDKit also produced two. The corrected profile returned
zero on Linux and one on macOS. All returned optimized candidates had
unfinished status, and no returned geometry contradicted or made declared
stereo unevaluable according to the strict external verifier. This does not
identify whether the Linux corrected embedding returned no raw candidate or
the subsequent internal stereo screen removed it; the additional diagnostic
is intended to distinguish those cases. Results are recorded in
`validation/results/macrocycle-stage-diagnostic-v11-linux-2026-11-22.json`
(SHA-256 `61baeacf39df21d9e3b4523bbfdb4e805f10b1e503d3b1d36e110cf5910e7575`)
and `validation/results/macrocycle-stage-diagnostic-v11-macos-2026-11-22.json`
(SHA-256 `a23f0333e923421b4b44a1d65d17d38ca633951a07a98a4f5a9d4531a7e804`).

A further [Rust stage diagnostic](https://github.com/kent-tokyo/chematic/actions/runs/44992153299)
checked out the exact source of the v11 wheels. On Linux the corrected
embedding itself returned zero candidates, before MMFF or the stereo gate. On
macOS it returned one; all 14 declared tetrahedral constraints were satisfied
before and after MMFF, with unfinished optimization status. The raw outputs are
`validation/results/macrocycle-raw-corrected-v11-linux-2026-11-22.json`
(SHA-256 `45134bea73df4f4402d0166ca7b4586aca5fb1f7f7c645695fd8fe6bd96d7116`)
and `validation/results/macrocycle-raw-corrected-v11-macos-2026-11-22.json`
(SHA-256 `399d544029b479561a81e659f38505145254e619ddf3046893b80a613e454b4a`).
The earlier diagnostic run did not validate the candidate source tree against
the wheel, and the next run masked a compiler failure behind a shell pipe;
neither provides results. The linked run fixed both problems.

An opt-in fallback API now retries the compatibility gradient only when the
corrected, stereo-checked ensemble returns no usable candidates. It applies the
same internal stereo gate to that retry; the existing corrected API is
unchanged. The [Linux source-built candidate regression](https://github.com/kent-tokyo/chematic/actions/runs/44992673854)
returned zero corrected candidates and two fallback candidates for this case.
Both were unfinished; the independent strict verifier found no contradicted or
unevaluable declarations. The result is
`validation/results/macrocycle-fallback-development-linux-2026-11-22.json`,
SHA-256 `25ced86eeb64345f34ebb7c7ba4ce982c97c6b22667e17807d79d7f504398322`.
The wheel SHA-256 was
`deaa345879874ef2d6b935c39d512807c175cbeaa2e203025f35993b65c51d38`,
built from merge commit `f8c48c8b44321e14f4648b94ec03035ad5ae099d`.
This is one exposed development case; a 32-record development repeat is
scheduled. The previously sealed holdout protocol does not describe this
fallback profile and remains unscored. The 32-record repeat uses up to 20
corrected attempts plus 20 compatibility attempts when fallback is needed,
while RDKit uses 20 attempts. Its fixed-count results therefore do not compare
equal attempt budgets; a separately timed comparison is also required.

The same candidate's macOS arm passed its Rust and Python checks. On the same
exposed 208-atom case the corrected and fallback calls both returned the same
one unfinished candidate, with no contradicted or unevaluable declaration in
the strict verifier. Its result is
`validation/results/macrocycle-fallback-development-macos-2026-11-22.json`
(SHA-256 `cd5ebc205d25eaf64ccb0258e0882813a52c6ad6c72441c7271bc89ba08f0a9a`);
the macOS wheel SHA-256 was
`637154a74bd4eb7000a284c5a6fe2420833dc87d1c9e7de2fb9b50b8569f7b2e`.
Both wheels came from merge commit `f8c48c8b44321e14f4648b94ec03035ad5ae099d`.

The [Linux fallback repeat on 32 exposed development references](https://github.com/kent-tokyo/chematic/actions/runs/44993494986)
completed with all 32 records evaluated. Recovery remained 2 / 4
(chematic / RDKit). Chematic returned 605 candidates versus RDKit's 607;
neither engine had contradicted or unevaluable declared stereo, and neither
had a record with no candidate. Compared with the earlier corrected-profile
run on this host, chematic's RMSD arrays changed only for the 208-atom case,
where its count rose from zero to two. Their reference RMSDs were 6.80 and
6.35 Å, so this did not recover that reference; both were unfinished. All 32
RDKit RMSD arrays matched the earlier run. The fixed candidate limit is 20,
but fallback can spend 20 additional embedding attempts; the result is not an
equal-work comparison. The raw result is
`validation/results/macrocycle-pdb-linux-development-fallback-32-strict-2026-11-22.json.gz`,
uncompressed SHA-256
`6d8590234425b725893b598c93b9c670950aa8934ad2d717120332bb614d5331`;
recovery and strict-stereo summaries are alongside it. This exposed sample is
not independent acceptance evidence.

The corresponding [macOS fallback development repeat](https://github.com/kent-tokyo/chematic/actions/runs/44993927884)
also evaluated all 32 records. Recovery was 3 / 4 (chematic / RDKit), with
605 / 606 candidates and no contradicted or unevaluable declarations in
either engine. Both returned at least one candidate for each record. Its raw
result is
`validation/results/macrocycle-pdb-macos-development-fallback-32-strict-2026-11-22.json.gz`,
uncompressed SHA-256
`c201d5551e0a781f72970ba568dce47c0de369002d61705fffff728e574ac9c5`;
recovery and strict-stereo summaries are alongside it. The Linux and macOS
hosts can produce different coordinates; this is an exposed development
comparison, not a platform-independent recovery guarantee.

In a [separate Linux 2-second-deadline development comparison](https://github.com/kent-tokyo/chematic/actions/runs/44995216291),
recovery was 2 / 4 (chematic / RDKit) across the same 32 records. The median
number of candidates completed before the deadline was 2 / 1. The engines
returned 185 / 164 candidates, with zero contradicted or unevaluable declared
stereo, and seven / eleven records had none before the deadline. Calls were not
interrupted: median actual elapsed time was 2.257 / 2.251 seconds and maxima
were 5.122 / 7.193 seconds; 27 / 29 records overran. Thus this compares
acceptance deadlines, not hard wall-clock limits. Returned unfinished
candidates were eligible for the main recovery measure; restricting recovery
to converged candidates gave 2 / 2. The exposed sample is not independent
acceptance evidence. The raw result is
`validation/results/macrocycle-pdb-linux-fallback-time-2s-2026-11-22.json`,
SHA-256 `3b19f32fbebd868b1d6b609e72e2933d10a594ec34c54edd3ca5e134d6b7aacd`,
with recovery and strict-stereo summaries alongside it.

The [macOS 2-second-deadline repeat](https://github.com/kent-tokyo/chematic/actions/runs/44995558202)
recovered 3 / 4 references, or 2 / 2 using only converged candidates. The
median number of accepted candidates was 2 / 2, with 196 / 172 total.
Neither engine had contradicted or unevaluable declared stereo; each had nine
records with no candidate before the deadline. Calls were not interrupted:
median actual elapsed time was 2.256 / 2.313 seconds and maxima were
5.110 / 8.151 seconds; 27 / 29 records overran. These 32 exposed records are
not independent acceptance evidence. The raw result is
`validation/results/macrocycle-pdb-macos-fallback-time-2s-2026-11-22.json`,
SHA-256 `ebbc4e42084ab81a66e9aa82b9ab276e8f01b3fb442832dc827592422d681a5b`,
with recovery and strict-stereo summaries alongside it.

For the exposed development record `8bi2__Ligand_peptide_QR7_A_701`, which
RDKit alone recovered in the earlier macOS 32-record fallback run, a
shared-start diagnostic used source `63c053ce218dff9418971d570fea899e716f1630`
and the macOS CPython 3.13 wheel with SHA-256
`a23175f9ef80635d7d568023fd2bfadca68abbf9e850ad4c2106305152dbb9fc`.
With 20 initial structures per engine and 200 MMFF94 iterations, neither
optimizer reached 1 Å from a chematic start, while each reached it once from
an RDKit start. These counts are RMSD-only. The per-conformer statuses,
energies and strict stereo checks are in
`validation/results/macrocycle-8bi2-cross-start-macos-63c053ce-2026-11-22.json`.
With the same seed and candidate count, disabling macrocycle torsions alone
or macrocycle 1–4 bounds alone produced no chematic-start geometry within
1 Å after either optimizer. The minimum RMSDs after chematic minimization
were 1.156 Å with both enabled, 1.183 Å without macrocycle torsions and
1.029 Å without macrocycle 1–4 bounds. These one-record exploratory results
are stored as `macrocycle-8bi2-no-torsions-macos-63c053ce-2026-11-22.json`
and `macrocycle-8bi2-no-14-macos-63c053ce-2026-11-22.json` in
`validation/results/` and do not justify a default change.

The 32-record macOS development run was repeated with the `63c053ce` wheel
and corrected profile. Recovery was 3 / 4 for chematic / RDKit; all 605 / 606
returned candidates had evaluable, consistent declared stereo. Candidate
counts and RMSD arrays for both engines matched the earlier macOS fallback
run. This exposed sample is not holdout evidence. The raw result is
`validation/results/macrocycle-pdb-macos-development-63c053ce-corrected-32-2026-11-22.json.gz`,
uncompressed SHA-256
`c01ce49a4e83b911cfd23e459d5cbcb0d4f35d68a18239ea01f7fde776611a10`.

The experimental checked compatibility-gradient wheel
`0af4d69edbcc363e036e550d27a721c4f781d2ea00d3e1ce7cffde0a472e7507`
was also evaluated on the same 32 macOS development records. It recovered
3 / 4 for chematic / RDKit, with the same per-record recovery outcomes as
the corrected baseline. All 605 / 606 candidates had consistent and
evaluable declared stereo. The scorer's `corrected` option is routed to the
experimental method by its wrapper. This sample is exposed and cannot
establish macrocycle noninferiority. Raw output is
`validation/results/macrocycle-pdb-macos-development-compatible-checked-32-2026-11-22.json.gz`,
uncompressed SHA-256
`1f6182c019c41a2e04e00e10d4f33916fd0cbd150814dc93a8715f6ae6ba0022`.

As an exposed development ablation, all seven of the 32 records with
corrected-profile minimum RMSD above 1 Å and at most 1.6 Å were tested with
either macrocycle torsions or macrocycle 1–4 bounds disabled. Each setting
used 20 chematic starting structures and 200 MMFF iterations. Neither
setting produced a 1 Å recovery that passed the strict stereo and geometry
screens for any of the seven; all tested bond and angle checks had reference
parameters. Raw per-candidate results are
`validation/results/macrocycle-near-miss-ablation-macos-2026-11-22.json.gz`,
uncompressed SHA-256
`1ba1ed8818cc5dc08671ca5ca54aad2de222fca24b7e5337891d919dc4b3ecae`.
This selected analysis does not establish the effect on the full population.

The experimental mixed-gradient wheel
`55ff3b6f14d01d7da5af6b00e35ded4ef2df03d22bfd7a4d803f2dc36f51f49b`
also recovered 3 / 4 on the same 32 macOS development records. Recovery was
unchanged for each record relative to corrected; 604 / 606 returned
chematic / RDKit candidates had consistent and evaluable declared stereo.
This exposed sample does not establish macrocycle noninferiority. Raw data are
`validation/results/macrocycle-pdb-macos-development-mixed-checked-32-2026-11-22.json.gz`,
uncompressed SHA-256
`0694f6ea4747e5044bd71eafa14b3ec78b2aa072fd3c89dcb2c7a225e683b3ff`.

With seed 43 and otherwise identical mixed-profile settings on the same 32
macOS records, chematic / RDKit recovered 4 / 6. Relative to chematic seed
42, two recoveries were gained and one lost; all 605 chematic candidates
had consistent and evaluable declared stereo. This selected sample cannot
establish that seed 43 improves macrocycle sampling. The remaining 331
development records were selected without inspecting their outcomes to test
both seeds with the same 20-attempt budget. The seed-43 32-record data are
`validation/results/macrocycle-pdb-macos-development-mixed-seed43-32-2026-11-22.json.gz`,
uncompressed SHA-256
`18f4ee07ea329fc088914df373389e54cb42491fab5ea0bb85d80d679f871140`.
The remaining-development split is
`validation/results/macrocycle-development-remaining331-split-2026-11-22.json.gz`,
uncompressed SHA-256
`7d82fc45028f1220d120a12947babdf7903a11350b7673a9adce5e95ad152caa`.
