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
violations among returned chematic candidates. Report the raw recovery and
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
`validation/results/macrocycle-uff-support-v11-2026-11-22.json.gz`. The audit
used the Linux v11 wheel and pinned structural files; macOS support has not
been separately established.

For a possible fixed-count holdout run, the 1,243 names are partitioned by
SHA-256(name) order into 16 round-robin shards before any conformer outcome
is inspected. Each derived split is a subset of the already audited holdout;
its derived audit retains the parent split and audit hashes. The manifests pin
the source, wheel, scorer, geometry checker, sharding program and fixed-count
summarizer. Their uncompressed SHA-256 values are
`163d98d3600d0db2c2e55b6a88a1f13a0b8d98f9125d1799c3c12aa9958f80a2`
for Linux and
`303b28e515421cce8ea95528839d396b1123c768057db3033613c682558794d2`
for macOS. They are stored in
`validation/results/macrocycle-shard-manifest-v11-linux-2026-11-22.json.gz`
and `validation/results/macrocycle-shard-manifest-v11-macos-2026-11-22.json.gz`.
The sharding utility rejects missing, duplicated or unexpected records,
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
`bf7a7314bae5dbda5f940972a3d6f5ed17e815b0d8c6851b2335d9f29b320603`,
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
