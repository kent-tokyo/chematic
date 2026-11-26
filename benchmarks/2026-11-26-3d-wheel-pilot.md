# 3D wheel and MMFF-unavailable pilot

These macOS arm64 results are development diagnostics on previously exposed
COD data, not an independent holdout. macOS 27.0.1 arm64, Python 3.13.6 and
NumPy 2.5.3 were used, with one thread per
engine, at most 20 attempts and a 0.2 second worker deadline per molecule.
The process startup and geometry/RMSD scoring were outside the deadline.
The candidate was built from main `7bd29b3` with the initial fallback changes,
before the later rejection-diagnostics refactor, using release optimization,
LTO disabled and 16 codegen units.

The published chematic 1.2.5 wheel SHA-256 is
`556fa68c8121e4ed49184b403cfd38a88acc213689d58b9e1a67ed9efce80e6c`;
the RDKit 2026.09.2 wheel SHA-256 is
`3860ba53505f27cb24d05383027f78e3f82aefd6caf193c217a94d876d91b023`.
The candidate shared library SHA-256 is
`9e5a44d1c56040b19fee9df3015a54cf7296265ed083ef04ffeb7ca3306d9031`.
The COD SDF SHA-256 is
`6ad24cf292a62561e55d3556395b841440f9d0d8588c7f2b85c8fb7e4675af74`.
The fixed-count selection input SHA-256 is
`72d2ccabd7bba29105d0cc6fda9894a9826a075397a09665df8330ad9ffe4211`.
The earlier support audit SHA-256 is
`f716ec708915a0357f9cf00eba2551a9687d485d31c52a8932c781c256bfc848`.

The first pilot took 20 records in SHA-256(name) order from the earlier COD
fixed-count set. On the published wheels, 13 chematic and 10 RDKit records
met the 1 Å, stereo and supported-geometry criteria; 10 and 8 respectively
also had converged MMFF. On a separate run the candidate and RDKit returned
the same respective counts. These are small, previously exposed samples and
do not establish relative performance. Candidate and wheel build settings
also differ.

The second pilot selected 20 records in the same order from 157 earlier
COD support-audit records for which chematic MMFF was unavailable at audit
time. Their SDF indices and names came from that audit, and their SMILES were
regenerated from the corresponding SDF structures with RDKit 2026.09.2.
The published wheel produced 19 chematic generation errors and one
record with a returned conformer. The candidate returned conformers for
19 records; 18 contained status -1 conformers with unavailable MMFF.
Sixteen of those 18 met the 1 Å, stereo and supported-geometry criteria.
The remaining input returned no candidate before its deadline. RDKit met
the same criteria for 17 records in the candidate run. This screen only
covers the bond and angle parameters reported in each record, and one of
the 20 inputs lacked full coverage. The fallback conformers were not MMFF
optimized.

The complete set of 157 support-audit records was subsequently run with the
same worker deadline and settings. The published chematic wheel returned an
accepted conformer for 1 record, had 151 generation errors and hit the
deadline on 6 records. The candidate returned an accepted conformer for 150
records, including 149 with MMFF-unavailable status, with no generation
errors. Among these 157, the full supported-geometry and stereo screen plus
1 Å RMSD was met for 1 published-wheel record and 99 candidate records. Of
the candidate recoveries, 98 included an MMFF-unavailable conformer. RDKit
recovered 93 in each run. Neither engine had a converged-MMFF recovery in
this selection. Thirteen records lacked complete geometry parameter coverage
and therefore could not pass this screen. Deadlines were reached for 114
candidate runs and 140 RDKit runs were forcibly stopped at the deadline;
returned candidates completed before the
deadline. The JSON records expose the timeout and the lower-bound attempted
count. Each number uses all 157 records as the denominator.

The 240 molecules with chematic generation errors in the earlier 3,579-row
Linux comparison were also rerun on this macOS host. The old result file's
SHA-256 is `2dc5f09bd074f827aeb9acd40612f131bf2124ffd5611db3c72e777ada5cfd70`.
It recorded 191 boron and 49 selenium type-assignment errors. The selection
manifest is `chematic-pr1282-240-input.json`; its source SDF hash is
`6e2f93fb05c8ff49aff24577aa62d3f04f6531b505e12a9f711c9b79cf3ba008`.
The published wheel had 236 setup errors and four deadlines with no accepted
coordinates. The candidate returned MMFF-unavailable coordinates for 229
molecules and recorded no generation errors. Eight of the remaining 11 ran
out of time without coordinates and three completed the attempt budget with
no accepted candidate. All status -1 conformers included an error message.
The candidate met the 1 Å, stereo and full supported-geometry criteria for
125 of 240 molecules, and RDKit met them for 114 of 240 in each run. Neither
engine had a converged-MMFF recovery. Twenty-five inputs lacked complete
geometry parameter coverage and could not pass the full screen. Candidate
and RDKit were forcibly stopped at the deadline on 178 and 223 inputs, respectively, including
those that returned earlier candidates. The longest candidate elapsed time
including worker termination was 0.2042 seconds. These results describe the
fixed corpus and host; they do not establish a hard upper bound on all hosts.

The three inputs that exhausted 20 attempts without a returned candidate
were indices 44, 7180 and 8931. Their compatibility embedding produced
coordinates in a separate diagnostic run, but the checked fallback rejected
all embedded candidates in a four-attempt regression test because declared
stereo was unevaluable by `verify_stereo`. Two inputs declare square-planar
stereo; the third has a chiral selenium center. This is a verifier coverage
limit, so silently accepting these coordinates as stereo-correct would not be
justified. The benchmark's independent stereo scorer also checks only a subset
of stereo specifications. It must not be taken as validation of square-planar
or all selenium stereochemistry.

An additional three-input diagnostic rerun used the later diagnostics API,
an optimized shared library with SHA-256
`dcb61d43c61ddf068ec08c6ea1f97ce08e1f59b485d630b39129fb3bb5ff0c52`,
and benchmark script SHA-256
`4e2032636d6c66a9f641605a6b8b1a1805f8fc2f7df4a1d976e76fe1096813f0`.
It recorded 20, 20 and 19 unevaluable-stereo rejections for indices 7180,
44 and 8931 respectively, with no stereo violations reported. The last
worker was stopped at the deadline, so its rejection count is a lower bound.

A subsequent change added verification for a three-single-bond, positively
charged selenium center with a stereo lone pair. The previously rejected
index 8931 then returned all 20 embedded candidates in a source-build
diagnostic without a time budget; RDKit's 3D CIP check agreed on each.
The square-planar cases remain unevaluable. Regression tests also compare
both selenium enantiomers and an alternate atom order with RDKit and check
one RDKit geometry against its mirror. This focused check does not establish
coverage of other selenium coordination or independent holdout performance.

For the remaining square-planar examples, the two experimentally referenced
trans-neighbor angles were approximately 169.2° and 174.7°. In a separate
20-attempt compatibility-embedding diagnostic, the corresponding angles
ranged from 111.4° to 122.3° and 94.9° to 137.4°, respectively. None reached
150°. These examples need square-planar constraints in coordinate generation
as well as a verifier; merely accepting unevaluable stereo would return the
wrong coordination geometry. These angle observations are development
diagnostics, not a general square-planar acceptance threshold.

Per-record files and the 157-record selection are in
`validation/results/2026-11-26-3d-wheel-pilot/`. Each output records the
hash of the benchmark script used for that run. The 157- and 240-record runs
used the script snapshot `platinum_time_budget_full_run.py` in that directory.
The candidate runs used
`--profile mixed-fallback`; with one attempt per call the mixed profile uses
the corrected torsion gradient. Pilot timing excludes worker startup and the
older pilot output did not yet record worker termination cleanup time.
