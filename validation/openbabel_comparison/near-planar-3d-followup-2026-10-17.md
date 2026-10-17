# Nearly planar 3D centers: source-gated replay (2026-10-17)

The [second browser holdout](stereo-extended-holdouts-2026-10-17.md) exposed
two false tetrahedral assignments in CheMatic's original 3D CTAB reader.
Records 136 and 188 have three distinct heavy-atom neighbors around carbon
and a small nonzero out-of-plane displacement, but the RDKit 2025.09.3
source/readback gate leaves their tetrahedral stereo unspecified. Open Babel
3.2.1 assigned both centers too. The earlier report correctly preserves the
pre-fix 198/200 result; this is a separate post-fix replay.

For three measured bond vectors from the center, the reader now requires
`abs(det(v0, v1, v2)) / (|v0| |v1| |v2|) >= 0.1` before using the sign to
infer a tetrahedral tag. This dimensionless test is independent
of coordinate scale. In the two pinned 200-record browser SDF sets, the five
unspecified candidate carbons had scores at most 0.0252, while the 339
assigned candidate carbons had scores at least 0.4300. The older 300-record
stress sets had 204 assigned candidate carbons with a minimum of 0.3601.
These observations motivate the threshold; they do not prove it is optimal
for all molecular geometries.

## Source-readback results

All 700 source-gated molecules were read in both CTAB dialects. Exact browser
SDFs and manifests are committed with the
[first browser holdout](unmarked-3d-browser-holdout-2026-10-17.md) and
[second browser holdout](stereo-extended-holdouts-2026-10-17.md). The 300
stress records were regenerated with the seeds and exclusion rule in the
[original 3D comparison](unmarked-3d-stereo-2026-10-16.md); their four SDF
SHA-256 values matched that document exactly. The source-built CheMatic
`sdf_v3000_accuracy_dump` binary had SHA-256
`ac086d25561cc48155ee04067374e1cb8840371657003631f172756c4beb9e36`.
The [per-row replay](results/near-planar-3d-replay-2026-10-17.json.gz)
includes graph, stereochemistry, atom-map, and Open Babel rewrite checks.

| Input | CheMatic V2000 | Open Babel V2000 | CheMatic V3000 | Open Babel V3000 |
| --- | ---: | ---: | ---: | ---: |
| Browser holdout 1 (200) | 200 | 197 | 200 | 197 |
| Browser holdout 2 (200) | 200 | 198 | 200 | 198 |
| Stress seed 20261016 (100) | 100 | 95 | 100 | 94 |
| Stress seed 20261017 (200) | 200 | 195 | 200 | 194 |
| **Total (700)** | **700** | **685** | **700** | **683** |

Counts are RDKit-normalized, atom-map-insensitive isomeric-SMILES agreement
with each source molecule after RDKit itself successfully reads the same
unmarked 3D CTAB. CheMatic did not lose an earlier passing source record.
The two near-planar records are pinned as V2000 and V3000 regression fixtures
in `crates/chematic-mol/tests/fixtures/`.

The comparison covers these readers, input distributions, and normalization
rules. It does not establish a toolkit-wide ranking or superiority to RDKit,
which supplied the source and readback gate. COSMolKit is not rerun in this
follow-up; its result on the second 200-record holdout remains in the earlier
report.
