# Preview geometry diagnostics: development comparison

The opt-in timed preview now reports the number of bonds, worst bond length,
fraction of bonds more than 50% away from the library's approximate
covalent-radius bond reference, and number of non-bonded atom pairs under
1.2 Å. The pair count includes hydrogens, 1–3 pairs and interfragment pairs.
These are observations, not a gate or a claim of physically valid geometry.

For an exploratory comparison, the diagnostic and the separate RDKit UFF
geometry screen were applied to the same 90 previously exposed,
element-stratified development molecules, both on their embedded coordinates
and on a separately saved deterministic ±0.35 Å Cartesian perturbation.
Inputs and coordinates are in
`validation/results/2026-11-27-forcefield-hetero/`; these are neither a new
holdout nor necessarily realistic coordinates. Six molecules had incomplete
UFF bond or angle parameter coverage and were excluded from the following
cross-tabulation. For this comparison only, an alert means either a bond
violation or a clash; the RDKit screen's angle violations are not included
because the preview diagnostic does not measure angles.

| Coordinates, 84 fully covered molecules | Both alert | Only preview alert | Only UFF screen alert | Neither alerts |
| --- | ---: | ---: | ---: | ---: |
| Embedded | 0 | 0 | 0 | 84 |
| Perturbed | 61 | 4 | 5 | 14 |

The two checks use different bond references, thresholds and atom-pair
coverage. A disagreement does not establish which check is correct. In
particular, the preview values cannot replace an independent geometry screen,
assess all chemistry or establish conformer accuracy. The non-gating API also
returns the coordinates with `geometry_error` if its diagnostic itself fails.
All 180 coordinate records were processed; the six incompletely covered
molecules in each coordinate condition remain in the per-record results.

The tested local candidate wheel SHA-256 was
`4d23fadd0b8c1258b157447206bc6a211b066d1873a3ed390ca952ae156b41e3`;
its extension SHA-256 was
`ed3083dfeb7b03ae9949a19de3eeb5064b67bca62e82c36cff6643c00dae829e`.
This is a source-built development wheel, not a published release. The host
was macOS 27.0.1 arm64, Python 3.13.6, and RDKit 2026.09.2. Inputs,
per-record outputs, scorer hash and runtime metadata are in
`validation/results/2026-11-27-preview-geometry/development.json`; the runner
is `scripts/benchmarks/preview_geometry_development.py`.

## Timed API development check

The same candidate wheel was also exercised with the 20 previously exposed
flexible development molecules and the existing timed-preview protocol. Each
molecule was run once with previews disabled and once with previews enabled,
alternating order, with a 0.2-second deadline including Python startup. All
40 calls timed out. The no-preview mode returned at least one completed
candidate for 14/20 inputs; preview mode returned a completed or in-progress
candidate for 15/20. Three previews were returned, all with diagnostics and
no diagnostic alerts or diagnostic errors. All three also passed the
independent strict-stereo and coarse UFF geometry screen; their heavy-atom
reference RMSDs were 1.285, 1.718 and 4.102 Å. These small-sample counts do
not establish speed or accuracy improvement and should not be compared as a
controlled speed measurement with earlier wheel runs.

The 56 completed or preview coordinates from this run were used to recompute
the recorded independent screens; all matched on the same host. The timed
per-record coordinates, diagnostics, screens and provenance are in
`validation/results/2026-11-27-preview-geometry/timed-development.json`; the
runner is `scripts/benchmarks/timed_preview_development.py`.
