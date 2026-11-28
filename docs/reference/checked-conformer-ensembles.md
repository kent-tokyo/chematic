# Checked conformer ensembles

Experimental Python methods for generating and MMFF94-minimizing ETKDG
conformers with additional checks are available in v1.2.4 and later. They are
opt-in and do not establish output equivalence to RDKit.

```python
import chematic

mol = chematic.from_smiles("C[C@H](O)F").add_hydrogens()
conformers = mol.etkdg_embed_multiple_confs_mmff_mixed_checked(
    20, random_seed=42, max_mmff_iters=200, repair_mmff_iters=400
)

for status, energy, coordinates in conformers:
    if status == 0:
        print(energy, coordinates)
```

The returned coordinates are in atom order and include the explicit
hydrogens. Each tuple contains the minimization status, MMFF94 energy in
kcal/mol, and coordinates in angstroms. Status zero requires a finite energy
and coordinates and a final MMFF gradient with no component greater than
0.01 kcal/mol/angstrom in magnitude. A nonzero status means the candidate
was not confirmed converged; it can still be returned. Some attempts may
fail or be discarded, so the returned list can be shorter than the requested
count or empty. Invalid options or force-field setup failures raise
`ValueError` in the Python binding.

`etkdg_embed_multiple_confs_mixed_checked_with_fallback` uses the same
arguments as the mixed checked method and returns `(status, energy,
coordinates, force_field_error)`. When MMFF setup fails, stereo-checked
embedded coordinates are returned with status `-1`, `None` energy and an
error message. Status `-1` with a `None` error instead reports a minimizer
line-search failure. Status `0` and `1` have the same meaning as above and their
error field is `None`. Coordinates returned because MMFF was unavailable have
not undergone MMFF optimization or an independent geometry-quality screen.

`etkdg_embed_multiple_confs_mixed_checked_diagnostics` has the same options
and returns `(conformers, embedded, force_field_error, stereo_violated,
stereo_unevaluable, invalid_coordinates)`. When MMFF setup fails, the last
three values count embedded candidates rejected for each reason; a candidate
with both violated and unevaluable stereo is counted as violated. When MMFF
setup succeeds, the counts classify discarded candidates whose original
embedding failed the declared-stereo or coordinate check after the optimized
geometry also failed that check. An optimized geometry replaced by a safe
original embedding is not counted as discarded. A non-finite fallback energy
is not included in these counts, so they are not necessarily exhaustive. If no
candidate is embedded, the counts are `None`. An empty conformer list with an
`embedded` count greater than zero can therefore be distinguished from a
failure to embed any candidate.

For a wall-clock limit around a whole Python call, use the opt-in process
wrapper:

```python
result = chematic.generate_conformers_timed("CB(C)C", 20, timeout_seconds=5.0)
if result["status"] == "completed":
    print(result["conformers"], result["force_field_error"])
```

The wrapper parses the SMILES and adds explicit hydrogens in a fresh Python
process. Its returned coordinates follow that molecule's atom order. A
`timeout` result means the worker was killed and, by default, contains no partial
ensemble; `error` reports an input or calculation exception and
`worker_error` reports a failed worker. A completed result includes the same
diagnostic counts and per-conformer status and optional energy as the
diagnostic API above. The deadline includes process startup; process creation,
scheduling and cleanup can make the observed elapsed time exceed it. This
wrapper needs a Python executable that can import the same chematic package;
the worker checks the package path and reports an error if it differs.

To retain completed attempts when the deadline expires, opt in explicitly:

```python
result = chematic.generate_conformers_timed(
    "CCCCCCCC", 20, timeout_seconds=1.0, return_partial=True
)
print(result["status"], result["attempts_completed"], result["conformers"])
```

This mode runs one checked attempt at a time and checkpoints each completed
attempt in a temporary directory. On timeout, the interrupted attempt is
excluded and `attempts_completed` counts only finished attempts, including
those that produced no accepted conformer. A timeout before the first
checkpoint returns an empty ensemble. Each attempt uses the corrected ETKDG
gradient: this is a different search profile from the mixed batch, so the
same seed and count do not imply the same ensemble. Each completed attempt is
appended to a temporary journal; if termination interrupts a write, its
incomplete final record is ignored. Atom order, individual
conformer statuses and force-field error semantics are unchanged. The
checkpoint has filesystem and serialization overhead; no general speed or
accuracy improvement is implied. Worker errors are reported separately and retain the last complete checkpoint
when one is available. If a complete journal record is malformed, the valid
prefix is retained and `checkpoint_error` identifies the first bad line; later
records are not used. An unreadable checkpoint is also reported as
`checkpoint_error`. In either case the original worker status and error are
preserved. If the first record is malformed, no checkpoint data is returned.

In this mode, `attempt_diagnostics` contains one record per completed attempt,
also on timeout when a checkpoint is available. Each record includes its
attempt number, embedded and returned counts, the force-field error and the
three rejection counts. `conformer_status_counts` distinguishes `converged`,
`unfinished`, `minimizer_failed`, `force_field_unavailable`, and `other`
returned conformers. `unfinished` includes any status 1 result, including a
safe embedded fallback or a failed final gradient check; it does not imply that
the iteration limit alone caused it. Rejection counts may be `None`, and do not
include every possible reason for discarding a candidate. The interrupted
attempt is not a completed diagnostic record; its separately requested preview
may appear in `in_progress`.

For each completed attempt, `unclassified_rejections` is the number of embedded
candidates that were neither returned nor included in the three rejection
counts or the optional `regular_geometry_rejected` count. The top-level field
sums that number over completed attempts. If a
rejection count is unknown and a remaining discarded candidate could belong to
it, the corresponding value is `None`, as is the aggregate. A zero is still
reported when no discarded candidates remain. The number does not identify the
cause: a non-finite fallback energy, for example, can contribute to it.
Impossible candidate totals are reported as a checkpoint error; on timeout or
worker failure, earlier valid records remain available.

Set `conformer_require_regular_geometry=True` to screen each returned conformer
after the ordinary checked-ensemble processing. This applies to both normal
and `return_partial=True` modes, including retained force-field fallback
coordinates. Rows that fail the conservative declared-stereo and tetrahedral
geometry screen are omitted. `regular_geometry_rejected` reports the number of
omitted rows separately from the backend's rejection counters. In partial mode
it appears in each completed `attempt_diagnostics` record and as an aggregate;
an interrupted attempt has not yet been counted. On timeout without partial
results the counter is unavailable. A zero counter does not establish that
any conformer was generated, that all stereo classes were checked, or that the
remaining geometry and force-field energies are physically valid. This screen
can also withhold plausible structures. The option is independent of the
preview screen below and is included in the worker's deadline.

`regular_geometry_rejection_reasons` splits that count into disjoint categories:
`stereo_violated`, `distorted_tetrahedral_geometry`, and `stereo_unevaluable`.
When a rejected conformer has multiple issues, the first applicable category
in that order is used; the categories are counts of conformers, not individual
stereo elements. These counts describe the additional post-processing screen
and do not replace or overlap the backend rejection counters. Coordinates
rejected by this screen are not returned. Invalid coordinates or an inconsistent
diagnostic cause a worker error instead of being guessed into a category.
Per-attempt and aggregate reason counts are retained in partial checkpoints.
When recovering an older checkpoint that records rejections without their
breakdown, the aggregate breakdown is `None` rather than guessed; a zero
rejection count still has a known all-zero breakdown.

```python
result = chematic.generate_conformers_timed(
    "C[C@H](O)F", 10, timeout_seconds=5.0, return_partial=True,
    conformer_require_regular_geometry=True,
)
print(result.get("regular_geometry_rejected"))
```

In partial mode, `regular_geometry_rejected` counts completed attempts only;
zero may mean no attempt completed. If the check itself fails, the worker
reports an error and previously checkpointed attempts remain available.
With partial recovery, checkpoints missing the selected screen's rejection
count are rejected instead of silently treated as zero.

To retain a provisional embedding even when its MMFF optimization is
interrupted, set `return_preview=True` together with `return_partial=True`.
On timeout or worker failure, `result.get("in_progress")` may contain an
`attempt` number and `coordinates`. These coordinates passed the library's
declared-stereo check but have not completed MMFF optimization; they have no
energy or convergence status and are not added to `conformers`. Unevaluable
or violated declared stereo is not published as a preview. This check does
not cover every possible stereochemical specification or establish general
geometric quality. The `in_progress.geometry` object reports `bond_count`,
`worst_bond_length` in Å, `bond_violation_rate_50pct`, and
`gross_clash_count`. Bond references are approximate covalent-radius estimates;
the violation rate counts bonds whose lengths differ from those estimates by
more than 50%. Gross clashes count non-bonded atom pairs below 1.2 Å,
including hydrogens, 1–3 pairs and interfragment pairs. For a molecule with
no bonds, worst length and violation rate are `None`. These indicators are
not a pass/fail gate and may miss or flag structures differently from other
geometry screens, especially for unusual elements or bonding. If diagnostic
calculation fails, `geometry` is `None` and `geometry_error` describes the
failure; the stereo-checked preview remains available. The optimizer reuses
the embedded coordinates rather than embedding again. This mode also prepares
distance bounds and torsion knowledge once per worker and reuses them across
attempts; preparation is included in the deadline. Completed attempts replace
their provisional record.

Set `preview_require_regular_geometry=True` together with `return_preview=True`
and `return_partial=True` to also require the conservative tetrahedral geometry
screen before publishing a preview. A retained preview then includes
`stereo_screen: "regular_geometry"` in `in_progress`. A preview may be absent
because embedding failed, the stereo/geometry check did not pass, or the worker
was interrupted before its checkpoint was written; absence does not identify
which occurred. The screen can withhold plausible structures and does not prove
chemical validity. See [3D stereo diagnostics](stereo-3d-diagnostics.md) for
its limits. This option affects previews only; use
`conformer_require_regular_geometry=True` separately for completed conformers.

The mixed method divides the requested attempts between the compatibility
and corrected ETKDG torsion gradients, using the first half for compatibility
and the remainder for corrected. The two modes share the same attempt seed
schedule and MMFF94 minimization. For a compatibility-gradient ensemble
without that split, use
`mol.etkdg_embed_multiple_confs_mmff_compatible_checked(...)` with the same
arguments. Both methods apply the source's declared-stereo checks. If a
minimized structure fails a check, the implementation can retain the original
embedded structure with a nonzero status, or omit the candidate. These checks
do not establish correctness for every possible stereochemical specification.
For square-planar declarations, the separate [angle diagnostic](square-planar-geometry.md)
measures available ligand angles without changing this acceptance policy.
The verifier also supports a [limited neutral selenium class](neutral-selenium-stereo.md).

`repair_mmff_iters` defaults to zero. A positive value allows that many extra
MMFF iterations when attempting to repair a near-linear three-coordinate
nitrogen; a repair is retained only if its checks pass. The example uses 400
to match the ongoing evaluation, not as a universally optimal setting.
The [evaluation protocol](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-22-platinum-holdout-protocol.md)
and [macrocycle evaluation](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-22-macrocycle-pdb-evaluation.md)
record the tested candidates, datasets, screens and results. The 0.2-second
acceptance deadline used in timed evaluations does not interrupt an attempt.
