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
error message. Status `0` and `1` have the same meaning as above and their
error field is `None`. These embedded coordinates have not undergone MMFF
optimization or an independent geometry-quality screen.

`etkdg_embed_multiple_confs_mixed_checked_diagnostics` has the same options
and returns `(conformers, embedded, force_field_error, stereo_violated,
stereo_unevaluable, invalid_coordinates)`. When MMFF setup fails, the last
three values count embedded candidates rejected for each reason; a candidate
with both violated and unevaluable stereo is counted as violated. When MMFF
setup succeeds, or no candidate is embedded, these counts are `None` because
the optimized path does not classify all its rejections. An empty conformer
list with an `embedded` count greater than zero can therefore be distinguished
from a failure to embed any candidate.

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
when one is available. An unreadable checkpoint is reported as
`checkpoint_error` while preserving the original worker status and error.

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

`repair_mmff_iters` defaults to zero. A positive value allows that many extra
MMFF iterations when attempting to repair a near-linear three-coordinate
nitrogen; a repair is retained only if its checks pass. The example uses 400
to match the ongoing evaluation, not as a universally optimal setting.
The [evaluation protocol](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-22-platinum-holdout-protocol.md)
and [macrocycle evaluation](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-22-macrocycle-pdb-evaluation.md)
record the tested candidates, datasets, screens and results. The 0.2-second
acceptance deadline used in timed evaluations does not interrupt an attempt.
