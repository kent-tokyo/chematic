# Element-stratified common-coordinate force-field development comparison

This development probe extends the common-coordinate checks to boron,
selenium and other elements. Before scoring, the selector chose the lowest
SHA-256 of record name within each category from the existing COD development
partition, using category precedence B, Se, P/S, then other. It excluded
molecules above 30 heavy atoms and those with an atom-ring basis member of at
least 12 atoms. Caps were 20, 20, 30 and 30 respectively; the available
selenium category had only 10 molecules. The resulting 90 records include all
20 from the earlier ordinary MMFF development comparison and may overlap other
development experiments. They are not an independent holdout or a random
population sample. Selection, source hash, SMILES and category labels are
recorded in `validation/results/2026-11-27-forcefield-hetero/input.json`.

RDKit 2026.09.2 generated explicit-hydrogen ETKDGv3 coordinates with seed 42.
The first run evaluated those same coordinates in both engines. A second run
applied independently sampled uniform displacements in [-0.35, 0.35] Å to
each Cartesian coordinate, with a deterministic SHA-256-derived per-record
seed, and assigned the resulting coordinates to both force fields. These
stressed coordinates need not be chemically realistic and are not evidence of
conformer generation accuracy. Atom order and bond graph were checked first.

Both MMFF94 and MMFF94s supported the same 60 of 90 records in both engines;
all 20 boron and all 10 selenium records were unsupported by both. Both UFF
implementations accepted all 90. No record raised a processing error. Among
the 60 common MMFF records, numeric atom types matched, and maximum absolute
partial-charge difference was 2.23e-16 elementary charge. The largest
componentwise differences in energy and gradient were:

| Coordinates | Force field | Maximum energy difference (kcal/mol) | Maximum gradient difference (kcal/mol/Å) |
| --- | --- | ---: | ---: |
| Embedded | MMFF94 or MMFF94s | 1.28e-13 | 3.59e-7 |
| Embedded | UFF | 3.41e-13 | 3.46e-7 |
| Perturbed | MMFF94 or MMFF94s | 1.37e-12 | 1.10e-11 |
| Perturbed | UFF | 1.37e-12 | 1.82e-12 |

MMFF bond, angle, stretch-bend, out-of-plane, torsion, van der Waals and
electrostatic terms were also compared separately. Six evenly spaced gradient
components per supported molecule were checked against central finite
differences with step 1e-5 Å; the maximum sampled errors under perturbed
coordinates were 2.03e-5 for MMFF and 2.17e-5 kcal/mol/Å for UFF. These
finite-difference checks sample only six components and do not prove that all
gradients are correct. No general accuracy or parameter-coverage claim follows
from this small, selectively sampled development set.

The chematic candidate wheel SHA-256 was
`f9a04b8dfaf3d64cc173ee884cc4835191a545db1d014e80ab970c0677ddbdeb`
with extension SHA-256
`17f0d46b52b0eb6e8bb30edae084f0534404f01a0483549029502c5923a2403f`;
it is a local source build, not a published wheel. The host ran macOS 27.0.1
arm64 and Python 3.13.6; the reference RDKit wheel SHA-256 was
`3860ba53505f27cb24d05383027f78e3f82aefd6caf193c217a94d876d91b023`.
The scorer SHA-256 was
`90d4cd0b4a1075fb842beff46144638af5f8f8730edffa68b7cc03963777b67c`.
Exact coordinates, individual errors and support failures are stored alongside
the input in `validation/results/2026-11-27-forcefield-hetero/`. The selector
is `scripts/validation/select_force_field_development.py`; the scorer is
`scripts/benchmarks/platinum_force_fields.py`.
