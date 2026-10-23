# GRO to LAMMPS atomic geometry with explicit units

`gro_to_lammps_atomic` converts one GRO frame to a LAMMPS `atom_style atomic`
data document. Caller-supplied atom types follow source row order; masses are
a dense table in grams/mole, with type `i+1` using `masses[i]`.

```python
from pathlib import Path
import chematic

report = chematic.gro_to_lammps_atomic(
    Path("geometry.gro").read_text(),
    atom_types=[1, 2], masses=[22.99, 35.45], units="real",
)
Path("geometry.data").write_text(report["data_text"])
```

The LAMMPS input script must independently specify `units real` and
`atom_style atomic` before `read_data geometry.data`. Data files do not declare
the simulation unit system. Rust callers use
`chematic_mol::gro_lammps::gro_to_lammps_atomic` with `GroLammpsUnits::Real`
or `Metal`; the report includes typed `LammpsData` and the original `GroFrame`.

| Quantity | GRO | LAMMPS real | LAMMPS metal |
| --- | --- | --- | --- |
| Position, cell edges, tilts | nm | Å, multiply by 10 | Å, multiply by 10 |
| Velocity | nm/ps | Å/fs, multiply by 0.01 | Å/ps, multiply by 10 |
| Caller mass table | g/mol | unchanged | unchanged |

Both positions and the supplied periodic cell are scaled together. Box sizes
are never recomputed from the atom extent. Orthogonal and restricted triclinic
cells are accepted: `v1=(ax,0,0)`, `v2=(bx,by,0)`, `v3=(cx,cy,cz)` with positive
diagonals. LAMMPS bounds have origin zero, high bounds `10*(ax,by,cz)`, and
tilts `10*(bx,cx,cy)`. These are data-file edge parameters, not a dump file's
enclosing bounds. Other orientations, singular/negative cells and unrepresentable
scaled values raise errors; no rotation or lattice reconstruction is performed.

Coordinates retain their source order and values after scaling, including
negative/outside-cell coordinates. The adapter neither wraps them nor infers
periodic image flags. LAMMPS can subsequently remap coordinates according to its
own boundary settings; this conversion does not reconstruct trajectory history.
Velocities must be present on every atom or none. Mixed presence is refused.

Output IDs are sequential starting at one, so wrapped or repeated GRO atom
numbers cannot collide. `atom_identity` maps each generated ID to the zero-based
source row, original atom number/name and residue number/name. `source_text`
retains the exact GRO input, including CRLF. The title and residue identity are
sidecar data; residue numbers are not interpreted as LAMMPS molecule IDs.

No elements, bonds, charges, force-field coefficients, angles or molecular
topology are inferred from atom names or distances. This is an atomic geometry
adapter, not a complete force-field or molecular simulation setup. It supports
only `real` and `metal`; other unit styles fail explicitly. Limits are 16 MiB
source input, 100,000 atoms, 1024-byte lines and 1..4096 positive finite masses.
Every supplied type must index the mass table. Empty frames are supported with
a nonempty mass table and a valid cell.

## Provenance and validation

The conversion follows the official [GROMACS GRO format](https://manual.gromacs.org/current/reference-manual/file-formats.html#gro),
[LAMMPS restricted triclinic geometry](https://docs.lammps.org/Howto_triclinic.html)
and [LAMMPS unit systems](https://docs.lammps.org/units.html). No implementation
from another conversion toolkit was copied.

[Open Babel #2669](https://github.com/openbabel/openbabel/issues/2669) reports
GRO atom coordinates being scaled while its box is replaced. An original
three-atom probe reproduced this on Open Babel 3.2.1: atom coordinates became
Å correctly, but the supplied 2×3×4 nm box became bounds -0.5..9.5 Å on every
axis. The expected bounds are 0..20, 0..30 and 0..40 Å. chematic had separate
GRO and LAMMPS readers/writers, without this explicit geometry adapter.

Original fixtures/probes are in `validation/openbabel_issues/gro-lammps/`.
Dated outputs, source/artifact hashes and logs are under
`validation/results/gro-lammps-units-2026-10-17/`. Tests cover unit scales,
true triclinic tilts, source-ID collisions, CRLF, velocity presence, invalid
cells/types/masses and resource limits. Native LAMMPS parse/write round trips
and independent analytic geometry checks validate the output; no LAMMPS engine
execution or complete simulation-readiness claim is included.
Local validation passed 868 Rust tests (five new), 75 Python tests (32 new)
and 96 analytic geometry cases across both units, velocity policies, cell
profiles and all six source-row permutations. Strict Clippy, formatting,
binding/stub checks, refactor budgets and strict MkDocs checks passed.
