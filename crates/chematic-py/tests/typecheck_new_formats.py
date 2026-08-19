"""
mypy --strict smoke script for the 7 v0.17.0 file-format Python bindings
(mmCIF, PQR, ORCA input/output, QCSchema, Gaussian Cube, OpenDX, LAMMPS
data/dump), including `values_3d` (added in this same PR).

NOT a pytest test (deliberately not named `test_*.py` so pytest's default
collection in this directory doesn't pick it up) -- it's a standalone
type-check target. Run:

    maturin build --release -m crates/chematic-py/Cargo.toml -o /tmp/chematic-wheel-check
    python -m venv /tmp/chematic-typecheck-venv
    /tmp/chematic-typecheck-venv/bin/pip install /tmp/chematic-wheel-check/chematic-*.whl mypy numpy
    cd /tmp && /tmp/chematic-typecheck-venv/bin/mypy \\
        <absolute path to this file> --strict

Run with an absolute path from outside the repo (as above) so mypy resolves
`chematic` from the venv's installed wheel, not this repo's source tree --
the whole point of this check is proving type info is genuinely shipped in
the distributed package, not just present in the source tree.

This intentionally exercises only the 7 new-format bindings' stub entries,
not the full `chematic` API surface (see this PR's scope notes).
"""
from typing import Any

import numpy as np
import numpy.typing as npt

import chematic
from chematic import VolumetricGrid, LammpsDumpFrame

MMCIF_TEXT = (
    "data_TEST\n"
    "loop_\n"
    "_atom_site.group_PDB\n"
    "_atom_site.id\n"
    "_atom_site.type_symbol\n"
    "_atom_site.label_atom_id\n"
    "_atom_site.label_comp_id\n"
    "_atom_site.label_asym_id\n"
    "_atom_site.Cartn_x\n"
    "_atom_site.Cartn_y\n"
    "_atom_site.Cartn_z\n"
    "ATOM   1  O  O1  HOH A 1.000 2.000 3.000\n"
)

# --- mmCIF ---
mmcif_result: dict[str, Any] = chematic.parse_mmcif(MMCIF_TEXT)
mmcif_mol: chematic.Mol = mmcif_result["mol"]
mmcif_text: str = chematic.write_mmcif([{"element": "C", "x": 0.0, "y": 0.0, "z": 0.0}])

# --- PQR ---
pqr_result: dict[str, Any] = chematic.parse_pqr(
    "ATOM      1  N   ALA     1     -0.966   1.523   1.412 -0.400  1.500\n"
)
pqr_text: str = chematic.write_pqr(pqr_result["atoms"])
inferred: "str | None" = chematic.infer_element("ATOM", "ALA", "CA")

# --- ORCA ---
orca_input_result: dict[str, Any] = chematic.parse_orca_input("* xyz 0 1\nO 0 0 0\n*\n")
orca_input_text: str = chematic.write_orca_input(orca_input_result)
orca_output_result: dict[str, Any] = chematic.parse_orca_output(
    "****ORCA TERMINATED NORMALLY****\n"
)

# --- QCSchema ---
qc_result: dict[str, Any] = chematic.parse_qcschema_molecule(
    '{"schema_name":"qcschema_molecule","schema_version":1,"symbols":["H"],'
    '"geometry":[0.0,0.0,0.0],"molecular_charge":0.0,"molecular_multiplicity":2}'
)
qc_text: str = chematic.write_qcschema_molecule(qc_result)

# --- Gaussian Cube / OpenDX (VolumetricGrid) ---
grid: VolumetricGrid = VolumetricGrid(
    origin=(0.0, 0.0, 0.0),
    axes=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    shape=(2, 2, 2),
    values=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
)
values_flat: npt.NDArray[np.float64] = grid.values
values_3d: npt.NDArray[np.float64] = grid.values_3d
axes_arr: npt.NDArray[np.float64] = grid.axes
grid_shape: "tuple[int, int, int]" = grid.shape
grid_units: str = grid.units
one_value: "float | None" = grid.get(0, 0, 0)
cube_text: str = grid.to_cube()
grid_mol, grid_coords = grid.to_molecule()
grid_mol_typed: chematic.Mol = grid_mol
grid_coords_typed: "list[list[float]]" = grid_coords

LAMMPS_DATA_TEXT = (
    "comment\n"
    "\n"
    "1 atoms\n"
    "1 atom types\n"
    "\n"
    "0.0 10.0 xlo xhi\n"
    "0.0 10.0 ylo yhi\n"
    "0.0 10.0 zlo zhi\n"
    "\n"
    "Atoms # atomic\n"
    "\n"
    "1 1 5.0 5.0 5.0\n"
)

# --- LAMMPS data/dump ---
lammps_data_result: dict[str, Any] = chematic.parse_lammps_data(LAMMPS_DATA_TEXT, "atomic")
lammps_data_text: str = chematic.write_lammps_data(lammps_data_result)

dump_frame: LammpsDumpFrame = LammpsDumpFrame(
    timestep=0,
    box_bounds={"lo": (0.0, 0.0, 0.0), "hi": (10.0, 10.0, 10.0), "tilt": None},
    column_names=["id", "x", "y", "z"],
    rows=[[1.0, 1.0, 1.0, 1.0]],
)
dump_timestep: int = dump_frame.timestep
dump_positions: "npt.NDArray[np.float64] | None" = dump_frame.cartesian_positions()
dump_column: "list[float] | None" = dump_frame.column("x")


def main() -> None:
    print("typecheck_new_formats: all bindings type-check under mypy --strict")


if __name__ == "__main__":
    main()
