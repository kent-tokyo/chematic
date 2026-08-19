"""
Materials-science file formats quickstart (Gaussian Cube, mmCIF, LAMMPS dump)
==============================================================================
Shows 3 of the v0.17.0 file-format bindings in realistic use: parsing a
Gaussian Cube scalar-field grid and inspecting it as a 3-D numpy array,
parsing an mmCIF structure file for its atom count/fields, and resolving
real Cartesian positions from a LAMMPS dump frame (including the triclinic
box case).

Run:
    python examples/materials_formats_quickstart.py

Dependencies:
    pip install chematic numpy
"""
import chematic
from chematic import VolumetricGrid, LammpsDumpFrame

# A tiny hand-authored Gaussian Cube file: a 2x2x2 grid of density values
# around a single oxygen atom (not a real DFT calculation -- just enough
# structure to demonstrate the API).
CUBE_TEXT = (
    "Water density\n"
    "Generated for chematic examples\n"
    "1    0.000000    0.000000    0.000000\n"
    "2    1.000000    0.000000    0.000000\n"
    "2    0.000000    1.000000    0.000000\n"
    "2    0.000000    0.000000    1.000000\n"
    "8    8.000000    0.500000    0.500000    0.500000\n"
    "0.0 1.0 2.0 3.0\n"
    "4.0 5.0 6.0 7.0\n"
)

# A tiny mmCIF structure: 2 atoms, no bond table (mmCIF's _atom_site
# category never carries connectivity).
MMCIF_TEXT = (
    "data_EXAMPLE\n"
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
    "ATOM   2  H  H1  HOH A 1.500 2.500 3.500\n"
)

# One frame of a LAMMPS dump/trajectory, orthogonal box, plain x/y/z columns.
LAMMPS_DUMP_TEXT = (
    "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n"
    "0.0 10.0\n0.0 10.0\n0.0 10.0\nITEM: ATOMS id x y z\n"
    "1 1.0 2.0 3.0\n2 4.0 5.0 6.0\n"
)


def cube_example() -> None:
    print("=== Gaussian Cube ===")
    grid = VolumetricGrid.from_cube(CUBE_TEXT)
    print(f"shape={grid.shape} units={grid.units} atoms={len(grid.atoms)}")
    # values_3d reshapes the flat `values` array to `grid.shape`, so
    # values_3d[i, j, k] == grid.get(i, j, k).
    v3d = grid.values_3d
    print(f"values_3d.shape={v3d.shape}, values_3d[1, 0, 0]={v3d[1, 0, 0]}")
    assert v3d[1, 0, 0] == grid.get(1, 0, 0)
    print()


def mmcif_example() -> None:
    print("=== mmCIF ===")
    result = chematic.parse_mmcif(MMCIF_TEXT)
    print(f"atom_count={len(result['atoms'])}")
    for atom in result["atoms"]:
        print(f"  {atom['element']:>2} ({atom['x']:.3f}, {atom['y']:.3f}, {atom['z']:.3f})")
    print()


def lammps_dump_example() -> None:
    print("=== LAMMPS dump ===")
    frame = chematic.parse_lammps_dump_frame(LAMMPS_DUMP_TEXT)
    print(f"timestep={frame.timestep} num_atoms={frame.num_atoms}")
    positions = frame.cartesian_positions()
    for i, (x, y, z) in enumerate(positions):
        print(f"  atom {i}: ({x}, {y}, {z})")
    print()


def main() -> None:
    cube_example()
    mmcif_example()
    lammps_dump_example()


if __name__ == "__main__":
    main()
