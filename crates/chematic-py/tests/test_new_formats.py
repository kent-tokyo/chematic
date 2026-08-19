"""Tests for the v0.18.0 Python bindings of the 7 file-format modules
chematic-mol gained in v0.17.0 (mmCIF, PQR, ORCA input/output, QCSchema,
Gaussian Cube, OpenDX, LAMMPS data/dump) -- previously zero Python exposure.

Fixture strings are hand-authored/adapted from the corresponding Rust-side
fixtures in crates/chematic-mol/src/{mmcif,pqr,orca,cube,opendx,lammps_data,
lammps_dump}.rs (already-verified-valid inputs for each format's grammar),
kept inline and minimal here per this test suite's existing convention.
"""
import pytest

import chematic
from chematic import VolumetricGrid, LammpsDumpFrame


# ---------------------------------------------------------------------------
# mmCIF
# ---------------------------------------------------------------------------

MMCIF_FIXTURE = """\
data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM   1  O  O1  HOH A 1.000 2.000 3.000
ATOM   2  H  H1  HOH A 1.500 2.500 3.500
"""


def test_parse_mmcif_basic():
    result = chematic.parse_mmcif(MMCIF_FIXTURE)
    assert isinstance(result["mol"], chematic.Mol)
    assert len(result["atoms"]) == 2
    assert result["atoms"][0]["element"] == "O"
    assert result["atoms"][0]["res_name"] == "HOH"
    assert result["atoms"][0]["x"] == pytest.approx(1.0)
    assert result["cell"] is None


def test_mmcif_round_trip():
    result = chematic.parse_mmcif(MMCIF_FIXTURE)
    text = chematic.write_mmcif(result["atoms"])
    result2 = chematic.parse_mmcif(text)
    assert len(result2["atoms"]) == len(result["atoms"])
    assert result2["atoms"][0]["element"] == result["atoms"][0]["element"]
    assert result2["atoms"][1]["x"] == pytest.approx(result["atoms"][1]["x"])


def test_mmcif_unknown_element_raises_value_error():
    with pytest.raises(ValueError):
        chematic.write_mmcif([{"element": "Zz", "x": 0.0, "y": 0.0, "z": 0.0}])


def test_mmcif_malformed_input_raises_value_error():
    with pytest.raises(ValueError):
        chematic.parse_mmcif("not an mmcif file at all")


# ---------------------------------------------------------------------------
# PQR
# ---------------------------------------------------------------------------

PQR_FIXTURE = """\
ATOM      1  N   ALA     1     -0.966   1.523   1.412 -0.400  1.500
ATOM      2  CA  ALA     1      0.257   0.679   1.911  0.100  1.700
HETATM    3  ZN  ZN      2      5.000   5.000   5.000  2.000  1.090
"""


def test_parse_pqr_basic():
    result = chematic.parse_pqr(PQR_FIXTURE)
    assert len(result["atoms"]) == 3
    assert result["atoms"][1]["element"] == "C"  # "CA" -> C, not Ca
    assert result["atoms"][2]["element"] == "Zn"
    assert result["coords"][0] == pytest.approx([-0.966, 1.523, 1.412])


def test_pqr_round_trip():
    result = chematic.parse_pqr(PQR_FIXTURE)
    text = chematic.write_pqr(result["atoms"])
    result2 = chematic.parse_pqr(text)
    assert len(result2["atoms"]) == 3
    assert result2["atoms"][2]["element"] == "Zn"


def test_pqr_element_inferred_when_omitted():
    text = chematic.write_pqr([
        {"atom_name": "CA", "res_name": "ALA", "res_seq": 1, "x": 0.0, "y": 0.0, "z": 0.0},
    ])
    result = chematic.parse_pqr(text)
    assert result["atoms"][0]["element"] == "C"


def test_infer_element():
    assert chematic.infer_element("ATOM", "ALA", "CA") == "C"
    assert chematic.infer_element("HETATM", "ZN", "ZN") == "Zn"
    assert chematic.infer_element("ATOM", "X", "123") is None


def test_pqr_malformed_input_raises_value_error():
    with pytest.raises(ValueError):
        chematic.parse_pqr("ATOM      1  N   ALA     1     -0.966   1.523\n")


# ---------------------------------------------------------------------------
# ORCA input / output
# ---------------------------------------------------------------------------

ORCA_INPUT_FIXTURE = "! B3LYP def2-SVP\n! Opt Freq\n\n* xyz 0 1\nO 0 0 0\n*\n"

ORCA_OUTPUT_FIXTURE = (
    "\n Total Charge           Charge          ....    0\n"
    " Multiplicity           Mult            ....    1\n\n"
    "-------------------------   --------------------\n"
    "FINAL SINGLE POINT ENERGY       -76.025678123456\n"
    "-------------------------   --------------------\n\n"
    "****ORCA TERMINATED NORMALLY****\n"
)


def test_parse_orca_input():
    result = chematic.parse_orca_input(ORCA_INPUT_FIXTURE)
    assert "B3LYP" in result["keywords"]
    assert "Opt" in result["keywords"]
    assert result["coords"]["kind"] == "xyz"
    assert result["coords"]["charge"] == 0
    assert result["coords"]["multiplicity"] == 1
    assert result["coords"]["atoms"][0]["element"] == "O"


def test_orca_input_round_trip():
    result = chematic.parse_orca_input(ORCA_INPUT_FIXTURE)
    text = chematic.write_orca_input(result)
    result2 = chematic.parse_orca_input(text)
    assert result2["keywords"] == result["keywords"]
    assert result2["coords"]["atoms"][0]["element"] == "O"


def test_orca_input_malformed_raises_value_error():
    with pytest.raises(ValueError):
        chematic.parse_orca_input("* xyz notanumber 1\nC 0 0 0\n*\n")


def test_parse_orca_output():
    result = chematic.parse_orca_output(ORCA_OUTPUT_FIXTURE)
    assert result["charge"] == 0
    assert result["multiplicity"] == 1
    assert result["final_energy_hartree"] == pytest.approx(-76.025678123456)
    assert result["termination"]["kind"] == "normal"
    assert result["optimization_convergence"] == "not_requested"


def test_orca_output_non_finite_energy_raises_value_error():
    with pytest.raises(ValueError):
        chematic.parse_orca_output(
            "-------------------------   --------------------\n"
            "FINAL SINGLE POINT ENERGY       NaN\n"
            "-------------------------   --------------------\n"
        )


# ---------------------------------------------------------------------------
# QCSchema
# ---------------------------------------------------------------------------

# r(O-H) = 0.9584 Angstrom equilibrium water geometry, converted to Bohr
# independently (matches crates/chematic-mol/tests/qcschema.rs's own fixture).
WATER_BOHR_GEOMETRY = [
    0.0, 0.0, 0.0,
    0.0, 1.431544636036637, 1.109419726497757,
    0.0, -1.431544636036637, 1.109419726497757,
]

WATER_QCSCHEMA_MOLECULE = {
    "schema_name": "qcschema_molecule",
    "schema_version": 1,
    "symbols": ["O", "H", "H"],
    "geometry": WATER_BOHR_GEOMETRY,
    "molecular_charge": 0.0,
    "molecular_multiplicity": 1,
    "connectivity": [[0, 1, 1.0], [0, 2, 1.0]],
}


def test_parse_qcschema_molecule():
    import json
    result = chematic.parse_qcschema_molecule(json.dumps(WATER_QCSCHEMA_MOLECULE))
    assert result["symbols"] == ["O", "H", "H"]
    assert isinstance(result["mol"], chematic.Mol)
    assert result["mol"].formula == "H2O"
    assert len(result["coords"]) == 3


def test_qcschema_molecule_round_trip():
    import json
    result = chematic.parse_qcschema_molecule(json.dumps(WATER_QCSCHEMA_MOLECULE))
    text = chematic.write_qcschema_molecule(result)
    result2 = chematic.parse_qcschema_molecule(text)
    assert result2["symbols"] == result["symbols"]
    assert result2["geometry"] == pytest.approx(result["geometry"])


def test_chematic_to_qc_molecule_and_back():
    # from_xyz gives explicit H atoms in the graph (unlike from_smiles("O"),
    # whose H count is implicit) so `coords` matches `mol`'s atom count.
    # Checked via heavy_atoms (1 oxygen), not `.formula` -- `.formula`
    # additionally counts *implicit* H filled in for unsatisfied valence,
    # which depends on from_xyz's distance-based bond perception (an
    # unrelated concern from what this test exercises: QCSchema round-trip
    # fidelity of the atoms actually present).
    xyz = "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.76 0.59\nH 0.0 -0.76 0.59\n"
    mol, coords = chematic.from_xyz(xyz)
    qc = chematic.chematic_to_qc_molecule(mol, coords)
    assert qc["symbols"] == ["O", "H", "H"]
    mol2, coords2, charge, mult = chematic.qc_molecule_to_chematic(qc)
    assert mol2.heavy_atoms == 1
    assert len(coords2) == 3
    assert charge == 0.0
    assert mult == 1


def test_qcschema_molecule_malformed_raises_value_error():
    with pytest.raises(ValueError):
        chematic.parse_qcschema_molecule('{"symbols": ["H"], "geometry": [1.0]}')


def test_parse_atomic_input():
    import json
    doc = {
        "molecule": WATER_QCSCHEMA_MOLECULE,
        "driver": "energy",
        "model": {"method": "hf", "basis": "sto-3g"},
    }
    result = chematic.parse_atomic_input(json.dumps(doc))
    assert result["driver"] == "energy"
    assert result["model"]["method"] == "hf"
    assert isinstance(result["mol"], chematic.Mol)


def test_atomic_input_round_trip():
    import json
    doc = {
        "molecule": WATER_QCSCHEMA_MOLECULE,
        "driver": "energy",
        "model": {"method": "hf"},
        "keywords": {"scf_type": "df"},
    }
    result = chematic.parse_atomic_input(json.dumps(doc))
    text = chematic.write_atomic_input(result)
    result2 = chematic.parse_atomic_input(text)
    assert result2["driver"] == "energy"
    assert result2["keywords"] == {"scf_type": "df"}


def test_atomic_input_malformed_raises_value_error():
    with pytest.raises(ValueError):
        chematic.parse_atomic_input("not json at all")


def test_parse_atomic_result():
    import json
    doc = {
        "molecule": WATER_QCSCHEMA_MOLECULE,
        "driver": "energy",
        "model": {"method": "hf"},
        "provenance": {"creator": "test"},
        "properties": {},
        "return_result": -76.02,
        "success": True,
    }
    result = chematic.parse_atomic_result(json.dumps(doc))
    assert result["success"] is True
    assert result["return_result"] == pytest.approx(-76.02)


def test_atomic_result_round_trip():
    import json
    doc = {
        "molecule": WATER_QCSCHEMA_MOLECULE,
        "driver": "energy",
        "model": {"method": "hf"},
        "provenance": {"creator": "test"},
        "properties": {},
        "return_result": -76.02,
        "success": True,
    }
    result = chematic.parse_atomic_result(json.dumps(doc))
    text = chematic.write_atomic_result(result)
    result2 = chematic.parse_atomic_result(text)
    assert result2["success"] is True
    assert result2["return_result"] == pytest.approx(-76.02)


def test_atomic_result_success_inconsistency_raises_value_error():
    import json
    doc = {
        "molecule": WATER_QCSCHEMA_MOLECULE,
        "driver": "energy",
        "model": {"method": "hf"},
        "provenance": {"creator": "test"},
        "success": True,
        # success=True but no return_result -- schema violation.
    }
    with pytest.raises(ValueError):
        chematic.parse_atomic_result(json.dumps(doc))


# ---------------------------------------------------------------------------
# Gaussian Cube
# ---------------------------------------------------------------------------

CUBE_FIXTURE = (
    "Water density\n"
    "Generated for chematic tests\n"
    "1    0.000000    0.000000    0.000000\n"
    "2    1.000000    0.000000    0.000000\n"
    "2    0.000000    1.000000    0.000000\n"
    "2    0.000000    0.000000    1.000000\n"
    "8    8.000000    0.500000    0.500000    0.500000\n"
    "0.0 1.0 2.0 3.0\n"
    "4.0 5.0 6.0 7.0\n"
)


def test_parse_cube():
    grid = VolumetricGrid.from_cube(CUBE_FIXTURE)
    assert grid.shape == (2, 2, 2)
    assert grid.units == "bohr"
    assert list(grid.values) == pytest.approx([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    assert len(grid.atoms) == 1
    assert grid.atoms[0][0] == "O"
    assert grid.get(1, 0, 0) == pytest.approx(4.0)


def test_cube_round_trip():
    grid = VolumetricGrid.from_cube(CUBE_FIXTURE)
    text = grid.to_cube()
    grid2 = VolumetricGrid.from_cube(text)
    assert grid2.shape == grid.shape
    assert list(grid2.values) == pytest.approx(list(grid.values))


def test_cube_malformed_raises_value_error():
    with pytest.raises(ValueError):
        VolumetricGrid.from_cube("not a cube file\nat all\n")


# ---------------------------------------------------------------------------
# OpenDX
# ---------------------------------------------------------------------------

OPENDX_FIXTURE = (
    "# Data from APBS\n"
    "# POTENTIAL (kT/e)\n"
    "object 1 class gridpositions counts 2 2 2\n"
    "origin -1.0 -1.0 -1.0\n"
    "delta 0.5 0.0 0.0\n"
    "delta 0.0 0.5 0.0\n"
    "delta 0.0 0.0 0.5\n"
    "object 2 class gridconnections counts 2 2 2\n"
    "object 3 class array type double rank 0 items 8 data follows\n"
    "0.0 1.0 2.0\n"
    "3.0 4.0 5.0\n"
    "6.0 7.0\n"
    'attribute "dep" string "positions"\n'
    'object "regular positions regular connections" class field\n'
    'component "positions" value 1\n'
    'component "connections" value 2\n'
    'component "data" value 3\n'
)


def test_parse_opendx():
    grid = VolumetricGrid.from_opendx(OPENDX_FIXTURE)
    assert grid.shape == (2, 2, 2)
    assert grid.units == "angstrom"
    assert grid.atoms == []
    assert list(grid.values) == pytest.approx([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])


def test_opendx_round_trip():
    grid = VolumetricGrid.from_opendx(OPENDX_FIXTURE)
    text = grid.to_opendx()
    grid2 = VolumetricGrid.from_opendx(text)
    assert grid2.origin == pytest.approx(grid.origin)
    assert list(grid2.values) == pytest.approx(list(grid.values))


def test_opendx_bohr_grid_fails_closed_but_lossy_succeeds():
    # PR #335 regression guard: a Bohr-unit grid must be refused by the
    # strict writer and only accepted through the explicitly-named lossy one.
    # Non-zero origin so the Bohr->Angstrom conversion is actually visible
    # (an all-zero origin would pass the old, vacuous assertion either way).
    grid = VolumetricGrid(
        origin=(1.0, 0.0, 0.0),
        axes=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        shape=(2, 2, 2),
        values=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        units="bohr",
    )
    with pytest.raises(ValueError):
        grid.to_opendx()
    text = grid.to_opendx_lossy()
    grid2 = VolumetricGrid.from_opendx(text)
    # Bohr -> Angstrom conversion applied to origin/axes (lengths)...
    assert grid2.origin[0] == pytest.approx(0.529177210903)
    # ...never to values.
    assert list(grid2.values) == pytest.approx(list(grid.values))


def test_opendx_grid_with_atoms_fails_closed_even_lossy():
    grid = VolumetricGrid(
        origin=(0.0, 0.0, 0.0),
        axes=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        shape=(1, 1, 1),
        values=[0.0],
        atoms=[("C", 6.0, (0.0, 0.0, 0.0))],
        units="angstrom",
    )
    with pytest.raises(ValueError):
        grid.to_opendx()
    with pytest.raises(ValueError):
        grid.to_opendx_lossy()


def test_volumetric_grid_to_molecule():
    grid = VolumetricGrid.from_cube(CUBE_FIXTURE)
    mol, coords = grid.to_molecule()
    assert isinstance(mol, chematic.Mol)
    assert len(coords) == 1


def test_values_3d_axis_order():
    # A non-cubic (2,3,4) shape so a transposed reshape can't coincidentally
    # pass: values = 0..23 in flat (k-fastest) order, per checked_index(i,j,k)
    # = i*shape[1]*shape[2] + j*shape[2] + k (chematic_mol::VolumetricGrid).
    values = list(range(24))
    grid = VolumetricGrid(
        origin=(0.0, 0.0, 0.0),
        axes=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        shape=(2, 3, 4),
        values=[float(v) for v in values],
    )
    v3d = grid.values_3d
    assert v3d.shape == (2, 3, 4)
    # Every (i, j, k) must agree with checked_index/get -- the strongest
    # oracle available, but shares its formula with the reshape under test,
    # so it alone can't catch a self-consistently-wrong reshape.
    for i in range(2):
        for j in range(3):
            for k in range(4):
                assert v3d[i, j, k] == pytest.approx(grid.get(i, j, k))
    # Independent, hand-computed spot checks a transpose/wrong reshape would
    # break: with k fastest, index(i,j,k) = i*12 + j*4 + k.
    assert v3d[0, 1, 0] == pytest.approx(4.0)   # index 0*12+1*4+0 = 4
    assert v3d[1, 0, 0] == pytest.approx(12.0)  # index 1*12+0*4+0 = 12
    assert v3d[0, 0, 1] == pytest.approx(1.0)
    assert v3d[1, 2, 3] == pytest.approx(23.0)  # last element


def test_values_3d_matches_flat_values_reshaped():
    grid = VolumetricGrid.from_cube(CUBE_FIXTURE)
    flat = list(grid.values)
    v3d = grid.values_3d
    assert v3d.shape == grid.shape
    nx, ny, nz = grid.shape
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                idx = grid.checked_index(i, j, k)
                assert v3d[i, j, k] == pytest.approx(flat[idx])


# ---------------------------------------------------------------------------
# LAMMPS data
# ---------------------------------------------------------------------------

LAMMPS_DATA_FIXTURE = (
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


def test_parse_lammps_data():
    data = chematic.parse_lammps_data(LAMMPS_DATA_FIXTURE, "atomic")
    assert data["counts"]["atoms"] == 1
    assert data["atom_style"] == "atomic"
    assert data["box"]["lo"] == pytest.approx((0.0, 0.0, 0.0))
    assert data["box"]["tilt"] is None
    assert len(data["atoms"]) == 1
    assert data["atoms"][0]["x"] == pytest.approx(5.0)
    assert data["atoms"][0]["charge"] is None


def test_lammps_data_round_trip():
    data = chematic.parse_lammps_data(LAMMPS_DATA_FIXTURE, "atomic")
    text = chematic.write_lammps_data(data)
    data2 = chematic.parse_lammps_data(text, "atomic")
    assert data2["atoms"][0]["x"] == pytest.approx(data["atoms"][0]["x"])
    assert data2["box"]["hi"] == pytest.approx(data["box"]["hi"])


def test_lammps_data_unsupported_atom_style_raises_value_error():
    with pytest.raises(ValueError):
        chematic.parse_lammps_data(LAMMPS_DATA_FIXTURE, "sphere")


def test_lammps_data_malformed_raises_value_error():
    with pytest.raises(ValueError):
        chematic.parse_lammps_data("not a lammps data file\n", "atomic")


# ---------------------------------------------------------------------------
# LAMMPS dump
# ---------------------------------------------------------------------------

LAMMPS_DUMP_FIXTURE = (
    "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n"
    "0.0 10.0\n0.0 10.0\n0.0 10.0\nITEM: ATOMS id type x y z\n1 1 5.0 5.0 5.0\n"
)


def test_parse_lammps_dump_frame():
    frame = chematic.parse_lammps_dump_frame(LAMMPS_DUMP_FIXTURE)
    assert isinstance(frame, LammpsDumpFrame)
    assert frame.timestep == 0
    assert frame.num_atoms == 1
    assert frame.column_names == ["id", "type", "x", "y", "z"]
    assert frame.column("x") == pytest.approx([5.0])
    positions = frame.cartesian_positions()
    assert positions is not None
    assert list(positions[0]) == pytest.approx([5.0, 5.0, 5.0])


def test_parse_lammps_dump_all():
    two_frames = LAMMPS_DUMP_FIXTURE + LAMMPS_DUMP_FIXTURE.replace("TIMESTEP\n0", "TIMESTEP\n1")
    frames = chematic.parse_lammps_dump_all(two_frames)
    assert len(frames) == 2
    assert frames[0].timestep == 0
    assert frames[1].timestep == 1


def test_lammps_dump_round_trip():
    frame = chematic.parse_lammps_dump_frame(LAMMPS_DUMP_FIXTURE)
    text = chematic.write_lammps_dump_frame(frame)
    frame2 = chematic.parse_lammps_dump_frame(text)
    assert frame2.timestep == frame.timestep
    assert frame2.column("x") == pytest.approx(frame.column("x"))


def test_write_lammps_trajectory():
    frame = chematic.parse_lammps_dump_frame(LAMMPS_DUMP_FIXTURE)
    text = chematic.write_lammps_trajectory([frame, frame])
    frames = chematic.parse_lammps_dump_all(text)
    assert len(frames) == 2


def test_lammps_dump_malformed_raises_value_error():
    with pytest.raises(ValueError):
        chematic.parse_lammps_dump_frame("ITEM: NOT TIMESTEP\n0\n")


def test_box_bounds_round_trip_utilities():
    box = chematic.box_bounds_to_true((0.0, 0.0, 0.0), (10.0, 10.0, 10.0), (1.0, 0.5, 0.0))
    lo, hi = chematic.true_to_box_bounds(box)
    assert lo == pytest.approx((0.0, 0.0, 0.0))


# ---------------------------------------------------------------------------
# Cross-language parity fixtures (Binding Quality Pack, v0.18.0)
#
# Same 4 fixtures (verbatim CUBE_FIXTURE/OPENDX_FIXTURE/MMCIF_FIXTURE above,
# plus the LAMMPS triclinic frame below) and the same hardcoded expected
# values as crates/chematic-mol/tests/format_binding_parity.rs (Rust) and
# crates/chematic-wasm/tests/format_parity.test.mjs (WASM) -- each computed
# once, independently, not by trusting another binding's output. See the
# Rust file's module doc comment for why this independently-hardcoded
# approach (rather than a shared fixture file) still proves 3-way parity.
# ---------------------------------------------------------------------------


def test_parity_cube_fixture_summary():
    grid = VolumetricGrid.from_cube(CUBE_FIXTURE)
    assert grid.shape == (2, 2, 2)
    assert grid.units == "bohr"
    assert grid.point_count() == 8
    assert len(grid.atoms) == 1
    # First/last values (reversed-flatten tripwire) plus two interior
    # values only a correctly k-fastest-ordered flatten reproduces.
    values = list(grid.values)
    assert values[0] == pytest.approx(0.0)
    assert values[-1] == pytest.approx(7.0)
    assert grid.get(0, 1, 0) == pytest.approx(2.0)
    assert grid.get(1, 0, 0) == pytest.approx(4.0)


def test_parity_opendx_fixture_summary():
    grid = VolumetricGrid.from_opendx(OPENDX_FIXTURE)
    assert grid.shape == (2, 2, 2)
    assert grid.units == "angstrom"
    assert grid.point_count() == 8
    assert len(grid.atoms) == 0
    # axes: a Bohr<->Angstrom unit-conversion bug would show up here directly.
    # pytest.approx doesn't support nested sequences, so flatten first.
    flat_axes = [v for row in grid.axes for v in row]
    assert flat_axes == pytest.approx([0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5])
    assert grid.get(0, 1, 0) == pytest.approx(2.0)
    assert grid.get(1, 0, 0) == pytest.approx(4.0)


def test_parity_mmcif_fixture_summary():
    result = chematic.parse_mmcif(MMCIF_FIXTURE)
    assert len(result["atoms"]) == 2
    # MMCIF_FIXTURE has no occupancy column at all, so this exercises both
    # the coordinate field-mapping AND the spec-mandated occupancy default
    # (1.0) in one pair of assertions.
    assert result["atoms"][0]["element"] == "O"
    assert result["atoms"][0]["x"] == pytest.approx(1.0)
    assert result["atoms"][0]["y"] == pytest.approx(2.0)
    assert result["atoms"][0]["z"] == pytest.approx(3.0)
    assert result["atoms"][0]["occupancy"] == pytest.approx(1.0)


def test_parity_lammps_dump_triclinic_fixture_summary():
    # Built directly via the LammpsDumpFrame constructor (box_bounds is
    # already the resolved TRUE box, matching chematic_mol::LammpsDumpFrame
    # exactly) rather than parsed from hand-written dump-file text: a dump
    # file's "ITEM: BOX BOUNDS" line carries the *bound* box, not the true
    # box, and hand-deriving that conversion for a fixture file would be
    # exactly the kind of step this parity check exists to make impossible
    # to get wrong. Same triclinic box/tilt/xs values as
    # chematic_mol::lammps_dump's own triclinic_frame() test fixture.
    frame = LammpsDumpFrame(
        timestep=2000,
        box_bounds={"lo": (0.0, 0.0, 0.0), "hi": (10.0, 10.0, 10.0), "tilt": (2.0, 1.0, 0.5)},
        column_names=["id", "xs", "ys", "zs"],
        rows=[[1.0, 0.5, 0.5, 0.5]],
        boundary_flags=("pp", "ff", "ss"),
    )
    assert frame.num_atoms == 1
    positions = frame.cartesian_positions()
    assert positions is not None
    # Hand-computed (see crates/chematic-mol/tests/format_binding_parity.rs
    # for the derivation): x=6.5, y=5.25, z=5.0.
    assert list(positions[0]) == pytest.approx([6.5, 5.25, 5.0])
