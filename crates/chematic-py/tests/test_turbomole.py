"""Synthetic MIT OR Apache-2.0 TURBOMOLE lattice preservation contracts."""
import itertools
import math
import pytest
import chematic

TEXT = "$periodic 3\n$lattice angs\n4 1 0\n-1 5 2\n1 0 -6\n$coord frac\n0.5 0.25 -0.5 si f\n0 0 0 cl\n$end\n"

def emit(g):
    return chematic.write_turbomole(g["elements"], g["positions_angstrom"], lattice_angstrom=g["lattice_angstrom"], frozen=g["frozen"])

def test_fractional_coordinate_and_cell_roundtrip():
    g = chematic.parse_turbomole(TEXT)
    assert g["elements"] == ["Si", "Cl"]
    assert g["positions_angstrom"][0] == pytest.approx([1.25, 1.75, 3.5])
    assert g["frozen"] == [True, False]
    assert g["periodicity"] == 3
    assert g["source_text"] == TEXT
    back = chematic.parse_turbomole(emit(g))
    assert back["elements"] == g["elements"]
    assert back["frozen"] == g["frozen"]
    for key in ("positions_angstrom", "lattice_angstrom"):
        for a, b in zip(back[key], g[key]):
            assert a == pytest.approx(b, abs=1e-12)

@pytest.mark.parametrize("order", list(itertools.permutations(range(3))))
def test_all_group_orders(order):
    groups = ["$periodic 3\n", "$lattice angs\n2 0 0\n0 3 0\n0 0 4\n", "$coord angs\n1 2 3 c f\n"]
    g = chematic.parse_turbomole("".join(groups[i] for i in order) + "$end\n")
    assert g["positions_angstrom"] == [[1.0, 2.0, 3.0]]
    assert g["periodicity"] == 3

def test_cell_parameter_conversion():
    text = "$periodic 3\n$cell angs\n4 5 6 90 60 90\n$coord frac\n0.5 0.25 -0.5 c\n$end\n"
    g = chematic.parse_turbomole(text)
    assert g["positions_angstrom"][0] == pytest.approx([0.5, 1.25, -1.5*math.sqrt(3)])

@pytest.mark.parametrize("text", [
    "$coord\n0 0 0 c\n", "$coord\n$end", "$end", "$coord frac\n0 0 0 c\n$end",
    "$coord file=coord\n0 0 0 c\n$end", "$coord\n0 0 0 c\n$periodic 2\n$end",
    "$coord\n0 0 0 c\n$periodic 3\n$end", "$coord\n0 0 0 c\n$eht charge=1\n$end",
    "$coord\n0 0 0 c x\n$end", "$coord\nNaN 0 0 c\n$end", "$coord\n0 0 0 13c\n$end",
    "$coord\n0 0 0 c\n$end\nextra", "$coord\n0 0 0 🧪\n$end",
    TEXT.replace("$periodic 3\n", ""), TEXT.replace("1 0 -6", "4 1 0"),
    TEXT.replace("$end", "$cell angs\n4 5 6 90 90 90\n$end"),
])
def test_invalid_document_fails(text):
    with pytest.raises(ValueError):
        chematic.parse_turbomole(text)

@pytest.mark.parametrize("kwargs", [{"max_input_bytes": 2}, {"max_atoms": 1}, {"max_line_len": 5}])
def test_limits(kwargs):
    with pytest.raises(ValueError):
        chematic.parse_turbomole(TEXT, **kwargs)

@pytest.mark.parametrize("elements,positions,kwargs", [
    (["c"], [[0,0,0]], {}), (["13C"], [[0,0,0]], {}),
    (["C"], [], {}), (["C"], [[0,0,0]], {"frozen": []}),
    (["C"], [[float("nan"),0,0]], {}), ([], [], {}),
    (["C"], [[0,0,0]], {"lattice_angstrom": [[0,0,0],[0,1,0],[0,0,1]]}),
    (["C"], [[float.fromhex("0x1.fffffffffffffp+1023"),0,0]], {}),
])
def test_writer_invalid_values(elements, positions, kwargs):
    with pytest.raises(ValueError):
        chematic.write_turbomole(elements, positions, **kwargs)

def test_nonperiodic_bohr_no_inferred_hydrogens():
    g = chematic.parse_turbomole("$coord\n1D0 0 0 n\n$end\n")
    assert g["positions_angstrom"][0][0] == pytest.approx(0.529177210544, abs=1e-14)
    assert g["elements"] == ["N"]
    assert g["lattice_angstrom"] is None
    assert g["periodicity"] == 0

@pytest.mark.parametrize("coord_unit,cell_unit", [("angs", "bohr"), ("bohr", "angs")])
def test_coordinate_and_lattice_units_are_independent(coord_unit, cell_unit):
    text = f"$periodic 3\n$lattice {cell_unit}\n2 0 0\n0 3 0\n0 0 4\n$coord {coord_unit}\n1 2 3 c\n$end\n"
    g = chematic.parse_turbomole(text)
    coord_factor = 1.0 if coord_unit == "angs" else 0.529177210544
    cell_factor = 1.0 if cell_unit == "angs" else 0.529177210544
    assert g["positions_angstrom"][0] == pytest.approx([coord_factor, 2*coord_factor, 3*coord_factor])
    assert g["lattice_angstrom"][1][1] == pytest.approx(3*cell_factor)
