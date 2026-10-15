"""Original MIT OR Apache-2.0 EXYZ contract cases; no third-party data."""
import itertools
import math
import pytest
import chematic

TEXT = "2\n title  %VIRTUAL %PBC \nSi 0.12345678901234567 -2 3\nCl -1 0 2 VIRTUAL\n\nVector1 4 1 0\nVector2 -1 5 2\nVector3 1 0 -6\nOffset 0.5 -1 2\n"

def emit(g):
    return chematic.write_exyz(g["elements"], g["positions_angstrom"], header=g["header"],
        lattice_angstrom=g["lattice_angstrom"], offset_angstrom=g["offset_angstrom"], virtual_atoms=g["virtual_atoms"])

def test_structured_roundtrip_and_source_metadata():
    import hashlib
    import json
    from pathlib import Path
    root = Path(__file__).resolve().parents[3]
    fixture = json.loads((root / "validation/openbabel_issues/exyz/contract-fixtures.json").read_text(encoding="utf-8"))["records"][0]
    assert fixture["payload"] == TEXT
    assert fixture["sha256"] == hashlib.sha256(TEXT.encode()).hexdigest()
    g = chematic.parse_exyz(TEXT)
    for key, value in fixture["expected"].items():
        if key not in {"units", "tolerance"}:
            assert g[key] == value
    assert g["source_text"] == TEXT
    assert g["header"] == " title  %VIRTUAL %PBC "
    assert g["elements"] == ["Si", "Cl"]
    assert g["positions_angstrom"] == [[0.12345678901234567,-2,3],[-1,0,2]]
    assert g["virtual_atoms"] == [False, True]
    assert g["lattice_angstrom"] == [[4,1,0],[-1,5,2],[1,0,-6]]
    assert g["offset_angstrom"] == [0.5,-1,2]
    back = chematic.parse_exyz(emit(g))
    assert {k:v for k,v in back.items() if k != "source_text"} == {k:v for k,v in g.items() if k != "source_text"}
    g["positions_angstrom"][0][0] = 0.25
    assert chematic.parse_exyz(emit(g))["positions_angstrom"][0][0] == 0.25

@pytest.mark.parametrize("periodic,virtual", list(itertools.product([False, True], repeat=2)))
def test_writer_adds_required_flags(periodic, virtual):
    out = chematic.write_exyz(["C"], [[1,2,3]], header="sample",
        lattice_angstrom=[[2,0,0],[0,3,0],[0,0,-4]] if periodic else None, virtual_atoms=[virtual])
    g = chematic.parse_exyz(out)
    assert ("%PBC" in g["header"]) == periodic
    assert ("%VIRTUAL" in g["header"]) == virtual
    assert g["virtual_atoms"] == [virtual]
    assert g["offset_angstrom"] == ([0,0,0] if periodic else None)

@pytest.mark.parametrize("values", [[-0.0,5e-324,-5e-324], [1e200,-1e-200,0.12345678901234567]])
def test_full_double_precision(values):
    g = chematic.parse_exyz(chematic.write_exyz(["Og"], [values]))
    assert g["positions_angstrom"][0] == values
    assert math.copysign(1, g["positions_angstrom"][0][0]) == math.copysign(1, values[0])

@pytest.mark.parametrize("text", [
    "", "0\n\n", "2\n\nC 0 0 0\n", "99999999999999999999999\n\n",
    "1\n\nC 0 0 0 VIRTUAL", "1\n%VIRTUAL\nC 0 0 0 EXTRA",
    "1\nProperties=species:S:1:pos:R:3\nC 0 0 0", "1\n%NO\nC 0 0 0",
    "1\n%VIRTUAL %VIRTUAL\nC 0 0 0", "1\n\n13C 0 0 0", "1\n\nc 0 0 0",
    "1\n\nC NaN 0 0", "1\n\nC 1e999 0 0", "1\n\nC 0 0 0\n1\n\nC 0 0 0",
    TEXT.replace("Offset 0.5 -1 2\n", ""), TEXT.replace("Vector2", "Vector1"),
    TEXT.replace("Vector3 1 0 -6", "Vector3 4 1 0"), TEXT.replace("\n\nVector1", "\nVector1"),
    TEXT.replace("%PBC", ""), TEXT.replace("Offset 0.5 -1 2", "Offset inf -1 2"),
])
def test_reader_rejects_invalid_data(text):
    with pytest.raises(ValueError):
        chematic.parse_exyz(text)

@pytest.mark.parametrize("kwargs", [{"max_atoms":1}, {"max_input_bytes":2}, {"max_line_len":5}])
def test_limits(kwargs):
    with pytest.raises(ValueError):
        chematic.parse_exyz(TEXT, **kwargs)

@pytest.mark.parametrize("elements,positions,kwargs", [
    ([],[],{}), (["C"],[],{}), (["c"],[[0,0,0]],{}), (["D"],[[0,0,0]],{}),
    (["C"],[[float("nan"),0,0]],{}), (["C"],[[0,0,0]],{"virtual_atoms":[]}),
    (["C"],[[0,0,0]],{"offset_angstrom":[1,2,3]}),
    (["C"],[[0,0,0]],{"header":"bad\nheader"}),
    (["C"],[[0,0,0]],{"header":"%PBC"}),
    (["C"],[[0,0,0]],{"lattice_angstrom":[[1,0,0],[0,1,0],[1,1,0]]}),
])
def test_writer_rejects_invalid_data(elements,positions,kwargs):
    with pytest.raises(ValueError):
        chematic.write_exyz(elements, positions, **kwargs)

def test_plain_xyz_subset_and_separate_ase_api():
    assert chematic.parse_exyz("1\r\n\r\nH 0 0 0")["elements"] == ["H"]
    text = '1\nLattice="1 0 0 0 1 0 0 0 1" Properties=species:S:1:pos:R:3\nC 0 0 0\n'
    with pytest.raises(ValueError, match="ASE ExtXYZ"):
        chematic.parse_exyz(text)
    assert chematic.from_extxyz(text)["coords"] == [[0.0,0.0,0.0]]

def test_shipped_type_stubs_declare_public_functions():
    import ast
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "python/chematic/__init__.pyi"
    names = {n.name for n in ast.parse(path.read_text(encoding="utf-8")).body if isinstance(n, ast.FunctionDef)}
    assert {"parse_exyz", "write_exyz"} <= names

def test_tabs_in_header_and_atom_rows_preserve_header():
    text = "1\nlabel\t%VIRTUAL\nC\t1\t2\t3\tVIRTUAL\n"
    g = chematic.parse_exyz(text)
    assert g["header"] == "label\t%VIRTUAL"
    assert chematic.parse_exyz(emit(g))["virtual_atoms"] == [True]
