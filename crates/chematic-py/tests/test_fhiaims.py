"""Original MIT OR Apache-2.0 FHI-aims geometry contracts."""
import ast
from pathlib import Path
import pytest
import chematic

TEXT = "# original\nlattice_vector 4 1 0\nconstrain_relaxation x\nlattice_vector -1 5 2\nlattice_vector 1 0 -6\natom_frac 0.5 0.25 -0.5 Si\nconstrain_relaxation y\ninitial_charge 0.25\ninitial_moment -1.5\natom 1 2 3 Cl\n"
FIELDS=("coordinate_kinds","lattice_angstrom","constraints","initial_charges","initial_moments","lattice_constraints")
def emit(g):
    return chematic.write_fhiaims(g["elements"],g["coordinates"],**{k:g[k] for k in FIELDS})

def test_fractional_geometry_and_metadata_roundtrip():
    import hashlib
    import json
    p=Path(__file__).resolve().parents[3]/"validation/openbabel_issues/fhiaims/contract-fixtures.json"
    fixture=json.loads(p.read_text(encoding="utf-8"))["records"][0]
    assert fixture["payload"]==TEXT
    assert fixture["sha256"]==hashlib.sha256(TEXT.encode()).hexdigest()
    g=chematic.parse_fhiaims(TEXT)
    for key,value in fixture["expected"].items():
        if key not in {"units","tolerance"}:assert g[key]==value
    assert g["positions_angstrom"]==[[1.25,1.75,3.5],[1,2,3]]
    assert g["coordinates"][0]==[0.5,0.25,-0.5]
    assert g["coordinate_kinds"]==["fractional","cartesian"]
    assert g["initial_charges"]==[0.25,None]
    assert g["initial_moments"]==[-1.5,None]
    assert g["constraints"]==[[False,True,False],[False,False,False]]
    assert g["lattice_constraints"][0]==[True,False,False]
    assert g["source_text"]==TEXT
    back=chematic.parse_fhiaims(emit(g))
    for key in g.keys()-{"source_text"}:assert back[key]==g[key]
    g["coordinates"][0]=[0,0,0]
    assert chematic.parse_fhiaims(emit(g))["positions_angstrom"][0]==[0,0,0]

@pytest.mark.parametrize("dimension",range(4))
def test_partial_periodicity(dimension):
    lattice=[[2,1,0],[0,3,1],[0,0,-4]][:dimension]
    g=chematic.parse_fhiaims(chematic.write_fhiaims(["C"],[[1,2,3]],lattice_angstrom=lattice))
    assert g["periodicity"]==dimension and g["lattice_angstrom"]==lattice

@pytest.mark.parametrize("mask",range(8))
def test_all_atom_and_lattice_constraint_masks(mask):
    row=[bool(mask&1),bool(mask&2),bool(mask&4)]
    g=chematic.parse_fhiaims(chematic.write_fhiaims(["Fe"],[[0,0,0]],constraints=[row],lattice_angstrom=[[1,0,0]],lattice_constraints=[row],initial_charges=[0.0],initial_moments=[None]))
    assert g["constraints"]==[row] and g["lattice_constraints"]==[row]
    assert g["initial_charges"]==[0.0] and g["initial_moments"]==[None]

@pytest.mark.parametrize("text",["", "atom 0 0 C", "atom 0 0 0 Fe_custom", "atom_frac 0 0 0 C", "atom 0 0 0 C\nvelocity 1 2 3", "initial_charge 1\natom 0 0 0 C", "atom NaN 0 0 C", "atom 0 0 0 C\ninitial_moment inf", "atom 0 0 0 C\ninitial_charge 0\ninitial_charge 1", "atom 0 0 0 C\nconstrain_relaxation x y", "lattice_vector 0 0 0\natom 0 0 0 C"])
def test_invalid_reader_inputs(text):
    with pytest.raises(ValueError):chematic.parse_fhiaims(text)

@pytest.mark.parametrize("kwargs",[{"max_atoms":1},{"max_input_bytes":2},{"max_line_len":5}])
def test_reader_limits(kwargs):
    with pytest.raises(ValueError):chematic.parse_fhiaims(TEXT,**kwargs)

@pytest.mark.parametrize("elements,coordinates,kwargs",[
    ([],[],{}),(["C"],[],{}),(["c"],[[0,0,0]],{}),(["C"],[[float("inf"),0,0]],{}),
    (["C"],[[0,0,0]],{"coordinate_kinds":["bad"]}),(["C"],[[0,0,0]],{"coordinate_kinds":["fractional"]}),
    (["C"],[[0,0,0]],{"constraints":[]}),(["C"],[[0,0,0]],{"initial_charges":[]}),
    (["C"],[[0,0,0]],{"initial_moments":[float("nan")]}),
    (["C"],[[0,0,0]],{"lattice_angstrom":[[1,0,0]],"lattice_constraints":[]}),
])
def test_invalid_writer_inputs(elements,coordinates,kwargs):
    with pytest.raises(ValueError):chematic.write_fhiaims(elements,coordinates,**kwargs)

def test_fortran_exponents_crlf_comments_and_stub_declarations():
    g=chematic.parse_fhiaims("atom\t1D0 0 0 C # original\r\ninitial_charge 0D0\r\n")
    assert g["coordinates"]==[[1,0,0]] and g["initial_charges"]==[0.0]
    p=Path(__file__).resolve().parents[1]/"python/chematic/__init__.pyi"
    names={n.name for n in ast.parse(p.read_text(encoding="utf-8")).body if isinstance(n,ast.FunctionDef)}
    assert {"parse_fhiaims","write_fhiaims"}<=names
