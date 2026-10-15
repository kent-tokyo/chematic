import json
import pytest
import chematic


def test_default_carbon_and_explicit_zero_hydrogens():
    expected = {None: "CH4", -1: "CH4", 0: "C", 2: "CH2", 4: "CH4"}
    for hydrogen, formula in expected.items():
        atom = dict(x=0, y=0)
        if hydrogen is not None: atom["h"] = hydrogen
        record = chematic.parse_cdjson(json.dumps(dict(a=[atom])))
        assert record["molecules"][0]["mol"].formula == formula


def test_all_molecules_ids_coordinates_and_exact_source():
    text = '{\r\n"m":[{"i":"mol","a":[{"i":"atom","l":"N","c":1,"m":15,"h":4,"x":1.2345678901234567,"y":-2,"z":3}]},{"a":[]}]}\r\n'
    d = chematic.parse_cdjson(text)
    assert d["source_text"] == text
    assert len(d["molecules"]) == 2
    m = d["molecules"][0]
    assert m["coordinates"] == [[1.2345678901234567, -2, 3]]
    assert m["id"] == "mol" and m["atom_ids"] == ["atom"]
    assert m["mol"].atom_count == 1


@pytest.mark.parametrize("field", [{"q": {}}, {"s2": {}}, {"r": 1}, {"p": 1}, {"l": "8"}, {"unknown": 1}])
def test_unsupported_atom_fields(field):
    with pytest.raises(ValueError, match="unsupported"):
        chematic.parse_cdjson(json.dumps(dict(a=[dict(x=0, y=0, **field)])))


@pytest.mark.parametrize("text", [
    '{"a":[{"x":0,"x":1,"y":0}]}', '{"a":[{"x":0,"\\u0078":1,"y":0}]}',
    '{"a":[{"x":0}]}', '{"m":[],"s":[{"t":"Line"}]}',
    '{"a":[{"x":0,"y":0,"h":null}]}',
    '{"a":[{"x":0,"y":0},{"x":1,"y":0}],"b":[{"b":0,"e":2}]}',
])
def test_malformed_or_unsupported_documents(text):
    with pytest.raises(ValueError): chematic.parse_cdjson(text)


@pytest.mark.parametrize("kwargs", [dict(max_input_bytes=1), dict(max_atoms=0), dict(max_molecules=0)])
def test_tighter_limits(kwargs):
    with pytest.raises(ValueError, match="resource limit"):
        chematic.parse_cdjson('{"a":[{"x":0,"y":0}]}', **kwargs)


@pytest.mark.parametrize("smiles", ["C", "CC", "C(=O)CO", "[C]", "[CH2]", "[NH4+]", "[13CH3]O", "[2H]O[2H]"])
def test_graph_writer_roundtrips_formula_and_coordinates(smiles):
    mol = chematic.from_smiles(smiles)
    coords = [[i*1.2345678901234567, -2, 3] for i in range(mol.atom_count)]
    text = chematic.write_cdjson_graph(mol, coords)
    m = chematic.parse_cdjson(text)["molecules"][0]
    assert m["mol"].formula == mol.formula
    assert m["coordinates"] == coords
    assert all(isinstance(a["l"], str) for a in json.loads(text)["m"][0]["a"])


@pytest.mark.parametrize("smiles", ["c1ccccc1", "F/C=C\\F", "[C@H](F)(Cl)Br", "[CH3:7]"])
def test_writer_refuses_unsupported_chemistry(smiles):
    mol = chematic.from_smiles(smiles)
    with pytest.raises(ValueError):
        chematic.write_cdjson_graph(mol, [[0, 0, 0]]*mol.atom_count)


def test_writer_checks_coordinates_and_limits_are_unsigned():
    mol = chematic.from_smiles("C")
    for coords in [[], [[float("inf"), 0, 0]]]:
        with pytest.raises(ValueError): chematic.write_cdjson_graph(mol, coords)
    with pytest.raises(OverflowError): chematic.parse_cdjson("{}", max_atoms=-1)


def test_coordinate_length_checked_before_numeric_extraction():
    mol = chematic.from_smiles("C")
    with pytest.raises(ValueError, match="coordinate count"):
        chematic.write_cdjson_graph(mol, [["bad", 0, 0], ["bad", 0, 0]])
    with pytest.raises(TypeError):
        chematic.write_cdjson_graph(mol, ((0, 0, 0),))
