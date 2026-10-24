"""V2000 source geometry is retained instead of being lost at the bridge."""
import json
from pathlib import Path
import struct

import pytest
import chematic

PACKET = json.loads((Path(__file__).parents[3] / "validation/v2000_coordinate_projection/cases.json").read_text())


def bits(values):
    return [struct.pack(">d", float(v)) for row in values for v in row]


@pytest.mark.parametrize("case", PACKET["cases"], ids=lambda c: c["id"])
@pytest.mark.parametrize("variant", PACKET["variants"])
@pytest.mark.parametrize("format", [".MOL", "sdf"])
def test_bridge_retains_all_source_xyz_and_declared_graph_fields(case, variant, format):
    text = case["input"]
    if variant == "crlf": text = text.replace("\n", "\r\n")
    if variant == "no-final-newline": text = text.rstrip("\n")
    mol, name, coords = chematic.from_mol_block_with_coords_3d(text.encode())
    assert name == case["id"]
    expected = case["expected"]
    assert bits(coords) == bits(expected["coordinates"])
    assert len(coords) == mol.atom_count
    legacy, old_name, xy = chematic.from_mol_block_with_coords(text)
    assert old_name == name
    assert legacy.smiles == mol.smiles
    assert bits(xy) == bits([row[:2] for row in coords])
    document = text + ("\n$$$$\n" if format == "sdf" else "")
    projection = case.get("cjson_projection", {"accepted": True})
    if not projection["accepted"]:
        with pytest.raises(ValueError, match=projection["error"]):
            chematic.convert_format(document, format, "cjson")
        return
    value = json.loads(chematic.convert_format(document, format, "cjson"))
    atoms, bonds = value["atoms"], value["bonds"]
    flat = atoms["coords"]["3d"]
    assert bits([flat[i:i+3] for i in range(0,len(flat),3)]) == bits(expected["coordinates"])
    assert atoms["elements"]["number"] == expected["numbers"]
    assert atoms.get("formalCharges", [0]*len(coords)) == expected["formal_charges"]
    assert atoms.get("isotopes", [0]*len(coords)) == expected["isotopes"]
    pairs = bonds["connections"]["index"]
    assert [[pairs[2*i],pairs[2*i+1],order] for i,order in enumerate(bonds["order"])] == expected["bonds"]


@pytest.mark.parametrize("axis", range(3))
@pytest.mark.parametrize("raw", ["", "NaN", "inf", "-inf", "1e999", "not-a-num"])
def test_bad_fields_refuse_instead_of_inventing_zero(axis, raw):
    lines = PACKET["cases"][0]["input"].splitlines()
    row = lines[4]
    lines[4] = row[:axis*10] + raw.rjust(10) + row[axis*10+10:]
    source = "\n".join(lines)
    with pytest.raises(ValueError, match=f"{'xyz'[axis]} coordinate"):
        chematic.from_mol_block_with_coords_3d(source)
    with pytest.raises(ValueError):
        chematic.convert_format(source, "mol", "cjson")


def test_overrides_keep_precedence_and_coordinate_outputs_receive_source_values():
    source = PACKET["cases"][0]["input"]
    override = [[9.125,8.25,7.5],[6,5,4]]
    out = json.loads(chematic.convert_format(source,"mol","cjson",coords=override))
    assert out["atoms"]["coords"]["3d"] == sum(override,[])
    out = json.loads(chematic.convert_format(source,"mol","cjson",coords=[]))
    assert "coords" not in out["atoms"]
    for target in ["pdb","xyz","pdbqt"]:
        assert "3.375" in chematic.convert_format(source,"mol",target)
    with pytest.raises(ValueError,match="records exceeds limit 1"):
        chematic.convert_format(source+"$$$$\n"+source+"$$$$\n","sdf","cjson")


def test_size_and_utf8_limits():
    with pytest.raises(ValueError,match="input bytes"):
        chematic.from_mol_block_with_coords_3d(" "*(16*1024*1024+1))
    with pytest.raises(UnicodeDecodeError):
        chematic.from_mol_block_with_coords_3d(b"\xff")
