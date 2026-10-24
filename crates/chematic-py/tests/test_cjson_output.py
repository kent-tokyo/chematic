"""CJSON graph projection must refuse unsupported in-memory meaning."""
import json
from pathlib import Path
import struct

import pytest
import chematic

PACKET = json.loads((Path(__file__).resolve().parents[3] / "validation/cjson_output/cases.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("case", PACKET["cases"], ids=lambda c: c["id"])
def test_cml_common_bridge_checks_bond_representability(case):
    for source in [case["input"], case["input"].replace("><", ">\n<")]:
        mol, coords = chematic.from_cml_with_3d(source)
        if case["accepted"]:
            direct = json.loads(mol.to_cjson_checked(coords))
            assert json.loads(chematic.convert_format(source, "cml", "cjson")) == direct
            assert direct["bonds"]["order"] == [case["expected"]["bonds"][0][2]]
        else:
            for call in [lambda: mol.to_cjson_checked(coords),
                         lambda: chematic.convert_format(source, "cml", "cjson")]:
                with pytest.raises(ValueError, match="cannot preserve bond order " + case["error_order"]):
                    call()
            # The compatibility writer remains explicitly lossy.
            assert json.loads(mol.to_cjson(coords))["bonds"]["order"] == [1]


@pytest.mark.parametrize("source", ["F/C=C/F", "N[C@@](F)(Cl)Br", "[C]", "[CH4]",
                                   "c1ccccc1", "[CH4:7]", "*C", "[Pt@SP1](Cl)(Br)(I)F"])
def test_direct_and_common_smiles_reject_unrepresented_properties(source):
    mol = chematic.from_smiles(source)
    with pytest.raises(ValueError, match="cannot preserve"):
        mol.to_cjson_checked()
    with pytest.raises(ValueError, match="cannot preserve"):
        chematic.convert_format(source, "smiles", "cjson")
    assert json.loads(mol.to_cjson())["chemicalJson"] == 1


@pytest.mark.parametrize("rows", [[[1, 2]], [[1, 2, 3, 4]], [[0, 0, 0], [0, 0, 0]],
                                [[float("nan"), 0, 0]], [[0, float("inf"), 0]], [[0, 0, -float("inf")]]])
def test_coordinates_are_never_padded_truncated_omitted_or_null(rows):
    mol = chematic.from_smiles("C")
    with pytest.raises(ValueError, match="CJSON:"):
        mol.to_cjson_checked(rows)
    with pytest.raises(ValueError, match="CJSON:"):
        chematic.convert_format("C", "smiles", "cjson", coords=rows)


def test_supported_geometry_and_absence_are_distinct():
    mol = chematic.from_smiles("C")
    coords = [[-0.0, 5e-324, -1.4175225513503917]]
    data = json.loads(mol.to_cjson_checked(coords))
    assert [struct.pack(">d", x) for x in data["atoms"]["coords"]["3d"]] == [struct.pack(">d", x) for x in coords[0]]
    assert "coords" not in json.loads(mol.to_cjson_checked())["atoms"]
