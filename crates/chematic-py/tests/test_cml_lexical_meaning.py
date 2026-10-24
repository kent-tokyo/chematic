"""Formatting and opaque XML must not change the authored graph/H projection."""
import json
from pathlib import Path

import chematic
import pytest

PACKET = json.loads((Path(__file__).parents[3] /
                     "validation/cml_lexical_meaning/cases.json").read_text(encoding="utf-8"))


def fields(mol):
    value = json.loads(mol.to_moljson())
    atoms = value["atoms"]
    ids = {a["id"]: i for i, a in enumerate(atoms)}
    return dict(numbers=[{"C": 6, "N": 7, "O": 8}[a["element"]] for a in atoms],
                charges=[a.get("charge", 0) for a in atoms],
                isotopes=[a.get("isotope", 0) for a in atoms],
                hydrogens=[a["hydrogens"] for a in atoms],
                bonds=sorted([min(ids[b["source_id"]], ids[b["target_id"]]),
                              max(ids[b["source_id"]], ids[b["target_id"]]), b["order"]]
                             for b in value["bonds"]))


@pytest.mark.parametrize("case", PACKET["cases"], ids=lambda c: c["id"])
@pytest.mark.parametrize("encoded", [False, True])
def test_public_readers_and_common_conversion_preserve_graph(case, encoded):
    source = case["input"].encode() if encoded else case["input"]
    expected = {k: v for k, v in case["expected"].items() if k not in ("coords2d", "atom_maps")}
    for reader in (chematic.from_cml, chematic.from_cml_strict):
        assert fields(reader(source)) == expected
    mol, xyz = chematic.from_cml_with_3d(source)
    assert fields(mol) == expected
    assert xyz is None
    # MolJSON retains input indexing, unlike canonical SMILES.
    converted = chematic.convert_format(case["input"], "cml", "moljson")
    assert fields(chematic.from_moljson(converted)) == expected
    report = chematic.convert_with_report(case["input"], "cml", "cml", policy="preserve_source")
    assert report["output"] == case["input"]
    assert report["candidate_graph_verified"] is False


def test_opaque_markup_cannot_make_an_empty_molecule_nonempty():
    for source in ("<molecule><!--\n<atom elementType='C'/>\n--></molecule>",
                   "<molecule><![CDATA[\n<atom elementType='C'/>\n]]></molecule>"):
        assert chematic.from_cml(source).atom_count == 0
        with pytest.raises(ValueError):
            chematic.from_cml_strict(source)
