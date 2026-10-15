"""Single-component CCD is distinct from atom_site mmCIF."""
import json
from pathlib import Path

import chematic
import pytest

WATER = (Path(__file__).resolve().parents[3] / "validation/openbabel_issues/ccd/water.cif").read_text(encoding="utf-8")


def test_record_and_explicit_graph_keep_template_fields():
    record = json.loads(chematic.parse_ccd_json(WATER))
    assert record["component_id"] == "SYN"
    assert (len(record["atoms"]), len(record["bonds"])) == (3, 2)
    atom = record["atoms"][0]
    assert atom["atom_id"] == "O#1"
    assert atom["model_coordinates"] == [1, 2, 3]
    assert atom["ideal_coordinates"] == [0, 0, 0]
    assert atom["extra_fields"]["_chem_comp_atom.alt_atom_id"] == "O1"
    assert chematic.from_ccd(WATER).formula == "H2O"
    with pytest.raises(ValueError):
        chematic.parse_mmcif(WATER)


@pytest.mark.parametrize("option,value", [("max_input_bytes", 1), ("max_atoms", 2), ("max_bonds", 1)])
def test_typed_resource_errors(option, value):
    with pytest.raises(ValueError) as error:
        chematic.parse_ccd_json(WATER, **{option: value})
    assert error.value.category == "resource_limit"
    assert error.value.format == "ccd"


def test_stereo_is_recorded_and_graph_requires_loss_consent():
    text = "data_T\n_chem_comp.id T\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_stereo_config\nT C C 0 R\nT F F 0 N\nT CL Cl 0 N\nT BR Br 0 N\nT H H 0 N\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nT C F SING\nT C CL SING\nT C BR SING\nT C H SING\n"
    assert json.loads(chematic.parse_ccd_json(text))["atoms"][0]["stereo_config"] == "R"
    with pytest.raises(ValueError) as error:
        chematic.from_ccd(text)
    assert error.value.code == "ccd_stereo_loss"
    assert chematic.from_ccd(text, allow_stereo_loss=True).formula == "CHBrClF"


@pytest.mark.parametrize("text", ["", WATER + "\ndata_OTHER\n", WATER.replace("SING N N", "DELO N N"), WATER.replace("'O#1' H2 SING", "BAD H2 SING"), WATER.replace("N N 1 2 3", "N N NaN 2 3")])
def test_invalid_inputs_have_ccd_diagnostics(text):
    with pytest.raises(ValueError) as error:
        chematic.parse_ccd_json(text)
    assert error.value.category == "malformed"
    assert error.value.format == "ccd"


def test_crlf_and_null_coordinate_family():
    assert json.loads(chematic.parse_ccd_json(WATER.replace("\n", "\r\n"))) == json.loads(chematic.parse_ccd_json(WATER))
    record = json.loads(chematic.parse_ccd_json(WATER.replace("N N 1 2 3", "N N ? ? ?")))
    assert record["atoms"][0]["model_coordinates"] is None
    assert record["atoms"][0]["ideal_coordinates"] == [0, 0, 0]


def test_unknown_charge_is_retained_and_projection_refuses_imputation():
    text = WATER.replace("'O#1' O 0", "'O#1' O ?")
    assert json.loads(chematic.parse_ccd_json(text))["atoms"][0]["charge"] is None
    with pytest.raises(ValueError) as error:
        chematic.from_ccd(text, allow_stereo_loss=True)
    assert error.value.code == "ccd_unknown_charge"


def test_long_repeated_tags_hit_the_expansion_budget():
    tag = "_chem_comp_atom.extra_" + "x" * 20_000
    text = "data_A\n_chem_comp.id A\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n" + tag + "\n"
    text += "".join(f"A C{i} C v\n" for i in range(2000))
    with pytest.raises(ValueError, match="expanded loop bytes") as error:
        chematic.parse_ccd_json(text)
    assert error.value.category == "resource_limit"
