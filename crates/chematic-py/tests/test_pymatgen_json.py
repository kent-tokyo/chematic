"""Source and explicit-loss contracts for Open Babel #2675."""

from copy import deepcopy
from dataclasses import FrozenInstanceError
import json

import numpy as np
import pytest
from chematic import PeriodicStructure, parse_pymatgen_structures


def structure():
    return {
        "@module": "pymatgen.core.structure",
        "@class": "Structure",
        "lattice": {
            "matrix": [[4.0, 0.5, 0.0], [0.0, 5.0, 0.75], [0.25, 0.0, 6.0]],
            "pbc": [True, True, True],
        },
        "sites": [
            {
                "species": [
                    {"element": "Fe", "occu": 0.6},
                    {"element": "Ni", "occu": 0.3},
                ],
                "abc": [-0.25, 0.5, 1.25],
                "label": "metal α",
                "properties": {},
            },
            {
                "species": [{"element": "O", "occu": 1.0}],
                "abc": [0.25, 0.0, 0.125],
                "label": "oxygen",
            },
        ],
        "properties": {},
    }


def parse(data):
    return parse_pymatgen_structures(json.dumps(data))[0]


@pytest.mark.parametrize("encoding", ["str", "bytes"])
def test_source_order_disorder_labels_unwrapped_coordinates(encoding):
    data = structure()
    text = json.dumps(data, indent=2, ensure_ascii=False) + "\n"
    value = text if encoding == "str" else text.encode()
    record = parse_pymatgen_structures(value)[0]
    assert record.source_text == text and record.record_index == 0
    assert record.structure_dict() == data
    assert record.projection_losses == ()
    native = record.to_periodic_structure()
    assert isinstance(native, PeriodicStructure)
    assert native.sites[0].species == [("Fe", 0.6), ("Ni", 0.3)]
    assert native.sites[0].label == "metal α"
    np.testing.assert_array_equal(
        native.fractional_positions(), [s["abc"] for s in data["sites"]]
    )
    np.testing.assert_allclose(
        native.cartesian_positions(),
        np.array([s["abc"] for s in data["sites"]]) @ data["lattice"]["matrix"],
        rtol=0,
        atol=1e-15,
    )
    copied = record.structure_dict()
    copied["sites"][0]["label"] = "changed"
    assert record.structure_dict() == data
    with pytest.raises(FrozenInstanceError):
        record.record_index = 7


def test_complete_metadata_retention_and_explicit_projection_loss():
    data = structure()
    data["charge"] = -2.0
    data["properties"] = {"temperature": 298.15, "nested": {"value": [1, None, "μ"]}}
    data["custom"] = {"@module": "never.import.this", "@class": "Arbitrary"}
    data["sites"][0]["species"][0].update(oxidation_state=2.0, spin=3.0)
    data["sites"][0]["properties"] = {"magmom": [0.0, 1.0, 0.0]}
    data["sites"][1]["custom_id"] = "stable-id"
    data["lattice"]["custom"] = [1, 2]
    record = parse(data)
    assert record.structure_dict() == data
    assert set(record.projection_losses) == {
        "structure.charge",
        "structure.properties",
        "structure.custom",
        "lattice.custom",
        "sites[0].properties",
        "sites[0].species[0].oxidation_state",
        "sites[0].species[0].spin",
        "sites[1].custom_id",
    }
    with pytest.raises(ValueError, match="loses metadata"):
        record.to_periodic_structure()
    assert record.to_periodic_structure(allow_metadata_loss=True).site_count() == 2
    with pytest.raises(TypeError):
        record.to_periodic_structure(allow_metadata_loss=1)


def test_response_records_and_envelope_metadata_are_retained():
    a = structure()
    b = structure()
    b["sites"].reverse()
    doc = {
        "data": [
            {"structure": a, "material_id": "synthetic-a"},
            {"structure": b, "material_id": "synthetic-b"},
        ],
        "meta": {"total_doc": 2},
    }
    source = json.dumps(doc)
    records = parse_pymatgen_structures(source)
    assert len(records) == 2
    for i, record in enumerate(records):
        assert record.source_text == source and record.record_index == i
        assert record.structure_dict() == doc["data"][i]["structure"]
        assert set(record.projection_losses) == {
            "response",
            f"data[{i}].material_id",
        }
        with pytest.raises(ValueError):
            record.to_periodic_structure()
    assert (
        records[1].to_periodic_structure(allow_metadata_loss=True).sites[0].label
        == "oxygen"
    )


@pytest.mark.parametrize("class_name", [None, "Structure", "IStructure"])
def test_optional_class_pbc_and_coordinate_fields(class_name):
    data = structure()
    if class_name is None:
        data.pop("@module")
        data.pop("@class")
    else:
        data["@class"] = class_name
    data["lattice"].pop("pbc")
    data.pop("properties")
    record = parse(data)
    assert record.to_periodic_structure().site_count() == 2


def test_same_element_distinct_oxidation_entries_are_not_merged():
    data = structure()
    data["sites"][0]["species"] = [
        {"element": "Fe", "occu": 0.4, "oxidation_state": 2},
        {"element": "Fe", "occu": 0.6, "oxidation_state": 3},
    ]
    record = parse(data)
    assert record.structure_dict()["sites"][0]["species"] == data["sites"][0]["species"]
    assert record.to_periodic_structure(allow_metadata_loss=True).sites[0].species == [
        ("Fe", 0.4),
        ("Fe", 0.6),
    ]


def test_zero_occupancy_entries_are_retained():
    data = structure()
    data["sites"][0]["species"].append({"element": "Zn", "occu": 0.0})
    assert parse(data).to_periodic_structure().sites[0].species[-1] == ("Zn", 0.0)


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: d.update({"@class": "Molecule"}),
        lambda d: d.update({"@module": "os"}),
        lambda d: d.update({"lattice": []}),
        lambda d: d["lattice"].update(matrix=[[1, 0, 0]]),
        lambda d: d["lattice"].update(matrix=[[0.0, 0.0, 0.0]] * 3),
        lambda d: d["lattice"].update(
            matrix=[[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
        ),
        lambda d: d["lattice"].update(
            matrix=[[True, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        ),
        lambda d: d["lattice"].update(pbc=[True, True, False]),
        lambda d: d["lattice"].update(pbc=[1, 1, 1]),
        lambda d: d["lattice"].update(pbc=[True, True]),
        lambda d: d["lattice"].update(a=99.0),
        lambda d: d["lattice"].update(alpha=12.0),
        lambda d: d["lattice"].update(volume=9.0),
        lambda d: d["lattice"].update(a="4"),
        lambda d: d.update(sites=[]),
        lambda d: d.update(sites={}),
        lambda d: d["sites"].__setitem__(0, None),
        lambda d: d["sites"][0].update(abc=[1.0, 2.0]),
        lambda d: d["sites"][0].update(abc=["1", 2.0, 3.0]),
        lambda d: d["sites"][0].update(abc=[True, 2.0, 3.0]),
        lambda d: d["sites"][0].update(abc=[1e101, 2.0, 3.0]),
        lambda d: d["sites"][0].update(xyz=[1.0, 2.0, 3.0]),
        lambda d: d["sites"][0].update(label=7),
        lambda d: d["sites"][0].update(label="\ud800"),
        lambda d: d["sites"][0].update(properties=[]),
        lambda d: d["sites"][0].update(species=[]),
        lambda d: d["sites"][0].update(species={}),
        lambda d: d["sites"][0].update(species=[None]),
        lambda d: d["sites"][0]["species"][0].update(element="Xx"),
        lambda d: d["sites"][0]["species"][0].update(element="\ud800"),
        lambda d: d["sites"][0]["species"][0].update(element=26),
        lambda d: d["sites"][0]["species"][0].update(occu=-0.1),
        lambda d: d["sites"][0]["species"][0].update(occu=1.1),
        lambda d: d["sites"][0]["species"][0].update(occu=True),
        lambda d: d["sites"][0]["species"][0].pop("occu"),
        lambda d: d.update(properties=[]),
        lambda d: d.update(charge=True),
        lambda d: d.update(data=[]),
    ],
)
def test_malformed_structures_raise(mutation):
    data = structure()
    mutation(data)
    with pytest.raises(ValueError):
        parse(data)


@pytest.mark.parametrize(
    "text",
    [
        "",
        "[]",
        "null",
        "1",
        "{",
        "{}",
        '{"data":[]}',
        '{"data":{}}',
        '{"data":[{}]}',
        '{"data":[null]}',
        '{"data":[{"structure":[]}]}',
        '{"a":1,"a":2}',
        '{"a":{"b":1,"b":2}}',
        '{"a":NaN}',
        '{"a":Infinity}',
        '{"a":1e999}',
        '{"a":1e-999}',
        '{"a":-1e-999}',
    ],
)
def test_invalid_json_and_responses(text):
    with pytest.raises(ValueError):
        parse_pymatgen_structures(text)


@pytest.mark.parametrize(
    "value", [None, {}, [], 1, bytearray(b"{}"), memoryview(b"{}")]
)
def test_wrong_input_types(value):
    with pytest.raises(TypeError):
        parse_pymatgen_structures(value)


@pytest.mark.parametrize("value", [b"\xff", "\ud800"])
def test_invalid_utf8(value):
    with pytest.raises(ValueError):
        parse_pymatgen_structures(value)


def test_consistent_cartesian_and_redundant_lattice_fields():
    data = structure()
    native = parse(data).to_periodic_structure()
    for site, xyz in zip(data["sites"], native.cartesian_positions().tolist()):
        site["xyz"] = xyz
    for name, val in zip(("a", "b", "c"), native.lattice.lengths):
        data["lattice"][name] = val
    for name, val in zip(("alpha", "beta", "gamma"), native.lattice.angles_degrees):
        data["lattice"][name] = val
    data["lattice"]["volume"] = native.lattice.volume
    assert parse(data).projection_losses == ()


def test_nesting_input_and_collection_limits(monkeypatch):
    from chematic import _pymatgen_json as module

    with pytest.raises(ValueError, match="nesting"):
        parse_pymatgen_structures("[" * 65 + "0" + "]" * 65)
    with pytest.raises(ValueError, match="16 MiB"):
        parse_pymatgen_structures(" " * (16 * 1024 * 1024 + 1))
    with pytest.raises(ValueError, match="16 MiB"):
        parse_pymatgen_structures("α" * (8 * 1024 * 1024 + 1))
    with pytest.raises(ValueError, match="256"):
        parse_pymatgen_structures(
            json.dumps({"data": [{"structure": structure()}] * 257})
        )
    data = structure()
    data["sites"][0]["species"] *= 33
    with pytest.raises(ValueError, match="1..32"):
        parse(data)
    monkeypatch.setattr(module, "MAX_SITES", 3)
    with pytest.raises(ValueError, match="total site"):
        parse_pymatgen_structures(
            json.dumps({"data": [{"structure": structure()}] * 2})
        )
    monkeypatch.setattr(module, "MAX_SPECIES", 2)
    with pytest.raises(ValueError, match="total site"):
        parse(structure())


def test_later_bad_record_aborts_instead_of_skipping():
    good = structure()
    bad = deepcopy(good)
    bad["sites"][0]["species"][0]["element"] = "Xx"
    with pytest.raises(ValueError):
        parse_pymatgen_structures(
            json.dumps({"data": [{"structure": good}, {"structure": bad}]})
        )


def test_quoted_delimiters_do_not_count_as_nesting_and_zero_underflow_is_zero():
    data = structure()
    data["properties"] = {"quoted": '[[["\\' * 70, "zero": 0.0}
    text = json.dumps(data).replace('"zero": 0.0', '"zero": 0e-999')
    assert (
        parse_pymatgen_structures(text)[0].structure_dict()["properties"]["zero"] == 0.0
    )
