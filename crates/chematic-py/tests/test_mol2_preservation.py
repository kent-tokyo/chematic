"""Shared source fixture with Rust, CLI, and Node/WASM preservation tests."""

import json
from pathlib import Path

import chematic
import pytest


SOURCE = (Path(__file__).parents[3] / "validation/mol2_preservation/amide.mol2").read_text()


def test_mol2_record_metadata_and_source():
    record = json.loads(chematic.parse_mol2_record_json(SOURCE))
    assert record["schema"] == "chematic.mol2-record.v1"
    assert record["name"] == "acetamide"
    assert record["charge_type"] == "USER_CHARGES"
    assert record["atoms"][2]["source_id"] == 90
    assert record["atoms"][2]["atom_type"] == "O.2"
    assert record["atoms"][2]["partial_charge"] == -0.70
    assert record["formal_charges"] == [0, 0, 0, 0]
    assert record["bonds"][2]["bond_type"] == "am"
    assert record["coordinates"][0] == [0.123456789, 0.0, 0.0]
    assert record["opaque_sections"] == ["SUBSTRUCTURE", "VENDOR_DATA"]
    assert record["graph_projection"]["opaque_sections"] == "not_interpreted"
    assert record["source"] == SOURCE


@pytest.mark.parametrize("source", [SOURCE, SOURCE.replace("\n", "\r\n")])
def test_mol2_self_conversion_retains_source(source):
    assert chematic.roundtrip_mol2_record(source) == source
    assert chematic.convert_format(source, ".mol2", "mol2") == source


@pytest.mark.parametrize("bad", [
    SOURCE.replace("ACM -0.70", "ACM NaN"),
    SOURCE.replace("12 7 3 am", "12 7 3 vendor"),
    SOURCE.replace("4 3 1 0 0", "4 4 1 0 0"),
    SOURCE + SOURCE,
])
def test_invalid_records_are_rejected_across_entry_points(bad):
    for reader in (chematic.parse_mol2_record_json, chematic.roundtrip_mol2_record):
        with pytest.raises(ValueError):
            reader(bad)
    with pytest.raises(ValueError):
        chematic.convert_format(bad, "mol2", "mol2")


@pytest.mark.parametrize("options", [
    {"coords": [[0, 0, 0]] * 4}, {"charges": [0.0] * 4},
    {"name": "override"}, {"comment": "override"},
])
def test_source_preserving_conversion_does_not_silently_ignore_overrides(options):
    with pytest.raises(ValueError, match="does not accept overrides"):
        chematic.convert_format(SOURCE, "mol2", "mol2", **options)
