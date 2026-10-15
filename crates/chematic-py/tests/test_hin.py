"""Original HIN fixtures: aromatic bond semantics and immutable source retention."""
import json
from pathlib import Path

import chematic
import pytest

SOURCE = (Path(__file__).parents[3] / "validation/openbabel_issues/hin/aromatic-ring.hin").read_text(encoding="utf-8")


def test_openbabel3000_aromatic_hin_to_mdl_is_type4():
    mol = chematic.from_hin(SOURCE)
    block = mol.to_mol_block()
    assert all(row[6:9].strip() == "4" for row in block.splitlines()[10:16])
    assert chematic.from_mol_block(block).to_mol_block().splitlines()[10:16] == block.splitlines()[10:16]


def test_typed_record_keeps_source_ids_partial_charges_and_precision():
    record = json.loads(chematic.parse_hin_record_json(SOURCE))
    assert record["schema"] == "chematic.hin-record.v1"
    assert record["molecule_id"] == 7
    assert record["title"] == "aromatic ring — source fields"
    assert [a["source_id"] for a in record["atoms"]] == [1, 2, 3, 4, 5, 6]
    assert record["atoms"][0]["coordinates"] == [1.2345678901234567, 0.0, 0.125]
    assert record["atoms"][0]["partial_charge"] == -0.75
    assert all(b["bond_type"] == "a" for b in record["bonds"])
    assert record["graph_projection"]["formal_charges"] == "default_zero_not_inferred_from_partial_charges"
    assert record["source"] == SOURCE


@pytest.mark.parametrize("source", [SOURCE, SOURCE.replace("\n", "\r\n"), SOURCE.rstrip()])
def test_exact_source_retention(source):
    assert chematic.roundtrip_hin_record(source) == source


@pytest.mark.parametrize("bad", [
    "", SOURCE + SOURCE, SOURCE.replace("endmol 7", "endmol 8"),
    SOURCE.replace("-0.75", "NaN"), SOURCE.replace("0.125 2", "inf 2"),
    SOURCE.replace("6 a 2 a", "6 s 2 a"),
    SOURCE.replace("6 a 2 a", "6 a 99 a"),
    SOURCE.replace("6 a 2 a", "6 a 6 a"),
    SOURCE.replace("6 a 2 a", "6 q 2 a"),
    SOURCE.replace("endmol 7", "res 1 ALA 1\nendmol 7"),
    pytest.param(" " * ((1 << 20) + 1), id="input-too-large"),
])
def test_invalid_or_unsupported_input_refuses(bad):
    for reader in [chematic.from_hin, chematic.parse_hin_record_json, chematic.roundtrip_hin_record]:
        with pytest.raises(ValueError):
            reader(bad)
