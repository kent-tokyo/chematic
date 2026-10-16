"""All-record, source-preserving MOL2/RMSD contracts (Open Babel #2671)."""
import json
from pathlib import Path

import chematic
import pytest

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "validation/openbabel_issues/mol2-rmsd"
REFERENCES = (FIXTURES / "references.mol2").read_text()
PROBE = (FIXTURES / "probe.mol2").read_text()


def records(text, **kwargs):
    return json.loads(chematic.parse_mol2_records_json(text, **kwargs))["records"]


def rmsds(reference=REFERENCES, probe=PROBE, **kwargs):
    return json.loads(chematic.mol2_reference_rmsds_json(reference, probe, **kwargs))


def test_every_reference_has_its_own_rmsd_identity_and_position():
    report = rmsds()
    assert report["schema"] == "chematic.mol2-reference-rmsds.v1"
    assert report["heavy_atoms_only"] is True
    assert report["alignment"] == "none"
    assert report["correspondence"] == "source_id_name_type_and_indexed_graph"
    rows = report["results"]
    assert [r["reference_index"] for r in rows] == list(range(5))
    assert [r["reference_name"] for r in rows] == [f"reference_{i}" for i in range(5)]
    assert [r["probe_name"] for r in rows] == ["probe"] * 5
    assert [r["heavy_atom_count"] for r in rows] == [3] * 5
    assert [r["rmsd"] for r in rows] == [0, 1, 2, 3, 5]


def test_alignment_is_explicit_and_applied_to_every_record():
    report = rmsds(align=True)
    assert report["alignment"] == "proper_rotation_translation"
    assert [r["rmsd"] for r in report["results"]] == pytest.approx([0]*5, abs=1e-12)


def test_records_preserve_source_charges_types_and_opaque_sections():
    parsed = records(REFERENCES)
    assert len(parsed) == 5
    assert "".join(r["source"] for r in parsed) == REFERENCES
    assert parsed[4]["coordinates"][0] == [3, 4, 0]
    for record in parsed:
        assert record["atoms"][0]["source_id"] == 10
        assert record["atoms"][2]["partial_charge"] == -0.3
        assert record["atoms"][2]["atom_type"] == "O.3"
        assert record["opaque_sections"] == ["SUBSTRUCTURE", "VENDOR_DATA"]


def test_preamble_crlf_and_final_line_are_not_normalized():
    text = "# preamble\r\n" + REFERENCES.rstrip().replace("\n", "\r\n")
    assert "".join(r["source"] for r in records(text)) == text


def test_single_record_apis_keep_their_previous_contract():
    with pytest.raises(ValueError):
        chematic.from_mol2(REFERENCES)
    with pytest.raises(ValueError):
        chematic.parse_mol2_record_json(REFERENCES)


@pytest.mark.parametrize("kwargs", [
    {"max_input_bytes": 1}, {"max_records": 4}, {"max_record_bytes": 1},
    {"max_atoms": 3}, {"max_bonds": 2},
])
def test_inspection_limits_raise_instead_of_returning_partial_results(kwargs):
    with pytest.raises(ValueError, match="limit"):
        records(REFERENCES, **kwargs)


@pytest.mark.parametrize("kwargs", [
    {"max_input_bytes": 1}, {"max_records": 4}, {"max_record_bytes": 1},
    {"max_atoms": 3}, {"max_bonds": 2},
])
def test_rmsd_limits_cover_reference_and_probe(kwargs):
    with pytest.raises(ValueError, match="limit"):
        rmsds(**kwargs)


def test_byte_and_record_limits_accept_the_exact_boundary():
    assert len(records(PROBE, max_input_bytes=len(PROBE.encode()),
                       max_record_bytes=len(PROBE.encode()), max_records=1)) == 1
    assert rmsds(PROBE, max_records=1, max_input_bytes=len(PROBE.encode()),
                 max_record_bytes=len(PROBE.encode()))["results"][0]["rmsd"] == 0


@pytest.mark.parametrize("changed", [
    PROBE.replace("10 C1", "10 other"),
    PROBE.replace("1 10 20 1", "1 10 30 1"),
    PROBE.replace("2 20 30 1", "2 20 30 2"),
    PROBE.replace("C.3", "C.2"),
    PROBE + "@<TRIPOS>UNITY_ATOM_ATTR\n10 1\ncharge 1\n",
])
def test_identity_mismatch_at_second_record_is_not_silently_skipped(changed):
    with pytest.raises(ValueError, match="reference 1.*identity/order"):
        rmsds(PROBE + changed)


def test_malformed_middle_record_has_an_indexed_error():
    bad = PROBE.replace("0.000000", "NaN")
    with pytest.raises(ValueError, match="record 1"):
        records(PROBE + bad + PROBE)
    with pytest.raises(ValueError, match="record 1"):
        rmsds(PROBE + bad + PROBE)


def test_empty_input_is_empty_inspection_but_not_a_valid_rmsd_reference():
    assert records("") == []
    with pytest.raises(ValueError, match="no records"):
        rmsds("")
    with pytest.raises(ValueError):
        records("not MOL2")


def test_probe_must_be_exactly_one_record():
    with pytest.raises(ValueError, match="probe"):
        rmsds(PROBE, REFERENCES)


@pytest.mark.parametrize("value", [None, b"MOL2", 1, [], Path("test.mol2")])
def test_text_inputs_do_not_implicitly_convert_other_types(value):
    with pytest.raises(TypeError):
        records(value)
    with pytest.raises(TypeError):
        rmsds(value)
    with pytest.raises(TypeError):
        rmsds(PROBE, value)


def test_hydrogen_motion_is_ignored_but_hydrogen_identity_is_checked():
    changed = PROBE.replace("0.000000 0.000000 1.000000 H", "100 200 300 H")
    assert rmsds(changed)["results"][0]["rmsd"] == 0
    with pytest.raises(ValueError, match="identity/order"):
        rmsds(PROBE.replace("40 H4", "40 H5"))


@pytest.mark.parametrize("align", [False, True])
def test_nonfinite_computation_never_becomes_json_null(align):
    with pytest.raises(ValueError, match="not finite"):
        rmsds(PROBE.replace("0.000000", "1e308"), align=align)
