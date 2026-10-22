"""Indexed formal charge is independent of ATOM partial charge columns."""
import hashlib
import json
from pathlib import Path
import struct

import pytest
import chematic

PACKET = json.loads((Path(__file__).parents[3] / 'validation/cjson_mol2_formal_charges/cases.json').read_text(encoding="utf-8"))


@pytest.mark.parametrize('case', PACKET['cases'], ids=lambda c:c['id'])
@pytest.mark.parametrize('newline', ['\n', '\r\n'])
def test_charged_atoms_preserve_graph_partial_charges_and_coordinates(case, newline):
    assert hashlib.sha256(case['input'].encode()).hexdigest() == case['sha256']
    source = case['input'].replace('\n', newline)
    report = chematic.convert_with_report(source, 'cjson', 'mol2', policy='allow_loss')
    assert report['status'] == 'completed_with_loss'
    assert report['candidate_graph_verified'] is True
    statuses = {o['field']:o['status'] for o in report['observations']}
    assert statuses['formal_charges'] == statuses['partial_charges'] == 'preserved'
    after = json.loads(chematic.parse_mol2_record_json(report['output']))
    expected = case['expected']
    assert after['formal_charges'] == expected['formal_charges']
    actual = [a['partial_charge'] for a in after['atoms']] + [v for point in after['coordinates'] for v in point]
    wanted = expected['partial_charges'] + [v for point in expected['coords'] for v in point]
    assert [struct.pack('>d',v) for v in actual] == [struct.pack('>d',v) for v in wanted]
    if any(expected['formal_charges']):
        assert report['output'].index('@<TRIPOS>UNITY_ATOM_ATTR') < report['output'].index('@<TRIPOS>BOND')
    else:
        assert '@<TRIPOS>UNITY_ATOM_ATTR' not in report['output']
    strict = chematic.convert_with_report(source, 'cjson', 'mol2')
    assert strict['status'] == 'refused' and strict['output'] is None
    assert strict['candidate_graph_verified'] is True


def test_missing_partial_dataset_does_not_fill_charge_columns():
    value = json.loads(PACKET['cases'][0]['input'])
    del value['partialCharges']
    report = chematic.convert_with_report(json.dumps(value), 'cjson', 'mol2', policy='allow_loss')
    after = json.loads(chematic.parse_mol2_record_json(report['output']))
    assert after['formal_charges'] == [-1]
    assert after['charge_type'] == 'NO_CHARGES'
    assert after['atoms'][0]['partial_charge'] is None


def test_isotope_loss_still_blocks_verified_mapping():
    value = json.loads(PACKET['cases'][0]['input'])
    value['atoms']['isotopes'] = [18]
    report = chematic.convert_with_report(json.dumps(value), 'cjson', 'mol2', policy='allow_loss')
    assert report['candidate_graph_verified'] is False
    assert report['candidate_atom_mapping'] == []
    statuses = {o['field']:o['status'] for o in report['observations']}
    assert statuses['isotopes'] == 'dropped'
    assert statuses['formal_charges'] == 'unverified'
    assert json.loads(chematic.parse_mol2_record_json(report['output']))['formal_charges'] == [-1]
