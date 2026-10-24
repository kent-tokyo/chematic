"""The local atom-record MolJSON dialect retains non-node H authority."""
import json
from pathlib import Path

import chematic
import pytest

PACKET = json.loads((Path(__file__).resolve().parents[3] / 'validation/moljson_hydrogens/cases.json').read_text(encoding="utf-8"))


@pytest.mark.parametrize('case', PACKET['cases'], ids=lambda c: c['id'])
def test_hydrogen_authority_survives_direct_and_common_routes(case):
    mol = chematic.from_moljson(case['input'])
    reference = chematic.from_smiles(case['expected_smiles'])
    expected = reference.smiles
    assert chematic.from_moljson(reference.to_moljson()).smiles == expected
    for _ in range(3):
        assert mol.smiles == expected
        encoded = mol.to_moljson()
        assert [a['hydrogens'] for a in json.loads(encoded)['atoms']] == case['expected_hydrogens']
        assert chematic.convert_format(encoded, 'moljson', 'smiles') == expected
        mol = chematic.from_moljson(encoded)


@pytest.mark.parametrize('value', PACKET['invalid_hydrogens'])
def test_present_invalid_hydrogens_are_rejected(value):
    text = json.dumps(dict(atoms=[dict(id='a1', element='C', hydrogens=value)], bonds=[]))
    with pytest.raises(ValueError, match='hydrogens'):
        chematic.from_moljson(text)
    with pytest.raises(ValueError, match='hydrogens'):
        chematic.convert_format(text, 'moljson', 'smiles')
