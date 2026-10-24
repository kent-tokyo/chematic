"""Original unixyz profile: actual native extension, editable token records."""
import copy
import json
from pathlib import Path
import pytest
import chematic

PACKET = json.loads((Path(__file__).resolve().parents[3] / 'validation/unixyz_preservation/cases.json').read_text(encoding="utf-8"))

@pytest.mark.parametrize('case', PACKET['valid'], ids=lambda c: c['id'])
def test_roundtrip(case):
    record = chematic.parse_unixyz(case['text'].encode())
    assert record == case['expected']
    for _ in range(3):
        record = chematic.parse_unixyz(chematic.write_unixyz(record))
        assert record == case['expected']

@pytest.mark.parametrize('case', PACKET['refused'], ids=lambda c: c['id'])
def test_refused(case):
    with pytest.raises(ValueError, match='^UniChem XYZ line '):
        chematic.parse_unixyz(case['text'])

def test_edits_fields_limits_and_strict_boundary():
    record = copy.deepcopy(PACKET['valid'][0]['expected'])
    record['title'] = 'edited'
    record['atoms'][0]['position_angstrom'][0] = '-1.2345678901234567'
    assert chematic.parse_unixyz(chematic.write_unixyz(record)) == record
    for key in ['max_atoms', 'max_input_bytes', 'max_line_len']:
        with pytest.raises(ValueError):
            chematic.write_unixyz(record, **{key: 0})
        with pytest.raises(ValueError):
            chematic.parse_unixyz(PACKET['valid'][0]['text'], **{key: 0})
    for field, value in [('charge', 1), ('unknown', None)]:
        bad = copy.deepcopy(record)
        bad['atoms'][0][field] = value
        with pytest.raises(ValueError, match='unknown field'):
            chematic.write_unixyz(bad)
    for value in [0.1, None, True, 'NaN', '1e-999']:
        bad = copy.deepcopy(record)
        bad['atoms'][0]['position_angstrom'][0] = value
        with pytest.raises(ValueError):
            chematic.write_unixyz(bad)
    with pytest.raises(TypeError):
        chematic.parse_unixyz(object())
    with pytest.raises(UnicodeDecodeError):
        chematic.parse_unixyz(b'\xff')
