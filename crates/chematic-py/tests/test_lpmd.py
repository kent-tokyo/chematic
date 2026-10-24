"""Original synthetic LPMD contract: all frames and current edited fields."""
import copy
import json
from pathlib import Path

import pytest
import chematic

PACKET = json.loads((Path(__file__).resolve().parents[3] / 'validation/lpmd_preservation/cases.json').read_text(encoding='utf-8'))

@pytest.mark.parametrize('case', PACKET['valid'], ids=lambda c: c['id'])
def test_roundtrip(case):
    record = chematic.parse_lpmd(case['text'].encode('utf-8'))
    for _ in range(3):
        assert record == case['record']
        assert isinstance(chematic.lpmd_positions_angstrom(record), list)
        record = chematic.parse_lpmd(chematic.write_lpmd(record))
    assert record == case['record']

@pytest.mark.parametrize('case', PACKET['refused'], ids=lambda c: c['id'])
def test_refused(case):
    with pytest.raises(ValueError, match='LPMD line '):
        chematic.parse_lpmd(case['text'])

@pytest.mark.parametrize('case', PACKET['projection_refused'], ids=lambda c: c['id'])
def test_unsafe_projection_keeps_raw_record(case):
    record = chematic.parse_lpmd(case['text'])
    assert chematic.parse_lpmd(chematic.write_lpmd(record)) == case['record']
    with pytest.raises(ValueError, match='LPMD line '):
        chematic.lpmd_positions_angstrom(record)

def test_edits_limits_and_json_schema():
    record = copy.deepcopy(PACKET['valid'][0]['record'])
    record['frames'][0]['cell'][0][0] = '10'
    record['frames'][0]['atoms'][0][1] = '.5'
    assert chematic.lpmd_positions_angstrom(chematic.parse_lpmd(chematic.write_lpmd(record)))[0][0] == [5.,2.,8.]
    for kwargs in ({'max_atoms':1},{'max_frames':0},{'max_line_len':10},{'max_input_bytes':10}):
        with pytest.raises(ValueError): chematic.parse_lpmd(PACKET['valid'][0]['text'], **kwargs)
        with pytest.raises(ValueError): chematic.write_lpmd(record, **kwargs)
    for mutation in (lambda r:r.update(unknown=1),lambda r:r['frames'][0].update(extra=1),lambda r:r['frames'][0]['atoms'][0].__setitem__(1,.5)):
        bad=copy.deepcopy(record); mutation(bad)
        with pytest.raises(ValueError): chematic.write_lpmd(bad)
    with pytest.raises(UnicodeDecodeError): chematic.parse_lpmd(b'\xff')
    with pytest.raises(TypeError): chematic.parse_lpmd(1)
