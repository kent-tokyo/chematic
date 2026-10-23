import json
from pathlib import Path
import struct
import pytest
import chematic

PACKET = json.loads((Path(__file__).resolve().parents[3] / 'validation/cml_3d/cases.json').read_text(encoding="utf-8"))
def variants(case):
    text = case['input']
    return [text, text.replace('\n', '\r\n'), ''.join(line.strip() for line in text.splitlines())]
def bits(values):
    return [struct.pack('>d', value) for row in values for value in row]

@pytest.mark.parametrize('case', PACKET['valid'], ids=lambda c:c['id'])
def test_cartesian_geometry_in_direct_and_common_api(case):
    for text in variants(case):
        for source in [text, text.encode()]:
            mol, coords = chematic.from_cml_with_3d(source)
            assert bits(coords) == bits(case['coords'])
            assert all(isinstance(row, list) for row in coords)
            assert mol.atom_count == len(case['numbers'])
        data = json.loads(chematic.convert_format(text, '.CML', 'cjson'))
        assert bits([data['atoms']['coords']['3d']]) == bits(case['coords'])
        assert data['atoms']['elements']['number'] == case['numbers']
        assert data['bonds']['connections']['index'] == case['bond_indices']
        assert data['bonds']['order'] == case['bond_orders']

@pytest.mark.parametrize('case', PACKET['absent'], ids=lambda c:c['id'])
def test_missing_geometry_is_not_generated(case):
    for text in variants(case):
        assert chematic.from_cml_with_3d(text)[1] is None
        assert 'coords' not in json.loads(chematic.convert_format(text,'cml','cjson'))['atoms']

@pytest.mark.parametrize('case', PACKET['invalid'], ids=lambda c:c['id'])
def test_invalid_input_never_returns_candidate(case):
    for text in variants(case):
        with pytest.raises(ValueError): chematic.from_cml_with_3d(text)
        with pytest.raises(ValueError): chematic.convert_format(text,'cml','cjson')

def test_override_and_utf8_type_boundaries():
    source = PACKET['valid'][3]['input']
    coords = [[-1.25,2.5,3.75]]
    assert json.loads(chematic.convert_format(source,'cml','cjson',coords=coords))['atoms']['coords']['3d'] == coords[0]
    with pytest.raises(UnicodeDecodeError): chematic.from_cml_with_3d(b'\xff')
    with pytest.raises(ValueError,match='input bytes'): chematic.from_cml_with_3d(' '*(16*1024*1024+1))

def test_graph_charge_and_isotope_are_not_changed():
    data = json.loads(chematic.convert_format(PACKET['valid'][0]['input'],'cml','cjson'))
    assert data['atoms']['formalCharges'] == [0,0,1]
    assert data['atoms']['isotopes'] == [0,13,0]
    # Existing graph-only API and its return type remain available.
    assert chematic.from_cml(PACKET['valid'][0]['input']).atom_count == 3
