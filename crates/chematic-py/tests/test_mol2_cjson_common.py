import json
from pathlib import Path
import struct
import pytest
import chematic
PACKET=json.loads((Path(__file__).parents[3]/'validation/mol2_cjson_common/cases.json').read_text(encoding='utf-8'))
def bits(values):return [struct.pack('>d',float(v)) for v in values]
@pytest.mark.parametrize('case',PACKET['cases'],ids=lambda c:c['id'])
@pytest.mark.parametrize('layout',PACKET['variants'])
def test_selected_projection_or_refusal(case,layout):
    text=case['input'].replace('\n','\r\n') if layout=='crlf' else case['input'].rstrip('\n') if layout=='no-final-newline' else case['input']
    calls=[lambda:chematic.project_mol2_to_cjson(text.encode()),lambda:chematic.convert_format(text,'.MOL2','cjson')]
    for call in calls:
        if case['outcome']!='preserved':
            with pytest.raises(ValueError):call()
            continue
        value=json.loads(call());e=case['expected'];a=value['atoms']
        assert a['elements']['number']==e['numbers']
        assert a.get('formalCharges',[0]*len(e['numbers']))==e['formal_charges']
        assert value['bonds']['connections']['index']==[0,1]
        assert value['bonds']['order']==[1]
        assert bits(a['coords']['3d'])==bits(v for row in e['coordinates'] for v in row)
        assert value['properties']['name']==e['name']
        if e['method'] is None:assert 'partialCharges' not in value
        else:
            assert list(value['partialCharges'])==[e['method']]
            assert bits(value['partialCharges'][e['method']])==bits(e['partial_charges'])
@pytest.mark.parametrize('kwargs',[dict(coords=[]),dict(charges=[]),dict(name='new'),dict(comment='new')])
def test_overrides_fail_instead_of_discarding_source_fields(kwargs):
    with pytest.raises(ValueError,match='does not accept overrides'):
        chematic.convert_format(PACKET['cases'][0]['input'],'mol2','cjson',**kwargs)
def test_common_smiles_mol2_cjson_path_remains_usable():
    source=chematic.convert_format('CCO','smiles','mol2')
    value=json.loads(chematic.convert_format(source,'mol2','cjson'))
    assert 'partialCharges' not in value
    assert value['atoms']['coords']['3d']==[0.0]*9
    assert chematic.convert_format(json.dumps(value),'cjson','smiles')
def test_declared_but_missing_values_are_visible_in_explicit_loss_report():
    source=next(c['input'] for c in PACKET['cases'] if c['id']=='all_missing')
    report=chematic.convert_with_report(source,'mol2','cjson',policy='allow_loss')
    assert report['output'] is not None
    assert any(o['field']=='partial_charges' and o['status']=='unverified' for o in report['observations'])
def test_utf8_and_input_bytes_are_bounded():
    with pytest.raises(UnicodeDecodeError):chematic.project_mol2_to_cjson(b'\xff')
    with pytest.raises(ValueError):chematic.project_mol2_to_cjson(' '*(16*1024*1024+1))
