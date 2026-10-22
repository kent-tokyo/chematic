import json
from pathlib import Path
import struct
import pytest
import chematic

PACKET=json.loads((Path(__file__).resolve().parents[3]/'validation/cjson_checked/cases.json').read_text(encoding="utf-8"))
def compact(text):
    result=[];quoted=False;escaped=False
    for c in text:
        if quoted:
            result.append(c)
            if escaped:escaped=False
            elif c=='\\':escaped=True
            elif c=='"':quoted=False
        elif c=='"':quoted=True;result.append(c)
        elif c not in ' \t\r\n':result.append(c)
    return ''.join(result)
def variants(c):return [c['input'],c['input'].replace('\n','\r\n'),compact(c['input'])]
def common_bridge_graph(text):
    cml=chematic.convert_format(text,'.CJSON','cml')
    return json.loads(chematic.convert_format(cml,'cml','cjson'))
def checked_graph(data,expected):
    atoms=data['atoms'];n=len(expected['numbers']);bonds=data['bonds']
    assert atoms['elements']['number']==expected['numbers']
    assert atoms.get('formalCharges',[0]*n)==expected['charges']
    assert atoms.get('isotopes',[0]*n)==expected['isotopes']
    assert bonds['connections']['index']==expected['indices']
    assert bonds['order']==expected['orders']

@pytest.mark.parametrize('case',PACKET['valid']+PACKET['compatibility'],ids=lambda c:c['id'])
def test_checked_graph_and_common_bridge(case):
    for text in variants(case):
        mol,coords=chematic.from_cjson_checked(text)
        checked_graph(json.loads(mol.to_cjson(coords)),case['expected'])
        checked_graph(common_bridge_graph(text),case['expected'])
        assert chematic.convert_format(text,'.CJSON','cjson')==text

@pytest.mark.parametrize('case',PACKET['invalid'],ids=lambda c:c['id'])
def test_invalid_input_is_refused_before_any_common_output(case):
    for text in variants(case):
        with pytest.raises(ValueError):chematic.from_cjson_checked(text)
        for target in ['smiles','mol','mol_v3000','mol2','cml','cjson','moljson']:
            with pytest.raises(ValueError):chematic.convert_format(text,'cjson',target)

def test_legacy_reader_and_missing_defaults_remain_available():
    text=PACKET['invalid'][0]['input']
    assert chematic.from_cjson(text)[0].atom_count==2
    mol,coords=chematic.from_cjson_checked('{"atoms":{"elements":{"number":[6]}}}')
    assert mol.atom_count==1 and coords==[]

def test_numeric_coordinate_and_byte_limits():
    text='{"atoms":{"elements":{"number":[6]},"coords":{"3d":[-1.4175225513503917,-0.0,5e-324]}}}'
    coords=chematic.from_cjson_checked(text)[1][0]
    assert [struct.pack('>d',v) for v in coords]==[struct.pack('>d',v) for v in [-1.4175225513503917,-0.0,5e-324]]
    with pytest.raises(ValueError,match='input bytes'):chematic.from_cjson_checked(' '*(16*1024*1024+1))
